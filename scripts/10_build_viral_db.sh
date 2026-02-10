#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
source "${SCRIPT_DIR}/lib/common.sh"

usage() {
  cat <<'USAGE'
Uso:
  bash scripts/10_build_viral_db.sh --target <key> [--query "<ncbi query>" | --taxid <id>]

Saída:
  db/<target>/
    sequences.fasta
    blastdb/
    metadata.json
USAGE
}

TARGET=""
QUERY=""
TAXID=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --target) TARGET="${2:-}"; shift 2 ;;
    --query) QUERY="${2:-}"; shift 2 ;;
    --taxid) TAXID="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERRO] argumento inválido: $1"; usage; exit 2 ;;
  esac
done

[[ -n "$TARGET" ]] || log_error "Parâmetro obrigatório ausente: --target"
if [[ -n "$QUERY" && -n "$TAXID" ]]; then
  log_error "Use apenas um entre --query e --taxid"
fi

command -v esearch >/dev/null 2>&1 || log_error "EDirect não encontrado (esearch)."
command -v efetch >/dev/null 2>&1 || log_error "EDirect não encontrado (efetch)."
command -v makeblastdb >/dev/null 2>&1 || log_error "makeblastdb não encontrado (blast+)."

TARGETS_FILE="$REPO_ROOT/config/targets.json"
[[ -f "$TARGETS_FILE" ]] || log_error "Catálogo de alvos não encontrado: config/targets.json"

TARGET_JSON="$(python3 - "$TARGETS_FILE" "$TARGET" <<'PY'
import json,sys
path,key=sys.argv[1],sys.argv[2]
with open(path,encoding='utf-8') as fh:
    data=json.load(fh)
for item in data:
    if item.get('key')==key:
        print(json.dumps(item,ensure_ascii=False))
        break
PY
)"
[[ -n "$TARGET_JSON" ]] || log_error "Target '$TARGET' não encontrado em config/targets.json"

DEFAULT_QUERY="$(python3 - "$TARGET_JSON" <<'PY'
import json,sys
obj=json.loads(sys.argv[1])
print(obj.get('query_default') or '')
PY
)"
DEFAULT_TAXID="$(python3 - "$TARGET_JSON" <<'PY'
import json,sys
obj=json.loads(sys.argv[1])
print(obj.get('taxid') or '')
PY
)"

if [[ -z "$QUERY" && -n "$TAXID" ]]; then
  QUERY="txid${TAXID}[Organism:exp]"
fi
if [[ -z "$QUERY" && -z "$TAXID" && -n "$DEFAULT_TAXID" ]]; then
  TAXID="$DEFAULT_TAXID"
  QUERY="txid${TAXID}[Organism:exp]"
fi
if [[ -z "$QUERY" ]]; then
  QUERY="$DEFAULT_QUERY"
fi
[[ -n "$QUERY" ]] || log_error "Nenhuma query definida para o alvo '$TARGET'."

OUT_DIR="$REPO_ROOT/db/$TARGET"
FASTA="$OUT_DIR/sequences.fasta"
BLAST_DIR="$OUT_DIR/blastdb"
BLAST_PREFIX="$BLAST_DIR/$TARGET"
META="$OUT_DIR/metadata.json"

mkdir -p "$BLAST_DIR"

log_info "[DB] target=$TARGET"
log_info "[DB] query=$QUERY"
log_info "[DB] baixando sequências do NCBI..."

esearch -db nucleotide -query "$QUERY" | efetch -format fasta > "$FASTA"

SEQ_COUNT="$(grep -c '^>' "$FASTA" || true)"
if [[ "$SEQ_COUNT" -eq 0 ]]; then
  log_error "NCBI query retornou 0 sequências"
fi

log_info "[DB] sequências recuperadas: $SEQ_COUNT"
log_info "[DB] construindo BLAST DB..."
makeblastdb -in "$FASTA" -dbtype nucl -out "$BLAST_PREFIX" -parse_seqids >/dev/null

EDIRECT_VERSION="$(esearch -version 2>/dev/null | head -n1 || true)"
[[ -n "$EDIRECT_VERSION" ]] || EDIRECT_VERSION="unknown"

python3 - "$META" "$TARGET_JSON" "$QUERY" "$TAXID" "$SEQ_COUNT" "$EDIRECT_VERSION" <<'PY'
import json,sys
meta_path,target_json,query,taxid,seq_count,version = sys.argv[1:7]
target=json.loads(target_json)
metadata={
    "target": target.get("key"),
    "display_name": target.get("display_name"),
    "query": query,
    "taxid": taxid or None,
    "date": __import__('datetime').datetime.utcnow().isoformat(timespec='seconds') + 'Z',
    "sequence_count": int(seq_count),
    "source": "NCBI Nucleotide (EDirect)",
    "version": version,
    "target_catalog": target,
}
with open(meta_path,'w',encoding='utf-8') as fh:
    json.dump(metadata,fh,indent=2,ensure_ascii=False)
PY

log_info "[DB] pronto: $OUT_DIR"
