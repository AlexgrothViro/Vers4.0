# Picornavirus Pipeline — WSL Bundle

## Windows (mais fácil)
1) Instale o WSL2 + Ubuntu (se ainda não tiver)
2) Extraia o bundle em uma pasta (ex.: Downloads)
3) Dê duplo clique em: `bundle/run.bat`

Exemplos:
- Smoke test:
  - `bundle/run.bat smoke-test`
- Rodar pipeline:
  - `bundle/run.bat pipeline SAMPLE=81554_S150`

## Dentro do WSL (Linux)
- `bash bundle/run.sh smoke-test`
- `bash bundle/run.sh pipeline SAMPLE=81554_S150`

Na primeira execução ele baixa o micromamba e cria o ambiente local em `.bundle/`.
