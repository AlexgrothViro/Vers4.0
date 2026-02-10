# UX Review (Dashboard + Pipeline Scripts)

| Falha | Tipo | Severidade | Impacto | Recomendação prática |
|---|---|---|---|---|
| DB hardcoded em Teschovirus (query fixa) | robustez | alta | Bloqueia reutilização para outros vírus e dificulta manutenção. | Usar catálogo versionado + script genérico por `target/query/taxid` (implementado). |
| Entrada de amostra exigia caminho manual completo | UX | alta | Alto erro operacional (paths inválidos, formato Windows). | Oferecer seleção por `sample_id` e upload/import R1/R2 ou `.zip` (implementado). |
| Falhas com mensagens técnicas/ambíguas | UX | média | Usuário não técnico perde tempo para diagnosticar. | Padronizar erros amigáveis: `R2 não encontrado`, `arquivo vazio`, `NCBI query retornou 0 sequências` (implementado). |
| Sem metadados rastreáveis do DB gerado | robustez | média | Perda de reprodutibilidade e auditoria do banco usado. | Salvar `metadata.json` com query/taxid, data, contagem e fonte/versão (implementado). |
| Pipeline UI com pouca orientação de importação | UX | média | Dúvida entre colar caminho e upload, gerando retrabalho. | Separar fluxos em formulários dedicados (caminho vs upload) e manter logs em tempo real (implementado). |
| Ausência de launcher Windows único para UX | UX | média | Fricção para iniciar WSL + dashboard no contexto do usuário final. | Adicionar `bundle/start_platform.bat` validando WSL/distro e abrindo navegador (implementado). |
