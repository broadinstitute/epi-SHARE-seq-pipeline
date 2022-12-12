# epi-SHARE-seq-pipeline: release branch
Epigenomics Program pipeline to analyze SHARE-seq data.

- Primary workflow is [share-seq.wdl](share_seq.wdl)
- The workflows directory contains sub workflows for preprocessing and for running specifically workflows separately (e.g. ATAC, RNA, DORCS)
- The tasks directory contains the tasks called from the workflows
- The src directory contains python and R code called within tasks
- The dockerfiles directory contains dockerfiles needed for running the workflows
