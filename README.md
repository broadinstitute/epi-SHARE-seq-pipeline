# epi-SHARE-seq-pipeline
Epigenomics Program pipeline to analyze SHARE-seq data.

Pipeline specifications can be found [here](https://docs.google.com/document/d/1J-NWpDLkEGLsLjVe6h6-Rx4nxzTdgy1TJZvuMnYiiyg/edit?usp=sharing).

Pipeline main page on [dockstore](https://dockstore.org/workflows/github.com/broadinstitute/epi-SHARE-seq-pipeline/SHARE-seq:release?tab=info).

- Primary workflow is [share-seq.wdl](share_seq.wdl)
- The workflows directory contains sub workflows for preprocessing and for running specifically workflows separately (e.g. ATAC, RNA, DORCS)
- The tasks directory contains the tasks called from the workflows
- The src directory contains python and R code called within tasks
- The dockerfiles directory contains dockerfiles needed for running the workflows

