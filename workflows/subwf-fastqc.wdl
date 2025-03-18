version 1.0

import "../tasks/task_fastqc.wdl" as fastqc

workflow wf_fastqc {
    call fastqc.run_fastqc as qc

    output {
        File multiqc_report = qc.multiqc_report
        File multiqc_tar = qc.multiqc_tar
    }
}