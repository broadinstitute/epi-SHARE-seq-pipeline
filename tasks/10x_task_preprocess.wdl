version 1.0

# TASK
# 10x_task_preprocess

task preprocess_tenx {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: preprocess 10x ATAC data.'
    }

    input {
        # This task takes in input the 3 fastqs coming out from cellranger mkfastqs and preprocess them.
        File fastq_R1 # Pair 1 reads
        File fastq_R3 # Pair 2 reads
        File fastq_R2 # Barcode fastq
        File? whitelist # Barcode whitelist (chemistry specific)
        Int? barcode_dist = 2
        Float? threshold_pct_barcode_matching = 0.60
        String chemistry
        String? prefix
        Int? cpus = 16
        Float? disk_factor = 8.0
        Float? memory_factor = 0.15
        #String docker_image = "us.gcr.io/buenrostro-share-seq/10x_task_preprocessing"
        String? docker_image = "polumechanos/10x_task_preprocess"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fastq_R1, "G") + size(fastq_R2, "G") + size(fastq_R3, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 5.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    # auto-detect barcode complementation outfiles
    String barcode_complementation_qc = "${default="10x" prefix}.atac.preprocess.complementation.qc.txt"
    String barcode_complementation_out = "${default="10x" prefix}.atac.preprocess.complementation.out.txt"

    # barcode correction and filtering outfiles
    String barcode_correction_qc = "${default="10x" prefix}.atac.preprocess.barcode.correction.qc.txt"
    String cleaned_fastq_R1 = "${default="10x" prefix}.atac.preprocess.cleaned.R1.fastq.gz"
    String cleaned_fastq_R2 = "${default="10x" prefix}.atac.preprocess.cleaned.R2.fastq.gz"

    # read trimming outfiles
    String final_fastq_R1 = "${default="10x" prefix}.atac.preprocess.cleaned.trimmed.R1.fastq.gz"
    String final_fastq_R2 = "${default="10x" prefix}.atac.preprocess.cleaned.trimmed.R2.fastq.gz"
    String trimming_log_json = "${default="10x" prefix}.atac.preprocess.trimming.log.json"
    String trimming_log_html = "${default="10x" prefix}.atac.preprocess.trimming.log.html"
    String trimming_stats = "${default="10x" prefix}.atac.preprocess.trimming.adapter.stats.txt"

    String barcode_conversion_dict = "barcode_conversion_dict.csv"

    String monitor_log = 'monitor_10x_preprocessing.log.txt'

    command <<<
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &

        # Strip read description
        zcat ~{fastq_R1} | sed 's/ .*//' | gzip > stripped_R1.fastq.gz
        zcat ~{fastq_R3} | sed 's/ .*//' | gzip > stripped_R2.fastq.gz
        zcat ~{fastq_R2} | sed 's/ .*//' | gzip > stripped_barcode.fastq.gz

        if [[ '~{whitelist}' == *.gz ]]; then
            gunzip -c ~{whitelist} > whitelist.txt
        else
            ln -s ~{whitelist} whitelist.txt
        fi

        # auto-detect barcode complementation
        # python3 barcode_revcomp_detect.py barcode_fastq chemistry whitelist qc_out out threshold

        python3 $(which barcode_revcomp_detect.py) stripped_barcode.fastq.gz ~{chemistry} whitelist.txt ~{barcode_complementation_qc} ~{barcode_complementation_out} ~{threshold_pct_barcode_matching}

        # barcode correction and filtering
        # python3 match_barcodes.py

        python3 $(which match_barcodes.py) stripped_R1.fastq.gz stripped_R2.fastq.gz stripped_barcode.fastq.gz ~{chemistry} ~{barcode_dist} ~{barcode_complementation_out} whitelist.txt ~{cleaned_fastq_R1} ~{cleaned_fastq_R2} ~{barcode_correction_qc} ~{cpus}

        # Cleaned old files
        rm stripped_R1.fastq.gz stripped_R2.fastq.gz stripped_barcode.fastq.gz
    >>>

    output {
        File fastq_R1_preprocessed = cleaned_fastq_R1
        File fastq_R2_preprocessed = cleaned_fastq_R2
        File tenx_barcode_complementation_qc = barcode_complementation_qc
        File tenx_barcode_correction_qc = barcode_correction_qc
        File? tenx_barcode_conversion_dict = barcode_conversion_dict
        #File tenx_trimming_log_json = trimming_log_json
        #File trimming_log_html = trimming_log_html
        #File tenx_trimming_stats = trimming_stats
    }

    runtime {
        cpu: cpus
        docker: "${docker_image}"
        disks: "local-disk ${disk_gb} ${disk_type}"
        maxRetries: 1
        memory: "${mem_gb} GB"
        memory_retry_multiplier: 2
    }

    parameter_meta {
        fastq_R1: {
                description: 'Pairs 1 fastq',
                help: 'Pairs 1 fastq',
            }
        fastq_R2: {
                description: 'Barcode fastq',
                help: 'Barcode fastq',
            }
        fastq_R3: {
                description: 'Pairs 2 fastq',
                help: 'Pairs 2 fastq',
            }
    }

}
