version 1.0

# TASK
# atac-qc-tss-enrichment

task qc_atac_tss_enrichment {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: ATAC qc tss enrichment'
    }

    input {
        # This function takes in input the raw and filtered bams
        # and compute some alignment metrics along with the TSS
        # enrichment plot.
        File? fragments
        File? fragments_index
        File? tss
        String? genome_name
        String? prefix

        # Runtime
        Int? cpus = 1
        Float? disk_factor = 10.0
        Float? memory_factor = 0.3
        String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_qc_atac"
        String singularity_image = "docker://us.gcr.io/buenrostro-share-seq/share_task_qc_atac"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fragments, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 24.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    String monitor_log = "atac_qc_monitor.log"


    command<<<
        set -e

        # I am not writing to a file anymore because Google keeps track of it automatically.
        bash $(which monitor_script.sh) 1>&2 &

        ln -s ~{fragments} in.fragments.tsv.gz
        ln -s ~{fragments_index} in.fragments.tsv.gz.tbi

        # TSS enrichment stats
        echo '------ START: Compute TSS enrichment ------' 1>&2
        time python3 $(which qc_atac_compute_tss_enrichment.py) \
            -e 2000 \
            --tss ~{tss} \
            --prefix "~{prefix}.atac.qc.~{genome_name}" \
            in.fragments.tsv.gz
    >>>

    output {
        File atac_qc_tss_enrichment_barcode_stats = "${prefix}.atac.qc.${genome_name}.tss_enrichment_barcode_stats.tsv"
        File atac_qc_tss_enrichment_plot = "${prefix}.atac.qc.${genome_name}.tss_enrichment_bulk.png"
        File atac_qc_tss_enrichment_score_bulk = "${prefix}.atac.qc.${genome_name}.tss_score_bulk.txt"
    }

    runtime {
        cpu: cpus
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker: "${docker_image}"
        singularity: "${singularity_image}"
        maxRetries: 1
        memory: "${mem_gb} GB"
        memory_retry_multiplier: 2
    }

    parameter_meta {
        fragments: {
                description: 'Fragment File',
                help: 'Fragment File',
                example: 'fragments.tsv.gz'
            }
        tss: {
                description: 'TSS bed file',
                help: 'List of TSS in bed format used for the enrichment plot.',
                example: 'refseq.tss.bed'
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.',
                examples: ['hg38', 'mm10', 'both']
            }
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by bowtie2',
                examples: '4'
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz',
                example: ['put link to gcr or dockerhub']
            }
    }
}
