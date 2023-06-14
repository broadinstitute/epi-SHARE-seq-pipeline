version 1.0

# TASK
# atac-qc-reads-in-peaks
task qc_atac_reads_in_peaks {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: ATAC qc reads in peaks'
    }

    input {
        # This function takes in input the raw and filtered bams
        # and compute some alignment metrics along with the TSS
        # enrichment plot.
        File? fragments
        File? fragments_index
        File? peaks
        String? genome_name
        String? prefix

        # Runtime
        Int? cpus = 8
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

        # Fragments in peaks
        # "~{prefix}.reads.in.peak.tsv"
        echo '------ START: Compute fragments in peaks------' 1>&2
        time python3 $(which qc_atac_compute_reads_in_peaks.py) \
            --peaks ~{peaks} \
            --prefix "~{prefix}.atac.qc.~{genome_name}" \
            in.fragments.tsv.gz
    >>>

    output {
        File atac_qc_fragments_in_peaks = "${prefix}.atac.qc.${genome_name}.reads.in.peak.tsv"
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
        peaks: {
                description: 'Peaks bed file',
                help: 'List of peaks in bed format used for the enrichment plot.',
                example: 'ccREs.bed'
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
