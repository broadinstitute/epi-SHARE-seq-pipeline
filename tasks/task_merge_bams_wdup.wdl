version 1.0

# TASK
# task-merge-filtered-bams-wdups

task atac_merge_filtered_bams_wdups {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: merge the individual bams together'
    }

    input {
        # This task takes in input the preprocessed ATAC fastqs and align them to the genome.
        Array[File] wdups_bams
        String genome_name
        String prefix = "sample-share"
        String? barcode_tag = "CB" # "XC"
        Int? multimappers # = 5
        Int? cpus = 16
        Float? disk_factor = 8.0
        Float? memory_factor = 0.15
        String? docker_image = "us.gcr.io/buenrostro-share-seq/share_task_merge_bams"
        String? singularity_image = "docker://us.gcr.io/buenrostro-share-seq/share_task_merge_bams"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(wdups_bams, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 16.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    # Determining memory for samtools.
    Float samtools_memory_gb = 0.8 * mem_gb # Samtools has overheads so reducing the memory to 80% of the total.

    # Number of threads to beable to use 4GB of memory per thread seems to be the fastest way
    Int samtools_threads_ = floor(samtools_memory_gb / 4)
    Int samtools_threads =  if samtools_threads_ == 0 then 1 else samtools_threads_

    Int sambamba_threads = floor(cpus/2)

    # Now that we know how many threads we can use to assure 4GB of memory per thread
    # we assign any remaining memory to the threads.
    Int samtools_memory_per_thread_ = floor(samtools_memory_gb * 1024 / samtools_threads) # Computing the memory per thread for samtools in MB.
    Int samtools_memory_per_thread = if samtools_memory_per_thread_ < 768 then 768 else samtools_memory_per_thread_

    # Tim parameters
    Int machine_mem_mb = 18150
    Int cpu = 1
    Int compression_level = 5
    # default to 500GiB of space
    Int disk = 500
    Int command_mem_mb = machine_mem_mb - 500

    # Define tmp file name
    String unsorted_bam = "${prefix}.atac.merge.${genome_name}.bam"

    # Define the output names
    String merged_bam = "${prefix}.atac.filtered.wdup.k${multimappers}.${genome_name}.sorted.bam"
    String merged_bai = "${prefix}.atac.filtered.wdup.k${multimappers}.${genome_name}.sorted.bam.bai"
    String stats_log = '${default="share-seq" prefix}.atac.filtered.wdup.k${multimappers}.${genome_name}.stats.log.txt'
    String duplicate_stats = '${default="share-seq" prefix}.atac.filtered.wdup.qc.duplicate.stats.${genome_name}.tsv'

    String monitor_log = "atac_merge_monitor.log"

    command <<<
        set -e

        bash $(which monitor_script.sh) 2>&1 &

        # Trying picard
        
        java -Dsamjdk.compression_level=~{compression_level} -Xms~{command_mem_mb}m -Xmx~{command_mem_mb}m -jar /usr/local/bin/picard.jar \
        MergeSamFiles \
        USE_THREADING=true \
        ASSUME_SORTED="queryname" \
        SORT_ORDER="queryname" \
        INPUT=~{sep=' INPUT=' wdups_bams} \
        OUTPUT=~{merged_bam}

        sambamba index -t ~{cpus} ~{merged_bam}
        
        # Duplicates per barcode
        echo '------ START: Compute duplication per barcode ------' 1>&2
        time python3 $(which qc_atac_count_duplicates_per_barcode.py) \
            -o ~{duplicate_stats} \
            --bc_tag ~{barcode_tag} \
            ~{merged_bam}

    >>>

    output {
        File atac_wdup_alignment = merged_bam
        File atac_wdup_alignment_index = merged_bai
        File atac_qc_duplicate_stats = duplicate_stats
    }

    runtime {
        cpu: cpu
        docker: "${docker_image}"
        singularity: "${singularity_image}"
        disks: "local-disk ${disk} HDD"
        disk: disk + " GB" # TES
        #disks: "local-disk ${disk_gb} ${disk_type}"
        maxRetries:1
        memory: "${machine_mem_mb} MiB"
        #memory: "${mem_gb} GB"
        memory_retry_multiplier: 2
    }

    parameter_meta {
        wdups_bams: {
                description: 'Individuals bams from the scatter alignment task',
                help: 'Individuals bams from the scatter alignment task',
                example: 'align.raw.L1.bam',
            }
        cpus: {
                description: 'Number of cpus.',
                help: 'Set the number of cpus used by bowtie2',
                default: 16
            }
        disk_factor: {
                description: 'Multiplication factor to determine disk required for task align.',
                help: 'This factor will be multiplied to the size of FASTQs to determine required disk of instance (GCP/AWS) or job (HPCs).',
                default: 8.0
            }
        memory_factor: {
                description: 'Multiplication factor to determine memory required for task align.',
                help: 'This factor will be multiplied to the size of FASTQs to determine required memory of instance (GCP/AWS) or job (HPCs).',
                default: 0.15
        }
        prefix: {
                description: 'Prefix for output files.',
                help: 'Prefix that will be used to name the output files',
                examples: 'my-experiment'
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for the alignment step.',
                example: ["us.gcr.io/buenrostro-share-seq/share_task_bowtie2"]
            }
    }


}
