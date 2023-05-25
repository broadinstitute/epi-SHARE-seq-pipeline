version 1.0

# TASK
# SHARE-atac-merge_bams

task share_atac_merge_bams {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: merge the individual bams together'
    }

    input {
        # This task takes in input the preprocessed ATAC fastqs and align them to the genome.
        Array[File] bams
        Array[File] logs
        String genome_name
        String prefix = "sample-share"
        Int? multimappers # = 5
        Int? cpus = 16
        Float? disk_factor = 8.0
        Float? memory_factor = 0.15
        String? docker_image = "us.gcr.io/buenrostro-share-seq/share_task_merge_bams"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(bams, "G")

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


    # Define tmp file name
    String unsorted_bam = "${prefix}.atac.merge.${genome_name}.bam"

    # Define the output names
    String merged_bam = "${prefix}.atac.merged.k${multimappers}.${genome_name}.sorted.bam"
    String merged_bai = "${prefix}.atac.merged.k${multimappers}.${genome_name}.sorted.bam.bai"
    String alignment_log = "${prefix}.atac.merged.k${multimappers}.${genome_name}.log"

    String monitor_log = "atac_merge_monitor.log"

    command <<<
        set -e

        bash $(which monitor_script.sh)  2>&1 &

        sambamba merge -t ~{cpus} ~{unsorted_bam} ~{sep=" " bams}

        sambamba sort -t ~{samtools_threads} -m ~{samtools_memory_per_thread}M -o ~{merged_bam} unsorted
        
        sambamba index -t ~{cpus} ~{merged_bam}

        awk '{ sum[FNR%15]+=$1 } END { for(idx in sum) if (idx>0) {print sum[idx]}}' ~{sep=" " logs} > ~{alignment_log}

    >>>

    output {
        File atac_merged_alignment = merged_bam
        File atac_merged_alignment_index = merged_bai
        File atac_merged_alignment_log = alignment_log
        File atac_merge_monitor_log = monitor_log
    }

    runtime {
        cpu: cpus
        docker: "${docker_image}"
        disks: "local-disk ${disk_gb} ${disk_type}"
        maxRetries:1
        memory: "${mem_gb} GB"
        memory_retry_multiplier: 2
    }

    parameter_meta {
        bams: {
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
