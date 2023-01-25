version 1.0

# TASK
# SHARE-atac-bowtie2

task share_atac_align {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: align ATAC task using bowtie2'
    }

    input {
        # This task takes in input the preprocessed ATAC fastqs and align them to the genome.
        Array[File] fastq_R1
        Array[File] fastq_R2
        Boolean encode_mode = false
        Int? cpus = 16
        Int? multimappers = 1
        Float? disk_factor = 8.0
        Float? memory_factor = 0.15
        String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_bowtie2"
        File genome_index_tar       # This is a tar.gz folder with all the index files.
        String genome_name          # GRCh38, mm10
        String prefix = "sample"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fastq_R1, "G") + size(fastq_R2, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 5.0 + size(genome_index_tar, "G") + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    # Determining memory for samtools.
    Float samtools_memory_gb = 0.8 * mem_gb # Samtools has overheads so reducing the memory to 80% of the total.

    # Number of threads to beable to use 4GB of memory per thread seems to be the fastest way
    Int samtools_threads_ = floor(samtools_memory_gb / 4)
    Int samtools_threads =  if samtools_threads_ == 0 then 1 else samtools_threads_

    # Now that we know how many threads we can use to assure 4GB of memory per thread
    # we assign any remaining memory to the threads.
    Int samtools_memory_per_thread_ = floor(samtools_memory_gb * 1024 / samtools_threads) # Computing the memory per thread for samtools in MB.
    Int samtools_memory_per_thread = if samtools_memory_per_thread_ < 768 then 768 else samtools_memory_per_thread_


    # Define tmp file name
    String unsorted_bam = "${prefix}.atac.align.${genome_name}.bam"

    # Define the output names
    String sorted_bam = "${prefix}.atac.align.k${multimappers}.${genome_name}.sorted.bam"
    String sorted_bai = "${prefix}.atac.align.k${multimappers}.${genome_name}.sorted.bam.bai"
    String alignment_log = "${prefix}.atac.align.k${multimappers}.${genome_name}.log"

    String monitor_log = "atac_align_monitor.log"

    command <<<
        set -e

        bash $(which monitor_script.sh) > ~{monitor_log} 2>&1 &

        tar zxvf ~{genome_index_tar} --no-same-owner -C ./
        genome_prefix=$(basename $(find . -type f -name "*.rev.1.bt2") .rev.1.bt2)

        # Aligning and adding the cell barcode to the CB tag and the barcodes plus pkr in the XC tag.
        bowtie2 -X2000 \
            -p ~{cpus} \
            -k ~{multimappers} \
            --rg-id ~{prefix + "."}atac \
            --rg "SM:None" \
            --rg "LB:None" \
            --rg "PL:Illumina" \
            ~{if encode_mode then "--sam-append-comment" else ""} \
            -x $genome_prefix \
            -1 ~{sep="," fastq_R1} \
            -2 ~{sep="," fastq_R2} 2> ~{alignment_log} | \
            samtools view \
                -b \
                -@ ~{samtools_threads} \
                - \
                -o ~{unsorted_bam}


        if [[ ~{encode_mode} ]] then
            samtools sort \
                -@ ~{samtools_threads} \
                -m ~{samtools_memory_per_thread}M \
                ~{unsorted_bam} \
                -o ~{sorted_bam}
        else
            # Splitting the read name to ge the cell barcode and adding it to the CB tag in the BAM file.
            samtools view ~{unsorted_bam} | \
            awk '{if ($0 ~ /^@/) {print $0} else {split($1,a,"[,_]"); print($0 "\tCB:Z:" a[2]a[3]a[4] "\tXC:Z:" a[2]a[3]a[4] "," a[5]);}}' | \
            samtools sort \
                -@ ~{samtools_threads} \
                -m ~{samtools_memory_per_thread}M \
                - \
                -o ~{sorted_bam}
        fi

        samtools index -@ ~{cpus} ~{sorted_bam}

    >>>

    output {
        File? atac_alignment = sorted_bam
        File? atac_alignment_index = sorted_bai
        File? atac_alignment_log = alignment_log
        File? atac_alignment_monitor_log = monitor_log
    }

    runtime {
        cpu : cpus
        memory : "${mem_gb} GB"
        memory_retry_multiplier: 2
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker : "${docker_image}"
        maxRetries:1
    }

    parameter_meta {
        fastq_R1: {
                description: 'Read1 fastq.',
                help: 'Processed fastq for read1.',
                example: 'input.atac.R1.fq.gz',
            }
        fastq_R2: {
                description: 'Read2 fastq.',
                help: 'Processed fastq for read2.',
                example: 'input.atac.R2.fq.gz'
            }
        multimappers: {
                    description: 'Specifiy the numbers of multimappers allowed.',
                    help: 'This is the integer that will be passed to the -k parameter of bowtie2',
                    example: [5]
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
        genome_index_tar: {
                description: 'Bowtie2 indexes',
                help: 'Index files for bowtie2 to use during alignment in tar format. No subfolders.',
                examples: ['GRCh38.tar.gz', 'mm10.tar.gz']
            }
        genome_name: {
                description: 'Reference name.',
                help: 'The name of the reference genome used by the aligner. This is appended to the output file name.',
                examples: ['GRCh38', 'mm10']
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
