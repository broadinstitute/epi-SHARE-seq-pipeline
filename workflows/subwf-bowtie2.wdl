version 1.0

workflow bowtie2 {
    input {
        Array[File] read1
        Array[File] read2
        File whitelist
        String prefix
        String chemistry
        String pkr
        String genome_name
        String genome_index_tar
        Int? align_multimappers
        Boolean correct_fastqs
    }

    if ( correct_fastqs ) {
        scatter (read_pair in zip(read1, read2)) {
            call correct_fastq {
                input:
                    fastq_R1 = read_pair.left,
                    fastq_R2 = read_pair.right,
                    whitelist = whitelist,
                    sample_type = "ATAC",
                    pkr = pkr,
                    prefix = prefix
            }
        }
    }

    scatter (read_pair in zip(select_first([correct_fastq.corrected_fastq_R1, read1]), select_first([correct_fastq.corrected_fastq_R2, read2]))) {
        call trim {
            input:
                fastq_R1 = read_pair.left,
                fastq_R2 = read_pair.right,
                chemistry = chemistry
        }
    }

    scatter(read_pair in zip(trim.fastq_R1_trimmed, trim.fastq_R2_trimmed)) {
        call bowtie2 {
            input:
                fastq_R1 = [read_pair.left],
                fastq_R2 = [read_pair.right],
                chemistry= chemistry,
                genome_name = genome_name,
                genome_index_tar = genome_index_tar,
                multimappers = align_multimappers,
                prefix = prefix
        }
    }
    
    call merge_bams {
        input:
            bams = bowtie2.atac_alignment,
            logs = bowtie2.atac_alignment_log,
            multimappers = align_multimappers,
            genome_name = genome_name,
            prefix = prefix
    }

    output{
    	File atac_mapped_and_unmapped_bam = merge_bams.atac_merged_alignment
    }
}


task correct_fastq {
    meta {
        version: 'v0.1'
        author: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Correct FASTQs task'
    }

    input {
        File fastq_R1
        File fastq_R2
        File whitelist
        String sample_type
        String? pkr
        String? prefix
        Boolean paired_rna = false

        Int? cpus = 16
        Float? disk_factor = 8.0
        Float? memory_factor = 2
        String? docker_image = "us.gcr.io/buenrostro-share-seq/share_task_correct_fastq:dev"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fastq_R1, "G") + size(fastq_R2, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 16.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    String corrected_fastq_R1 = basename(fastq_R1, ".fastq.gz") + "_corrected.fastq"
    String corrected_fastq_R2 = basename(fastq_R2, ".fastq.gz") + "_corrected.fastq"
    String corrected_fastq_barcode = basename(fastq_R1, ".fastq.gz") + "_corrected_barcode.fastq"
    String monitor_log = "correct_fastqs_monitor.log"

    command <<<
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &

        # Perform barcode error correction on FASTQs
        python3 $(which correct_fastq.py) \
            ~{fastq_R1} \
            ~{fastq_R2} \
            ~{corrected_fastq_R1} \
            ~{corrected_fastq_R2} \
            ~{corrected_fastq_barcode} \
            ~{whitelist} \
            ~{sample_type} \
            ~{prefix} \
            ~{pkr} \
            ~{if paired_rna then "--paired_rna" else ""}

        pigz -p ~{cpus} *.fastq
    >>>

    output {
        File corrected_fastq_R1 = "~{corrected_fastq_R1}.gz"
        File corrected_fastq_R2 = "~{corrected_fastq_R2}.gz"
        File corrected_fastq_barcode= "~{corrected_fastq_barcode}.gz"
        File barcode_qc = "~{prefix}_barcode_qc.txt"
	    File monitor_log = "~{monitor_log}"
    }

    runtime {
        cpu : cpus
        memory : "~{mem_gb} GB"
        disks: "local-disk ~{disk_gb} ~{disk_type}"
        docker : "~{docker_image}"
    }
}


task trim {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: trim ATAC FASTQs.'
    }

    input {
        File fastq_R1 # Pair 1 reads
        File fastq_R2 # Pair 2 reads
        String chemistry

        Int? cpus = 16
        Float? disk_factor = 8.0
        Float? memory_factor = 0.15
        String? docker_image = "us.gcr.io/buenrostro-share-seq/share_task_trim_fastqs_atac:dev"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fastq_R1, "G") + size(fastq_R2, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 16.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    # Read trimming outfiles
    String fastq_R1_trimmed = basename(fastq_R1, ".fastq.gz") + "_trimmed.fastq"
    String fastq_R2_trimmed = basename(fastq_R2, ".fastq.gz") + "_trimmed.fastq"
    String trimming_log_json = basename(fastq_R1, "R1.fastq.gz") + ".atac.preprocess.trimming.log.json"
    String trimming_log_html = basename(fastq_R1, "R1.fastq.gz") + ".atac.preprocess.trimming.log.html"
    String trimming_stats = basename(fastq_R1, "R1.fastq.gz") + ".atac.preprocess.trimming.adapter.stats.txt"
    String monitor_log = 'trim_fastqs_atac_monitor.log'

    command <<<
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &

        # Use trim_fastq script for SHARE ATAC trimming
        if [ '~{chemistry}' == 'shareseq' ]; then
            python3 $(which trim_fastq.py) ~{fastq_R1} ~{fastq_R2} ~{fastq_R1_trimmed} ~{fastq_R2_trimmed} ~{trimming_stats}

        # Use fastp for 10X ATAC trimming
        else
            fastp -i ~{fastq_R1} -I ~{fastq_R2} -o ~{fastq_R1_trimmed} -O ~{fastq_R2_trimmed} -h ~{trimming_log_html} -j ~{trimming_log_json} -G -Q -L -w ~{cpus} 2> ~{trimming_stats}
      
        fi
  
        pigz -p ~{cpus} *.fastq
    >>>

    output {
        File fastq_R1_trimmed = fastq_R1_trimmed + ".gz"
        File fastq_R2_trimmed = fastq_R2_trimmed + ".gz"
        File? tenx_trimming_log_json = trimming_log_json
        File? tenx_trimming_log_html = trimming_log_html
        File trimming_stats = trimming_stats
        File trim_fastqs_atac_monitor = monitor_log
    }

    runtime {
        cpu: cpus
        docker: "${docker_image}"
        disks: "local-disk ${disk_gb} ${disk_type}"
        memory: "${mem_gb} GB"
    }

    parameter_meta {
        fastq_R1: {
                description: 'Pairs 1 fastq',
                help: 'Pairs 1 fastq',
            }
        fastq_R2: {
                description: 'Pairs 2 fastq',
                help: 'Pairs 2 fastq',
            }
    }
}

task bowtie2 {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: align ATAC task using bowtie2'
    }

    input {
        # This task takes in input the preprocessed ATAC fastqs and align them to the genome.
        Array[File] fastq_R1
        Array[File] fastq_R2
        Int? multimappers # = 5
        File genome_index_tar       # This is a tar.gz folder with all the index files.
        String chemistry
        String genome_name          # GRCh38, mm10
        String prefix = "sample-share"
        Int? cpus = 16
        Float? disk_factor = 8.0
        Float? memory_factor = 0.15
        String? docker_image = "us.gcr.io/buenrostro-share-seq/share_task_bowtie2:v1.0.0"
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
    String sorted_bam = "${prefix}.atac.align.k${multimappers}.${genome_name}.read.name.sorted.bam"
    String sorted_bai = "${prefix}.atac.align.k${multimappers}.${genome_name}.read.name.sorted.bam.bai"
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
            ~{if defined(multimappers) then "-k ~{multimappers}" else ""} --rg-id ~{prefix + "."}atac \
            --rg "SM:None" \
            --rg "LB:None" \
            --rg "PL:Illumina" \
            ~{if "~{chemistry}" != "shareseq" then "--sam-append-comment" else ""} \
            -x $genome_prefix \
            -1 ~{sep="," fastq_R1} \
            -2 ~{sep="," fastq_R2} 2> ~{alignment_log} | \
            samtools view \
                -b \
                -S \
                -@ ~{samtools_threads} \
                - \
                -o ~{unsorted_bam}


        if [ '~{chemistry}' != 'shareseq' ]; then
            samtools sort \
                -@ ~{samtools_threads} \
                -m ~{samtools_memory_per_thread}M \
                -n \
                ~{unsorted_bam} \
                -o ~{sorted_bam}
        else
            # Splitting the read name to ge the cell barcode and adding it to the CB tag in the BAM file.
            samtools view -h ~{unsorted_bam} | \
            awk '{if ($0 ~ /^@/) {print $0} else {n=split($1,a,"[,_]"); if ( n > 4 ) {print($0 "\tCB:Z:" a[2]a[3]a[4] "\tXC:Z:" a[2]a[3]a[4] "_" a[5]);}else{print($0 "\tCB:Z:" a[2]a[3]a[4] "\tXC:Z:" a[2]a[3]a[4])}}}' | \
            samtools sort \
                -@ ~{samtools_threads} \
                -m ~{samtools_memory_per_thread}M \
                - \
                -o ~{sorted_bam}
        fi

        samtools index -@ ~{cpus} ~{sorted_bam}

    >>>

    output {
        File atac_alignment = sorted_bam
        File atac_alignment_index = sorted_bai
        File atac_alignment_log = alignment_log
        File atac_alignment_monitor_log = monitor_log
    }

    runtime {
        cpu: cpus
        docker: "${docker_image}"
        disks: "local-disk ${disk_gb} ${disk_type}"
        memory: "${mem_gb} GB"
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


task merge_bams {
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
        String? docker_image = "us.gcr.io/buenrostro-share-seq/share_task_merge_bams:dev"
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
    String merged_bam = "${prefix}.atac.merged.k${multimappers}.${genome_name}.sorted.bam"
    String merged_bai = "${prefix}.atac.merged.k${multimappers}.${genome_name}.sorted.bam.bai"
    String alignment_log = "${prefix}.atac.merged.k${multimappers}.${genome_name}.log"

    String monitor_log = "atac_merge_monitor.log"

    command <<<
        set -e

        bash $(which monitor_script.sh) 2>&1 &

        #sambamba merge -t ~{cpus} ~{unsorted_bam} ~{sep=" " bams}

        #sambamba sort -t ~{cpus} -m ~{command_mem_mb}M -o ~{merged_bam} ~{unsorted_bam}
      
        #sambamba index -t ~{cpus} ~{merged_bam}

        # Trying picard
        
        java -Dsamjdk.compression_level=~{compression_level} -Xms~{command_mem_mb}m -Xmx~{command_mem_mb}m -jar /usr/local/bin/picard.jar \
        MergeSamFiles \
        USE_THREADING=true \
        SORT_ORDER="coordinate" \
        INPUT=~{sep=' INPUT=' bams} \
        OUTPUT=~{merged_bam}

        sambamba index -t ~{cpus} ~{merged_bam}

        sed 's/^[[:space:]]*//g' ~{sep=" " logs} | cut -f 1 -d ' ' | awk '{ sum[FNR%15]+=$1 } END {n_total=length(sum);for (idx=1; idx <= n_total; idx++){print sum[idx]}}' > ~{alignment_log}

    >>>

    output {
        File atac_merged_alignment = merged_bam
        File atac_merged_alignment_index = merged_bai
        File atac_merged_alignment_log = alignment_log
    }

    runtime {
        cpu: cpu
        docker: "${docker_image}"
        disks: "local-disk ${disk} HDD"
        disk: disk + " GB" # TES
        memory: "${machine_mem_mb} MiB"
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