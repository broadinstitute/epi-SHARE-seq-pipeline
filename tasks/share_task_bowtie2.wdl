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
        Int? cpus = 16
        File fastq_R1
        File fastq_R2
        File genome_index       # This is a tar.gz folder with all the index files.
        String genome_name      # GRCh38, mm10
        String? docker_image = "polumechanos/share_task_bowtie2"
        String? prefix
    }

    Float input_file_size_gb = size(fastq_R1, "G")
    # This is almost fixed for either mouse or human genome
    Int mem_gb = 16
    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    Int disk_gb = 50

    # Define tmp file name
    String unsorted_bam = "${default="share-seq" prefix}.atac.align.${genome_name}.bam"
    # Define the output names
    String sorted_bam = "${default="share-seq" prefix}.atac.align.${genome_name}.sorted.bam"
    String sorted_bai = "${default="share-seq" prefix}.atac.align.${genome_name}.sorted.bai"
    String alignment_log = "${default="share-seq" prefix}.atac.align.${genome_name}.log"

    command {
        set -e

        tar zxvf ${genome_index} --no-same-owner -C ./
        genome_prefix=$(basename $(find . -type f -name "*.rev.1.bt2") .rev.1.bt2)


        bowtie2 -X2000 \
            -p ${cpus} \
            --rg-id ${prefix + "."}atac \
            -x $genome_prefix \
            -1 ${fastq_R1} \
            -2 ${fastq_R2} 2> ${alignment_log} |\
            samtools view \
                -bS \
                -@ ${cpus} \
                - \
                -o ${unsorted_bam}


        samtools sort \
            -@ ${cpus} \
            -m ${mem_gb}G \
            ${unsorted_bam} > ${sorted_bam}
        samtools index -@ ${cpus} ${sorted_bam}

    }

    output {
        File atac_alignment = sorted_bam
        File atac_alignment_index = sorted_bai
        File atac_alignment_log = alignment_log
    }

    runtime {
        cpu : cpus
        memory : mem_gb+'G'
        disks : 'local-disk ${disk_gb} SSD'
        docker : "${docker_image}"
    }

    parameter_meta {
        fastq_R1: {
                description: 'Read1 fastq',
                help: 'Processed fastq for read1.',
                example: 'processed.atac.R1.fq.gz',
            }
        fastq_R2: {
                description: 'Read2 fastq',
                help: 'Processed fastq for read2.',
                example: 'processed.atac.R2.fq.gz'
            }
        genome_index: {
                description: 'Bowtie2 indexes',
                help: 'Index files for bowtie2 to use during alignment.',
                examples: ['hg19.tar.gz']
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.',
                examples: ['GRCh38', 'mm10']
            }
        prefix: {
                description: 'Prefix for output files',
                help: 'Prefix that will be used to name the output files',
                examples: 'MyExperiment'
            }
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by bowtie2',
                default: 16
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz',
                example: ['put link to gcr or dockerhub']
            }
    }


}
