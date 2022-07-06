version 1.0

# TASK
# SHARE-atac-STAR

task share_rna_align {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: align RNA task'
    }

    input {
        # This function takes in input the pre-processed fastq and align it to the genome
        # using STAR.

        Array[File] fastq_R1
        Array[File]? fastq_R2
        File genome_index_tar
        String genome_name
        String? prefix
        String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_star"
        Int cpus = 16
    }
    #Float input_file_size_gb = size(input[0], "G")
    Int samtools_cpus = 6
    Int samtools_mem_gb = 8
    Int mem_gb = 64
    Int disk_gb = 250
    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)

    # Define the output names
    String sorted_bam = "${default="share-seq" prefix}.rna.align.${genome_name}.sorted.bam"
    String sorted_bai = "${default="share-seq" prefix}.rna.align.${genome_name}.sorted.bam.bai"
    String alignment_log = "${default="share-seq" prefix}.rna.align.${genome_name}.log"

    command {
        set -e
        # Untar the genome
        tar xvzf ${genome_index_tar} --no-same-owner -C ./

        mkdir out

        $(which STAR) \
            --runThreadN ${cpus} \
            --chimOutType WithinBAM \
            --genomeDir ./ \
            --readFilesIn ${sep=',' fastq_R1} ${sep=',' fastq_R2}  \
            --outFileNamePrefix out/${default="share-seq" prefix}.rna.align.${genome_name}. \
            --outFilterMultimapNmax 20 \
            --outFilterScoreMinOverLread 0.3 \
            --outFilterMatchNminOverLread 0.3 \
            --outSAMattributes NH HI AS nM MD \
            --limitOutSJcollapsed 2000000 \
            --outSAMtype BAM Unsorted \
            --limitIObufferSize 400000000 \
            --outReadsUnmapped Fastx \
            --readFilesCommand zcat

        $(which samtools) sort \
            -@ ${samtools_cpus} \
            -m ${samtools_mem_gb}G \
            -o out/${sorted_bam} \
            out/${default="share-seq" prefix}.rna.align.${genome_name}.Aligned.out.bam

        $(which samtools) index \
            -@ ${cpus} \
            out/${sorted_bam}
    }

    output {
        File rna_alignment = "out/${sorted_bam}"
        File rna_alignment_index = "out/${sorted_bai}"
        File rna_alignment_log = glob('out/*.Log.final.out')[0]
    }

    runtime {
        cpu : cpus
        memory : mem_gb+'G'
        disks : 'local-disk ${disk_gb} SSD'
        maxRetries: 0
        docker: docker_image
    }

    parameter_meta {
        fastq_R1: {
                description: 'Read1 fastq',
                help: 'Processed fastq for read1.',
                example: 'processed.atac.R1.fq.gz'
            }
        genome_index_tar: {
                description: 'STAR indexes',
                help: 'Index files for STAR to use during alignment in tar.gz.',
                example: ['']
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.',
                example: ['hg38', 'mm10', 'both']
            }
        prefix: {
                description: 'Prefix for output files',
                help: 'Prefix that will be used to name the output files',
                example: 'MyExperiment'
            }
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by bowtie2',
                example: '4'
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: STAR',
                example: ['put link to gcr or dockerhub']
            }
    }
}
