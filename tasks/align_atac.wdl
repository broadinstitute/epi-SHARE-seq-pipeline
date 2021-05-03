task align_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: align ATAC task'
    }
    
    input {
        # This function takes in input the raw fastqs from Novaseq or
        # Nextseq and perform trimming, adapter removal, appending
        # cell barcode to the read name, and splitting ATAC and RNA.
        
        File fastq_R1
        File fastq_R2
        File genome_index
        String genome_name
        String? prefix
        String docker_image
        Int cpus= 4
        
    }
    
    Float input_file_size_gb = size(input[0], "G")
    Float mem_gb = 8.0
    Int disk_gb = round(20.0 + 4 * input_file_size_gb)

    command {
        set -e
        
        bowtie2 -X2000 \
            -p $cores \
            --rg-id $Name \
            -x ${genome_index} \
            -1 fastq_R1 \
			-2 fastq_R2 |\
        samtools view \
            -bS
            -@ ${cpus} \
            - \
            -o atac.align.${genome_name}.bam) \
            2> atac.bowtie2.align.${genome_name}.log
        
    }
    
    output {
        File atac_bowtie2_align= atac.align.${genome_name}.bam
        File log= atac.bowtie2.align.${genome_name}.log
    }

    runtime {
        cpu : ${cpus}
        memory : '${mem_gb} GB'
        disks : 'local-disk ${disk_gb} SSD'
        preemptible: 0 
        docker: docker_image
    }
    
    parameter_meta {
        fastq_R1: {
                description: 'Read1 fastq',
                help: 'Processed fastq for read1.',
                example: 'processed.atac.R1.fq.gz'
            }
        fastq_R2: {
                description: 'Read2 fastq',
                help: 'Processed fastq for read2.',
                example: 'processed.atac.R2.fq.gz'
            }
        genome_index: {
                description: 'Bowtie2 indexes',
                help: 'Index files for bowtie2 to use during alignemnt.'
                examples: ['']
            },
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.'
                examples: ['hg38', 'mm10', 'both']
            },
        prefix: {
                description: 'Prefix for output files',
                help: 'Prefix that will be used to name the output files'
                examples: 'MyExperiment'
            },
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by bowtie2'
                examples: '4'
            },
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz'
                example: ['put link to gcr or dockerhub']
            }
    }


}
