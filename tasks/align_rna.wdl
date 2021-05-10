task align_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: align RNA task'
    }
    
    input {
        # This function takes in input the raw fastqs from Novaseq or
        # Nextseq and perform trimming, adapter removal, appending
        # cell barcode to the read name, and splitting ATAC and RNA.
        
        File fastq_R1
        File genome_index_tar
        String genome_name
        String? prefix
        String docker_image
        Int cpus= 4
        
    }
    
    Float input_file_size_gb = size(input[0], "G")
    Float mem_gb = 40.0
    Int disk_gb = round(20.0 + 4 * input_file_size_gb)

    command {
        set -e
        
        # Untar the genome
        tar xvzf ${genome_index_tar} --no-same-owner -C ./index
        
        # Find the prefix
        PREFIX_GENOME=`ls ./index/*.sa | sed 's/.sa//'`
        

        STAR --chimOutType WithinBAM \
             --runThreadN ${cpus} \
             --genomeDir $PREFIX_GENOME \
             --readFilesIn ${fastq_R1}  \
             --outFileNamePrefix out/${prefix}.${genome_name}. \
             --outFilterMultimapNmax 20 \
             --outFilterMismatchNoverLmax 0.06 \
             --limitOutSJcollapsed 2000000 \
             --outSAMtype BAM Unsorted \
             --limitIObufferSize 400000000 \
             --outReadsUnmapped Fastx \
             --readFilesCommand zcat
             
    }
    
    output {
        File rna_align_bam= glob('out/*.bam')[0]
        File rna_align_bai= glob('out/*.bai')[0]
        File rna_align_log= glob('out/*.Log.final.out')[0]
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
            },
        genome_index_tar: {
                description: 'STAR indexes',
                help: 'Index files for STAR to use during alignment in tar.gz.'
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
                help: 'Docker image for preprocessing step. Dependencies: STAR'
                example: ['put link to gcr or dockerhub']
            }
    }


}
