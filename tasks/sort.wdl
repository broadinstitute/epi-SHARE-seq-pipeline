task sort {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: sort bam task'
    }
    
    input {
        # This function takes in input the raw fastqs from Novaseq or
        # Nextseq and perform trimming, adapter removal, appending
        # cell barcode to the read name, and splitting ATAC and RNA.
        
        File bam
        String genome_name
        String type
        String docker_image
        Int cpus= 4
        
    }
    
    Float input_file_size_gb = size(bam, "G")
    Float mem_gb = 8.0
    Int disk_gb = round(20.0 + 4 * input_file_size_gb)

    command {
        set -e
        
#        java -Xmx24g -Djava.io.tmpdir=$dir/temp/ -jar  picard SortSam \
#            SO=coordinate \
#            I=$Name.$Species.bam \
#            O=$Name.$Species.st.bam \
#            VALIDATION_STRINGENCY=SILENT \
#            TMP_DIR=$dir/temp/ 2>>$dir/Run.log
#        samtools index -@ $cores $Name.$Species.st.bam
        
        samtools sort \
            -@ ${cpus} \
            -m ${mem_gb}G \
            ${bam} > ${type}.align.${genome_name}.bam
        samtools index -@ ${cpus} ${type}.align.${genome_name}.bam
        
    }
    
    output {
        File bam_sorted= align.${genome_name}.bam
        File bam_sorted_index= align.${genome_name}.bam.bai
    }

    runtime {
        cpu : ${cpus}
        memory : '${mem_gb} GB'
        disks : 'local-disk ${disk_gb} SSD'
        preemptible: 0 
        docker: docker_image
    }
    
    parameter_meta {
        bam: {
                description: 'bam file',
                help: 'Aligned reads in bam format',
                example: 'aligned.hg38.bam'
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.'
                examples: ['hg38', 'mm10', 'both']
            },
        type: {
                description: 'Assay analyzed',
                help: 'The name of the assay being sorted.'
                examples: ['ATAC', 'RNA']
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
