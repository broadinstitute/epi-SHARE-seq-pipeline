task update_rgid_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: update RGID task'
    }
    
    input {
        # This function takes in input the raw fastqs from Novaseq or
        # Nextseq and perform trimming, adapter removal, appending
        # cell barcode to the read name, and splitting ATAC and RNA.
        
        Bool multimapper= false
        File bam
        String genome_name
        String? prefix= "rna"
        String docker_image
        Int cpus= 4
        
    }
    
    Float input_file_size_gb = size(input[0], "G")
    Float mem_gb = 40.0
    Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    
    String updated_bam= "${prefix}.${genome_name}.rigid.reheader.unique.st.bam"
    String updated_bam_index= "${prefix}.${genome_name}.rigid.reheader.unique.st.bam.bai"

    command {
        set -e
        
       
        # Update RGID, remove low quality reads and unwanted chrs
        # If keepig multimappers, keep primary aligned reads only,
        # otherwise filter by quality score (-q 30)
        
        samtools view -h -@ ${cpus} ${bam} | \
            sed 's/chrMT/chrM/g' | \
            awk -v OFS='\t' '{$1=substr($1,1,length($1)-34)""substr($1,length($1)-22,23)""substr($1,length($1)-34,11); print $0}') | \
            samtools view -@ $cores -bS ${if multimapper then "-F 256" else "-q 30"} > ${updated_bam}
            
        samtools index -@ ${cpus} ${updated_bam}
		
    }
    
    output {
        File rna_rgid_updated_bam= ${updated_bam}
        File rna_rgid_updated_bai= ${updated_bam_index}
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
                description: 'Alignment bam file',
                help: 'Aligned reads in bam format.',
                example: 'hg38.aligned.bam'
            },
        multimapper: {
                description: 'Multimappers flag',
                help: 'Flag to set if you want to keep the multimapping reads.'
                default: false
                examples: [true, false]
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
                help: 'Docker image for preprocessing step. Dependencies: samtools'
                example: ['put link to gcr or dockerhub']
            }
    }


}
