task assign_features_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: assign features rna task'
    }
    
    input {
        # This function takes in input the bam file produced by the STAR
        # aligner run on a mixed index (e.g. mouse + human) and split
        # the reads into the two genomes
        
        Boolean multimapper= false
        Boolean intron= false
        File bam
        File gtf
        String gene_naming= "gene_name"
        String genome_name
        String? prefix= "rna.featureCounts"
        String docker_image
        Int cpus= 4
        
    }
    
    Float input_file_size_gb = size(input[0], "G")
    Float mem_gb = 40.0
    Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    
    command {
        set -e
        
        # Count reads in exons
        # If multimappers are selected use '-Q 0 -M' options.
        # For unique mappers use '-Q 30'
        featureCounts -T ${cpus} \
            -Q ${if multimapper then "0 -M " alse "30"} \
            -a ${gtf} \
            -t exon \
            -g ${gene_naming} \
            -o ${prefix}.${genome_name}.feature.count.txt \
            -R BAM \
            ${bam} >> featureCount.log
        
        temp_filename=${bam}.featureCounts.bam
        
        # Extract reads that assigned to genes
        if [[ ${intron} == "true" ]]; then        
            featureCounts -T ${cpus} \
            -Q ${if multimapper then "0 -M " alse "30"} \
            -a ${gtf} \
            -t gene \
            -g ${gene_naming} \
            -o ${prefix}.${genome_name}.feature.count.txt \
            -R BAM \
            ${bam}.featureCounts.bam >> featureCount.log
            
            temp_filename=$temp_filename.featureCounts.bam
        fi
        
        samtools sort -@ ${cpus} -m 2G -o ${prefix}.${genome_name}.wdup.bam $temp_filename
        samtools index -@ ${cpus} ${prefix}.${genome_name}.wdup.bam
        
    }
    
    output {
        File assigned_features_rna_bam= "${prefix}.${genome_name}.wdup.bam"
        File assigned_features_rna_bam= "${prefix}.${genome_name}.wdup.bam.bai"
        File featureCounts_log= "featureCount.log"
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
        gtf: {
                description: 'GTF file',
                help: 'Genes definitions in GTF format.',
                example: 'hg38.refseq.gtf'
            },
        multimapper: {
                description: 'Multimappers flag',
                help: 'Flag to set if you want to keep the multimapping reads.'
                default: false
                examples: [true, false]
            },
        intron: {
                description: 'Introns flag',
                help: 'Flag to set if you want to include reads overlapping introns.'
                default: false
                examples: [true, false]
            },
        genome_name: {
                description: 'Reference name',
                help: 'The name genome reference used to align.'
                examples: ['hg38', 'mm10', 'hg19', 'mm9']
            },
        gene_naming: {
                description: 'Gene nomenclature',
                help: 'Choose if you want to use the official gene symbols (gene_name) or ensemble gene names (gene_id).'
                default: 'gene_name'
                examples: ['gene_name', 'gene_id']
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
