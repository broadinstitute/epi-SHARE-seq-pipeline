version 1.0

# TASK
# SHARE-rna-feature-count

task feature_counts_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: assign features rna task'
    }

    input {
        # This function takes in input the bam file produced by the STAR
        # aligner run on a mixed index (e.g. mouse + human) and split
        # the reads into the two genomes

        Boolean multimapper
        Boolean intron
        File bam
        File gtf
        String gene_naming = "gene_name"
        String genome_name
        String? prefix
        String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_count_rna"
        Int cpus= 6
    }

    #Float input_file_size_gb = size(input[0], "G")
    Int mem_gb = 64
    Int disk_gb = 200
    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)


    String out_bam = "${default="share-seq" prefix}.rna.featurecounts.alignment.wdup.${if multimapper then "multi" else "unique"}.${if intron then "intron" else "exon"}.${genome_name}.bam"
    String out_bai = "${default="share-seq" prefix}.rna.featurecounts.alignment.wdup.${if multimapper then "multi" else "unique"}.${if intron then "intron" else "exon"}.${genome_name}.bam.bai"
    String featurecount_log = "${default="share-seq" prefix}.rna.featurecounts.alignment.wdup.${if multimapper then "multi" else "unique"}.${if intron then "intron" else "exon"}.${genome_name}.featurecount.log"


    command {
        set -e

        ln -s ${bam} temp_input.bam

        # Count reads in exons
        # If multimappers are selected use '-Q 0 -M' options.
        # For unique mappers use '-Q 30'
        featureCounts -T ${cpus} \
            -Q ${if multimapper then "0 -M " else "30"} \
            -a ${gtf} \
            -t exon \
            -g ${gene_naming} \
            -o ${default="share-seq" prefix}.rna.featurecounts.alignment.wdup.${if multimapper then "multi" else "unique"}.${if intron then "intron" else "exon"}.${genome_name}.featurecount.txt \
            -R BAM \
            temp_input.bam >> ${featurecount_log}

        temp_filename="temp_input.bam.featureCounts.bam"

        # Extract reads that assigned to genes
        if [[ ${intron} == "true" ]]; then
            featureCounts -T ${cpus} \
            -Q ${if multimapper then "0 -M " else "30"} \
            -a ${gtf} \
            -t gene \
            -g ${gene_naming} \
            -o ${default="share-seq" prefix}.rna.featurecounts.alignment.wdup.${if multimapper then "multi" else "unique"}.${if intron then "intron" else "exon"}.${genome_name}.featurecount.txt \
            -R BAM \
            $temp_filename >> ${featurecount_log}

            temp_filename="$temp_filename.featureCounts.bam"
        fi

        samtools sort -@ ${cpus} -m 8G -o ${out_bam} $temp_filename
        samtools index -@ ${cpus} ${out_bam}

    }

    output {
        File rna_featurecount_alignment = out_bam
        File rna_featurecount_alignment_index = out_bai
        #File rna_featurecount_log = "featureCount.log"
        File rna_featurecount_txt = glob("*.featurecount.txt")[0]
        File rna_featurecount_summary = glob("*.featurecount.txt.summary")[0]
    }

    runtime {
        cpu : cpus
        memory : mem_gb+'G'
        disks : 'local-disk ${disk_gb} SSD'
        maxRetries : 0
        docker: docker_image
    }

    parameter_meta {
        bam: {
                description: 'Alignment bam file',
                help: 'Aligned reads in bam format.',
                example: 'hg38.aligned.bam'
            }
        gtf: {
                description: 'GTF file',
                help: 'Genes definitions in GTF format.',
                example: 'hg38.refseq.gtf'
            }
        multimapper: {
                description: 'Multimappers flag',
                help: 'Flag to set if you want to keep the multimapping reads.',
                default: false,
                example: [true, false]
            }
        intron: {
                description: 'Introns flag',
                help: 'Flag to set if you want to include reads overlapping introns.',
                default: false,
                example: [true, false]
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name genome reference used to align.',
                examples: ['hg38', 'mm10', 'hg19', 'mm9'],
            }
        gene_naming: {
                description: 'Gene nomenclature',
                help: 'Choose if you want to use the official gene symbols (gene_name) or ensemble gene names (gene_id).',
                default: 'gene_name',
                examples: ['gene_name', 'gene_id']
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
                help: 'Docker image for preprocessing step. Dependencies: samtools',
                example: ['put link to gcr or dockerhub']
            }
    }
}
