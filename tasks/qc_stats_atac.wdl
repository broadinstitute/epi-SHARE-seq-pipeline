version 1.0

task qc_stats_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: ATAC qc statistics task'
    }

    input {
        # This function takes in input the raw and filtered bams
        # and compute some alignment metrics along with the TSS
        # enrichment plot.

        Int cpus= 4
        File raw_bam
        File raw_bam_index
        File filtered_bam
        File filtered_bam_index
        File tss
        String genome_name
        String? prefix
        String docker_image


    }

    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    Int disk_gb = 50
    Float input_file_size_gb = size(raw_bam, "G")
    Int mem_gb = 16

    String stats_log = '${prefix + '.'}atac.stats.${genome_name}.log.txt'
    String hist_log = '${prefix + '.'}atac.hist.${genome_name}.log.txt'
    String hist_log_pdf = '${prefix + '.'}atac.hist.${genome_name}.log.pdf'
    String tss_pileup_prefix = '${prefix + '.'}atac.tss.pileup.${genome_name}.log'
    String tss_pileup_out = '${prefix + '.'}atac.tss.pileup.${genome_name}.log.png'


    command {
        set -e

        mv ${raw_bam} in.raw.bam
        mv ${raw_bam_index} in.raw.bai
        mv ${filtered_bam} in.filtered.bam
        mv ${filtered_bam_index} in.filtered.bai


        echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Raw" > ${stats_log}
        samtools idxstats in.raw.bam >> ${stats_log}

        echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Filtered" >> ${stats_log}
        samtools idxstats in.filtered.bam >> ${stats_log}

        echo '' > ${hist_log}
        java -jar $(which picard.jar) CollectInsertSizeMetrics \
            VALIDATION_STRINGENCY=SILENT \
            I=in.raw.bam \
            O=${hist_log} \
            H=${hist_log_pdf} \
            W=1000  2>> picard_run.log

        # make TSS pileup fig # original code has a 'set +e' why?
        # the pyMakeVplot is missing
        python $(which make-tss-pileup-jbd.py) \
            -a in.filtered.bam \
            -b ${tss} \
            -e 2000 \
            -p ends \
            -v \
            -u \
            -o ${tss_pileup_prefix}

    }

    output {
        File final_stats = stats_log
        File final_hist_stats_pdf = hist_log_pdf
        File final_hist_stats = stats_log
        File tss_pileup = tss_pileup_out
    }

    runtime {
        cpu : cpus
        memory : mem_gb+'G'
        disks : 'local-disk ${disk_gb} SSD'
        maxRetries : 0
        docker: docker_image
    }

    parameter_meta {
        raw_bam: {
                description: 'Unfiltered bam',
                help: 'Not filtered alignment bam file.',
                example: 'aligned.hg38.bam'
            }
        raw_bam: {
                description: 'Filtered bam',
                help: 'Filtered alignment bam file. Typically, no duplicates and quality filtered.',
                example: 'aligned.hg38.rmdup.filtered.bam'
            }
        tss: {
                description: 'TSS bed file',
                help: 'List of TSS in bed format used for the enrichment plot.',
                example: 'refseq.tss.bed'
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.',
                examples: ['hg38', 'mm10', 'both']
            }
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by bowtie2',
                examples: '4'
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz',
                example: ['put link to gcr or dockerhub']
            }
    }
}
