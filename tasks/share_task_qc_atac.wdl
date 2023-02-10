version 1.0

# TASK
# SHARE-atac-qc-atac

task qc_atac {
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
        File peaks
        File tss
        Int mapq_threshold = 30
        Int minimum_number_fragments = -1
        String? barcode_tag = "CB"
        String genome_name
        String? prefix
        #String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_qc_atac:dev"
        String docker_image = "polumechanos/share_task_qc_atac:dev"
    }

    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    Int disk_gb = 100
    Float input_file_size_gb = size(raw_bam, "G")
    Int mem_gb = 16

    String stats_log = '${default="share-seq" prefix}.atac.qc.stats.${genome_name}.log.txt'
    String hist_log = '${default="share-seq" prefix}.atac.qc.hist.${genome_name}.log.txt'
    String duplicate_stats = '${default="share-seq" prefix}.atac.qc.duplicate.stats.${genome_name}.tsv'
    # pdf string needed as required input to Picard CollectInsertSizeMetrics
    String hist_log_pdf = '${default="share-seq" prefix}.atac.qc.hist.${genome_name}.log.pdf'
    String hist_log_png = '${default="share-seq" prefix}.atac.qc.hist.${genome_name}.log.png'
    String tss_pileup_prefix = '${default="share-seq" prefix}.atac.qc.tss.pileup.${genome_name}.log'
    String tss_pileup_out = '${default="share-seq" prefix}.atac.qc.tss.pileup.${genome_name}.log.png'
    String samstats_raw_log = "${prefix}.atac.qc.${genome_name}.samstats.raw.log.txt"
    String samstats_raw_out = "${prefix}.atac.qc.${genome_name}.samstats.raw.txt"
    String samstats_filtered_log = "${prefix}.atac.qc.${genome_name}.samstats.filtered.log.txt"
    String samstats_filtered_out = "${prefix}.atac.qc.${genome_name}.samstats.filtered.txt"
    String pbc_stats = "${prefix}.atac.qc.${genome_name}.pbcstats.log"


    command {
        set -e

        ln -s ${raw_bam} in.raw.bam
        ln -s ${raw_bam_index} in.raw.bam.bai
        ln -s ${filtered_bam} in.filtered.bam
        ln -s ${filtered_bam_index} in.filtered.bam.bai

        # samstats raw
        # output of bowtie2
        samtools view -o - in.raw.bam | SAMstats --sorted_sam_file - --outf ${samstats_raw_out} > ${samstats_raw_log}

        # SAMstat final filtered file
        # final bam
        samtools view in.filtered.bam |  SAMstats --sorted_sam_file - --outf ${samstats_filtered_out}  > ${samstats_filtered_log}

        # library complexity
        # queryname_final_bam from filter
        samtools view ${queryname_final_bam} | python3 $(which pbc_stats.py) ${pbc_stats}

        # TSS enrichment stats
        python3 $(which qc-atac-tss-enrichment.py) \
            -e 2000 \
            --tss ${tss} \
            --bc_tag CB \
            --mapq_threshold {mapq_threshold} \
            --prefix "${prefix}.atac.qc.${genome_name}" \
            in.filtered.bam

        # Fragments in peaks
        python3 $(which qc-atac-compute-reads-in-peaks.py) \
            --peaks ${peaks} \
            --mapq_threshold {mapq_threshold} \
            --bc_tag CB \
            --prefix "${prefix}.atac.qc.${genome_name}" \
            in.filtered.bam

        # Duplicates per barcode
        python3 $(which count-duplicates-per-barcode.py) \
            -o ${duplicate_stats} \
            --bc_tag CB \
            --mapq_threshold {mapq_threshold} \
            --prefix "${prefix}.atac.qc.${genome_name}" \
            in.filtered.bam

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

        # Insert size plot bulk
        python3 $(which plot_insert_size_hist.py) ${hist_log} ${prefix} ${hist_log_png}


    }

    output {
        File atac_qc_samstats_raw = samstats_raw_out
        File atac_qc_samstats_filtered = samstats_filtered_out
        File atac_qc_pbc_stats = pbc_stats

        File atac_qc_final_stats = stats_log
        File atac_qc_final_hist_png = hist_log_png
        File atac_qc_final_hist = hist_log

        File atac_qc_tss_enrichment_barcode_stats = "${prefix}.atac.qc.${genome_name}.tss_enrichment_barcode_stats.tsv"
        File atac_qc_tss_enrichment_plot = "${prefix}.atac.qc.${genome_name}.tss_enrichment_bulk.png"
        File atac_qc_tss_enrichment_score_bulk = "${prefix}.atac.qc.${genome_name}.tss_score_bulk.txt"
        File atac_qc_duplicate_stats = duplicate_stats
        File atac_qc_fragments_in_peaks = "${prefix}.atac.qc.${genome_name}.fragments.in.peak.tsv"
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
