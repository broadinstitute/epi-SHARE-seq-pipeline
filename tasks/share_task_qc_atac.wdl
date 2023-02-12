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
        File? raw_bam
        File? raw_bam_index
        File? filtered_bam
        File? filtered_bam_index
        File? queryname_final_bam
        File? peaks
        File? tss
        Int? mapq_threshold = 30
        String? barcode_tag = "CB"
        String? genome_name
        String? prefix

        # Runtime
        Int? cpus = 2
        Float? disk_factor = 8.0
        Float? memory_factor = 0.15
        #String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_qc_atac:dev"
        String docker_image = "polumechanos/share_task_qc_atac:dev"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(raw_bam, "G") + size(filtered_bam, "G") + size(queryname_final_bam, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 6.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    # Determining memory for samtools.
    Float samtools_memory_gb = 0.8 * mem_gb # Samtools has overheads so reducing the memory to 80% of the total.

    # Number of threads to beable to use 4GB of memory per thread seems to be the fastest way
    Int samtools_threads_ = floor(samtools_memory_gb / 4)
    Int samtools_threads =  if samtools_threads_ == 0 then 1 else samtools_threads_

    # Now that we know how many threads we can use to assure 4GB of memory per thread
    # we assign any remaining memory to the threads.
    Int samtools_memory_per_thread_ = floor(samtools_memory_gb * 1024 / samtools_threads) # Computing the memory per thread for samtools in MB.
    Int samtools_memory_per_thread = if samtools_memory_per_thread_ < 768 then 768 else samtools_memory_per_thread_

    # Memory for picard
    Float picard_java_heap_factor = 0.9
    Int picard_java_memory = round(picard_java_heap_factor * mem_gb)

    String stats_log = '${default="share-seq" prefix}.atac.qc.stats.${genome_name}.log.txt'
    String duplicate_stats = '${default="share-seq" prefix}.atac.qc.duplicate.stats.${genome_name}.tsv'
    # pdf string needed as required input to Picard CollectInsertSizeMetrics
    String hist_log_pdf = '${default="share-seq" prefix}.atac.qc.hist.${genome_name}.log.pdf'
    String hist_log_png = '${default="share-seq" prefix}.atac.qc.hist.${genome_name}.log.png'
    String hist_log = '${default="share-seq" prefix}.atac.qc.hist.${genome_name}.log.txt'
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
            --bc_tag ${barcode_tag} \
            --mapq_threshold ${mapq_threshold} \
            --prefix "${prefix}.atac.qc.${genome_name}" \
            in.filtered.bam

        # Duplicates per barcode
        python3 $(which count-duplicates-per-barcode.py) \
            -o ${duplicate_stats} \
            --bc_tag ${barcode_tag} \
            --prefix "${prefix}.atac.qc.${genome_name}" \
            in.filtered.bam

        # Fragments in peaks
        # "${prefix}.fragments.in.peak.tsv"
        python3 $(which qc-atac-compute-reads-in-peaks.py) \
            --peaks ${peaks} \
            --mapq_threshold ${mapq_threshold} \
            --bc_tag ${barcode_tag} \
            --prefix "${prefix}.atac.qc.${genome_name}" \
            in.filtered.bam

        echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Raw" > ${stats_log}
        samtools idxstats -@ ${samtools_threads} in.raw.bam >> ${stats_log}

        echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Filtered" >> ${stats_log}
        samtools idxstats -@ ${samtools_threads} in.filtered.bam >> ${stats_log}

        echo '' > ${hist_log}
        java -Xmx${picard_java_memory}G -jar $(which picard.jar) CollectInsertSizeMetrics \
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
        cpu: cpus
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker: "${docker_image}"
        maxRetries: 1
        memory: "${mem_gb} GB"
        memory_retry_multiplier: 2
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
