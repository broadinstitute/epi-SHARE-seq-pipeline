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
        File? wdup_bam
        File? wdup_bam_index
        File? fragments
        File? fragments_index
        File? queryname_final_bam
        File? mito_metrics_bulk # From filter step
        File? mito_metrics_barcode # From filter step
        File? peaks
        File? tss
        File? barcode_conversion_dict

        Int? mapq_threshold = 30
        Int? fragment_cutoff = 10
        String? barcode_tag = "CB"
        String? genome_name
        String? prefix

        # Runtime
        Int? cpus = 8
        Float? disk_factor = 10.0
        Float? memory_factor = 0.3
        String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_qc_atac"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(raw_bam, "G") + size(filtered_bam, "G") + size(queryname_final_bam, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 24.0 + memory_factor * input_file_size_gb

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
    String final_barcode_metadata = '${default="share-seq" prefix}.atac.qc.${genome_name}.metadata.tsv'
    String fragment_barcode_rank_plot = "${default="share-seq" prefix}.atac.qc.${genome_name}.fragment.barcode.rank.plot.png"

    String monitor_log = "atac_qc_monitor.log"


    command<<<
        set -e

        # I am not writing to a file anymore because Google keeps track of it automatically.
        bash $(which monitor_script.sh) 1>&2 &

        cp ~{mito_metrics_barcode} mito_metrics

        ln -s ~{raw_bam} in.raw.bam
        ln -s ~{raw_bam_index} in.raw.bam.bai
        ln -s ~{filtered_bam} in.filtered.bam
        ln -s ~{filtered_bam_index} in.filtered.bam.bai
        ln -s ~{wdup_bam} in.wdup.bam
        ln -s ~{wdup_bam_index} in.wdup.bam.bai
        ln -s ~{fragments} in.fragments.tsv.gz
        ln -s ~{fragments_index} in.fragments.tsv.gz.tbi

        # samstats raw
        # output of bowtie2
        echo '------ SKIPPING: SAMstats for raw bam ------' 1>&2
        #time sambamba view -t ~{cpus} in.raw.bam | SAMstats --sorted_sam_file - --outf ~{samstats_raw_out} > ~{samstats_raw_log}

        # SAMstat final filtered file
        # final bam
        echo '------ SKIPPING: SAMstats for final bam ------' 1>&2
        #time sambamba view -t ~{cpus}  in.filtered.bam |  SAMstats --sorted_sam_file - --outf ~{samstats_filtered_out}  > ~{samstats_filtered_log}

        # library complexity
        # queryname_final_bam from filter
        echo '------ START: Compute library complexity ------' 1>&2
        time sambamba view -t ~{cpus} ~{queryname_final_bam} | python3 $(which pbc_stats.py) ~{pbc_stats}

        # TSS enrichment stats
        echo '------ START: Compute TSS enrichment ------' 1>&2
        time python3 $(which qc_atac_compute_tss_enrichment.py) \
            -e 2000 \
            --tss ~{tss} \
            --prefix "~{prefix}.atac.qc.~{genome_name}" \
            in.fragments.tsv.gz

        # Duplicates per barcode
        echo '------ START: Compute duplication per barcode ------' 1>&2
        time python3 $(which qc_atac_count_duplicates_per_barcode.py) \
            -o ~{duplicate_stats} \
            --bc_tag ~{barcode_tag} \
            in.wdup.bam

        if [ ~{if defined(barcode_conversion_dict) then "true" else "false"} == "true" ];then
            cp ~{duplicate_stats} tmp_duplicate_stats
            echo '------ START: Convert barcode duplicate metrics for 10x ------' 1>&2
            time awk -F ",|\t" -v OFS="\t" 'FNR==NR{map[$1]=$2; next}FNR==1{print "barcode","reads_unique","reads_duplicate","pct_duplicates"}FNR>1{print map[$1],$2,$3,$4}' ~{barcode_conversion_dict} tmp_duplicate_stats > ~{duplicate_stats}

            cp mito_metrics tmp_mito_metrics
            echo '------ START: Convert barcode mito metrics for 10x ------' 1>&2
            time awk -F ",|\t" -v OFS="\t" 'FNR==NR{map[$1]=$2; next}FNR==1{print "barcode","raw_reads_nonmito","raw_reads_mito"}FNR>1{print map[$1],$2,$3}' ~{barcode_conversion_dict} tmp_mito_metrics > mito_metrics
        fi

        # Fragments in peaks
        # "~{prefix}.reads.in.peak.tsv"
        echo '------ START: Compute fragments in peaks------' 1>&2
        time python3 $(which qc_atac_compute_reads_in_peaks.py) \
            --peaks ~{peaks} \
            --prefix "~{prefix}.atac.qc.~{genome_name}" \
            in.fragments.tsv.gz

        echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Raw" > ~{stats_log}
        echo '------ START: Samtools idxstats on raw ------' 1>&2
        time samtools idxstats -@ ~{cpus} in.raw.bam >> ~{stats_log}

        echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Filtered" >> ~{stats_log}
        echo '------ START: Samtools idxstats on final ------' 1>&2
        time samtools idxstats -@ ~{cpus} in.filtered.bam >> ~{stats_log}

        echo '' > ~{hist_log}
        echo '------ START: Picard CollectInsertSizeMetrics ------' 1>&2
        time java -Xmx~{picard_java_memory}G -jar $(which picard.jar) CollectInsertSizeMetrics \
            VALIDATION_STRINGENCY=SILENT \
            I=in.raw.bam \
            O=~{hist_log} \
            H=~{hist_log_pdf} \
            W=1000  2>> picard_run.log

        # Insert size plot bulk
        echo '------ START: Generate TSS enrichment plot for bulk ------' 1>&2
        time python3 $(which plot_insert_size_hist.py) ~{hist_log} ~{prefix} ~{hist_log_png}

        echo '------ START: Generate metadata ------' 1>&2
        time join -j 1  <(cat ~{prefix}.atac.qc.~{genome_name}.tss_enrichment_barcode_stats.tsv | (sed -u 1q;sort -k1,1)) <(cat ~{duplicate_stats} | (sed -u 1q;sort -k1,1)) | \
        join -j 1 - <(cat ~{prefix}.atac.qc.~{genome_name}.reads.in.peak.tsv | (sed -u 1q;sort -k1,1)) | \
        join -j 1 - <(cat mito_metrics| (sed -u 1q;sort -k1,1)) | \
        awk -v FS=" " -v OFS=" " 'NR==1{print $0,"pct_reads_promoter","pct_reads_peaks","pct_mito_reads"}NR>1{print $0,$4*100/$7,$10*100/$7,$13*100/($12+$13)}' | sed 's/ /\t/g'> ~{final_barcode_metadata}

        # Barcode rank plot
        echo '------ START: Generate barcod rank plot ------' 1>&2
        time Rscript $(which atac_qc_plots.R) ~{final_barcode_metadata} ~{fragment_cutoff} ~{fragment_barcode_rank_plot}
    >>>

    output {
        File? atac_qc_samstats_raw = samstats_raw_out
        File? atac_qc_samstats_filtered = samstats_filtered_out
        File atac_qc_pbc_stats = pbc_stats

        File atac_qc_final_stats = stats_log
        File atac_qc_final_hist_png = hist_log_png
        File atac_qc_final_hist = hist_log

        File atac_qc_tss_enrichment_barcode_stats = "${prefix}.atac.qc.${genome_name}.tss_enrichment_barcode_stats.tsv"
        File atac_qc_tss_enrichment_plot = "${prefix}.atac.qc.${genome_name}.tss_enrichment_bulk.png"
        File atac_qc_tss_enrichment_score_bulk = "${prefix}.atac.qc.${genome_name}.tss_score_bulk.txt"
        File atac_qc_duplicate_stats = duplicate_stats
        File atac_qc_fragments_in_peaks = "${prefix}.atac.qc.${genome_name}.reads.in.peak.tsv"

        File atac_qc_barcode_metadata = final_barcode_metadata

        File? atac_qc_barcode_rank_plot = fragment_barcode_rank_plot

        #File? atac_qc_monitor_log = monitor_log

        #added so that it can be passed to outer one to get arguments out
        Int? fragment_cutoff_out = fragment_cutoff
    
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
        fragment_cutoff: {
                description: 'Fragment cutoff',
                help: 'Cutoff for number of fragments required when making fragment barcode rank plot.',
                example: 10
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
