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
        File? fragments
        File? fragments_index
        File? barcode_summary
        File? peaks
        File? tss
        File? barcode_conversion_dict

        Int? fragment_cutoff = 10
        String? genome_name
        String? prefix
        String? subpool="none"

        # Runtime
        Int? cpus = 1
        Float? disk_factor = 10.0
        Float? memory_factor = 0.3
        #String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_qc_atac"
        String docker_image = "docker.io/polumechanos/share_task_qc_atac"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fragments, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 24.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(100.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    # pdf string needed as required input to Picard CollectInsertSizeMetrics
    String hist_log_pdf = '${default="share-seq" prefix}.atac.qc.hist.${genome_name}.log.pdf'
    String hist_log_png = '${default="share-seq" prefix}.atac.qc.hist.${genome_name}.log.png'
    String hist_log = '${default="share-seq" prefix}.atac.qc.hist.${genome_name}.log.txt'
    String tss_pileup_prefix = '${default="share-seq" prefix}.atac.qc.tss.pileup.${genome_name}.log'
    String tss_pileup_out = '${default="share-seq" prefix}.atac.qc.tss.pileup.${genome_name}.log.png'
    String final_barcode_metadata = '${default="share-seq" prefix}.atac.qc.${genome_name}.metadata.tsv'
    String fragment_barcode_rank_plot = "${default="share-seq" prefix}.atac.qc.${genome_name}.fragment.barcode.rank.plot.png"

    String monitor_log = "atac_qc_monitor.log"


    command<<<
        set -e

        # I am not writing to a file anymore because Google keeps track of it automatically.
        bash $(which monitor_script.sh) 1>&2 &

        ln -s ~{fragments} in.fragments.tsv.gz
        ln -s ~{fragments_index} in.fragments.tsv.gz.tbi

        if [ ~{if defined(barcode_conversion_dict) then "true" else "false"} == "true" ]; then
            echo '------ There is a conversion list ------' 1>&2
            if [ '~{subpool}' != "none" ]; then
                echo '------ There is a subpool ------' 1>&2
                awk -v subpool=~{subpool} -v OFS="\t" '{print $1"_"subpool,$2"_"subpool}' ~{barcode_conversion_dict} > temp_conversion
            else
                cp ~{barcode_conversion_dict} > temp_conversion
            fi
            awk -v FS='[,|\t]' -v OFS=',' 'FNR==NR{map[$2]=$1; next}FNR==1{print $0}FNR>1 && map[$1] {print map[$1],$2,$3,$4,$5}' temp_conversion ~{barcode_summary} > temp_summary
        else
            cp ~{barcode_summary} temp_summary
        fi
        echo '------ Filtering fragments ------' 1>&2
        time awk -v threshold=~{fragment_cutoff} -v FS='[,|\t]' 'NR==FNR && ($2-$3-$4-$5)>threshold {Arr[$1]++;next} Arr[$4] {print $0}' temp_summary <( zcat in.fragments.tsv.gz ) | bgzip -c > no-singleton.bed.gz

        tabix --zero-based --preset bed no-singleton.bed.gz

        # TSS enrichment stats
        echo '------ START: Compute TSS enrichment ------' 1>&2
        time python3 $(which qc_atac_compute_tss_enrichment.py) \
            -e 2000 \
            --tss ~{tss} \
            --prefix "~{prefix}.atac.qc.~{genome_name}" \
            no-singleton.bed.gz

        # Fragments in peaks
        # "~{prefix}.reads.in.peak.tsv"
        echo '------ START: Compute fragments in peaks------' 1>&2
        time python3 $(which qc_atac_compute_reads_in_peaks.py) \
            --peaks ~{peaks} \
            --prefix "~{prefix}.atac.qc.~{genome_name}" \
            no-singleton.bed.gz

        echo "insert_size" > ~{hist_log}
        time awk '{print $3-$2}' <(zcat in.fragments.tsv.gz ) | sort --parallel 4 -n | uniq -c | awk -v OFS="\t" '{print $2,$1}' >> ~{hist_log}

        # Insert size plot bulk
        echo '------ START: Generate TSS enrichment plot for bulk ------' 1>&2
        time python3 $(which plot_insert_size_hist.py) ~{hist_log} ~{prefix} ~{hist_log_png}

        echo '------ START: Generate metadata ------' 1>&2

        awk -v FS=',' -v OFS=" " 'NR==1{$1=$1;print $0,"unique","pct_dup","pct_unmapped";next}{$1=$1;if ($2-$3-$4-$5>0){print $0,($2-$3-$4-$5),$3/($2-$4-$5),($5+$4)/$2} else { print $0,0,0,0}}' ~{barcode_summary} > tmp-barcode-stats

        cut -f 1 ~{prefix}.atac.qc.~{genome_name}.reads.in.peak.tsv > barcodes_passing_threshold


        time join -j 1  <(cat ~{prefix}.atac.qc.~{genome_name}.tss_enrichment_barcode_stats.tsv | (sed -u 1q;sort -k1,1)) <(grep -wFf barcodes_passing_threshold tmp-barcode-stats | (sed -u 1q;sort -k1,1)) | \
        join -j 1 - <(cat ~{prefix}.atac.qc.~{genome_name}.reads.in.peak.tsv | (sed -u 1q;sort -k1,1)) | \
        awk -v FS=" " -v OFS=" " 'NR==1{print $0,"pct_reads_promoter","pct_reads_peaks","pct_mito_reads"}NR>1{print $0,$4*100/$7,$10*100/$7,$13*100/($12+$13)}' | sed 's/ /\t/g'> ~{final_barcode_metadata}

        # Barcode rank plot
        echo '------ START: Generate barcod rank plot ------' 1>&2
        time Rscript $(which atac_qc_plots.R) ~{final_barcode_metadata} ~{fragment_cutoff} ~{fragment_barcode_rank_plot}
    >>>

    output {
        File atac_qc_final_hist_png = hist_log_png
        File atac_qc_final_hist = hist_log

        File temp_frag = "no-singleton.bed.gz"

        File atac_qc_tss_enrichment_barcode_stats = "${prefix}.atac.qc.${genome_name}.tss_enrichment_barcode_stats.tsv"
        File atac_qc_tss_enrichment_plot = "${prefix}.atac.qc.${genome_name}.tss_enrichment_bulk.png"
        File atac_qc_tss_enrichment_score_bulk = "${prefix}.atac.qc.${genome_name}.tss_score_bulk.txt"
        File atac_qc_fragments_in_peaks = "${prefix}.atac.qc.${genome_name}.reads.in.peak.tsv"

        File atac_qc_barcode_metadata = final_barcode_metadata

        File? atac_qc_barcode_rank_plot = fragment_barcode_rank_plot

        #File? atac_qc_monitor_log = monitor_log
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
