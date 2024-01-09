version 1.0

# TASK
# qc-atac

task qc_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard IGVF pipeline: ATAC qc statistics task'
    }

    input {
        # This function takes in input the raw and filtered bams
        # and compute some alignment metrics along with the TSS
        # enrichment plot.
        File? fragments
        File? fragments_index
        File? barcode_summary
        File? tss
        File? barcode_conversion_dict

        Int? fragment_cutoff = 10
        String? genome_name
        String? prefix
        String? subpool="none"

        # Runtime
        Int? cpus = 8
        Float? disk_factor = 10.0
        Float? memory_factor = 20
        String docker_image = "us.gcr.io/buenrostro-share-seq/task_qc_atac:dev"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fragments, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 32.0 + memory_factor * input_file_size_gb

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

        echo '------ Number of barcodes BEFORE filtering------' 1>&2
        wc -l temp_summary

        echo '------ Filtering fragments ------' 1>&2
        time awk -v threshold=~{fragment_cutoff} -v FS='[,|\t]' 'NR==FNR && ($2-$3-$4-$5)>threshold {Arr[$1]++;next} Arr[$4] {print $0}' temp_summary <( zcat in.fragments.tsv.gz ) > no-singleton.bed

        bgzip -l 5 -@ ~{cpus} -c no-singleton.bed > no-singleton.bed.gz
        
        echo '------ Number of barcodes AFTER filtering------' 1>&2
        cat temp_summary | grep -v barcode | awk -v FS="," -v threshold=~{fragment_cutoff} '($2-$3-$4-$5)>threshold' | wc -l
        
        tabix --zero-based --preset bed no-singleton.bed.gz

        # TSS enrichment stats
        echo '------ START: Compute TSS enrichment ------' 1>&2
        time python3 $(which compute_tss_enrichment.py) \
            -e 2000 \
            -p ~{cpus} \
            --regions ~{tss} \
            --prefix "~{prefix}.atac.qc.~{genome_name}" \
            no-singleton.bed.gz

         

        # Insert size plot bulk
        echo '------ START: Generate Insert size plot ------' 1>&2

        echo "insert_size" > ~{hist_log}
        time awk '{print $3-$2}' <(zcat in.fragments.tsv.gz ) | sort --parallel 4 -n | uniq -c | awk -v OFS="\t" '{print $2,$1}' >> ~{hist_log}
        time python3 $(which plot_insert_size_hist.py) ~{hist_log} ~{prefix} ~{hist_log_png}

        echo '------ START: Generate metadata ------' 1>&2

        awk -v FS=',' -v OFS=" " 'NR==1{$1=$1;print $0,"unique","pct_dup","pct_unmapped";next}{$1=$1;if ($2-$3-$4-$5>0){print $0,($2-$3-$4-$5),$3/($2-$4-$5),($5+$4)/$2} else { print $0,0,0,0}}' temp_summary > tmp-barcode-stats

        cut -f 1 ~{prefix}.atac.qc.~{genome_name}.tss_enrichment_barcode_stats.tsv > barcodes_passing_threshold

        time join -j 1  <(cat ~{prefix}.atac.qc.~{genome_name}.tss_enrichment_barcode_stats.tsv | (sed -u 1q;sort -k1,1)) <(grep -wFf barcodes_passing_threshold tmp-barcode-stats | (sed -u 1q;sort -k1,1)) | \
        (sed -u 1q;sort -k1,1) | \
        awk -v FS=" " -v OFS=" " 'NR==1{print $0,"pct_reads_promoter"}NR>1{print $0,$4*100/$7}' | sed 's/ /\t/g' > ~{final_barcode_metadata}

        head ~{final_barcode_metadata}
        head ~{final_barcode_metadata} | awk '{print NF}'

        # Barcode rank plot
        echo '------ START: Generate barcod rank plot ------' 1>&2
        time Rscript $(which atac_qc_plots.R) ~{final_barcode_metadata} ~{fragment_cutoff} ~{fragment_barcode_rank_plot}
    >>>

    output {
        File atac_qc_final_hist_png = hist_log_png
        File atac_qc_final_hist = hist_log

        File temp_frag = "no-singleton.bed.gz"
        
        File temp_summary = "temp_summary"

        File atac_qc_tss_enrichment_barcode_stats = "${prefix}.atac.qc.${genome_name}.tss_enrichment_barcode_stats.tsv"
        File atac_qc_tss_enrichment_plot = "${prefix}.atac.qc.${genome_name}.tss_enrichment_bulk.png"
        File atac_qc_tss_enrichment_score_bulk = "${prefix}.atac.qc.${genome_name}.tss_score_bulk.txt"

        File atac_qc_barcode_metadata = final_barcode_metadata

        File? atac_qc_barcode_rank_plot = fragment_barcode_rank_plot
    }

    runtime {
        cpu: cpus
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker: "${docker_image}"
        memory: "${mem_gb} GB"
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
