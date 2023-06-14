version 1.0

# TASK
# atac-merge-qc

task atac_merge_qc {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: ATAC qc statistics'
    }

    input {
        # This function takes in input the raw and filtered bams
        # and compute some alignment metrics along with the TSS
        # enrichment plot.
        File? tss_enrichment_barcode_stats
        File? duplicate_stats
        File? reads_in_peaks
        File? mito_metrics_barcode # From filter step
        File? barcode_conversion_dict
        Int? fragment_cutoff = 10
        String? genome_name
        String? prefix

        # Runtime
        Int? cpus = 1
        Float? disk_factor = 0.08
        Float? memory_factor = 0.3
        String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_qc_atac"
        String singularity_image = "docker://us.gcr.io/buenrostro-share-seq/share_task_qc_atac"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(tss_enrichment_barcode_stats, "G") + size(duplicate_stats, "G") + size(reads_in_peaks, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 4.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(100.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    String final_barcode_metadata = '${default="share-seq" prefix}.atac.qc.${genome_name}.metadata.tsv'
    String fragment_barcode_rank_plot = "${default="share-seq" prefix}.atac.qc.${genome_name}.fragment.barcode.rank.plot.png"

    String monitor_log = "atac_qc_monitor.log"


    command<<<
        set -e

        # I am not writing to a file anymore because Google keeps track of it automatically.
        bash $(which monitor_script.sh) 1>&2 &

        cp ~{mito_metrics_barcode} mito_metrics
        cp ~{duplicate_stats} duplicate_stats


        if [ ~{if defined(barcode_conversion_dict) then "true" else "false"} == "true" ];then
            cp duplicate_stats tmp_duplicate_stats
            echo '------ START: Convert barcode duplicate metrics for 10x ------' 1>&2
            time awk -F ",|\t" -v OFS="\t" 'FNR==NR{map[$1]=$2; next}FNR==1{print "barcode","reads_unique","reads_duplicate","pct_duplicates"}FNR>1{print map[$1],$2,$3,$4}' ~{barcode_conversion_dict} tmp_duplicate_stats > duplicate_stats

            cp mito_metrics tmp_mito_metrics
            echo '------ START: Convert barcode mito metrics for 10x ------' 1>&2
            time awk -F ",|\t" -v OFS="\t" 'FNR==NR{map[$1]=$2; next}FNR==1{print "barcode","raw_reads_nonmito","raw_reads_mito"}FNR>1{print map[$1],$2,$3}' ~{barcode_conversion_dict} tmp_mito_metrics > mito_metrics
        fi

        echo '------ START: Generate metadata ------' 1>&2
        time join -j 1  <(cat ~{   tss_enrichment_barcode_stats} | (sed -u 1q;sort -k1,1)) <(cat duplicate_stats | (sed -u 1q;sort -k1,1)) | \
        join -j 1 - <(cat ~{reads_in_peaks} | (sed -u 1q;sort -k1,1)) | \
        join -j 1 - <(cat mito_metrics| (sed -u 1q;sort -k1,1)) | \
        awk -v FS=" " -v OFS=" " 'NR==1{print $0,"pct_reads_promoter","pct_reads_peaks","pct_mito_reads"}NR>1{print $0,$4*100/$7,$10*100/$7,$13*100/($12+$13)}' | sed 's/ /\t/g'> ~{final_barcode_metadata}

        # Barcode rank plot
        echo '------ START: Generate barcod rank plot ------' 1>&2
        time Rscript $(which atac_qc_plots.R) ~{final_barcode_metadata} ~{fragment_cutoff} ~{fragment_barcode_rank_plot}
    >>>

    output {
        File atac_qc_barcode_metadata = final_barcode_metadata
        File? atac_qc_barcode_rank_plot = fragment_barcode_rank_plot
    }

    runtime {
        cpu: cpus
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker: "${docker_image}"
        singularity: "${singularity_image}"
        maxRetries: 1
        memory: "${mem_gb} GB"
        memory_retry_multiplier: 2
    }

    parameter_meta {
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
