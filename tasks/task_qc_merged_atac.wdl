version 1.0

# TASK
# qc-merged-atac

task qc_merged_atac {
    meta {
        version: 'v0.1'
        author: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-seq pipeline: merged ATAC QC task'
    }

    input {
        Array[File] barcode_metadata
        File fragments
        File fragments_index
        File tss
        String? prefix
        String? genome_name
        Int? fragment_cutoff = 10

        # Runtime
        Float? disk_factor = 10.0
        Float? memory_factor = 2.0
        String docker_image = "us.gcr.io/buenrostro-share-seq/task_qc_atac:dev"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(barcode_metadata, "G") + size(fragments, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 4.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    String final_barcode_metadata = '${default='merged' prefix}.atac.qc.${genome_name}.metadata.tsv'
    String insert_size_hist = '${default='merged' prefix}.atac.qc.hist.${genome_name}.png'
    String fragment_barcode_rank_plot = '${default='merged' prefix}.atac.qc.${genome_name}.fragment.barcode.rank.plot.png'

    command <<<
        set -e

        bash $(which monitor_script.sh) 1>&2 &

        # TSS enrichment stats
        echo '------ START: Compute TSS enrichment ------' 1>&2
        time python3 $(which qc_atac_compute_tss_enrichment.py) \
            -e 2000 \
            --tss ~{tss} \
            --prefix "~{prefix}.atac.qc.~{genome_name}" \
            ~{fragments}

        # Merge barcode metadata
        echo '------ START: Merge barcode metadata files ------' 1>&2
        time python3 $(which merge_atac_barcode_metadata.py) merged_barcode_metadata ~{sep=' ' barcode_metadata}
        
        # Add TSS enrichment to barcode metadata
        join -j 1 -t $'\t' <(cat ~{prefix}.atac.qc.~{genome_name}.tss_enrichment_barcode_stats.tsv | (sed -u 1q; sort -k1,1)) <(cut -f 1,6-15 merged_barcode_metadata | (sed -u 1q; sort -k1,1)) > ~{final_barcode_metadata}

        # Insert size plot bulk
        gzip -dc ~{fragments} | awk '{print $3-$2}' > insert_sizes
        echo '------ START: Generate TSS enrichment plot for bulk ------' 1>&2
        time python3 $(which plot_insert_size_hist.py) insert_sizes ~{prefix} ~{insert_size_hist}

        # Barcode rank plot
        echo '------ START: Generate barcode rank plot ------' 1>&2
        time Rscript $(which atac_qc_plots.R) ~{final_barcode_metadata} ~{fragment_cutoff} ~{fragment_barcode_rank_plot}
    >>>

    output {
        File tss_enrichment_barcode_stats = "${prefix}.atac.qc.${genome_name}.tss_enrichment_barcode_stats.tsv"
        File tss_enrichment_plot = "${prefix}.atac.qc.${genome_name}.tss_enrichment_bulk.png"
        File enrichment_score_bulk = "${prefix}.atac.qc.${genome_name}.tss_score_bulk.txt"
        File barcode_metadata = final_barcode_metadata
        File insert_size_hist = insert_size_hist
        File? barcode_rank_plot = fragment_barcode_rank_plot
    }

    runtime {
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker: "${docker_image}"
        memory: "${mem_gb} GB"
    }

    parameter_meta {
        barcode_metadata: {
                description: 'ATAC barcode metadata files',
                help: 'Array of TSV files, each containing barcodes and associated statistics; one per entity to be merged.',
                example: ['first.atac.qc.metadata.tsv', 'second.atac.qc.metadata.tsv']
            }  
        fragments: {
                description: 'Merged fragment file',
                help: 'Merged ATAC fragment file.',
                example: 'merged.fragments.tsv.gz'
            }
        fragments_index: {
                description: 'Merged fragment file index',
                help: 'Index for merged ATAC fragment file.',
                example: 'merged.fragments.tsv.gz.tbi'
            }
        tss: {
                description: 'TSS bed file',
                help: 'List of TSSs in BED format.',
                example: 'refseq.tss.bed'
            }
        prefix: {
                description: 'Prefix for output files',
                help: 'Prefix that will be used to name the output files',
                example: 'MyExperiment'
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.',
                examples: ['hg38', 'mm10', 'both']
            }
        fragment_cutoff: {
                description: 'Fragment cutoff',
                help: 'Cutoff for number of fragments required when making fragment barcode rank plot.',
                example: 10
            }
    }    
}
