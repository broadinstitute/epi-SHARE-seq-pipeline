version 1.0

# Import the tasks called by the pipeline
import "../tasks/share_task_bowtie2.wdl" as share_task_align
import "../tasks/share_task_filter_atac.wdl" as share_task_filter
import "../tasks/share_task_qc_atac.wdl" as share_task_qc_atac
#import "../tasks/share_task_count_atac.wdl" as share_task_count
#import "../tasks/share_task_qc_library.wdl" as share_task_qc_library
#import "../tasks/share_task_archr.wdl" as share_task_archr
#import "../tasks/share_task_log_atac.wdl" as share_task_log_atac

workflow wf_atac {
    meta {
        version: 'v1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) and Sai Ma @ Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Sub-workflow to process the ATAC portion of SHARE-seq libraries.'
    }

    input {
        # ATAC sub-workflow inputs
        File chrom_sizes
        File tss_bed
        File peak_set
        Int? mapq_threshold = 30
        String? barcode_tag = "CB"
        String? prefix = "sample"
        String genome_name
        Int? cutoff
        Int? cpus = 4

        # Align-specific inputs
        ## Biological
        Array[File] read1
        Array[File] read2
        Int? align_multimappers = 1
        File genome_index_tar
        ## Runtime
        Int? align_cpus
        Float? align_disk_factor = 8.0
        Float? align_memory_factor = 0.15
        String? align_docker_image

        # Filter-specific inputs
        ## Biological
        Int? filter_minimum_fragments_cutoff = 1
        Int? filter_shift_plus = 4
        Int? filter_shift_minus = -4
        ## Runtime
        Int? filter_cpus = 16
        Float? filter_disk_factor = 8.0
        Float? filter_memory_factor = 0.15
        String? filter_docker_image

        # QC-specific inputs
        ## Biological
        File? raw_bam
        File? raw_bam_index
        File? filtered_bam
        File? filtered_bam_index
        Int? qc_fragment_cutoff
        String? barcode_tag_fragments
        ## Runtime
        Int? qc_cpus = 16
        Float? qc_disk_factor = 8.0
        Float? qc_memory_factor = 0.15
        String? qc_docker_image
    }

    call share_task_align.share_atac_align as align {
        input:
            fastq_R1 = read1,
            fastq_R2 = read2,
            genome_name = genome_name,
            genome_index_tar = genome_index_tar,
            multimappers = align_multimappers,
            prefix = prefix,
            disk_factor = align_disk_factor,
            memory_factor = align_memory_factor,
            cpus = align_cpus,
            docker_image = align_docker_image
    }

    call share_task_filter.share_atac_filter as filter {
        input:
            bam = align.atac_alignment,
            bam_index = align.atac_alignment_index,
            multimappers = align_multimappers,
            shift_plus = filter_shift_plus,
            shift_minus = filter_shift_minus,
            barcode_tag = barcode_tag,
            barcode_tag_fragments = select_first([barcode_tag_fragments,barcode_tag]),
            mapq_threshold = mapq_threshold,
            genome_name = genome_name,
            minimum_fragments_cutoff = filter_minimum_fragments_cutoff,
            prefix = prefix,
            cpus = filter_cpus,
            disk_factor = filter_disk_factor,
            docker_image = filter_docker_image,
            memory_factor = filter_memory_factor
    }

    call share_task_qc_atac.qc_atac as qc_atac{
        input:
            raw_bam = align.atac_alignment,
            raw_bam_index = align.atac_alignment_index,
            filtered_bam = filter.atac_filter_alignment_dedup,
            filtered_bam_index = filter.atac_filter_alignment_dedup_index,
            queryname_final_bam = filter.atac_filter_alignment_dedup_queryname,
            mito_metrics_bulk = filter.atac_filter_mito_metrics_bulk,
            mito_metrics_barcode = filter.atac_filter_mito_metrics_barcode,
            fragments = filter.atac_filter_fragments,
            fragments_index = filter.atac_filter_fragments_index,
            peaks = peak_set,
            tss = tss_bed,
            fragment_cutoff = qc_fragment_cutoff,
            mapq_threshold = mapq_threshold,
            barcode_tag = select_first([barcode_tag_fragments,barcode_tag]),
            genome_name = genome_name,
            prefix = prefix,
            cpus = qc_cpus,
            disk_factor = qc_disk_factor,
            docker_image = qc_docker_image,
            memory_factor = qc_memory_factor
    }

#    call share_task_log_atac.log_atac as log_atac {
#       input:
#           alignment_log = align.atac_alignment_log,
#           dups_log = qc_library.lib_size_log
#    }
#    call share_task_archr.archr as archr{
#        input:
#            atac_frag = count.atac_fragments_filtered,
#            genome = genome_name,
#            peak_set = peak_set,
#            prefix = prefix
#    }
#    call share_task_archr.archr as archr_strict{
#        input:
#            atac_frag = count.atac_fragments_filtered,
#            genome = genome_name,
#            peak_set = peak_set,
#            prefix = '${prefix}_strict',
#            min_tss = 5,
#            min_frags = 1000
#    }

    output {
        # Align
        File? share_atac_alignment_raw = align.atac_alignment
        File? share_atac_alignment_raw_index = align.atac_alignment_index
        File? share_atac_alignment_log = align.atac_alignment_log
        File? share_atac_alignment_monitor_log = align.atac_alignment_monitor_log

        # Filter
        File? share_atac_filter_alignment_dedup = filter.atac_filter_alignment_dedup
        File? share_atac_filter_alignment_dedup_index = filter.atac_filter_alignment_dedup_index
        File? share_atac_filter_alignment_wdup = filter.atac_filter_alignment_wdup
        File? share_atac_filter_alignment_wdup_index = filter.atac_filter_alignment_wdup_index
        File? share_atac_filter_fragments = filter.atac_filter_fragments
        File? share_atac_filter_fragments_index = filter.atac_filter_fragments_index
        File? share_atac_filter_monitor_log = filter.atac_filter_monitor_log
        File? share_atac_filter_mito_metrics_bulk = filter.atac_filter_mito_metrics_bulk
        File? share_atac_filter_mito_metrics_barcode = filter.atac_filter_mito_metrics_barcode

        # QC
        File share_atac_barcode_metadata = qc_atac.atac_qc_barcode_metadata
#        File share_atac_qc_library_counts = qc_library.lib_size_counts
#        File share_atac_qc_library_duplicates = qc_library.lib_size_log
#        File share_atac_qc_library_plot = qc_library.plot

#        File share_atac_qc_final = qc_atac.atac_final_stats
#        File share_atac_qc_hist_plot = qc_atac.atac_final_hist_png
#        File share_atac_qc_hist_txt = qc_atac.atac_final_hist
#        File share_atac_qc_tss_enrichment = qc_atac.atac_tss_pileup_png


#        File share_atac_archr_notebook_output = archr.notebook_output
#        File share_atac_archr_notebook_log = archr.notebook_log
#        File share_atac_archr_barcode_metadata = archr.archr_barcode_metadata
#        File? share_atac_archr_raw_tss_enrichment = archr.archr_raw_tss_by_uniq_frags_plot
#        File? share_atac_archr_filtered_tss_enrichment = archr.archr_filtered_tss_by_uniq_frags_plot
#        File? share_atac_archr_raw_fragment_size_plot = archr.archr_raw_frag_size_dist_plot
#        File? share_atac_archr_filtered_fragment_size_plot = archr.archr_filtered_frag_size_dist_plot

#        File? share_atac_archr_umap_doublets = archr.archr_umap_doublets
#        File? share_atac_archr_umap_cluster_plot = archr.archr_umap_cluster_plot
#        File? share_atac_archr_umap_num_frags_plot = archr.archr_umap_num_frags_plot
#        File? share_atac_archr_umap_tss_score_plot = archr.archr_umap_tss_score_plot
#        File? share_atac_archr_umap_frip_plot = archr.archr_umap_frip_plot

#        File? share_atac_archr_gene_heatmap_plot = archr.archr_heatmap_plot
#        File? share_atac_archr_arrow = archr.archr_arrow
#        File? share_atac_archr_obj = archr.archr_raw_obj
#        File? share_atac_archr_plots_zip = archr.plots_zip
#        Int share_atac_total_reads = log_atac.atac_total_reads
#        Int share_atac_aligned_uniquely = log_atac.atac_aligned_uniquely
#        Int share_atac_unaligned = log_atac.atac_unaligned
#        Int share_atac_feature_reads = log_atac.atac_feature_reads
#        Int share_atac_duplicate_reads = log_atac.atac_duplicate_reads

#        File share_atac_archr_strict_notebook_output = archr_strict.notebook_output
#        File share_atac_archr_strict_notebook_log = archr_strict.notebook_log

#        File? share_atac_archr_strict_raw_tss_enrichment = archr_strict.archr_raw_tss_by_uniq_frags_plot
#        File? share_atac_archr_strict_filtered_tss_enrichment = archr_strict.archr_filtered_tss_by_uniq_frags_plot
#        File? share_atac_archr_strict_raw_fragment_size_plot = archr_strict.archr_raw_frag_size_dist_plot
#        File? share_atac_archr_strict_filtered_fragment_size_plot = archr_strict.archr_filtered_frag_size_dist_plot

#        File? share_atac_archr_strict_umap_doublets = archr_strict.archr_umap_doublets
#        File? share_atac_archr_strict_umap_cluster_plot = archr_strict.archr_umap_cluster_plot
#        File? share_atac_archr_strict_umap_num_frags_plot = archr_strict.archr_umap_num_frags_plot
#        File? share_atac_archr_strict_umap_tss_score_plot = archr_strict.archr_umap_tss_score_plot
#        File? share_atac_archr_strict_umap_frip_plot = archr_strict.archr_umap_frip_plot

#        File? share_atac_archr_strict_gene_heatmap_plot = archr_strict.archr_heatmap_plot
#        File? share_atac_archr_strict_arrow = archr_strict.archr_arrow
#        File? share_atac_archr_strict_obj = archr_strict.archr_raw_obj
#        File? share_atac_archr_strict_plots_zip = archr_strict.plots_zip

    }
}
