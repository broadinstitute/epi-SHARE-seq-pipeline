version 1.0

# Import the tasks called by the pipeline
import "../tasks/share_task_correct_fastq.wdl" as share_task_correct_fastq
import "../tasks/share_task_trim_fastqs_atac.wdl" as share_task_trim
import "../tasks/share_task_bowtie2.wdl" as share_task_align
import "../tasks/share_task_merge_bams.wdl" as share_task_merge_bams
import "../tasks/share_task_filter_atac.wdl" as share_task_filter
import "../tasks/share_task_qc_atac.wdl" as share_task_qc_atac
import "../tasks/share_task_log_atac.wdl" as share_task_log_atac
import "../tasks/share_task_archr.wdl" as share_task_archr


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
        String? barcode_tag_fragments
        String chemistry
        String? prefix = "sample"
        String genome_name
        Int? cutoff
        String pipeline_modality = "full"
        Boolean trim_fastqs = true
        File? barcode_conversion_dict # For 10X multiome

        # Correct-specific inputs
        Boolean correct_barcodes = true
        File whitelist
        String? pkr
        # Runtime parameters
        Int? correct_cpus = 16
        Float? correct_disk_factor = 8.0
        Float? correct_memory_factor = 0.08
        String? correct_docker_image

        # Align-specific inputs
        Array[File] read1
        Array[File] read2
        Array[File] fastq_barcode = []
        Int? align_multimappers
        File genome_index_tar
        # Runtime parameters
        Int? align_cpus
        Float? align_disk_factor = 8.0
        Float? align_memory_factor = 0.15
        String? align_docker_image

        # Merge-specific inputs
        # Runtime parameters
        Int? merge_cpus
        Float? merge_disk_factor = 8.0
        Float? merge_memory_factor = 0.15
        String? merge_docker_image

        # Filter-specific inputs
        Int? filter_minimum_fragments_cutoff
        Int? filter_shift_plus = 4
        Int? filter_shift_minus = -4
        # Runtime parameters
        Int? filter_cpus = 16
        Float? filter_disk_factor = 8.0
        Float? filter_memory_factor = 0.15
        String? filter_docker_image

        # QC-specific inputs
        File? raw_bam
        File? raw_bam_index
        File? filtered_bam
        File? filtered_bam_index
        Int? qc_fragment_cutoff
        # Runtime parameters
        Int? qc_cpus = 16
        Float? qc_disk_factor = 8.0
        Float? qc_memory_factor = 0.15
        String? qc_docker_image

        # Trim-specific inputs
        # Runtime parameters
        Int? trim_cpus = 16
        Float? trim_disk_factor = 8.0
        Float? trim_memory_factor = 0.15
        String? trim_docker_image
        
        # ArchR-specific inputs
        # Runtime parameters
        Float? archr_disk_factor
        Float? archr_memory_factor 
        String? archr_docker_image
    }

    String barcode_tag_fragments_ = if chemistry=="shareseq" then select_first([barcode_tag_fragments, "XC"]) else select_first([barcode_tag_fragments, barcode_tag])

    File? none

    # Perform barcode error correction on FASTQs.
    if ( chemistry == "shareseq" && correct_barcodes ) {
        scatter (idx in range(length(read1))) {
            call share_task_correct_fastq.share_correct_fastq as correct {
                input:
                    fastq_R1 = read1[idx],
                    fastq_R2 = read2[idx],
                    fastq_barcode = if length(fastq_barcode) > 0 then fastq_barcode[idx] else none,
                    whitelist = whitelist,
                    sample_type = "ATAC",
                    pkr = pkr,
                    prefix = prefix,
                    cpus = correct_cpus,
                    disk_factor = correct_disk_factor,
                    memory_factor = correct_memory_factor,
                    docker_image = correct_docker_image
            }
        }
    }

    if ( trim_fastqs ){
        # Remove dovetail in the ATAC reads.
        scatter (read_pair in zip(select_first([correct.corrected_fastq_R1, read1]), select_first([correct.corrected_fastq_R2, read2]))) {
            call share_task_trim.share_trim_fastqs_atac as trim {
                input:
                    fastq_R1 = read_pair.left,
                    fastq_R2 = read_pair.right,
                    chemistry = chemistry,
                    cpus = trim_cpus,
                    disk_factor = trim_disk_factor,
                    memory_factor = trim_memory_factor,
                    docker_image = trim_docker_image
            }
        }
    }

    if (  "~{pipeline_modality}" != "no_align" ) {
        scatter(read_pair in zip(select_first([trim.fastq_R1_trimmed, correct.corrected_fastq_R1, read1]), select_first([trim.fastq_R2_trimmed, correct.corrected_fastq_R2, read2]))) {
            call share_task_align.share_atac_align as align {
                input:
                    fastq_R1 = [read_pair.left],
                    fastq_R2 = [read_pair.right],
                    chemistry= chemistry,
                    genome_name = genome_name,
                    genome_index_tar = genome_index_tar,
                    multimappers = align_multimappers,
                    prefix = prefix,
                    disk_factor = align_disk_factor,
                    memory_factor = align_memory_factor,
                    cpus = align_cpus,
                    docker_image = align_docker_image
            }
        }
        
        call share_task_merge_bams.share_atac_merge_bams as merge{
            input:
                bams = align.atac_alignment,
                logs = align.atac_alignment_log,
                multimappers = align_multimappers,
                genome_name = genome_name,
                prefix = prefix,
                cpus = merge_cpus,
                memory_factor = merge_memory_factor,
                disk_factor = align_disk_factor,
                docker_image = merge_docker_image
        }


        call share_task_filter.share_atac_filter as filter {
            input:
                bam = merge.atac_merged_alignment,
                bam_index = merge.atac_merged_alignment_index,
                multimappers = align_multimappers,
                shift_plus = filter_shift_plus,
                shift_minus = filter_shift_minus,
                barcode_tag = barcode_tag,
                barcode_tag_fragments = barcode_tag_fragments_,
                mapq_threshold = mapq_threshold,
                genome_name = genome_name,
                minimum_fragments_cutoff = filter_minimum_fragments_cutoff,
                prefix = prefix,
                barcode_conversion_dict = barcode_conversion_dict,
                cpus = filter_cpus,
                disk_factor = filter_disk_factor,
                docker_image = filter_docker_image,
                memory_factor = filter_memory_factor
        }

        call share_task_qc_atac.qc_atac as qc_atac{
            input:
                raw_bam = merge.atac_merged_alignment,
                raw_bam_index = merge.atac_merged_alignment_index,
                filtered_bam = filter.atac_filter_alignment_dedup,
                filtered_bam_index = filter.atac_filter_alignment_dedup_index,
                queryname_final_bam = filter.atac_filter_alignment_dedup_queryname,
                wdup_bam = filter.atac_filter_alignment_wdup,
                wdup_bam_index = filter.atac_filter_alignment_wdup_index,
                mito_metrics_bulk = filter.atac_filter_mito_metrics_bulk,
                mito_metrics_barcode = filter.atac_filter_mito_metrics_barcode,
                fragments = filter.atac_filter_fragments,
                fragments_index = filter.atac_filter_fragments_index,
                barcode_conversion_dict = barcode_conversion_dict,
                peaks = peak_set,
                tss = tss_bed,
                fragment_cutoff = qc_fragment_cutoff,
                mapq_threshold = mapq_threshold,
                barcode_tag = barcode_tag_fragments_,
                genome_name = genome_name,
                prefix = prefix,
                cpus = qc_cpus,
                disk_factor = qc_disk_factor,
                docker_image = qc_docker_image,
                memory_factor = qc_memory_factor
        }

        call share_task_log_atac.log_atac as log_atac {
        input:
            alignment_log = merge.atac_merged_alignment_log,
            dups_log = qc_atac.atac_qc_duplicate_stats
        }

        if (  "~{pipeline_modality}" == "full" ) {
            call share_task_archr.archr as archr{
                input:
                    atac_frag = filter.atac_filter_fragments,
                    genome = genome_name,
                    peak_set = peak_set,
                    prefix = prefix,
                    memory_factor = archr_memory_factor,
                    disk_factor = archr_disk_factor,
                    docker_image = archr_docker_image
            }
        }
    }

    output {
        # Correction/trimming
        Array[File]? atac_read1_processed = if defined(trim.fastq_R1_trimmed) then trim.fastq_R1_trimmed else correct.corrected_fastq_R1
        Array[File]? atac_read2_processed = if defined(trim.fastq_R1_trimmed) then trim.fastq_R2_trimmed else correct.corrected_fastq_R2

        # Align
        File? share_atac_alignment_raw = merge.atac_merged_alignment
        File? share_atac_alignment_raw_index = merge.atac_merged_alignment_index
        File? share_atac_alignment_log = merge.atac_merged_alignment_log

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
        File? share_atac_barcode_metadata = qc_atac.atac_qc_barcode_metadata
        File? share_atac_qc_final = qc_atac.atac_qc_final_stats
        File? share_atac_qc_hist_plot = qc_atac.atac_qc_final_hist_png
        File? share_atac_qc_hist_txt = qc_atac.atac_qc_final_hist
        File? share_atac_qc_tss_enrichment = qc_atac.atac_qc_tss_enrichment_plot
        File? share_atac_qc_barcode_rank_plot = qc_atac.atac_qc_barcode_rank_plot

        # Log
        Int? share_atac_total_reads = log_atac.atac_total_reads
        Int? share_atac_aligned_uniquely = log_atac.atac_aligned_uniquely
        Int? share_atac_unaligned = log_atac.atac_unaligned
        Int? share_atac_feature_reads = log_atac.atac_feature_reads
        Int? share_atac_duplicate_reads = log_atac.atac_duplicate_reads
        Float? share_atac_percent_duplicates = log_atac.atac_pct_dup

        # ArchR
        File? share_atac_archr_notebook_output = archr.notebook_output
        File? share_atac_archr_notebook_log = archr.notebook_log
        File? share_atac_archr_barcode_metadata = archr.archr_barcode_metadata
        File? share_atac_archr_raw_tss_enrichment = archr.archr_raw_tss_by_uniq_frags_plot
        File? share_atac_archr_filtered_tss_enrichment = archr.archr_filtered_tss_by_uniq_frags_plot
        File? share_atac_archr_raw_fragment_size_plot = archr.archr_raw_frag_size_dist_plot
        File? share_atac_archr_filtered_fragment_size_plot = archr.archr_filtered_frag_size_dist_plot

        File? share_atac_archr_umap_doublets = archr.archr_umap_doublets
        File? share_atac_archr_umap_cluster_plot = archr.archr_umap_cluster_plot
        File? share_atac_archr_umap_num_frags_plot = archr.archr_umap_num_frags_plot
        File? share_atac_archr_umap_tss_score_plot = archr.archr_umap_tss_score_plot
        File? share_atac_archr_umap_frip_plot = archr.archr_umap_frip_plot

        File? share_atac_archr_gene_heatmap_plot = archr.archr_heatmap_plot
        File? share_atac_archr_arrow = archr.archr_arrow
        File? share_atac_archr_obj = archr.archr_raw_obj
        File? share_atac_archr_plots_zip = archr.plots_zip
    }
}
