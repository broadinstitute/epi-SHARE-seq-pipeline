version 1.0

# Import the tasks called by the pipeline
import "../tasks/share_task_bowtie2.wdl" as share_task_align
import "../tasks/share_task_bam2bed.wdl" as share_task_bam2bed
import "../tasks/share_task_count_atac.wdl" as share_task_count
import "../tasks/share_task_qc_atac.wdl" as share_task_qc_atac
import "../tasks/share_task_qc_library.wdl" as share_task_qc_library
import "../tasks/share_task_archr.wdl" as share_task_archr
import "../tasks/share_task_log_atac.wdl" as share_task_log_atac

workflow wf_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) and Sai Ma @ Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Sub-workflow to process the ATAC portion of SHARE-seq libraries.'
    }

    input {
        # ATAC Sub-worflow inputs
        Array[File] read1
        Array[File] read2

        File chrom_sizes
        File idx_tar
        File tss_bed
        File peak_set

        # Parameters for bam2bed
        Int? bam2bed_cpus = 16
        Float? bam2bed_disk_factor = 8.0
        Float? bam2bed_memory_factor = 0.15

        String prefix = "shareseq-project"
        String genome_name
        Int? cutoff
        Int? cpus = 4
        String? docker
    }

    call share_task_align.share_atac_align as align {
        input:
            fastq_R1 = read1,
            fastq_R2 = read2,
            genome_name = genome_name,
            genome_index = idx_tar,
            prefix = prefix
    }

    call share_task_bam2bed.share_atac_bam2bed as bam2bed {
        input:
            cpus = bam2bed_cpus,
            bam = align.atac_alignment,
            bam_index = align.atac_alignment_index,
            disk_factor = bam2bed_disk_factor,
            memory_factor = bam2bed_memory_factor,
            genome_name = genome_name,
            chrom_sizes = chrom_sizes,
            prefix = prefix
    }

    call share_task_count.count_reads_atac as count {
        input:
            cutoff = cutoff,
            fragments_raw = bam2bed.atac_fragments_raw,
            genome_name = genome_name,
            prefix = prefix,
            cpus = cpus
    }

    call share_task_qc_library.qc_library as qc_library{
        input:
            raw_counts = count.atac_counts_unfiltered,
            filtered_counts = count.atac_counts_filtered,
            cutoff = cutoff,
            genome_name = genome_name,
            prefix = prefix,
            assay = "ATAC"
    }

    call share_task_qc_atac.qc_atac as qc_atac{
        input:
            raw_bam = align.atac_alignment,
            raw_bam_index = align.atac_alignment_index,
            filtered_bam = bam2bed.atac_alignment_filtered,
            filtered_bam_index = bam2bed.atac_alignment_filtered_index,
            tss = tss_bed,
            genome_name = genome_name,
            prefix = prefix,
            cpus = cpus
    }
    call share_task_log_atac.log_atac as log_atac {
       input:
           alignment_log = align.atac_alignment_log,
           dups_log = qc_library.lib_size_log
    }
    call share_task_archr.archr as archr{
        input:
            atac_frag = count.atac_fragments_filtered,
            genome = genome_name,
            peak_set = peak_set,
            prefix = prefix
    }
    call share_task_archr.archr as archr_strict{
        input:
            atac_frag = count.atac_fragments_filtered,
            genome = genome_name,
            peak_set = peak_set,
            prefix = '${prefix}_strict',
            min_tss = 5,
            min_frags = 1000
    }

    output {
        File share_atac_alignment_raw = align.atac_alignment
        File share_atac_alignment_raw_index = align.atac_alignment_index
        File share_atac_alignment_log = align.atac_alignment_log

        File share_atac_alignment_filtered = bam2bed.atac_alignment_filtered
        File share_atac_alignment_filtered_index = bam2bed.atac_alignment_filtered_index
        File share_atac_fragments_raw = bam2bed.atac_fragments_raw


        File share_atac_barcodes = count.atac_barcodes
        File share_atac_fragments_filtered = count.atac_fragments_filtered
        File share_atac_counts_raw = count.atac_counts_unfiltered
        File share_atac_counts_filtered = count.atac_counts_filtered

        File share_atac_qc_library_counts = qc_library.lib_size_counts
        File share_atac_qc_library_duplicates = qc_library.lib_size_log
        File share_atac_qc_library_plot = qc_library.plot

        File share_atac_qc_final = qc_atac.atac_final_stats
        File share_atac_qc_hist_plot = qc_atac.atac_final_hist_png
        File share_atac_qc_hist_txt = qc_atac.atac_final_hist
        File share_atac_qc_tss_enrichment = qc_atac.atac_tss_pileup_png

        File share_atac_archr_notebook_output = archr.notebook_output
        File share_atac_archr_notebook_log = archr.notebook_log
        File share_atac_archr_barcode_metadata = archr.archr_barcode_metadata
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
        Int share_atac_total_reads = log_atac.atac_total_reads
        Int share_atac_aligned_uniquely = log_atac.atac_aligned_uniquely
        Int share_atac_unaligned = log_atac.atac_unaligned
        Int share_atac_feature_reads = log_atac.atac_feature_reads
        Int share_atac_duplicate_reads = log_atac.atac_duplicate_reads

        File share_atac_archr_strict_notebook_output = archr_strict.notebook_output
        File share_atac_archr_strict_notebook_log = archr_strict.notebook_log

        File? share_atac_archr_strict_raw_tss_enrichment = archr_strict.archr_raw_tss_by_uniq_frags_plot
        File? share_atac_archr_strict_filtered_tss_enrichment = archr_strict.archr_filtered_tss_by_uniq_frags_plot
        File? share_atac_archr_strict_raw_fragment_size_plot = archr_strict.archr_raw_frag_size_dist_plot
        File? share_atac_archr_strict_filtered_fragment_size_plot = archr_strict.archr_filtered_frag_size_dist_plot

        File? share_atac_archr_strict_umap_doublets = archr_strict.archr_umap_doublets
        File? share_atac_archr_strict_umap_cluster_plot = archr_strict.archr_umap_cluster_plot
        File? share_atac_archr_strict_umap_num_frags_plot = archr_strict.archr_umap_num_frags_plot
        File? share_atac_archr_strict_umap_tss_score_plot = archr_strict.archr_umap_tss_score_plot
        File? share_atac_archr_strict_umap_frip_plot = archr_strict.archr_umap_frip_plot

        File? share_atac_archr_strict_gene_heatmap_plot = archr_strict.archr_heatmap_plot
        File? share_atac_archr_strict_arrow = archr_strict.archr_arrow
        File? share_atac_archr_strict_obj = archr_strict.archr_raw_obj
        File? share_atac_archr_strict_plots_zip = archr_strict.plots_zip

    }
}
