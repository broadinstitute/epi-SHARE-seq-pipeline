version 1.0

import "../tasks/share_task_correct_fastq.wdl" as share_task_correct_fastq
import "../tasks/share_task_starsolo.wdl" as share_task_starsolo
import "../tasks/share_task_generate_h5.wdl" as share_task_generate_h5
import "../tasks/share_task_qc_rna.wdl" as share_task_qc_rna
import "../tasks/share_task_log_rna.wdl" as share_task_log_rna
import "../tasks/share_task_seurat.wdl" as share_task_seurat

# Import the tasks called by the pipeline
workflow wf_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) and Sai Ma @ Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Sub-workflow to process the RNA portion of SHARE-seq libraries.'
    }

    input {
        # RNA Sub-workflow inputs

        # Correct
        File barcode_whitelist
        Int? correct_cpus = 1
        Float? correct_disk_factor = 8.0
        Float? correct_memory_factor = 0.15
        String? correct_docker_image

        # Align
        Array[File] read1
        Array[File] read2
        File idx_tar
        String chemistry
        String genome_name
        String prefix
        String? barcode_tag
        String? pkr
        Int? cpus = 16
        File? whitelist
        String? docker
        # QC
        Int? umi_cutoff
        Int? gene_cutoff

        # Seurat
        Boolean count_only = false

        #Seurat filtering parameters
        Int? rna_seurat_min_features
        Float? rna_seurat_percent_mt
        Int? rna_seurat_min_cells

        #Seurat UMAP
        Int? rna_seurat_umap_dim
        Float? rna_seurat_umap_resolution

        # Seurat runtime parameters
        Float? rna_seurat_disk_factor
        Float? rna_seurat_memory_factor

    }

    scatter (read_pair in zip(read1, read2)) {
        call share_task_correct_fastq.share_correct_fastq as correct {
            input:
                fastq_R1 = read_pair.left,
                fastq_R2 = read_pair.right,
                barcode_whitelist = barcode_whitelist,
                sample_type = "RNA",
                pkr = pkr,
                prefix = prefix,
                cpus = correct_cpus,
                disk_factor = correct_disk_factor,
                memory_factor = correct_memory_factor,
                docker_image = correct_docker_image
        }
    }

    call share_task_starsolo.share_rna_align as align {
        input:
            chemistry = chemistry,
            fastq_R1 = correct.corrected_fastq_R1,
            fastq_R2 = correct.corrected_fastq_R2,
            whitelist = whitelist,
            genome_name = genome_name,
            genome_index_tar = idx_tar,
            prefix = prefix,
            cpus = cpus
    }

    call share_task_generate_h5.generate_h5 as generate_h5 {
        input:
            tar = align.raw_tar,
            genome_name = genome_name,
            prefix = prefix,
            pkr = pkr
    }

    call share_task_qc_rna.qc_rna as qc_rna {
        input:
            bam = align.output_bam,
            umi_cutoff = umi_cutoff,
            gene_cutoff = gene_cutoff,
            pkr = pkr,
            barcode_tag = barcode_tag,
            genome_name = genome_name,
            prefix = prefix
    }

    call share_task_log_rna.log_rna as log_rna {
       input:
           alignment_log = align.log_final_out,
           dups_log = qc_rna.rna_duplicates_log
    }

    if (!count_only) {
        call share_task_seurat.seurat as seurat {
            input:
                rna_matrix = generate_h5.h5_matrix,
                genome_name = genome_name,
                min_features = rna_seurat_min_features,
                percent_mt = rna_seurat_percent_mt,
                min_cells = rna_seurat_min_cells,
                umap_dim = rna_seurat_umap_dim,
                umap_resolution = rna_seurat_umap_resolution,
                prefix = prefix,
                disk_factor = rna_seurat_disk_factor,
                memory_factor = rna_seurat_memory_factor
        }
    }

    output {
        File share_task_starsolo_output_bam = align.output_bam
        File share_rna_alignment_log = align.log_final_out
        File share_task_starsolo_log_out = align.log_out
        File share_task_starsolo_log_progress_out = align.log_progress_out
        File share_task_starsolo_output_sj = align.output_sj
        File share_task_starsolo_barcodes_stats = align.barcodes_stats
        File share_task_starsolo_features_stats = align.features_stats
        File share_task_starsolo_summary_csv = align.summary_csv
        File share_task_starsolo_umi_per_cell = align.umi_per_cell
        File share_task_starsolo_raw_tar = align.raw_tar

        File share_rna_h5 = generate_h5.h5_matrix

        File share_rna_barcode_metadata  = qc_rna.rna_barcode_metadata
        File share_rna_duplicates_log = qc_rna.rna_duplicates_log
        File? share_rna_umi_barcode_rank_plot = qc_rna.rna_umi_barcode_rank_plot
        File? share_rna_gene_barcode_rank_plot = qc_rna.rna_gene_barcode_rank_plot
        File? share_rna_gene_umi_scatter_plot = qc_rna.rna_gene_umi_scatter_plot

        File? share_rna_seurat_notebook_output = seurat.notebook_output
        File? share_rna_seurat_notebook_log = seurat.notebook_log
        File? share_rna_seurat_raw_violin_plot = seurat.seurat_raw_violin_plot
        File? share_rna_seurat_filtered_violin_plot = seurat.seurat_filtered_violin_plot
        File? share_rna_seurat_raw_qc_scatter_plot = seurat.seurat_raw_qc_scatter_plot
        File? share_rna_seurat_filtered_qc_scatter_plot = seurat.seurat_filtered_qc_scatter_plot
        File? share_rna_seurat_variable_genes_plot = seurat.seurat_variable_genes_plot
        File? share_rna_seurat_PCA_dim_loadings_plot = seurat.seurat_PCA_dim_loadings_plot
        File? share_rna_seurat_PCA_plot = seurat.seurat_PCA_plot
        File? share_rna_seurat_heatmap_plot = seurat.seurat_heatmap_plot
        File? share_rna_seurat_jackstraw_plot = seurat.seurat_jackstraw_plot
        File? share_rna_seurat_elbow_plot = seurat.seurat_elbow_plot
        File? share_rna_seurat_umap_cluster_plot = seurat.seurat_umap_cluster_plot
        File? share_rna_seurat_umap_rna_count_plot = seurat.seurat_umap_rna_count_plot
        File? share_rna_seurat_umap_gene_count_plot = seurat.seurat_umap_gene_count_plot
        File? share_rna_seurat_umap_mito_plot = seurat.seurat_umap_mito_plot
        File? share_rna_seurat_obj = seurat.seurat_filtered_obj
        File? share_rna_plots_zip = seurat.plots_zip

        Int share_rna_total_reads = log_rna.rna_total_reads
        Int share_rna_aligned_uniquely = log_rna.rna_aligned_uniquely
        Int share_rna_aligned_multimap = log_rna.rna_aligned_multimap
        Int share_rna_unaligned = log_rna.rna_unaligned
        Int share_rna_feature_reads = log_rna.rna_feature_reads
        Int share_rna_duplicate_reads = log_rna.rna_duplicate_reads
    }
}
