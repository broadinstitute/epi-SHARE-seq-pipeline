version 1.0

import "../tasks/share_task_seurat.wdl" as share_task_seurat
import "../tasks/share_task_starsolo.wdl" as share_task_starsolo

# Import the tasks called by the pipeline
workflow wf_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) and Sai Ma @ Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Sub-workflow to process the RNA portion of SHARE-seq libraries.'
    }

    input {
        # RNA Sub-workflow inputs

        # Align
        Array[File] read1
        Array[File] read2
        File idx_tar
        String prefix = "shareseq-project"
        String genome_name
        Int? cpus = 16
        String? docker
        File whitelist
        # Seurat
        Int umap_dim = 10
        Float umap_resolution = 0.5 
   }

    call share_task_starsolo.share_rna_align as align {
        input:
            fastq_R1 = read1,
            fastq_R2 = read2,
            genome_name = genome_name,
            genome_index_tar = idx_tar,
            prefix = prefix,
            whitelist = whitelist,
            cpus = cpus
    }

    call share_task_seurat.seurat as seurat{
        input:
            rna_matrix = align.raw_tar,
            genome_name = genome_name,
            umap_dim = umap_dim,
            umap_resolution = umap_resolution,
            prefix = prefix
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

        File share_rna_seurat_notebook_output = seurat.notebook_output
        File share_rna_seurat_notebook_log = seurat.notebook_log
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
        File? share_rna_seurat_raw_h5 = seurat.seurat_raw_matrix
        File? share_rna_seurat_seurat_barcode_metadata = seurat.seurat_barcode_metadata 
    }
}
