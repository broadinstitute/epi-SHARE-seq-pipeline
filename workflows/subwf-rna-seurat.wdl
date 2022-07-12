version 1.0


# Import the tasks called by the pipeline
import "../tasks/share_task_seurat.wdl" as share_task_seurat

workflow wf_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) and Sai Ma @ Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Sub-workflow to process the RNA portion of SHARE-seq libraries.'
    }

    input {
        # RNA Sub-worflow inputs

        # Align
        String prefix = "shareseq-project"
        String genome_name
        Int? cpus = 16
        String? docker
        File h5_matrix
        # Seurat
        Int umap_dim = 10
        Float umap_resolution = 0.5
    }

    call share_task_seurat.seurat as seurat{
        input:
            rna_matrix = h5_matrix,
            genome_name = genome_name,
            umap_dim = umap_dim,
            umap_resolution = umap_resolution,
            prefix = prefix,
            docker_image = docker
    }

    output {
        File share_rna_seurat_notebook_output = seurat.notebook_output
        File share_rna_seurat_notebook_log = seurat.notebook_log
        File share_rna_seurat_papermill_log = seurat.papermill_log
        File? share_rna_seurat_violin_plot = seurat.seurat_violin_plot
        File? share_rna_seurat_mitochondria_qc_plot = seurat.seurat_mitochondria_qc_plot
        File? share_rna_seurat_features_plot = seurat.seurat_features_plot
        File? share_rna_seurat_PCA_dim_loadings_plot = seurat.seurat_PCA_dim_loadings_plot
        File? share_rna_seurat_PCA_plot = seurat.seurat_PCA_plot
        File? share_rna_seurat_heatmap_plot = seurat.seurat_heatmap_plot
        File? share_rna_seurat_jackstraw_plot = seurat.seurat_jackstraw_plot
        File? share_rna_seurat_elbow_plot = seurat.seurat_elbow_plot
        File? share_rna_seurat_umap_plot = seurat.seurat_umap_plot
        File? share_rna_seurat_obj = seurat.seurat_obj
        File? share_rna_plots_zip = seurat.plots_zip
    }
}
