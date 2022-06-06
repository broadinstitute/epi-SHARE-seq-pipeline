version 1.0

# Import the tasks called by the pipeline
import "../tasks/share_task_archr.wdl" as share_task_archr


workflow wf_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) and Sai Ma @ Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Sub-workflow to process the ATAC portion of SHARE-seq libraries.'
    }

    input {
        # ATAC Sub-worflow inputs
        File atac_fragments_filtered
        String genome_name
        Int? cpus = 4
        String? docker
    }

    call share_task_archr.archr as archr{
        input:
            atac_frag = atac_fragments_filtered,
            genome = genome_name,
            min_tss = 4,
            min_frags = 100,
            doublet_k = 10,
            doublet_knn_method = "UMAP",
            lsi_method = 1,
            docker_image = docker
    }

    output {
        File share_atac_archr_notebook_output = archr.notebook_output
        File share_atac_notebook_log = archr.notebook_log
        File share_atac_papermill_log = archr.papermill_log
        File? share_atac_archr_gene_heatmap_plot = archr.archr_heatmap_plot
        File? share_atac_archr_tss_enrichment_raw = archr.archr_TSS_uniq_frags_plot
        File? share_atac_archr_tss_enrichment_filtered = archr.archr_TSS_uniq_frags_filtered_plot
        File? share_atac_archr_fragment_size_plot = archr.archr_fragment_size_dist_plot
        File? share_atac_archr_doublet_plot = archr.archr_doublet_plot
        File? share_atac_archr_umap_plot = archr.archr_umap_plot
        File? share_atac_archr_arrow = archr.archr_arrow
        File? share_atac_archr_obj = archr.archr_obj
        File? share_atac_archr_plots_zip = archr.plots_zip
    }
}
