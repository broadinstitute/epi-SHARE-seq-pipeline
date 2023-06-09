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
        String peak_set
        Int? cpus = 4
        String? docker
        String? prefix
        Int? min_tss = 4
        Int? min_frags = 100
        Float? archr_disk_factor = 8.0
        Float? archr_memory_factor = 4.0
    }

    call share_task_archr.archr as archr{
        input:
            atac_frag = atac_fragments_filtered,
            genome = genome_name,
            peak_set = peak_set,
            min_tss = min_tss,
            min_frags = min_frags,
            doublet_k = 10,
            doublet_knn_method = "UMAP",
            lsi_method = 1,
            docker_image = docker,
            prefix = prefix,
            disk_factor = archr_disk_factor,
            memory_factor = archr_memory_factor
    }

    output {
        File share_atac_archr_notebook_output = archr.notebook_output
        File share_atac_archr_notebook_log = archr.notebook_log
        
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
