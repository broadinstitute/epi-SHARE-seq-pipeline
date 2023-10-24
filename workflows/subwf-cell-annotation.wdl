version 1.0

# Import the tasks called by the pipeline
import "../tasks/share_task_cell_annotation.wdl" as share_task_cell_annotation

workflow wf_cell_annotation {
    meta {
        version: 'v0.1'
        author: 'Zhijian Li'
        affiliation: 'Broad Institute of MIT and Harvard'
        email: 'lizhijia@broadinstitute.org'
        description: 'SHARE-Seq pipeline: cell type annotation using RNA-seq data.'
    }

    input {
        # Sample name
        String? prefix="prefix"

        # Reference genome
        String genome

        # Reference data for cell annotation
        String reference_data_id
        String reference_data_name
        String reference_label

        # Set true if the reference data uses gene id as feature name. 
        # This is usually true for data downloaded from cellxgene server
        String? gene_id_to_symbol = "TRUE"
        
        String? save_seurat = "TRUE"

        # Query data
        File query_data
        
        # Proportion for down-sampling
        Float? downsample_frac = 1
    
        # Integration
        String? normalization_method = "SCT"
        String? anchors_reduction = "pcaproject"
    
        # Docker images
        String? docker_image="lzj1769/cell_annotation"
    
        # Runtime parameters
        Float? memory_factor = 5
        Float? disk_factor = 10
    }

    call share_task_cell_annotation.cell_annotation as cell_annotation{
        input:
            reference_data_id = reference_data_id,
            reference_data_name = reference_data_name,
            reference_label = reference_label,
            query_data = query_data,
            genome = genome,
            gene_id_to_symbol = gene_id_to_symbol,
            prefix = prefix,
            docker_image = docker_image,
            disk_factor = disk_factor,
            memory_factor = memory_factor,
            downsample_frac = downsample_frac,
            normalization_method = normalization_method,
            anchors_reduction = anchors_reduction,
            save_seurat = save_seurat

    }

    output {
        File share_cell_annotation_reference_h5ad = cell_annotation.reference_h5ad
        File share_cell_annotation_notebook_log = cell_annotation.notebook_log
        File share_cell_annotation_monitor_log = cell_annotation.monitor_log
        File share_cell_annotation_prediction = cell_annotation.prediction
        File share_cell_annotation_prediction_labels = cell_annotation.prediction_labels
        File share_cell_annotation_prediction_scores = cell_annotation.prediction_scores
        File share_cell_annotation_seurat_object = cell_annotation.seurat_object
    }
}
