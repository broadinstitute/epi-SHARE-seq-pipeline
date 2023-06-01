version 1.0

# Import the tasks called by the pipeline
import "../tasks/get_cellxgene_data.wdl" as get_cellxgene_data
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
        String reference_data_id
        String reference_data_name
        String reference_label

        File query_data
        String genome

        String? prefix="prefix"
        
        # Docker images
        String? docker_image_get_cellxgene_data="lzj1769/get_cellxgene_data"
        String? docker_image_cell_annotation="lzj1769/cell-annotation"
    
        # Runtime parameters
        Float? disk_factor = 0.1
        Float? memory_factor = 0.15
    }

    call get_cellxgene_data.get_cellxgene_data as get_cellxgene_data{
        input:
            reference_data_id = reference_data_id,
            reference_data_name = reference_data_name,
            docker_image = docker_image_get_cellxgene_data,
            disk_factor = disk_factor
    }

    call share_task_cell_annotation.cell_annotation as cell_annotation{
        input:
            reference_data_id = reference_data_id,
            reference_data_name = reference_data_name,
            reference_label = reference_label,
            query_data = query_data,
            genome = genome,
            prefix = prefix,
            docker_image = docker_image_cell_annotation,
            disk_factor = disk_factor,
            memory_factor = memory_factor
    }

    output {
        File share_rna_cell_annotation_notebook_output = cell_annotation.notebook_output
        File share_rna_cell_annotation_notebook_log = cell_annotation.notebook_log
        File share_rna_cell_annotation_monitor_log = cell_annotation.monitor_log
        File share_rna_cell_annotation_prediction = cell_annotation.prediction
        File share_rna_cell_annotation_predicted_labels_plot = cell_annotation.predicted_labels_plot
    }
}
