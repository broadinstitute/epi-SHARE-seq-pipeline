version 1.0

# Import the tasks called by the pipeline
import "../tasks/share_task_cell_annotation.wdl" as share_task_cell_annotation

workflow wf_rna {
    meta {
        version: 'v0.1'
        author: 'Zhijian Li'
        affiliation: 'Broad Institute of MIT and Harvard'
        email: 'lizhijia@broadinstitute.org'
        description: 'SHARE-Seq pipeline: cell type annotation using RNA-seq data.'
    }

    input {
        String reference_data_id
        String reference_label

        File query_data
        String genome

        String? prefix="prefix"

        #Seurat runtime parameters
        Float? disk_factor = 0.1
        Float? memory_factor = 0.15
                
        String? docker_image="lzj1769/cell-annotation"
    }

    call share_task_cell_annotation.cell_annotation as cell_annotation{
        input:
            reference_data_id = reference_data_id,
            reference_label = reference_label,
            query_data = query_data,
            genome = genome,
            prefix = prefix,
            docker_image = docker_image,
            disk_factor = disk_factor,
            memory_factor = memory_factor
    }

    output {
        File share_rna_cell_annotation_notebook_output = cell_annotation.notebook_output
        File share_rna_cell_annotation_notebook_log = cell_annotation.notebook_log
        File share_rna_cell_annotation_monitor_log = cell_annotation.monitor_log
        File share_rna_cell_annotation_prediction = cell_annotation.prediction
        File share_rna_cell_annotation_predicted_labels_plot = cell_annotation.predicted_labels_plot
        File share_rna_cell_annotation_predicted_scores_plot = cell_annotation.predicted_scores_plot
    }
}
