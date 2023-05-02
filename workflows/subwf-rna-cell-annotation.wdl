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
        # RNA Seurat inputs
        String prefix
        String genome_name
        String? docker_image="lzj1769/cell-annotation"
        File reference_data
        File query_data
        
        #Seurat runtime parameters
        Float? disk_factor
        Float? memory_factor
    }

    call share_task_cell_annotation.cell_annotation as cell_annotation{
        input:
            prefix = prefix,
            reference_data = reference_data,
            query_data = query_data,
            genome_name = genome_name,
            docker_image = docker_image,
            disk_factor = disk_factor,
            memory_factor = memory_factor
    }

    output {
        File share_rna_cell_annotation_notebook_output = cell_annotation.notebook_output
        File share_rna_cell_annotation_notebook_log = cell_annotation.notebook_log
    }
}
