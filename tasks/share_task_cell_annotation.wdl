version 1.0

task cell_annotation {
    meta {
        version: 'v0.1'
        author: 'Zhijian Li'
        affiliation: 'Broad Institute of MIT and Harvard'
        email: 'lizhijia@broadinstitute.org'
        description: 'SHARE-Seq pipeline: cell type annotation using RNA-seq data.'    
    }

    input {
        # 
        String genome

        # Input
        String reference_data_id
        String reference_data_name
        String reference_label
        File query_data


        # Output
        String prefix
        
        String output_filename = "${prefix}.rna.cell.annotation.notebook.${genome}.ipynb"
        String log_filename = "log/${prefix}.rna.cell.annotation.logfile.${genome}.txt"

        # Docker image
        String docker_image
        
        # Runtime parameter
        Float? disk_factor = 0.1
        Float? memory_factor = 0.15

        String? papermill = "TRUE"
    }
    
    # Determine the size of the input
    Float input_file_size_mb = size(query_data, "M")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 32.0 + memory_factor * input_file_size_mb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(disk_factor * input_file_size_mb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"
    
    String monitor_log = "rna_cell_annotation_monitor.log"

    #Plot filepaths
    String plots_filepath = '${prefix}.rna.cell.annotation.plots.${genome}'
    String predicted_labels_plot = '${plots_filepath}/${prefix}.rna.cell.annotation.predicted.labels.${genome}.png'
    String predicted_scores_plot = '${plots_filepath}/${prefix}.rna.cell.annotation.predicted.scores.${genome}.png'

    #Other filepaths
    String prediction = '${prefix}.rna.cell.annotation.prediction.${genome}.csv'

    command {
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &
        
        papermill $(which cell_annotation_notebook.ipynb) ${output_filename} \
        -p reference_data_id ${reference_data_id} \
        -p reference_label ${reference_label} \
        -p query_data ${query_data} \
        -p genome ${genome} \
        -p prefix ${prefix} \
        -p papermill ${papermill}
    }

    output {
        File notebook_output = output_filename
        File notebook_log = log_filename
        File monitor_log = monitor_log
        File prediction = prediction
        File predicted_labels_plot = predicted_labels_plot
        File predicted_scores_plot = predicted_scores_plot
    }

    runtime {
        memory : "${mem_gb} GB"
        memory_retry_multiplier: 2
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker : "${docker_image}"
        maxRetries:1
    }

    parameter_meta {
        reference_data_id: {
            description: 'Reference dataset id',
            help: 'The dataset id from cellxgene data base.',
            examples: ['3bbb6cf9-72b9-41be-b568-656de6eb18b5']
        }
        
        query_data: {
            description: 'Query data',
            help: 'scRNA-seq data used as query',
            examples: ['put link to gcr']
        }

        genome: {
            description: 'Reference name',
            help: 'Reference genome.',
            examples: ['hg38', 'mm10', 'hg19', 'mm9']
        }

        prefix: {
            description: 'Project name',
            help: 'String used to name your project and associated file names',
            example: "shareseq"
        }

        papermill: {
            description: 'Boolean papermill flag',
            help: 'Flag to notebook run in papermill mode',
            example: 'TRUE'
        }

        output_filename: {
            description: 'Output jupyter notebook name',
            help: 'The name assigned to output jupyter notebook',
            examples: 'output.ipynb'
        }

        docker_image: {
            description: 'Docker image.',
            help: 'Docker image for preprocessing step.',
            example: ['put link to gcr or dockerhub']
        }
        
        disk_factor: {
            description: 'Disk factor',
            help: 'Multiply this value to input .h5 file size (MB) to determine disk space (GB)',
            example: 16.0
        }
        
        memory_factor: {
            description: 'Memory factor',
            help: 'Multiply this value to input .h5 file size (MB) and add to default 32GB memory to determine RAM (GB)',
            example: 1.0
        }
    }
}
