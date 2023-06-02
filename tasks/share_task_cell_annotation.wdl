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
        # Sample or project name
        String prefix

        # Reference genome
        String genome

        # Reference data name and id
        String reference_data_id
        String reference_data_name
        String reference_label

        # Query data
        File query_data
 
        # Docker image
        String? docker_image
        
        # Runtime parameter
        Float? disk_factor
        Float? memory_factor

        String? gene_id_to_symbol = "TRUE"
    }
    
    # Determine the size of the input
    Float input_file_size_mb = size(query_data, "M")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 32.0 + memory_factor * input_file_size_mb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(disk_factor * input_file_size_mb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"
    
    #Plot filepaths
    String plots_filepath = '${prefix}.rna.cell.annotation.plots.${genome}'

    command {
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &
        
        # Download data from cellxgene
        python3 $(which get_cellxgene_data.py) ${reference_data_id} ${reference_data_name}


        # Perform cell annotation
        Rscript $(which cell_annotation.R) \
        ${prefix} \
        ${reference_data_name} \
        ${reference_label} \
        ${query_data} \
        ${genome} \
        ${gene_id_to_symbol}
        
        # papermill $(which cell_annotation_notebook.ipynb) ${output_filename} \
        # -p reference_data_name ${reference_data_name} \
        # -p reference_label ${reference_label} \
        # -p query_data ${query_data} \
        # -p genome ${genome} \
        # -p prefix ${prefix} \
        # -p papermill ${papermill}
    }

    output {
        File reference_h5ad = "${reference_data_name}"
        File monitor_log = "cell_annotation_monitor.log"
        File notebook_log = "log/${prefix}.cell.annotation.logfile.${genome}.txt"
        File prediction = '${prefix}.cell.annotation.prediction.${genome}.csv'
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
            help: 'The dataset id from cellxgene server.',
            examples: ['3bbb6cf9-72b9-41be-b568-656de6eb18b5']
        }

        reference_data_name: {
            description: 'Reference data',
            help: 'This file will be used as reference',
            examples: ['reference.h5ad']
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
