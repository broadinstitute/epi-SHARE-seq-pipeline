version 1.0

task get_cellxgene_data {
    meta {
        version: 'v0.1'
        author: 'Zhijian Li'
        affiliation: 'Broad Institute of MIT and Harvard'
        email: 'lizhijia@broadinstitute.org'
        description: 'SHARE-Seq pipeline: get data from cellxgene server.'    
    }

    input {
        #This tasks takes in an RNA matrix file, processes using Seurat and creates plots
        String reference_data_id
        String reference_data_name

        String docker_image = "lzj1769/get_cellxgene_data"
        Float? disk_factor = 0.1
    }
    
    # Determine the size of the input
    Float input_file_size_mb = size(query_data, "M")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 32.0

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(disk_factor * input_file_size_mb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    command {
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &
        
        # Download data from cellxgene
        python3 $(which get_cellxgene_data.py) \
        --dataset_id ${reference_data_id} --save_path ${reference_data_name}

    }

    output {
        File joint_calling_monitor = "monitoring.log"
        File joint_calling_log = "get_cellxgene_data.log"
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
    
        prefix: {
            description: 'Project name',
            help: 'String used to name your project and associated file names',
            example: "shareseq"
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
