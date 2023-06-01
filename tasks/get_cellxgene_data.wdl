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

        String? docker_image = "lzj1769/get_cellxgene_data"
    }
    
    # Determining memory size base on the size of the input files.
    Float mem_gb = 32.0

    # Determining disk size base on the size of the input files.
    Int disk_gb = 100.0

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    command {
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &
        
        # Download data from cellxgene
        python3 $(which get_cellxgene_data.py) \
        --dataset_id ${reference_data_id} --output_finame ${reference_data_name}

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

        reference_data_name: {
            description: 'Reference dataset name',
            help: 'String used to name the reference data.',
            examples: ['reference']
        }

        docker_image: {
            description: 'Docker image.',
            help: 'Docker image for preprocessing step.',
            example: ['put link to gcr or dockerhub']
        }
        
    }
}
