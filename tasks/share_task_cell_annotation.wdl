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
        String? prefix = "prefix"

        # Reference genome
        String genome

        # Reference data name and id
        String reference_data_id
        String reference_data_name
        String reference_label

        # Query data
        File query_data
 
        # Down sampling cells for reference data
        Float? downsample_frac
 
        String? gene_id_to_symbol 

        # Integration
        String? anchors_reduction = "pcaproject"

        #Output option
        String? save_seurat

        # Docker image
        String? docker_image
        
        # Runtime parameter
        Float? memory_factor
        Float? disk_factor
    }
    
    # Determine the size of the input
    Float input_file_size_gb = size(query_data, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 64.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"
    
    #Output files
    String reference_h5ad = "${reference_data_name}.h5ad"
    String monitor_log = "cell_annotation_monitor.log"
    String notebook_log = "log/${prefix}.cell.annotation.logfile.${genome}.txt"
    String prediction = "${prefix}.cell.annotation.prediction.${genome}.csv"
    String prediction_labels = "${prefix}.cell.annotation.labels.${genome}.png"
    String prediction_scores = "${prefix}.cell.annotation.scores.${genome}.pdf"
    String seurat_object = "${prefix}.cell.annotation.${genome}.rds"


    command {
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &
        
        # Download data from cellxgene
        python3 $(which get_cellxgene_data.py) \
        --id ${reference_data_id} \
        --out ${reference_data_name} \
        --downsample_frac ${downsample_frac} \
        --reference_label ${reference_label}
        
        # Convert h5ad to Seurat object
        Rscript $(which h5ad_to_seurat.R) \
        --prefix ${prefix} \
        --reference_data_name ${reference_data_name} \
        --genome ${genome}
        
        # Perform cell annotation
        Rscript $(which cell_annotation.R) \
        --prefix ${prefix} \
        --reference_data_name ${reference_data_name} \
        --reference_label ${reference_label} \
        --query_data ${query_data} \
        --genome ${genome} \
        --gene_id_to_symbol ${gene_id_to_symbol} \
        --anchors_reduction ${anchors_reduction} \
        --save_seurat ${save_seurat}

    }

    output {
        File reference_h5ad = "${reference_h5ad}"
        File monitor_log = "${monitor_log}"
        File notebook_log = "${notebook_log}"
        File prediction = '${prediction}'
        File prediction_labels = '${prediction_labels}'
        File prediction_scores = '${prediction_scores}'
        File seurat_object = '${seurat_object}'
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
