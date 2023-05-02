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
        #This tasks takes in an RNA matrix file, processes using Seurat and creates plots
        File reference_data
        File query_data

        String? genome_name = "hg38"

        String? normalization_method = "LogNormalize"
        Float? normalization_scale_factor = 10000

        String? variable_features_method = "vst"
        Int? variable_features_num = 2000

        String? weight_reduction = "pca"
        Int? n_dims = 30

        String prefix = "prefix"
        Int? threads = 8

        String papermill = "TRUE"
        
        String output_filename = "${prefix}.rna.cell.annotation.notebook.${genome_name}.ipynb"
        String log_filename = "log/${prefix}.rna.cell.annotation.logfile.${genome_name}.txt"

        #Int mem_gb = 128
        String docker_image = "lzj1769/cell-annotation"
        
        Float? disk_factor = 0.1
        Float? memory_factor = 0.15
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
    String plots_filepath = '${prefix}.rna.cell.annotation.plots.${genome_name}'

    #Other filepaths
    String filtered_seurat_h5 = '${prefix}.rna.seurat.filtered_matrix.${genome_name}.h5'
    String barcode_metadata = '${prefix}.rna.seurat.barcode_metadata.${genome_name}.tsv'
    String plots_zip_dir = '${plots_filepath}.zip'
    #String papermill_log_filename = 'papermill.logfile.txt'

    command {
    
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &
        
        papermill $(which cell_annotation_notebook.ipynb) ${output_filename} \
        -p reference_data ${reference_data} \
        -p query_data ${query_data} \
        -p genome ${genome_name} \
        -p normalization_method ${normalization_method} \
        -p normalization_scale_factor ${normalization_scale_factor} \
        -p variable_features_method ${variable_features_method} \
        -p variable_features_num ${variable_features_num} \
        -p weight_reduction ${weight_reduction} \
        -p n_dims ${n_dims} \
        -p prefix ${prefix} \
        -p papermill ${papermill}
    }


    output {
        File notebook_output = output_filename
        File notebook_log = log_filename
        File? plots_zip = plots_zip_dir
        File? cell_annotation_monitor_log = monitor_log
    }

    runtime {
        memory : "${mem_gb} GB"
        memory_retry_multiplier: 2
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker : "${docker_image}"
        maxRetries:1
    }

    parameter_meta {
        papermill: {
            description: 'Boolean papermill flag',
            help: 'Flag to notebook run in papermill mode',
            example: 'TRUE'
        }

        query_data: {
            description: 'Query data',
            help: 'The name genome_name reference used to align.',
            examples: ['hg38', 'mm10', 'hg19', 'mm9']
        }

        reference_data: {
            description: 'Reference data',
            help: 'The name genome_name reference used to align.',
            examples: ['hg38', 'mm10', 'hg19', 'mm9']
        }

        genome_name: {
            description: 'Reference name',
            help: 'The name genome_name reference used to align.',
            examples: ['hg38', 'mm10', 'hg19', 'mm9']
        }

        normalization_method: {
            description: 'Normalization method used in Seurat',
            help: 'Seurat normalization method used in Seurat::NormalizeData()',
            examples: ["LogNormalize","CLR","RC"]
        }

        normalization_scale_factor: {
            description: 'Scaling factor used in Seurat normalization',
            help: 'Scaling factor parameter used in Seurat::NormalizeData()',
            example: 10000
        }

        variable_features_method: {
            description: 'Method used to select variable features',
            help: 'Parameter used in Seurat::FindVariableFeatures()',
            example: "vst"
        }

        variable_features_num: {
            description: 'Number of variable features used to find',
            help: 'Parameter used in Seurat::FindVariableFeatures()',
            example: 2000
        }

        papermill: {
            description: 'Boolean papermill flag',
            help: 'Flag to notebook run in papermill mode',
            example: 'TRUE'
        }

        prefix: {
            description: 'Project name',
            help: 'String used to name your project and associated file names',
            example: "shareseq"
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
            example: 0.15
        }
    }
}
