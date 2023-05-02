version 1.0

# TASK
# SHARE-rna-generate-h5


task generate_h5 {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: RNA gene x cell matrix'
    }

    input {
        # This task computes the the gene x cell matrix.
        File tar
        String genome_name
        String? pkr
        String prefix
        Boolean? ensembl = false

        Float? disk_factor = 8.0
        Float? memory_factor = 2.0
        String? docker_image = "us.gcr.io/buenrostro-share-seq/share_task_generate_h5"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(tar, "G")

    # Determining memory size based on the size of the input files.
    Float mem_gb = 10.0 + memory_factor * input_file_size_gb

    # Determining disk size based on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type based on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    String h5 = "${default="share-seq" prefix}.${genome_name}.rna.h5"
    String monitor_log = "monitor.log"

    command {
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &

        # Untar
        tar xzvf ${tar}

        # Generate h5 file
        python3 $(which generate_h5_rna.py) \
            ./matrix.mtx.gz \
            ./features.tsv.gz \
            ./barcodes.tsv.gz \
            ${h5} \
            ${pkr} \
            $(if ensembl then "--ensembl" else "")
    }

    output {
        File h5_matrix = "${h5}"
    }

    runtime {
        memory : "${mem_gb} GB"
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker : "${docker_image}"
    }

    parameter_meta {
        tar: {
                description: 'STARsolo output tar.gz file',
                help: 'tar.gz file containing raw matrix, features, and barcodes file from STARsolo.',
                example: 'raw.tar.gz'
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name genome reference used to align.',
                example: ['hg38', 'mm10', 'hg19', 'mm9']
            }
        prefix: {
                description: 'Prefix for output files',
                help: 'Prefix that will be used to name the output files.',
                example: 'MyExperiment'
            }
	ensembl: {
                description: 'Use ENSEMBL gene IDs in h5 matrix',
                help: 'Boolean for if ENSEMBL gene IDs should be outputted in h5 matrix (rather than gene names).'
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install h5py scipy',
                example: ['put link to gcr or dockerhub']
            }
    }
}
