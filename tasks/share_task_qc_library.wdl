version 1.0

# TASK
# SHARE-atac-lib-size-rna


task qc_library {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: ATAC library size task'
    }

    input {
        # This task computs the the library size for the library.

        File raw_counts
        File filtered_counts
        Int cutoff
        Int? memory_gb = 16
        String genome_name
        String? prefix
        #String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_qc_library"
        String docker_image = "nchernia/share_task_qc_library:5"        
        String assay
    }

    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    Int disk_gb = 50
    Float input_file_size_gb = size(filtered_counts, "G")
    Int mem_gb = memory_gb


    command {
        set -e

        # TODO remove the hard coded file paths from R scripts
        # TODO create only one R script that uses the parameters to discriminate
        # Estimate lib size
        # both
        #Rscript $(which lib_size_sc_V5_species_mixing.R)./ '${prefix + '.'}atac.${genome_name}' ${cutoff} atac --save
        # hg38/mm10
        Rscript $(which lib_size_sc_V5_single_species.R) ${raw_counts} ${filtered_counts} ${cutoff} ${genome_name} ${assay} --save

    }

    output {
        File lib_size_counts = glob('*.libsize.counts.csv')[0]
        File lib_size_log = glob('*.dups.log')[0]
        File plot = glob('*.jpg')[0]
    }

    runtime {
        #cpu : cpus
        memory : mem_gb+'G'
        disks : 'local-disk ${disk_gb} SSD'
        maxRetries : 0
        docker: docker_image
    }

    parameter_meta {
        raw_counts: {
                description: 'Barcode count csv',
                help: 'Barcode counts from raw bam in csv format.',
                example: 'raw.counts.csv'
            }
        filtered_counts: {
            description: 'Barcode count csv',
            help: 'Barcode counts from filtered bam in csv format.',
            example: 'filtered.counts.csv'
        }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.',
                examples: ['hg38', 'mm10', 'both']
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz',
                example: ['put link to gcr or dockerhub']
            }
    }


}
