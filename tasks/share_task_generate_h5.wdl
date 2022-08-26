version 1.0

# TASK
# SHARE-atac-generate-h5


task generate_h5 {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: RNA gene cell matrix'
    }

    input {
        # This task computs the the gene by barcode matrix.

        File filtered_bed
        #String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_generate_h5"
        String docker_image = "swekhande/shareseq-prod:share-task-generate-h5"


    }

    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    Int disk_gb = 50
    Float input_file_size_gb = size(filtered_bed, "G")
    Int mem_gb = 32


    command {
        set -e
        Rscript $(which UMI_gene_perCell_plot_v2.R) ${filtered_bed} --save

    }

    output {
        File h5_matrix = glob('*.h5')[0]
        Array[File] umi_qc_plots = glob('*.pdf')
    }

    runtime {
        #cpu : cpus
        memory : mem_gb+'G'
        disks : 'local-disk ${disk_gb} SSD'
        maxRetries : 0
        docker: docker_image
    }

    parameter_meta {
            filtered_bed: {
            description: 'Barcode count csv',
            help: 'Barcode counts from filtered bam in csv format.',
            example: 'filtered.counts.csv'
        }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz',
                example: ['put link to gcr or dockerhub']
            }
    }
}
