version 1.0

# TASK
# SHARE-joint-cell-calling


task joint_cell_calling {
    meta {
        version: 'v0.1'
        author: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Joint cell calling plot'
    }

    input {
        # This task generates a plot of barcodes QC'd jointly by RNA and ATAC metrics.
        String prefix
        File? atac_barcode_metadata
        File? rna_barcode_metadata
	Int min_UMIs = 100
        Int min_genes = 200
        Int min_TSS = 4
        Int min_frags = 100
        Int disk_gb = 50
        Int mem_gb = 64
        String genome_name
        String docker_image = "mknudson/share_task_joint_cell_calling:test"
    }

    Int dsk_gb = disk_gb
    Int memory_gb = 64
    String barcode_metadata = '${default="share-seq" prefix}.joint.barcode.metadata.${genome_name}.csv'

    command {
        set -e

        bash $(which monitor_script.sh) > monitoring.log &

        python3 $(which joint_cell_plotting.py) ${prefix} ${rna_barcode_metadata} ${atac_barcode_metadata} ${min_UMIs} ${min_genes} ${min_TSS} ${min_frags} ${barcode_metadata}

        Rscript $(which joint_cell_plotting_density.R) ${prefix} ${barcode_metadata} ${min_UMIs} ${min_genes} ${min_TSS} ${min_frags}
    }

    output {
        File joint_calling_monitor = "monitoring.log"
        File joint_cell_plot = glob('*_joint_cell_plot.png')[0]
        File joint_cell_density_plot = glob('*_joint_cell_density_plot.png')[0]
    }

    runtime {
        #cpu : cpus
        memory : memory_gb+'G'
        disks : 'local-disk ${dsk_gb} SSD'
        maxRetries : 0
        docker: docker_image
    }

    parameter_meta {
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz',
                example: ['put link to gcr or dockerhub']
            }
    }
}
