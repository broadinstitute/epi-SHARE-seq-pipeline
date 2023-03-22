version 1.0

# TASK
# SHARE-joint-qc-plotting


task joint_qc_plotting {
    meta {
        version: 'v0.1'
        author: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Joint QC plot'
    }

    input {
        # This task generates a plot of barcodes QC'd jointly by RNA and ATAC metrics, as well as a
        # density plot of all barcodes passing at least one filter.
        File? atac_barcode_metadata
        File? rna_barcode_metadata
        Int remove_low_yielding_cells = 10
        Int min_umis = 100
        Int min_genes = 200
        Int min_tss = 4
        Int min_frags = 100

        Float? disk_factor = 8.0
        Float? memory_factor = 2.0

        String? prefix
        String genome_name

        String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_joint_qc"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(atac_barcode_metadata, "G") + size(rna_barcode_metadata, "G")

    # Determine memory size based on the size of the input files
    Float mem_gb = 5.0 + memory_factor * input_file_size_gb

    # Determine disk size based on the size of the input files
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    String joint_qc_plot = '${default="share-seq" prefix}.${genome_name}.joint.qc.plot.png'
    String joint_density_plot = '${default="share-seq" prefix}.${genome_name}.joint.density.plot.png'
    String joint_barcode_metadata = '${default="share-seq" prefix}.joint.barcode.metadata.${genome_name}.csv'

    command {
        set -e

        bash $(which monitor_script.sh) > monitoring.log &

        # Make joint qc plot
        python3 $(which joint_cell_plotting.py) ${default="share-seq" prefix} ${rna_barcode_metadata} ${atac_barcode_metadata} ${remove_low_yielding_cells} ${min_umis} ${min_genes} ${min_tss} ${min_frags} ${joint_qc_plot} ${joint_barcode_metadata}

        # Make joint density plot
        Rscript $(which joint_cell_plotting_density.R) ${default="share-seq" prefix} ${joint_barcode_metadata} ${joint_density_plot}
    }

    output {
        File joint_calling_monitor = "monitoring.log"
        File joint_calling_log = "joint_cell_plotting.log"
        File? joint_qc_plot = "${joint_qc_plot}"
        File? joint_density_plot = "${joint_density_plot}"
        File joint_barcode_metadata = "${joint_barcode_metadata}"
    }

    runtime {
        memory : "${mem_gb} GB"
        memory_retry_multiplier: 2
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker : docker_image
        maxRetries:1
    }

    parameter_meta {
        atac_barcode_metadata: {
                description: 'File containing ATAC barcode metrics.',
                help: 'tsv file with ATAC barcode (R1,R2,R3,PKR), fragments, TSS enrichment.',
                example: 'qc.atac.barcode.metadata.tsv'
            }
        rna_barcode_metadata: {
                description: 'File containing RNA barcode metrics.',
                help: 'tsv file with RNA barcode (R1,R2,R3,PKR), UMIs, genes.',
                example: 'qc.rna.barcode.metadata.tsv'
           }
        remove_low_yielding_cells: {
                description: 'UMI and fragments cutoff for plotting.',
                help: 'Minimum number of UMIs/fragments required for barcode to be plotted.',
                example: 10
           }
        min_umis: {
                description: 'UMI cutoff for RNA QC.',
                help: 'Minimum number of UMIs required for barcode to pass RNA QC.',
                example: 100
           }
        min_genes: {
                description: 'Gene cutoff for RNA QC.',
                help: 'Minimum number of genes required for barcode to pass RNA QC.',
                example: 200
           }
        min_tss: {
                description: 'TSS cutoff for ATAC QC.',
                help: 'Minimum TSS score required for barcode to pass ATAC QC.',
                example: 4
           }
        min_frags: {
                description: 'Fragments cutoff for ATAC QC.',
                help: 'Minimum number of fragments required for barcode to pass ATAC QC.',
                example: 100
           }
        prefix: {
                description: 'Prefix for output files',
                help: 'Prefix that will be used to name the output files',
                examples: 'MyExperiment'
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name genome reference used to align.',
                example: ['hg38', 'mm10', 'hg19', 'mm9']
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step.',
                example: ['put link to gcr or dockerhub']
            }
    }
}
