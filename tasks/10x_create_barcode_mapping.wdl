version 1.0

# TASK
# 10x_barcode_mapping

task mapping_tenx_barcodes {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: preprocess 10x ATAC data.'
    }

    input {
        # This task takes in input the 3 fastqs coming out from cellranger mkfastqs and preprocess them.
        File whitelist_atac # Barcode whitelist (chemistry specific)
        File whitelist_rna # Barcode whitelist (chemistry specific)

        Int? cpus = 16
        Float? disk_factor = 0.5
        Float? memory_factor = 0.15
        String? docker_image = "debian:bullseye-slim"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(whitelist_rna, "G") + size(whitelist_atac, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 5.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    String barcode_conversion_dict = "barcode_conversion_dict.csv"

    command <<<
        set -e

        if [ "$(zcat ~{whitelist_atac} | wc -l)" -eq "$(zcat ~{whitelist_rna} | wc -l)" ]; then
            zcat ~{whitelist_atac} | tr ACGTacgt TGCAtgca | rev | paste -d '\t' - <(zcat ~{whitelist_rna}) > ~{barcode_conversion_dict}
            paste -d '\t' <(zcat ~{whitelist_atac}) <(zcat ~{whitelist_rna}) >> ~{barcode_conversion_dict}
        fi
        # Fix for chromap.
        awk -v OFS="\t" '{print $2,$1}' barcode_conversion_dict.csv > temp
        mv temp ~{barcode_conversion_dict}
    >>>

    output {
        File? tenx_barcode_conversion_dict = barcode_conversion_dict
    }

    runtime {
        cpu: cpus
        docker: "${docker_image}"
        disks: "local-disk ${disk_gb} ${disk_type}"
        maxRetries: 1
        memory: "${mem_gb} GB"
        memory_retry_multiplier: 2
    }
}
