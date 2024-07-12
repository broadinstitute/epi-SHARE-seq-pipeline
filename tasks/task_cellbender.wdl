version 1.0

# TASK
# cellbender

task cellbender {
    meta {
        version: 'v0.1'
        source: 'Broad Institute Data Sciences Platform (copyright 2020); https://portal.firecloud.org/#methods/cellbender/remove-background/'
        modified_by: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: remove ambient RNA with CellBender'
    }

    input {
        String prefix
        File h5

        # Docker image for cellbender remove-background version 
        String docker_image = "us.gcr.io/broad-dsde-methods/cellbender:0.2.0"

        # Method configuration inputs
        Int? expected_cells = 5000
        Int? total_droplets_included
        String? model
        Int? low_count_threshold
        String? fpr  # in quotes: floats separated by whitespace: the output false positive rate(s)
        Int? epochs
        Int? z_dim
        String? z_layers  # in quotes: integers separated by whitespace
        Float? empty_drop_training_fraction
        String? blacklist_genes  # in quotes: integers separated by whitespace
        Float? learning_rate
        Boolean? exclude_antibody_capture = false

        # Hardware-related inputs
        String? hardware_zones = "us-east1-d us-east1-c us-central1-a us-central1-c us-west1-b"
        Int? hardware_disk_size_GB = 50
        Int? hardware_boot_disk_size_GB = 20
        Int? hardware_preemptible_tries = 0
        Int? hardware_cpu_count = 4
        Int? hardware_memory_GB = 15
        String? hardware_gpu_type = "nvidia-tesla-t4"
    }

    command {
        cellbender remove-background \
            --input "${h5}" \
            --output "${prefix}_cellbender.h5" \
            --cuda \
            ${"--expected-cells " + expected_cells} \
            ${"--total-droplets-included " + total_droplets_included} \
            ${"--fpr " + fpr} \
            ${"--model " + model} \
            ${"--low-count-threshold " + low_count_threshold} \
            ${"--epochs " + epochs} \
            ${"--z-dim " + z_dim} \
            ${"--z-layers " + z_layers} \
            ${"--empty-drop-training-fraction " + empty_drop_training_fraction} \
            ${"--blacklist-genes " + blacklist_genes} \
            ${"--learning-rate " + learning_rate} \
            ${true="--exclude-antibody-capture" false=" " exclude_antibody_capture}
    }

    output {
        File cellbender_log = "${prefix}_cellbender.log"
        File cellbender_pdf = "${prefix}_cellbender.pdf"
        File cellbender_csv = "${prefix}_cellbender_cell_barcodes.csv"
        File cellbender_h5 = "${prefix}_cellbender.h5"  # v2 creates a number of outputs depending on "fpr"
        File cellbender_filtered_h5 = "${prefix}_cellbender_filtered.h5"
    }

    runtime {
        docker: "${docker_image}"
        bootDiskSizeGb: hardware_boot_disk_size_GB
        disks: "local-disk ${hardware_disk_size_GB} HDD"
        memory: "${hardware_memory_GB}G"
        cpu: hardware_cpu_count
        zones: "${hardware_zones}"
        gpuCount: 1
        gpuType: "${hardware_gpu_type}"
        preemptible: hardware_preemptible_tries
        maxRetries: 0
    }

}

