version 1.0

import "../tasks/share_task_merge_h5.wdl" as share_task_merge_h5

workflow aggregate_counts {
    meta {
        version: 'v0.1'
        author: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-seq pipeline: sub-workflow to aggregate RNA and ATAC counts.'
    }

    input {
        # Common inputs
        String? prefix = 'merged'
        Array[String]? pkrs

        # RNA aggregation inputs
        Array[File] tars
        String? gene_naming = 'gene_name'
        Float? merge_h5_disk_factor
        Float? merge_h5_memory_factor
        String? merge_h5_docker_image

        # ATAC aggregation inputs
        Array[File] fragments

    }

    Boolean aggregate_rna = if length(tars) > 0 then true else false
    Boolean aggregate_atac = if length(fragments) > 0 then true else false

    if (aggregate_rna) {
        call share_task_merge_h5.share_merge_h5 as merge_h5 {
            input:
                tars = tars,
                pkrs = pkrs,
                prefix = prefix,
                gene_naming = gene_naming,
                disk_factor = merge_h5_disk_factor,
                memory_factor = merge_h5_memory_factor,
                docker_image = merge_h5_docker_image
        }
    }

    output {
        File? aggregated_h5 = merge_h5.h5_matrix
    }
}   
