version 1.0

import '../tasks/share_task_merge_h5.wdl' as share_task_merge_h5
import '../tasks/share_task_merge_fragments.wdl' as share_task_merge_fragments
import '../tasks/share_task_merge_barcode_metadata.wdl' as share_task_merge_barcode_metadata
import '../tasks/share_task_seurat.wdl' as share_task_seurat 
import '../tasks/share_task_archr.wdl' as share_task_archr
import '../tasks/share_task_joint_qc.wdl' as share_task_joint_qc
import './subwf-find-dorcs.wdl' as find_dorcs

workflow aggregate_counts {
    meta {
        version: 'v0.1'
        author: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-seq pipeline: sub-workflow to aggregate RNA and ATAC counts, and perform downstream processing.'
    }

    input {
        # Common inputs
        File genome_tsv
        String? genome_name
        String? prefix = 'merged'
        Boolean merge_only = false
        Boolean? concat_barcodes = false       
        
        # RNA aggregation inputs
        Array[File] tars
        Array[String]? pkrs
        String? gene_naming = 'gene_name'
        Float? merge_h5_disk_factor
        Float? merge_h5_memory_factor
        String? merge_h5_docker_image

        # ATAC aggregation inputs
        Array[File] fragments
        Float? merge_fragments_disk_factor
        String? merge_fragments_docker_image
        Int? merge_fragments_cpus

        # Merge barcode metadata inputs
        Array[File] rna_barcode_metadata
        Array[File] atac_barcode_metadata
        Float? merge_barcode_metadata_disk_factor
        Float? merge_barcode_metadata_memory_factor
        String? merge_barcode_metadata_docker_image

        # Seurat inputs 
        Int? seurat_min_features
        Float? seurat_percent_mt
        Int? seurat_min_cells
        Int? seurat_umap_dim
        Float? seurat_umap_resolution
        Float? seurat_disk_factor
        Float? seurat_memory_factor
        String? seurat_docker_image

        # ArchR inputs
        File? peak_set

        # Joint QC inputs
        Int? remove_low_yielding_cells = 10
        String? joint_qc_docker_image
    }

    Boolean aggregate_rna = if length(tars) > 0 then true else false
    Boolean aggregate_atac = if length(fragments) > 0 then true else false

    Map[String, File] annotations = read_map(genome_tsv)
    String genome_name_ =  select_first([genome_name, annotations["genome_name"]])
    File peak_set_ = select_first([peak_set, annotations["ccre"]])

    if (aggregate_rna) {
        call share_task_merge_h5.share_merge_h5 as merge_h5 {
            input:
                tars = tars,
                pkrs = pkrs,
                prefix = prefix,
                concat_barcodes = concat_barcodes,
                gene_naming = gene_naming,
                disk_factor = merge_h5_disk_factor,
                memory_factor = merge_h5_memory_factor,
                docker_image = merge_h5_docker_image
        }

        call share_task_merge_barcode_metadata.share_merge_barcode_metadata as merge_rna_barcode_metadata {
            input:
                barcode_metadata = rna_barcode_metadata,
                modality = 'RNA',
                prefix = prefix,
                concat_barcodes = concat_barcodes,
                disk_factor = merge_barcode_metadata_disk_factor,
                memory_factor = merge_barcode_metadata_memory_factor,
                docker_image = merge_barcode_metadata_docker_image
        }

        if (!merge_only) {
            call share_task_seurat.seurat as seurat {
                input:
                    rna_matrix = merge_h5.h5_matrix,
                    genome_name = genome_name_,
                    min_features = seurat_min_features,
                    percent_mt = seurat_percent_mt,
                    min_cells = seurat_min_cells,
                    umap_dim = seurat_umap_dim,
                    umap_resolution = seurat_umap_resolution,
                    prefix = prefix,
                    disk_factor = seurat_disk_factor,
                    memory_factor = seurat_memory_factor,
                    docker_image = seurat_docker_image
            }
        }
    }

    if (aggregate_atac) {
        call share_task_merge_fragments.share_merge_fragments as merge_fragments {
            input:
                fragments = fragments,
                prefix = prefix,
                disk_factor = merge_fragments_disk_factor,
                docker_image = merge_fragments_docker_image,
                cpus = merge_fragments_cpus
        }

        call share_task_merge_barcode_metadata.share_merge_barcode_metadata as merge_atac_barcode_metadata {
            input:
                barcode_metadata = atac_barcode_metadata,
                modality = 'ATAC',
                prefix = prefix,
                concat_barcodes = concat_barcodes,
                disk_factor = merge_barcode_metadata_disk_factor,
                memory_factor = merge_barcode_metadata_memory_factor,
                docker_image = merge_barcode_metadata_docker_image
        }

        if (!merge_only) {
            call share_task_archr.archr as archr {
                input:
                    atac_frag = merge_fragments.fragments,
                    genome = genome_name_,
                    peak_set = peak_set_,
                    prefix = prefix
            }
        }
    }

    if (aggregate_atac && aggregate_rna) {
        call share_task_joint_qc.joint_qc_plotting as joint_qc {
            input:
                atac_barcode_metadata = merge_atac_barcode_metadata.barcode_metadata,
                rna_barcode_metadata = merge_rna_barcode_metadata.barcode_metadata,
                genome_name = genome_name_,
                prefix = prefix,
                docker_image = joint_qc_docker_image
        }

        if (!merge_only) {
            call find_dorcs.wf_dorcs as dorcs {
                input:
                    rna_matrix = merge_h5.h5_matrix,
                    atac_fragments = merge_fragments.fragments,
                    peak_file = peak_set_,
                    genome = genome_name_,
                    prefix = prefix
            }
        }
    }

    output {
        File? aggregated_h5 = merge_h5.h5_matrix
        File? aggreated_rna_barcode_metadata = merge_rna_barcode_metadata.barcode_metadata
        
        File? aggregated_fragments = merge_fragments.fragments
        File? aggregated_atac_barcode_metadata = merge_atac_barcode_metadata.barcode_metadata

        File? share_rna_seurat_notebook_output = seurat.notebook_output
        File? share_rna_seurat_obj = seurat.seurat_filtered_obj
        File? share_rna_plots_zip = seurat.plots_zip

        File? share_atac_archr_notebook_output = archr.notebook_output
        File? share_atac_archr_arrow = archr.archr_arrow
        File? share_atac_archr_obj = archr.archr_raw_obj
        File? share_atac_archr_plots_zip = archr.plots_zip

        File? joint_barcode_metadata = joint_qc.joint_barcode_metadata

        File? dorcs_notebook_output = dorcs.dorcs_notebook_output
        File? dorcs_genes_summary = dorcs.dorcs_genes_summary
        File? dorcs_regions_summary = dorcs.dorcs_regions_summary
    }
}   
