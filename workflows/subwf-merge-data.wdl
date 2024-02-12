version 1.0

import "../tasks/task_merge_rna_counts.wdl" as task_merge_rna_counts
import "../tasks/task_merge_atac_fragments.wdl" as task_merge_atac_fragments
import "../tasks/task_qc_merged_atac.wdl" as task_qc_merged_atac
import "../tasks/task_seurat.wdl" as task_seurat 
import "../tasks/task_archr.wdl" as task_archr
import "../tasks/task_joint_qc.wdl" as task_joint_qc
import "./subwf-find-dorcs.wdl" as find_dorcs
import "../tasks/task_html_report.wdl" as html_report

workflow merge {
    meta {
        version: 'v0.1'
        author: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-seq pipeline: sub-workflow to aggregate RNA and ATAC counts, and perform downstream processing.'
    }

    input {
        # Common inputs
        File genome_tsv
        Array[String] dataset_names
        String? genome_name
        String? prefix = 'merged'
        
        # RNA merge counts inputs
        Array[File] tars
        Array[String]? subpool_names
        String? gene_naming = 'gene_name'
        Float? merge_counts_disk_factor
        Float? merge_counts_memory_factor
        String? merge_counts_docker_image

        # ATAC merge fragments inputs
        Array[File] fragments
        Float? merge_fragments_disk_factor
        String? merge_fragments_docker_image
        Int? merge_fragments_cpus

        # ATAC QC inputs
        Boolean run_qc_merged_atac = true
        Array[File] atac_barcode_metadata
        File? tss_bed
        Int? fragment_min_cutoff
        Int? hist_max_fragment
        Float? qc_merged_atac_disk_factor
        Float? qc_merged_atac_memory_factor
        String? qc_merged_atac_docker_image

        # Seurat inputs
        Boolean run_seurat = true 
        Int? seurat_min_features
        Float? seurat_percent_mt
        Int? seurat_min_cells
        Int? seurat_umap_dim
        Float? seurat_umap_resolution
        Float? seurat_disk_factor
        Float? seurat_memory_factor
        String? seurat_docker_image

        # ArchR inputs
        Boolean run_archr = true
        File? peak_set

        # Joint QC inputs
        Boolean run_joint_qc = true
        Int? remove_low_yielding_cells = 10
        String? joint_qc_docker_image

        # DORCs inputs
        Boolean run_dorcs = true

        # HTML report inputs
        Boolean run_html_report = true
    }

    Boolean aggregate_rna = if length(tars) > 0 then true else false
    Boolean aggregate_atac = if length(fragments) > 0 then true else false

    Map[String, File] annotations = read_map(genome_tsv)
    String genome_name_ =  select_first([genome_name, annotations["genome_name"]])
    File tss_bed_ = select_first([tss_bed, annotations["tss"]])
    File peak_set_ = select_first([peak_set, annotations["ccre"]])

    if (aggregate_rna) {
        call task_merge_rna_counts.merge_counts as merge_counts {
            input:
                tars = tars,
                subpool_names = subpool_names,
                prefix = prefix,
                dataset_names = dataset_names,
                gene_naming = gene_naming,
                disk_factor = merge_counts_disk_factor,
                memory_factor = merge_counts_memory_factor,
                docker_image = merge_counts_docker_image
        }

        if (run_seurat) {
            call task_seurat.seurat as seurat {
                input:
                    rna_matrix = merge_counts.h5_matrix,
                    dataset_barcodes = merge_counts.rna_dataset_barcodes,
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
        call task_merge_atac_fragments.merge_fragments as merge_fragments {
            input:
                fragments = fragments,
                prefix = prefix,
                disk_factor = merge_fragments_disk_factor,
                docker_image = merge_fragments_docker_image,
                cpus = merge_fragments_cpus
        }

        if (run_qc_merged_atac) {
            call task_qc_merged_atac.qc_merged_atac as qc_merged_atac {
                input:
                    barcode_metadata = atac_barcode_metadata,
                    fragments = merge_fragments.fragments,
                    fragments_index = merge_fragments.fragments_index,
                    tss = tss_bed_,
                    dataset_names = dataset_names,
                    prefix = prefix,
                    genome_name = genome_name_,
                    fragment_min_cutoff = fragment_min_cutoff,
                    hist_max_fragment = hist_max_fragment,
                    disk_factor = qc_merged_atac_disk_factor,
                    memory_factor = qc_merged_atac_memory_factor,
                    docker_image = qc_merged_atac_docker_image
            }
        }

        if (run_archr) {
            call task_archr.archr as archr {
                input:
                    atac_frag = merge_fragments.fragments,
                    genome = genome_name_,
                    peak_set = peak_set_,
                    prefix = prefix
            }
        }
    }

    if (aggregate_atac && aggregate_rna) {
        if (run_joint_qc && run_qc_merged_atac) {
            call task_joint_qc.joint_qc_plotting as joint_qc {
                input:
                    atac_barcode_metadata = qc_merged_atac.atac_barcode_metadata,
                    rna_barcode_metadata = merge_counts.rna_barcode_metadata,
                    genome_name = genome_name_,
                    prefix = prefix,
                    docker_image = joint_qc_docker_image
            }
        }

        if (run_dorcs) {
            call find_dorcs.wf_dorcs as dorcs {
                input:
                    rna_matrix = merge_counts.h5_matrix,
                    atac_fragments = merge_fragments.fragments,
                    peak_file = peak_set_,
                    genome = genome_name_,
                    prefix = prefix
            }
        }
    }

    if (run_html_report) {
        call html_report.html_report as html_report {
            input:
                prefix = prefix,
                image_files = [joint_qc.joint_qc_plot, joint_qc.joint_density_plot, seurat.seurat_raw_violin_plot, seurat.seurat_raw_qc_scatter_plot, seurat.seurat_filtered_violin_plot, seurat.seurat_filtered_qc_scatter_plot, seurat.seurat_variable_genes_plot, seurat.seurat_PCA_dim_loadings_plot, seurat.seurat_PCA_plot, seurat.seurat_heatmap_plot, seurat.seurat_jackstraw_plot, seurat.seurat_elbow_plot, seurat.seurat_umap_cluster_plot, seurat.seurat_umap_dataset_plot, seurat.seurat_umap_rna_count_plot, seurat.seurat_umap_gene_count_plot, seurat.seurat_umap_mito_plot, qc_merged_atac.fragment_barcode_rank_plot, qc_merged_atac.fragment_histogram, qc_merged_atac.insert_size_hist, archr.archr_raw_tss_by_uniq_frags_plot, archr.archr_filtered_tss_by_uniq_frags_plot, archr.archr_raw_frag_size_dist_plot, archr.archr_filtered_frag_size_dist_plot, archr.archr_umap_doublets, archr.archr_umap_cluster_plot, archr.archr_umap_doublets, archr.archr_umap_num_frags_plot, archr.archr_umap_tss_score_plot, archr.archr_umap_frip_plot, archr.archr_heatmap_plot, dorcs.j_plot],
                log_files = []
        }
    }

    output {
        File? merged_h5 = merge_counts.h5_matrix
        File? merged_rna_barcode_metadata = merge_counts.rna_barcode_metadata
 
        File? merged_fragments = merge_fragments.fragments

        File? merged_atac_barcode_metadata = qc_merged_atac.atac_barcode_metadata

        File? seurat_notebook_output = seurat.notebook_output
        File? seurat_obj = seurat.seurat_filtered_obj
        File? seurat_plots_zip = seurat.plots_zip

        File? archr_notebook_output = archr.notebook_output
        File? archr_arrow = archr.archr_arrow
        File? atac_archr_obj = archr.archr_raw_obj
        File? atac_archr_plots_zip = archr.plots_zip

        File? joint_barcode_metadata = joint_qc.joint_barcode_metadata

        File? dorcs_notebook_output = dorcs.dorcs_notebook_output
        File? dorcs_genes_summary = dorcs.dorcs_genes_summary
        File? dorcs_regions_summary = dorcs.dorcs_regions_summary

        File? html_report = html_report.html_report_file
    }
}   
