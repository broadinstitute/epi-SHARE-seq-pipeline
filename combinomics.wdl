version 1.0

# Import the sub-workflow for preprocessing the fastqs.
import "tasks/10x_task_preprocess.wdl" as preprocess_tenx
import "tasks/10x_create_barcode_mapping.wdl" as tenx_barcode_map
import "workflows/subwf-atac.wdl" as atac
import "workflows/subwf-rna.wdl" as rna
import "workflows/subwf-find-dorcs.wdl" as find_dorcs
import "tasks/task_joint_qc.wdl" as joint_qc
import "tasks/task_html_report.wdl" as html_report


# WDL workflow for SHARE-seq

workflow combinomics {

    input {
        # Common inputs

        Boolean trim_fastqs = true
        Boolean dorcs_flag = true
        String chemistry
        String prefix = "shareseq-project"
        String? subpool
        String pipeline_modality = "full" # "full": run everything; "count_only": stops after producing fragment file and count matrix; "no_align": correct and trim raw fastqs.

        File whitelists_tsv = 'gs://broad-buenrostro-pipeline-genome-annotations/whitelists/whitelists.tsv'
        File? whitelist
        File? whitelist_atac
        File? whitelist_rna

        # ATAC-specific inputs
        Array[File] read1_atac
        Array[File] read2_atac
        Array[File] fastq_barcode = []
        Boolean count_only = false
        File? chrom_sizes
        File? atac_genome_index_tar
        File? tss_bed
        String? barcode_tag = "CB"

        # ATAC - Filter
        ## Biological
        Int? atac_filter_minimum_fragments_cutoff = 1

        # RNA-specific inputs
        Array[File] read1_rna
        Array[File] read2_rna

        File? gtf
        File? idx_tar_rna

        String? gene_naming = "gene_name"

        # DORCs specific inputs
        File? peak_set

        # Joint qc
        Int remove_low_yielding_cells = 10

        File genome_tsv
        String? genome_name
    }

    Map[String, File] annotations = read_map(genome_tsv)
    String genome_name_ =  select_first([genome_name, annotations["genome_name"]])
    File peak_set_ = select_first([peak_set, annotations["ccre"]])
    File idx_tar_atac_ = select_first([atac_genome_index_tar, annotations["bowtie2_idx_tar"]])
    File chrom_sizes_ = select_first([chrom_sizes, annotations["chrsz"]])
    File tss_bed_ = select_first([tss_bed, annotations["tss"]])

    File idx_tar_rna_ = select_first([idx_tar_rna, annotations["star_idx_tar"]])
    File gtf_ = select_first([gtf, annotations["genesgtf"]])

    Boolean process_atac = if length(read1_atac)>0 then true else false
    Boolean process_rna = if length(read1_rna)>0 then true else false

    Map[String, File] whitelists = read_map(whitelists_tsv)
    File? whitelist_ = if chemistry=='10x_multiome' then whitelist else select_first([whitelist, whitelists[chemistry]])
    File? whitelist_rna_ = if chemistry=="10x_multiome" then select_first([whitelist_rna, whitelists["${chemistry}_rna"]]) else whitelist_rna
    File? whitelist_atac_ = if chemistry=="10x_multiome" then select_first([whitelist_atac, whitelists["${chemistry}_atac"]]) else whitelist_atac

    if ( chemistry != "shareseq" && process_atac) {
        scatter (idx in range(length(read1_atac))) {
            call preprocess_tenx.preprocess_tenx as preprocess_tenx{
                    input:
                        fastq_R1 = read1_atac[idx],
                        fastq_R3 = read2_atac[idx],
                        fastq_R2 = fastq_barcode[idx],
                        whitelist = select_first([whitelist_atac, whitelist_atac_]),
                        chemistry = chemistry,
                        prefix = prefix
            }
        }
        if ( chemistry == "10x_multiome" ){
            call tenx_barcode_map.mapping_tenx_barcodes as barcode_mapping{
                input:
                    whitelist_atac = select_first([whitelist_atac, whitelist_atac_]),
                    whitelist_rna = select_first([whitelist_rna, whitelist_rna_, whitelist_]),
            }
        }
    }

    if ( process_rna ) {
        if ( read1_rna[0] != "" ) {
            call rna.wf_rna as rna{
                input:
                    chemistry = chemistry,
                    read1 = read1_rna,
                    read2 = read2_rna,
                    whitelist = select_first([whitelist_rna, whitelist_rna_, whitelist, whitelist_]),
                    idx_tar = idx_tar_rna_,
                    prefix = prefix,
                    subpool = subpool,
                    genome_name = genome_name_,
                    pipeline_modality = pipeline_modality,
                    gene_naming = gene_naming
            }
        }
    }

    if ( process_atac ) {
        if ( read1_atac[0] != "" ) {
            call atac.wf_atac as atac{
                input:
                    read1 = select_first([preprocess_tenx.fastq_R1_preprocessed ,read1_atac]),
                    read2 = select_first([preprocess_tenx.fastq_R2_preprocessed ,read2_atac]),
                    chemistry = chemistry,
                    subpool = subpool,
                    whitelist = select_first([whitelist_atac, whitelist_atac_, whitelist, whitelist_]),
                    trim_fastqs = trim_fastqs,
                    chrom_sizes = chrom_sizes_,
                    genome_index_tar = idx_tar_atac_,
                    tss_bed = tss_bed_,
                    peak_set = peak_set_,
                    prefix = prefix,
                    genome_name = genome_name_,
                    barcode_conversion_dict = barcode_mapping.tenx_barcode_conversion_dict,
                    pipeline_modality = pipeline_modality
            }
        }
    }

    if ( process_atac && process_rna ) {
        if ( read1_atac[0] != "" && read1_rna[0] != "" ) {
            if ( pipeline_modality == "full" ) {
                if ( dorcs_flag ){
                    call find_dorcs.wf_dorcs as dorcs{
                        input:
                            rna_matrix = rna.rna_h5,
                            atac_fragments = atac.atac_fragments,
                            peak_file = peak_set_,
                            genome = genome_name_,
                            prefix = prefix
                    }
                }
            }
            
            if ( pipeline_modality != "no_align" ) {
                call joint_qc.joint_qc_plotting as joint_qc {
                    input:
                        atac_barcode_metadata = atac.atac_barcode_metadata,
                        rna_barcode_metadata = rna.rna_barcode_metadata,
                        prefix = prefix,
                        genome_name = genome_name_
                }
            }
        }
    }

    if ( pipeline_modality != "no_align" ) {
        call html_report.html_report as html_report {
            input:
                prefix = prefix,
                atac_metrics = atac.atac_qc_metrics_csv,

                rna_total_reads = rna.rna_total_reads,
                rna_aligned_uniquely = rna.rna_aligned_uniquely,
                rna_aligned_multimap = rna.rna_aligned_multimap,
                rna_unaligned = rna.rna_unaligned,
                rna_feature_reads = rna.rna_feature_reads,
                rna_duplicate_reads = rna.rna_duplicate_reads,
                rna_frig = rna.rna_frig,
                ## JPEG files to be encoded and appended to html
                # RNA plots
                image_files = [joint_qc.joint_qc_plot, joint_qc.joint_density_plot,
                               rna.rna_umi_barcode_rank_plot, rna.rna_gene_barcode_rank_plot, rna.rna_gene_umi_scatter_plot, rna.rna_seurat_raw_violin_plot, rna.rna_seurat_raw_qc_scatter_plot, rna.rna_seurat_filtered_violin_plot, rna.rna_seurat_filtered_qc_scatter_plot, rna.rna_seurat_variable_genes_plot, rna.rna_seurat_PCA_dim_loadings_plot, rna.rna_seurat_PCA_plot, rna.rna_seurat_heatmap_plot, rna.rna_seurat_jackstraw_plot, rna.rna_seurat_elbow_plot, rna.rna_seurat_umap_cluster_plot, rna.rna_seurat_umap_rna_count_plot, rna.rna_seurat_umap_gene_count_plot, rna.rna_seurat_umap_mito_plot,
                               atac.atac_qc_barcode_rank_plot, atac.atac_qc_hist_plot, atac.atac_qc_tss_enrichment, atac.atac_archr_raw_tss_enrichment, atac.atac_archr_filtered_tss_enrichment, atac.atac_archr_raw_fragment_size_plot, atac.atac_archr_filtered_fragment_size_plot, atac.atac_archr_umap_doublets, atac.atac_archr_umap_cluster_plot, atac.atac_archr_umap_doublets, atac.atac_archr_umap_num_frags_plot, atac.atac_archr_umap_tss_score_plot, atac.atac_archr_umap_frip_plot, atac.atac_archr_gene_heatmap_plot,
                               dorcs.j_plot],
                ## Links to files and logs to append to end of html
                log_files = [rna.rna_alignment_log,  rna.task_starsolo_barcodes_stats, rna.task_starsolo_features_stats, rna.task_starsolo_summary_csv, rna.task_starsolo_umi_per_cell, rna.task_starsolo_raw_tar,rna.rna_seurat_notebook_log, atac.atac_alignment_log, atac.atac_archr_notebook_log, dorcs.dorcs_notebook_log]
        }
    }

    output{
        # Fastq after correction/trimming
        Array[File]? atac_read1_processed = atac.atac_read1_processed
        Array[File]? atac_read2_processed = atac.atac_read2_processed

        Array[File]? rna_read1_processed = rna.rna_read1_processed
        Array[File]? rna_read2_processed = rna.rna_read2_processed

        # RNA outputs
        File? rna_final_bam = rna.task_starsolo_output_bam
        File? rna_starsolo_raw_tar = rna.task_starsolo_raw_tar
        File? rna_h5 = rna.rna_h5
        File? rna_barcode_metadata  = rna.rna_barcode_metadata
        File? rna_seurat_notebook_output = rna.rna_seurat_notebook_output
        File? rna_seurat_obj = rna.rna_seurat_obj

        # ATAC ouputs
        File? atac_fragments = atac.atac_fragments
        File? atac_fragments_index = atac.atac_fragments_index
        File? atac_barcode_metadata = atac.atac_barcode_metadata
        File? atac_archr_notebook_output = atac.atac_archr_notebook_output
        File? atac_archr_arrow = atac.atac_archr_arrow
        File? atac_track_bigwig = atac.atac_track_bigwig
        File? atac_track_bigwig_no_nucleosome = atac.atac_track_bigwig_no_nucleosome
        File? atac_track_bigwig_mono_nucleosome = atac.atac_track_bigwig_mono_nucleosome
        File? atac_track_bigwig_multi_nucleosome = atac.atac_track_bigwig_multi_nucleosome


        # DORCS output
        File? dorcs_notebook_output = dorcs.dorcs_notebook_output
        File? dorcs_genes_summary = dorcs.dorcs_genes_summary
        File? dorcs_regions_summary = dorcs.dorcs_regions_summary

        # Joint outputs
        File? joint_barcode_metadata = joint_qc.joint_barcode_metadata

        # Report
        File? html_summary = html_report.html_report_file
        File? csv_summary_file = html_report.csv_summary_file
    }

}

