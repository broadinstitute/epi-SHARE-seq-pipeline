version 1.0

# Import the sub-workflow for preprocessing the fastqs.
import "tasks/10x_task_preprocess.wdl" as preprocess_tenx
import "tasks/10x_create_barcode_mapping.wdl" as tenx_barcode_map
import "workflows/subwf-atac.wdl" as share_atac
import "workflows/subwf-rna-starsolo.wdl" as share_rna
import "workflows/subwf-find-dorcs.wdl" as find_dorcs
import "tasks/share_task_joint_qc.wdl" as joint_qc
import "tasks/share_task_html_report.wdl" as html_report
import "tasks/raise_exception.wdl" as raise_exception


# WDL workflow for SHARE-seq

workflow share {

    input {
        # Common inputs

        Boolean trim_fastqs = true
        Boolean dorcs_flag = true
        String chemistry
        String prefix = "shareseq-project"
        String? pkr=""
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

        #Int? cpus_atac
        #Int? cutoff_atac = 100
        #Int? atac_mapq_threshold = 30


        # ATAC - Align
        #Int? atac_align_multimappers

        # ATAC - Filter
        ## Biological
        Int? atac_filter_minimum_fragments_cutoff = 1
        #Int? atac_filter_shift_plus = 4
        #Int? atac_filter_shift_minus = -4

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

    if (chemistry != "shareseq" && chemistry != "10x_multiome" && chemistry != "10x_v2") {
      call raise_exception.raise_exception {
        input:
          msg = "Chemistry '~{chemistry}' not supported. Please use 'shareseq', '10x_multiome' or '10x_v2'"
      }
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
            call share_rna.wf_rna as rna{
                input:
                    chemistry = chemistry,
                    read1 = read1_rna,
                    read2 = read2_rna,
                    whitelist = select_first([whitelist_rna, whitelist_rna_, whitelist, whitelist_]),
                    idx_tar = idx_tar_rna_,
                    prefix = prefix,
                    pkr = pkr,
                    genome_name = genome_name_,
                    pipeline_modality = pipeline_modality,
                    gene_naming = gene_naming
            }
        }
    }

    if ( process_atac ) {
        if ( read1_atac[0] != "" ) {
            call share_atac.wf_atac as atac{
                input:
                    read1 = select_first([preprocess_tenx.fastq_R1_preprocessed ,read1_atac]),
                    read2 = select_first([preprocess_tenx.fastq_R2_preprocessed ,read2_atac]),
                    fastq_barcode = fastq_barcode,
                    chemistry = chemistry,
                    pkr = pkr,
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
                            rna_matrix = rna.share_rna_h5,
                            atac_fragments = atac.share_atac_filter_fragments,
                            peak_file = peak_set_,
                            genome = genome_name_,
                            prefix = prefix
                    }
                }
            }
            
            if ( pipeline_modality != "no_align" ) {
                call joint_qc.joint_qc_plotting as joint_qc {
                    input:
                        atac_barcode_metadata = atac.share_atac_barcode_metadata,
                        rna_barcode_metadata = rna.share_rna_barcode_metadata,
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
                atac_total_reads = atac.share_atac_total_reads,
                atac_aligned_uniquely = atac.share_atac_aligned_uniquely,
                atac_unaligned = atac.share_atac_unaligned,
                atac_feature_reads = atac.share_atac_feature_reads,
                atac_duplicate_reads = atac.share_atac_duplicate_reads,
                atac_percent_duplicates = atac.share_atac_percent_duplicates,
                rna_total_reads = rna.share_rna_total_reads,
                rna_aligned_uniquely = rna.share_rna_aligned_uniquely,
                rna_aligned_multimap = rna.share_rna_aligned_multimap,
                rna_unaligned = rna.share_rna_unaligned,
                rna_feature_reads = rna.share_rna_feature_reads,
                rna_duplicate_reads = rna.share_rna_duplicate_reads,

                ## JPEG files to be encoded and appended to html
                # RNA plots
                image_files = [joint_qc.joint_qc_plot, joint_qc.joint_density_plot, rna.share_rna_umi_barcode_rank_plot, rna.share_rna_gene_barcode_rank_plot, rna.share_rna_gene_umi_scatter_plot, rna.share_rna_seurat_raw_violin_plot, rna.share_rna_seurat_raw_qc_scatter_plot, rna.share_rna_seurat_filtered_violin_plot, rna.share_rna_seurat_filtered_qc_scatter_plot, rna.share_rna_seurat_variable_genes_plot, rna.share_rna_seurat_PCA_dim_loadings_plot, rna.share_rna_seurat_PCA_plot, rna.share_rna_seurat_heatmap_plot, rna.share_rna_seurat_jackstraw_plot, rna.share_rna_seurat_elbow_plot, rna.share_rna_seurat_umap_cluster_plot, rna.share_rna_seurat_umap_rna_count_plot, rna.share_rna_seurat_umap_gene_count_plot, rna.share_rna_seurat_umap_mito_plot, atac.share_atac_qc_barcode_rank_plot, atac.share_atac_qc_hist_plot, atac.share_atac_qc_tss_enrichment,  atac.share_atac_archr_raw_tss_enrichment, atac.share_atac_archr_filtered_tss_enrichment, atac.share_atac_archr_raw_fragment_size_plot, atac.share_atac_archr_filtered_fragment_size_plot, atac.share_atac_archr_umap_doublets, atac.share_atac_archr_umap_cluster_plot, atac.share_atac_archr_umap_doublets, atac.share_atac_archr_umap_num_frags_plot, atac.share_atac_archr_umap_tss_score_plot, atac.share_atac_archr_umap_frip_plot,atac.share_atac_archr_gene_heatmap_plot, dorcs.j_plot],

                ## Links to files and logs to append to end of html
                log_files = [rna.share_rna_alignment_log,  rna.share_task_starsolo_barcodes_stats, rna.share_task_starsolo_features_stats, rna.share_task_starsolo_summary_csv, rna.share_task_starsolo_umi_per_cell, rna.share_task_starsolo_raw_tar,rna.share_rna_seurat_notebook_log, atac.share_atac_alignment_log, atac.share_atac_archr_notebook_log, dorcs.dorcs_notebook_log]
        }
    }

    output{
        # Fastq after correction/trimming
        Array[File]? atac_read1_processed = atac.atac_read1_processed
        Array[File]? atac_read2_processed = atac.atac_read2_processed

        Array[File]? rna_read1_processed = rna.rna_read1_processed
        Array[File]? rna_read2_processed = rna.rna_read2_processed

        # RNA outputs
        File? share_rna_final_bam = rna.share_task_starsolo_output_bam
        File? share_rna_starsolo_raw_tar = rna.share_task_starsolo_raw_tar
        File? share_rna_h5 = rna.share_rna_h5
        File? share_rna_barcode_metadata  = rna.share_rna_barcode_metadata
        File? share_rna_seurat_notebook_output = rna.share_rna_seurat_notebook_output
        File? share_rna_seurat_obj = rna.share_rna_seurat_obj

        # ATAC ouputs
        File? share_atac_final_bam_dedup = atac.share_atac_filter_alignment_dedup
        File? share_atac_filter_fragments = atac.share_atac_filter_fragments
        File? share_atac_filter_fragments_index = atac.share_atac_filter_fragments_index
        File? share_atac_barcode_metadata = atac.share_atac_barcode_metadata
        File? share_atac_archr_notebook_output = atac.share_atac_archr_notebook_output
        File? share_atac_archr_arrow = atac.share_atac_archr_arrow

        # DORCS output
        File? dorcs_notebook_output = dorcs.dorcs_notebook_output
        File? dorcs_genes_summary = dorcs.dorcs_genes_summary
        File? dorcs_regions_summary = dorcs.dorcs_regions_summary

        # Joint outputs
        File? joint_barcode_metadata = joint_qc.joint_barcode_metadata

        # Report
        File? html_summary = html_report.html_report_file
    }

}

