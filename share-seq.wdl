version 1.0

# Import the sub-workflow for preprocessing the fastqs.
import "workflows/subwf-atac-single-organism.wdl" as share_atac
import "workflows/subwf-rna-starsolo.wdl" as share_rna
import "workflows/subwf-find-dorcs.wdl" as find_dorcs
import "tasks/share_task_joint_qc.wdl" as joint_qc
import "tasks/share_task_html_report.wdl" as html_report

# WDL workflow for SHARE-seq

workflow ShareSeq {

    input {
        # Common inputs
        String chemistry
        String prefix = "shareseq-project"
        String genome_name_input
        Int? cpus = 16

        # ATAC-specific inputs
        Array[File] read1_atac
        Array[File] read2_atac
        Boolean trim_fastqs = true
        File? chrom_sizes
        File? idx_tar_atac
        File? tss_bed
        Int? cpus_atac

        Int? cutoff_atac = 100

        # RNA-specific inputs
        Boolean? count_only = false
        Boolean? multimappers = false
        Boolean? include_multimappers = false
        Boolean? include_introns = true
        Array[File] read1_rna
        Array[File] read2_rna
        File whitelists_tsv = 'gs://broad-buenrostro-pipeline-genome-annotations/whitelists/whitelists.tsv'
        File? whitelist
        File? genes_annotation_bed
        File? gtf
        File? idx_tar_rna
        Int? umi_cutoff = 100
        Int? gene_cutoff = 100
        Int? cpus_rna
        String? gene_naming = "gene_name"

        # Seurat
        Int? rna_seurat_min_features = 200
        Int? rna_seurat_max_features = 2500 #currently not used in Seurat, but used in DORCs
        Float? rna_seurat_percent_mt = 5
        Int? rna_seurat_min_cells = 3
        Int? rna_seurat_umap_dim
        Float? rna_seurat_umap_resolution
        Float? rna_seurat_disk_factor
        Float? rna_seurat_memory_factor

        # DORCs specific inputs
        File? peak_set
        Int? cpus_dorcs
        String save_plots_to_dir = "TRUE"
        String? dorcs_output_filename

        # DORCs filter
        Int dorcGeneCutOff = 10
        Float fripCutOff = 0.3
        Float corrPVal = 0.05
        Int topNGene = 20

        # Joint qc
        Int remove_low_yielding_cells = 10

        # Regulatory region around TSS. Default is +/- 50Kb
        Int windowPadSize = 50000
        #Int bootstraps = 100

        String docker_image_dorcs = "us.gcr.io/buenrostro-share-seq/dorcs_task_find_dorcs"
        Int? mem_gb_dorcs

        File human_genome_tsv = "gs://broad-buenrostro-pipeline-genome-annotations/IGVF_human/GRCh38_genome_files_hg38.tsv"
        File mouse_genome_tsv = "gs://broad-buenrostro-pipeline-genome-annotations/mm10/mm10_genome_files_STARsolo.tsv"
    }

    String genome_name = if genome_name_input == "GRCh38" then "hg38" else genome_name_input

    Map[String, File] annotations = if genome_name == "mm10" then read_map(mouse_genome_tsv) else read_map(human_genome_tsv)
    File peak_set_ = select_first([peak_set, annotations["ccre"]])
    File idx_tar_atac_ = select_first([idx_tar_atac, annotations["bowtie2_idx_tar"]])
    File chrom_sizes_ = select_first([chrom_sizes, annotations["chrsz"]])
    File tss_bed_ = select_first([tss_bed, annotations["tss"]])

    File idx_tar_rna_ = select_first([idx_tar_rna, annotations["star_idx_tar"]])
    File gtf_ = select_first([gtf, annotations["genesgtf"]])
    File genes_annotation_bed_ = select_first([genes_annotation_bed, annotations["genesbed"]])

    Map[String, File] whitelists = read_map(whitelists_tsv)
    File? whitelist_ = if chemistry=='shareseq' then whitelist else select_first([whitelist, whitelists[chemistry]])

    Boolean process_atac = if length(read1_atac)>0 then true else false
    Boolean process_rna = if length(read1_rna)>0 then true else false

    if ( process_rna ) {
        if ( read1_rna[0] != "" ) {
            call share_rna.wf_rna as rna{
                input:
                    chemistry = chemistry,
                    read1 = read1_rna,
                    read2 = read2_rna,
                    whitelist = whitelist_,
                    idx_tar = idx_tar_rna_,
                    umi_cutoff = umi_cutoff,
                    gene_cutoff = gene_cutoff,
                    prefix = prefix,
                    genome_name = genome_name,
                    cpus = cpus_rna,
                    count_only = count_only,
                    rna_seurat_min_features = rna_seurat_min_features,
                    rna_seurat_percent_mt = rna_seurat_percent_mt,
                    rna_seurat_min_cells = rna_seurat_min_cells,
                    rna_seurat_umap_dim = rna_seurat_umap_dim,
                    rna_seurat_umap_resolution = rna_seurat_umap_resolution,
                    rna_seurat_disk_factor = rna_seurat_disk_factor,
                    rna_seurat_memory_factor = rna_seurat_memory_factor
            }
        }
    }

    if ( process_atac ) {
        if ( read1_atac[0] != "" ) {
            call share_atac.wf_atac as atac{
                input:
                    read1 = read1_atac,
                    read2 = read2_atac,
                    trim_fastqs = trim_fastqs,
                    chrom_sizes = chrom_sizes_,
                    idx_tar = idx_tar_atac_,
                    tss_bed = tss_bed_,
                    peak_set = peak_set_,
                    prefix = prefix,
                    genome_name = genome_name,
                    cutoff = cutoff_atac,
                    cpus = cpus_atac
            }
        }
    }

    if ( process_atac && process_rna ) {
        if ( read1_atac[0] != "" && read1_rna[0] != "" ) {
            call find_dorcs.wf_dorcs as dorcs{
                input:
                    rna_matrix = rna.share_rna_h5,
                    atac_fragments = atac.share_atac_fragments_filtered,
                    peak_file = peak_set_,

                    genome = genome_name,
                    n_cores = cpus_dorcs,
                    save_plots_to_dir = save_plots_to_dir,
                    output_filename = dorcs_output_filename,
                    prefix = prefix,

                    minFeature_RNA = rna_seurat_min_features,
                    maxFeature_RNA = rna_seurat_max_features,
                    percentMT_RNA = rna_seurat_percent_mt,
                    minCells_RNA = rna_seurat_min_cells,

                    dorcGeneCutOff = dorcGeneCutOff,
                    fripCutOff = fripCutOff,
                    corrPVal = corrPVal,
                    topNGene = topNGene,

                    windowPadSize = windowPadSize,
                    mem_gb = mem_gb_dorcs
            }
        }
        call joint_qc.joint_qc_plotting as joint_qc {
            input:
                atac_barcode_metadata = atac.share_atac_archr_barcode_metadata,
                rna_barcode_metadata = rna.share_rna_barcode_metadata,
                remove_low_yielding_cells = remove_low_yielding_cells,
                prefix = prefix,
                genome_name = genome_name
        }
    }

    call html_report.html_report as html_report {
        input:
            prefix = prefix,
            atac_total_reads = atac.share_atac_total_reads,
            atac_aligned_uniquely = atac.share_atac_aligned_uniquely,
            atac_unaligned = atac.share_atac_unaligned,
            atac_feature_reads = atac.share_atac_feature_reads,
            atac_duplicate_reads = atac.share_atac_duplicate_reads,
            rna_total_reads = rna.share_rna_total_reads,
            rna_aligned_uniquely = rna.share_rna_aligned_uniquely,
            rna_aligned_multimap = rna.share_rna_aligned_multimap,
            rna_unaligned = rna.share_rna_unaligned,
            rna_feature_reads = rna.share_rna_feature_reads,
            rna_duplicate_reads = rna.share_rna_duplicate_reads,

            ## JPEG files to be encoded and appended to html
            # RNA plots
            image_files = [joint_qc.joint_qc_plot, joint_qc.joint_density_plot, rna.share_rna_umi_barcode_rank_plot, rna.share_rna_gene_barcode_rank_plot, rna.share_rna_gene_umi_scatter_plot, rna.share_rna_seurat_raw_violin_plot, rna.share_rna_seurat_raw_qc_scatter_plot, rna.share_rna_seurat_filtered_violin_plot, rna.share_rna_seurat_filtered_qc_scatter_plot, rna.share_rna_seurat_variable_genes_plot, rna.share_rna_seurat_PCA_dim_loadings_plot, rna.share_rna_seurat_PCA_plot, rna.share_rna_seurat_heatmap_plot, rna.share_rna_seurat_jackstraw_plot, rna.share_rna_seurat_elbow_plot, rna.share_rna_seurat_umap_cluster_plot, rna.share_rna_seurat_umap_rna_count_plot, rna.share_rna_seurat_umap_gene_count_plot, rna.share_rna_seurat_umap_mito_plot, atac.share_atac_qc_library_plot, atac.share_atac_qc_hist_plot, atac.share_atac_qc_tss_enrichment, atac.share_atac_archr_gene_heatmap_plot, atac.share_atac_archr_raw_tss_enrichment, atac.share_atac_archr_filtered_tss_enrichment, atac.share_atac_archr_raw_fragment_size_plot, atac.share_atac_archr_filtered_fragment_size_plot, atac.share_atac_archr_umap_doublets, atac.share_atac_archr_umap_cluster_plot, atac.share_atac_archr_umap_doublets, atac.share_atac_archr_umap_num_frags_plot, atac.share_atac_archr_umap_tss_score_plot, atac.share_atac_archr_umap_frip_plot,atac.share_atac_archr_gene_heatmap_plot, atac.share_atac_archr_strict_raw_tss_enrichment, atac.share_atac_archr_strict_filtered_tss_enrichment, atac.share_atac_archr_strict_raw_fragment_size_plot, atac.share_atac_archr_strict_filtered_fragment_size_plot, atac.share_atac_archr_strict_umap_doublets, atac.share_atac_archr_strict_umap_cluster_plot, atac.share_atac_archr_umap_doublets, atac.share_atac_archr_strict_umap_num_frags_plot, atac.share_atac_archr_strict_umap_tss_score_plot, atac.share_atac_archr_strict_umap_frip_plot,atac.share_atac_archr_strict_gene_heatmap_plot, dorcs.j_plot],

            ## Links to files and logs to append to end of html
            log_files = [rna.share_rna_alignment_log,  rna.share_task_starsolo_barcodes_stats, rna.share_task_starsolo_features_stats, rna.share_task_starsolo_summary_csv, rna.share_task_starsolo_umi_per_cell, rna.share_task_starsolo_raw_tar,rna.share_rna_seurat_notebook_log, atac.share_atac_alignment_log, atac.share_atac_archr_notebook_log, dorcs.dorcs_notebook_log]
    }

    output{
        File? share_rna_output_bam = rna.share_task_starsolo_output_bam
        File? share_rna_alignment_log = rna.share_task_starsolo_log_out
        File? share_rna_summary_stats = rna.share_task_starsolo_summary_csv
        File? share_rna_barcodes_stats = rna.share_task_starsolo_barcodes_stats
        File? share_rna_features_stats = rna.share_task_starsolo_features_stats
        File? share_rna_summary_csv = rna.share_task_starsolo_summary_csv
        File? share_rna_umi_per_cell = rna.share_task_starsolo_umi_per_cell
        File? share_rna_raw_tar = rna.share_task_starsolo_raw_tar

        File? share_rna_h5 = rna.share_rna_h5

        File? share_rna_barcode_metadata  = rna.share_rna_barcode_metadata
        File? share_rna_duplicates_log = rna.share_rna_duplicates_log
        File? share_rna_umi_barcode_rank_plot  = rna.share_rna_umi_barcode_rank_plot
        File? share_rna_gene_barcode_rank_plot = rna.share_rna_gene_barcode_rank_plot
        File? share_rna_gene_umi_scatter_plot = rna.share_rna_gene_umi_scatter_plot

        File? share_rna_seurat_notebook_output = rna.share_rna_seurat_notebook_output
        File? share_rna_seurat_notebook_log = rna.share_rna_seurat_notebook_log
        File? share_rna_seurat_raw_violin_plot = rna.share_rna_seurat_raw_violin_plot
        File? share_rna_seurat_filtered_violin_plot = rna.share_rna_seurat_filtered_violin_plot
        File? share_rna_seurat_raw_qc_scatter_plot = rna.share_rna_seurat_raw_qc_scatter_plot
        File? share_rna_seurat_filtered_qc_scatter_plot = rna.share_rna_seurat_filtered_qc_scatter_plot
        File? share_rna_seurat_variable_genes_plot = rna.share_rna_seurat_variable_genes_plot
        File? share_rna_seurat_PCA_dim_loadings_plot = rna.share_rna_seurat_PCA_dim_loadings_plot
        File? share_rna_seurat_PCA_plot = rna.share_rna_seurat_PCA_plot
        File? share_rna_seurat_heatmap_plot = rna.share_rna_seurat_heatmap_plot
        File? share_rna_seurat_jackstraw_plot = rna.share_rna_seurat_jackstraw_plot
        File? share_rna_seurat_elbow_plot = rna.share_rna_seurat_elbow_plot
        File? share_rna_seurat_umap_cluster_plot = rna.share_rna_seurat_umap_cluster_plot
        File? share_rna_seurat_umap_rna_count_plot = rna.share_rna_seurat_umap_rna_count_plot
        File? share_rna_seurat_umap_gene_count_plot = rna.share_rna_seurat_umap_gene_count_plot
        File? share_rna_seurat_umap_mito_plot = rna.share_rna_seurat_umap_mito_plot
        File? share_rna_seurat_obj = rna.share_rna_seurat_obj
        File? share_rna_plots_zip = rna.share_rna_plots_zip

        File? share_atac_alignment_raw_index = atac.share_atac_alignment_raw_index
        File? share_atac_alignment_log = atac.share_atac_alignment_log

        File? share_atac_alignment_filtered = atac.share_atac_alignment_filtered
        File? share_atac_alignment_filtered_index = atac.share_atac_alignment_filtered_index
        File? share_atac_fragments_raw = atac.share_atac_fragments_raw


        File? share_atac_barcodes = atac.share_atac_barcodes
        File? share_atac_fragments_filtered = atac.share_atac_fragments_filtered
        File? share_atac_counts_raw = atac.share_atac_counts_raw
        File? share_atac_counts_filtered = atac.share_atac_counts_filtered

        File? share_atac_qc_library_counts = atac.share_atac_qc_library_counts
        File? share_atac_qc_library_duplicates = atac.share_atac_qc_library_duplicates
        File? share_atac_qc_library_plot = atac.share_atac_qc_library_plot

        File? share_atac_qc_final = atac.share_atac_qc_final
        File? share_atac_qc_hist_plot = atac.share_atac_qc_hist_plot
        File? share_atac_qc_hist_txt = atac.share_atac_qc_hist_txt
        File? share_atac_qc_tss_enrichment = atac.share_atac_qc_tss_enrichment

        File? share_atac_archr_notebook_output = atac.share_atac_archr_notebook_output
        File? share_atac_archr_notebook_log = atac.share_atac_archr_notebook_log
        File? share_atac_archr_gene_heatmap_plot = atac.share_atac_archr_gene_heatmap_plot
        File? share_atac_archr_raw_tss_enrichment = atac.share_atac_archr_raw_tss_enrichment
        File? share_atac_archr_filtered_tss_enrichment = atac.share_atac_archr_filtered_tss_enrichment
        File? share_atac_archr_filtered_fragment_size_plot = atac.share_atac_archr_filtered_fragment_size_plot
        File? share_atac_archr_umap_doublets = atac.share_atac_archr_umap_doublets
        File? share_atac_archr_umap_cluster_plot = atac.share_atac_archr_umap_cluster_plot
        File? share_atac_archr_arrow = atac.share_atac_archr_arrow
        File? share_atac_archr_obj = atac.share_atac_archr_obj
        File? share_atac_archr_plots_zip = atac.share_atac_archr_plots_zip

        File? dorcs_notebook_output = dorcs.dorcs_notebook_output
        File? dorcs_notebook_log = dorcs.dorcs_notebook_log
        File? seurat_violin_plot = dorcs.seurat_violin_plot
        File? j_plot = dorcs.j_plot
        File? plots_zip = dorcs.plots_zip
        File? dorcs_genes_summary = dorcs.dorcs_genes_summary
        File? dorcs_regions_summary = dorcs.dorcs_regions_summary

        File? joint_qc_plot = joint_qc.joint_qc_plot
        File? joint_density_plot = joint_qc.joint_density_plot
        File? joint_barcode_metadata = joint_qc.joint_barcode_metadata

        File? html_summary = html_report.html_report_file
    }

}

