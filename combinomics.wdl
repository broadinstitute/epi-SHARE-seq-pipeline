version 1.0

# Import the sub-workflow for preprocessing the fastqs.
import "tasks/10x_create_barcode_mapping.wdl" as tenx_barcode_map
import "workflows/subwf-atac.wdl" as atac
import "workflows/subwf-rna.wdl" as rna
import "workflows/subwf-find-dorcs.wdl" as find_dorcs
import "tasks/task_joint_qc.wdl" as joint_qc
import "tasks/task_html_report.wdl" as html_report
import "structs/atac_output_struct.wdl"

# WDL workflow for SHARE-seq

workflow combinomics {

    input {
        # Common inputs

        String chemistry
        String prefix = "combinomics"
        String? subpool
        String pipeline_modality = "full" # "full": run everything; "count_only": stops after producing fragment file and count matrix

        File whitelists_tsv = 'gs://broad-buenrostro-pipeline-genome-annotations/whitelists/whitelists.tsv'
        File? whitelist
        File? whitelist_atac
        File? whitelist_rna

        # ATAC-specific inputs
        Array[File] read1_atac
        Array[File] read2_atac
        Array[File] fastq_barcode
        File? chrom_sizes
        File? tss_bed

        # ATAC - Filter
        ## Biological

        # RNA-specific inputs
        Array[File] read1_rna
        Array[File] read2_rna

        File? gtf
        File? idx_tar_rna
        File? idx_tar_atac

        String? gene_naming = "gene_name"

        # Peaks for QC?
        # File? peak_set

        File genome_tsv
        String? genome_name
    }

    Map[String, File] annotations = read_map(genome_tsv)
    String genome_name_ =  select_first([genome_name, annotations["genome_name"]])
    # File peak_set_ = select_first([peak_set, annotations["ccre"]])
    #File idx_tar_atac_ = select_first([atac_genome_index_tar, annotations["bowtie2_idx_tar"]])
    File chrom_sizes_ = select_first([chrom_sizes, annotations["chrsz"]])
    File tss_bed_ = select_first([tss_bed, annotations["tss"]])

    File idx_tar_rna_ = select_first([idx_tar_rna, annotations["star_idx_tar"]])
    File idx_tar_atac_ = select_first([idx_tar_atac, annotations["chromap_idx_tar"]])
    File gtf_ = select_first([gtf, annotations["genesgtf"]])

    Boolean process_atac = if length(read1_atac)>0 then true else false
    Boolean process_rna = if length(read1_rna)>0 then true else false

    Map[String, File] whitelists = read_map(whitelists_tsv)
    File? whitelist_ = if (chemistry=="10x_multiome" || chemistry=="shareseq") then whitelist else select_first([whitelist, whitelists[chemistry]])
    File? whitelist_rna_ = if (chemistry=="10x_multiome" || chemistry=="shareseq") then select_first([whitelist_rna, whitelists["${chemistry}_rna"]]) else whitelist_rna
    File? whitelist_atac_ = if (chemistry=="10x_multiome" || chemistry=="shareseq") then select_first([whitelist_atac, whitelists["${chemistry}_atac"]]) else whitelist_atac

    if ( chemistry != "shareseq" && process_atac) {
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
                    read1 = read1_atac,
                    read2 = read2_atac,
                    fastq_barcode = fastq_barcode,
                    chemistry = chemistry,
                    subpool = subpool,
                    gtf = gtf_,
                    whitelist = select_first([whitelist_atac, whitelist_atac_, whitelist, whitelist_]),
                    chrom_sizes = chrom_sizes_,
                    reference_index_tar_gz = idx_tar_atac_,
                    tss_bed = tss_bed_,
                    prefix = prefix,
                    genome_name = genome_name_,
                    barcode_conversion_dict = barcode_mapping.tenx_barcode_conversion_dict
            }
        }
    }

    if ( process_atac && process_rna ) {
        if ( read1_atac[0] != "" && read1_rna[0] != "" ) {
            call joint_qc.joint_qc_plotting as joint_qc {
                input:
                    atac_barcode_metadata = atac.atac_qc_barcode_metrics,
                    rna_barcode_metadata = rna.rna_barcode_metadata,
                    prefix = prefix,
                    genome_name = genome_name_
            }
        }
    }

    call html_report.html_report as html_report {
        input:
            prefix = prefix,
            atac_metrics = atac.atac_qc_barcode_metrics,
            rna_metrics = rna.rna_qc_metrics,
            ## JPEG files to be encoded and appended to html
            # RNA plots
            image_files = [joint_qc.joint_qc_plot, joint_qc.joint_density_plot,
                            rna.rna_umi_barcode_rank_plot, rna.rna_gene_barcode_rank_plot, rna.rna_gene_umi_scatter_plot, rna.rna_umi_histogram, rna.rna_seurat_raw_violin_plot, rna.rna_seurat_raw_qc_scatter_plot, rna.rna_seurat_filtered_violin_plot, rna.rna_seurat_filtered_qc_scatter_plot, rna.rna_seurat_variable_genes_plot, rna.rna_seurat_PCA_dim_loadings_plot, rna.rna_seurat_PCA_plot, rna.rna_seurat_heatmap_plot, rna.rna_seurat_jackstraw_plot, rna.rna_seurat_elbow_plot, rna.rna_seurat_umap_cluster_plot, rna.rna_seurat_umap_rna_count_plot, rna.rna_seurat_umap_gene_count_plot, rna.rna_seurat_umap_mito_plot,
                            atac.atac_qc_tss_enrichment_library_plot, atac.atac_qc_fragment_size_distribution_plot, atac.atac_qc_atac_knee_plot, atac.atac_qc_fraction_of_duplicates_distribution_plot, atac.atac_qc_n_fragment_vs_tss_enrichment_plot, atac.atac_qc_atac_n_fragment_vs_tss_enrichment_filtered_plot, atac.atac_qc_umap_leiden_plot
                        ],
            ## Links to files and logs to append to end of html
            #rna.task_starsolo_umi_per_cell,
            log_files = [rna.rna_alignment_log,  rna.task_starsolo_barcodes_stats, rna.task_starsolo_features_stats, rna.task_starsolo_summary_csv,  rna.task_starsolo_mtx_unique_tar, rna.rna_seurat_notebook_log, atac.atac_align_log]
    }

    output{
        # RNA outputs
        File? rna_final_bam = rna.task_starsolo_output_bam
        File? rna_bam_index = rna.task_starsolo_output_bam_index
        File? rna_starsolo_raw_tar = rna.task_starsolo_mtx_unique_tar
        File? rna_align_output_folder_tar = rna.task_starsolo_output_folder_tar

        File? rna_h5 = rna.rna_h5
        File? rna_barcode_metadata  = rna.rna_barcode_metadata
        File? rna_seurat_notebook_output = rna.rna_seurat_notebook_output
        File? rna_seurat_obj = rna.rna_seurat_obj

        # ATAC ouputs
        File? atac_final_bam = atac.atac_bam
        File? atac_bam_index = atac.atac_bam_index
        File? atac_bam_log = atac.atac_bam_alignment_stats
        File? atac_fragments = atac.atac_fragment_file
        File? atac_fragments_index = atac.atac_fragment_file_index
        File? atac_barcode_alignment_stats = atac.atac_align_barcode_statistics
        File? atac_barcode_metrics = atac.atac_qc_barcode_metrics
        File? atac_h5ad = atac.atac_qc_snapatac2_h5ad

        # Joint outputs
        File? joint_barcode_metadata = joint_qc.joint_barcode_metadata

        # Report
        File? html_summary = html_report.html_report_file
        File? csv_summary_file = html_report.csv_summary_file

        Atac_outputs? atac_struct_output = atac.atac_struct_output
    }

}

