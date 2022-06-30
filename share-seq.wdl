version 1.0

# Import the sub-workflow for preprocessing the fastqs.
import "workflows/subwf-atac-single-organism.wdl" as share_atac
import "workflows/subwf-rna-single-organism.wdl" as share_rna
import "workflows/subwf-find-dorcs.wdl" as find_dorcs

# WDL workflow for SHARE-seq

workflow ShareSeq {

    input {
        # Common inputs
        String prefix = "shareseq-project"
        String genome_name
        Int? cpus = 16


        # ATAC specific inputs
        File chrom_sizes
        Array[File] read1_atac
        Array[File] read2_atac
        File idx_tar_atac
        File tss_bed
        Int? cpus_atac
        Int cutoff_atac

        # RNA specific inputs
        Boolean multimappers = false
        Boolean include_multimappers = false
        Boolean include_introns = false
        Array[File] read1_rna
        File genes_annotation_bed
        File gtf
        File idx_tar_rna
        Int? cpus_rna
        String? gene_naming = "gene_name"

        # Group UMI
        Boolean remove_single_umi = true
        String mode = "fast"
        Int cutoff_rna = 100

        # Lib_size QC
        Boolean qc = false


        # DORCs specific inputs
        File peak_set
        Int? cpus_dorcs = 4
        String save_plots_to_dir = "TRUE"
        String output_filename = "output.ipynb"

        # Seurat filters
        Int minFeature_RNA = 200
        Int maxFeature_RNA = 2500
        Float percentMT_RNA = 5
        Int minCells_RNA = 3

        # DORCs filter
        Int dorcGeneCutOff = 10
        Float fripCutOff = 0.3
        Float corrPVal = 0.05
        Int topNGene = 20

        # Regulatory region around TSS. Default is +/- 50Kb
        Int windowPadSize = 50000
        #Int bootstraps = 100

        String docker_image_dorcs = "us.gcr.io/buenrostro-share-seq/dorcs_task_find_dorcs"
        Int mem_gb_dorcs = 16
    }

    call share_rna.wf_rna as rna{
        input:
            read1 = read1_rna,
            idx_tar = idx_tar_rna,
            prefix = prefix,
            genome_name = genome_name,
            cpus = cpus_rna,
            # Update RGID
            multimappers = multimappers,
            # Assign features
            include_multimappers = include_multimappers,
            include_introns = include_introns,
            gtf = gtf,
            gene_naming = gene_naming,
            # Group UMI
            remove_single_umi = remove_single_umi,
            mode = mode,
            cutoff = cutoff_rna,
            # Lib_size QC
            qc = qc,
            genes_annotation_bed = genes_annotation_bed
    }
    call share_atac.wf_atac as atac{
        input:
            read1 = read1_atac,
            read2 = read2_atac,
            chrom_sizes = chrom_sizes,
            idx_tar = idx_tar_atac,
            tss_bed = tss_bed,
            prefix = prefix,
            genome_name = genome_name,
            cutoff = cutoff_atac,
            cpus = cpus_atac
    }
    call find_dorcs.wf_dorcs as dorcs{
        input:
            rna_matrix = rna.share_rna_h5_matrix,
            atac_fragments = atac.share_atac_fragments_filtered,
            peak_file = peak_set,

            genome = genome_name,
            n_cores = cpus_dorcs,
            save_plots_to_dir = save_plots_to_dir,
            output_filename = output_filename,

            minFeature_RNA = minFeature_RNA,
            maxFeature_RNA = maxFeature_RNA,
            percentMT_RNA = percentMT_RNA,
            minCells_RNA = minCells_RNA,

            dorcGeneCutOff = dorcGeneCutOff,
            fripCutOff = fripCutOff,
            corrPVal = corrPVal,
            topNGene = topNGene,

            # Regulatory region around TSS. Default is +/- 50Kb
            windowPadSize = windowPadSize,
            #bootstraps = bootstraps,
            mem_gb = mem_gb_dorcs
    }

    output{
        File share_rna_alignment_raw = rna.rna_alignment
        File share_rna_alignment_index = rna.rna_alignment_index
        File share_rna_alignment_log = rna.rna_alignment_log

        File share_rna_reheaded_alignment = rna.rna_reheaded_alignment
        File share_rna_reheaded_alignment_index = rna.rna_reheaded_alignment_index

        File share_rna_featurecount_alignment = rna.rna_featurecount_alignment
        File share_rna_featurecount_alignment_index = rna.rna_featurecount_alignment_index
        #File share_rna_featurecount_log = rna.rna_featurecount_log
        File share_rna_featurecount_txt = rna.rna_featurecount_txt
        File share_rna_featurecount_summary = rna.rna_featurecount_summary

        File share_rna_umi_barcodes = rna.rna_umi_barcodes_filtered
        File share_rna_umi_bed_filtered = rna.rna_umi_bed_filtered
        File share_rna_umi_bed_unfiltered = rna.rna_umi_bed_unfiltered
        File share_rna_umi_counts_filtered = rna.rna_umi_counts_filtered
        File share_rna_umi_counts_unfiltered = rna.rna_umi_counts_unfiltered

        File share_rna_qc_reads_distribution = rna.rna_qc_reads_distribution
        File share_rna_qc_reads_distribution2 = rna.rna_qc_reads_distribution2
        File share_rna_qc_reads_distribution_plot = qrnac_rna.rna_qc_reads_distribution_plot

        File share_rna_qc_library_counts = rna.lib_size_counts
        File share_rna_qc_library_duplicates = rna.lib_size_log
        Array[File] share_rna_qc_library_plots = rna.plots

        File share_rna_h5_matrix = rna.h5_matrix
        Array[File] share_rna_umi_qc_plots = rna.umi_qc_plots

        File share_rna_seurat_notebook_output = rna.notebook_output
        File share_rna_seurat_notebook_log = rna.notebook_log
        File share_rna_seurat_papermill_log = rna.papermill_log
        File? share_rna_seurat_violin_plot = rna.seurat_violin_plot
        File? share_rna_seurat_mitochondria_qc_plot = rna.seurat_mitochondria_qc_plot
        File? share_rna_seurat_features_plot = rna.seurat_features_plot
        File? share_rna_seurat_PCA_dim_loadings_plot = rna.seurat_PCA_dim_loadings_plot
        File? share_rna_seurat_PCA_plot = rna.seurat_PCA_plot
        File? share_rna_seurat_heatmap_plot = rna.seurat_heatmap_plot
        File? share_rna_seurat_jackstraw_plot = rna.seurat_jackstraw_plot
        File? share_rna_seurat_elbow_plot = rna.seurat_elbow_plot
        File? share_rna_seurat_umap_plot = rna.seurat_umap_plot
        File? share_rna_seurat_obj = rna.seurat_obj
        File? share_rna_plots_zip = rna.plots_zip

        File share_atac_alignment_raw = atac.atac_alignment
        File share_atac_alignment_raw_index = atac.atac_alignment_index
        File share_atac_alignment_log = atac.atac_alignment_log

        File share_atac_alignment_filtered = atac.atac_alignment_filtered
        File share_atac_alignment_filtered_index = atac.atac_alignment_filtered_index
        File share_atac_fragments_raw = atac.atac_fragments_raw


        File share_atac_barcodes = atac.atac_barcodes
        File share_atac_fragments_filtered = atac.atac_fragments_filtered
        File share_atac_counts_raw = atac.atac_counts_unfiltered
        File share_atac_counts_filtered = atac.atac_counts_filtered

        File share_atac_qc_library_counts = atac.lib_size_counts
        File share_atac_qc_library_duplicates = atac.lib_size_log
        Array[File] share_atac_qc_library_plots = atac.plots

        File share_atac_qc_final = atac.atac_final_stats
        File share_atac_qc_hist_plot = atac.atac_final_hist_pdf
        File share_atac_qc_hist_txt = atac.atac_final_hist
        File share_atac_qc_tss_enrichment = atac.atac_tss_pileup_png

        File share_atac_archr_notebook_output = atac.notebook_output
        File share_atac_archr_notebook_log = atac.notebook_log
        File share_atac_archr_papermill_log = atac.papermill_log
        File? share_atac_archr_gene_heatmap_plot = atac.archr_heatmap_plot
        File? share_atac_archr_tss_enrichment_raw = atac.archr_TSS_uniq_frags_plot
        File? share_atac_archr_tss_enrichment_filtered = atac.archr_TSS_uniq_frags_filtered_plot
        File? share_atac_archr_fragment_size_plot = atac.archr_fragment_size_dist_plot
        File? share_atac_archr_doublet_plot = atac.archr_doublet_plot
        File? share_atac_archr_umap_plot = atac.archr_umap_plot
        File? share_atac_archr_arrow = atac.archr_arrow
        File? share_atac_archr_obj = atac.archr_obj
        File? share_atac_archr_plots_zip = atac.plots_zip

        File notebook_output = dorcs.notebook_output
        File seurat_violin_plot = dorcs.seurat_violin_plot
        File j_plot = dorcs.j_plot
        File plots_zip = dorcs.plots_zip
        File dorcs_genes_summary = dorcs.dorcs_genes_summary
        File dorcs_regions_summary = dorcs.dorcs_regions_summary

    }

}


task merge_fastqs{
    input{
        Array[File] atac_read1
        Array[File] atac_read2
        Array[File] rna_read1
        String? prefix
    }
    command{
        cat ${sep=' ' atac_read1} > ${prefix + "."}merged.atac.R1.fq.gz
        cat ${sep=' ' atac_read2} > ${prefix + "."}merged.atac.R2.fq.gz
        cat ${sep=' ' rna_read1} > ${prefix + "."}merged.rna.R1.fq.gz
    }
    output{
        File merged_atac_fastq_R1 = glob('*.merged.atac.R1.fq.gz')[0]
        File merged_atac_fastq_R2 = glob('*.merged.atac.R2.fq.gz')[0]
        File merged_rna_fastq_R1 = glob('*.merged.rna.R1.fq.gz')[0]
    }
}


# Task to report errors to user.
# From https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/chip.wdl
task raise_exception {
  input {
    String msg
    Array[String]? vals
  }
  command {
    echo -e "\n* Error: ${msg}\n" >&2
    echo -e "* Vals: ${sep=',' vals}\n" >&2
    exit 2
  }
  output {
    String error_msg = '${msg}'
  }
  runtime {
    maxRetries : 0
    cpu : 1
    memory : '2 GB'
    time : 1
    disks : 'local-disk 10 SSD'
  }
}
