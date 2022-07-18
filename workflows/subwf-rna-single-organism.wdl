version 1.0


# Import the tasks called by the pipeline
import "../tasks/share_task_star.wdl" as share_task_align
import "../tasks/share_task_update_rgid.wdl" as share_task_update_rgid
import "../tasks/share_task_count_rna.wdl" as share_task_feature_counts
import "../tasks/share_task_group_umi.wdl" as share_task_group_umi
import "../tasks/share_task_qc_rna.wdl" as share_task_qc_rna
import "../tasks/share_task_qc_library.wdl" as share_task_qc_library
import "../tasks/share_task_generate_h5.wdl" as share_task_generate_h5
import "../tasks/share_task_seurat.wdl" as share_task_seurat

workflow wf_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) and Sai Ma @ Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Sub-workflow to process the RNA portion of SHARE-seq libraries.'
    }

    input {
        # RNA Sub-worflow inputs

        # Align
        Array[File] read1
        File idx_tar
        String prefix = "shareseq-project"
        String genome_name
        Int? cpus = 16
        String? docker
        # Update RGID
        Boolean multimappers = false
        # Assign features
        Boolean include_multimappers = false
        Boolean include_introns = false
        File gtf
        String gene_naming = "gene_name"
        # Group UMI
        Boolean remove_single_umi = true
        String mode = "fast"
        Int cutoff = 100
        # Lib_size QC
        Boolean qc = false
        File genes_annotation_bed
        # Seurat
        Int umap_dim = 10
        Float umap_resolution = 0.5
    }

    call share_task_align.share_rna_align as align {
        input:
            fastq_R1 = read1,
            genome_name = genome_name,
            genome_index_tar = idx_tar,
            prefix = prefix,
            cpus = cpus
    }

    call share_task_update_rgid.share_rna_update_rgid as update_rgid{
        input:
            bam = align.rna_alignment,
            multimapper = multimappers,
            genome_name = genome_name,
            prefix = prefix
    }

    call share_task_feature_counts.feature_counts_rna as count{
        input:
            multimapper = include_multimappers,
            intron = include_introns,
            bam = update_rgid.rna_reheaded_alignment,
            gtf = gtf,
            gene_naming = gene_naming,
            genome_name = genome_name,
            prefix = prefix
    }

    call share_task_group_umi.group_umi_rna as group_umi{
        input:
            bam = count.rna_featurecount_alignment,
            mode = mode,
            cutoff = cutoff,
            remove_single_umi = remove_single_umi,
            genome_name = genome_name,
            prefix = prefix
    }

    # TODO: the genes annotation needs to be passed just like the gtf
    # Check which one is the bam in input
    call share_task_qc_rna.qc_rna as qc_rna{
        input:
            qc = qc,
            bam = count.rna_featurecount_alignment,
            genes_annotations_bed = genes_annotation_bed,
            genome_name = genome_name,
            prefix = prefix
    }

    call share_task_qc_library.qc_library as qc_library {
        input:
            raw_counts = group_umi.rna_umi_counts_unfiltered,
            filtered_counts = group_umi.rna_umi_counts_filtered,
            cutoff = cutoff,
            genome_name = genome_name,
            prefix = prefix,
            assay = "RNA"
    }

    call share_task_generate_h5.generate_h5 as generate_h5 {
        input:
            filtered_bed = group_umi.rna_umi_bed_filtered
    }

    call share_task_seurat.seurat as seurat{
        input:
            rna_matrix = generate_h5.h5_matrix,
            genome_name = genome_name,
            umap_dim = umap_dim,
            umap_resolution = umap_resolution
    }

    output {
        File share_rna_alignment_raw = align.rna_alignment
        File share_rna_alignment_index = align.rna_alignment_index
        File share_rna_alignment_log = align.rna_alignment_log

        File share_rna_reheaded_alignment = update_rgid.rna_reheaded_alignment
        File share_rna_reheaded_alignment_index = update_rgid.rna_reheaded_alignment_index

        File share_rna_featurecount_alignment = count.rna_featurecount_alignment
        File share_rna_featurecount_alignment_index = count.rna_featurecount_alignment_index
        #File share_rna_featurecount_log = count.rna_featurecount_log
        File share_rna_featurecount_exon_txt = count.rna_featurecount_exon_txt
        File? share_rna_featurecount_intron_txt = count.rna_featurecount_intron_txt
        #File share_rna_featurecount_summary = count.rna_featurecount_summary

        File share_rna_umi_barcodes = group_umi.rna_umi_barcodes_filtered
        File share_rna_umi_bed_filtered = group_umi.rna_umi_bed_filtered
        File share_rna_umi_bed_unfiltered = group_umi.rna_umi_bed_unfiltered
        File share_rna_umi_counts_filtered = group_umi.rna_umi_counts_filtered
        File share_rna_umi_counts_unfiltered = group_umi.rna_umi_counts_unfiltered
        File share_rna_umi_rm_dup_log = group_umi.rna_umi_rm_dup_log

        File share_rna_qc_reads_distribution = qc_rna.rna_qc_reads_distribution
        File share_rna_qc_reads_distribution2 = qc_rna.rna_qc_reads_distribution2
        File share_rna_qc_reads_distribution_plot = qc_rna.rna_qc_reads_distribution_plot

        File share_rna_qc_library_counts = qc_library.lib_size_counts
        File share_rna_qc_library_duplicates = qc_library.lib_size_log
        Array[File] share_rna_qc_library_plots = qc_library.plots

        File share_rna_h5_matrix = generate_h5.h5_matrix
        Array[File] share_rna_umi_qc_plots = generate_h5.umi_qc_plots

        File share_rna_seurat_notebook_output = seurat.notebook_output
        File share_rna_seurat_notebook_log = seurat.notebook_log
        #File share_rna_seurat_papermill_log = seurat.papermill_log
        File? share_rna_seurat_filtered_violin_plot = seurat.seurat_filtered_violin_plot
        File? share_rna_seurat_filtered_qc_scatter_plot = seurat.seurat_filtered_qc_scatter_plot
        File? share_rna_seurat_variable_genes_plot = seurat.seurat_variable_genes_plot
        File? share_rna_seurat_PCA_dim_loadings_plot = seurat.seurat_PCA_dim_loadings_plot
        File? share_rna_seurat_PCA_plot = seurat.seurat_PCA_plot
        File? share_rna_seurat_heatmap_plot = seurat.seurat_heatmap_plot
        File? share_rna_seurat_jackstraw_plot = seurat.seurat_jackstraw_plot
        File? share_rna_seurat_elbow_plot = seurat.seurat_elbow_plot
        File? share_rna_seurat_umap_cluster_plot = seurat.seurat_umap_cluster_plot
        File? share_rna_seurat_obj = seurat.seurat_obj
        File? share_rna_plots_zip = seurat.plots_zip
    }
}
