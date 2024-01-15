version 1.0

import "../tasks/share_task_correct_fastq.wdl" as share_task_correct_fastq
import "../tasks/task_starsolo.wdl" as task_starsolo
import "../tasks/task_generate_h5.wdl" as task_generate_h5
import "../tasks/task_qc_rna.wdl" as task_qc_rna
import "../tasks/task_log_rna.wdl" as task_log_rna
import "../tasks/task_seurat.wdl" as task_seurat

# Import the tasks called by the pipeline
workflow wf_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) and Sai Ma @ Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Sub-workflow to process the RNA portion of SHARE-seq libraries.'
    }

    input {
        # RNA sub-workflow inputs
        String? subpool
        String prefix
        String genome_name
        String chemistry
        String pipeline_modality = "full"
        File whitelist
        Array[File] read1
        Array[File] read2

        # Do we want to use the multimappers assigned with EM method
        Boolean multimappers = false

        # Correct-specific inputs
        Boolean correct_barcodes = true
        Boolean? paired_rna = false
        # Runtime parameters
        Int? correct_cpus
        Float? correct_disk_factor
        Float? correct_memory_factor
        String? correct_docker_image

        # Align-specific inputs
        File idx_tar
        String? barcode_tag
        String soloUMIdedup = "1MM_All"
        String soloMultiMappers = "Unique EM"
        Int outFilterMultimapNmax = 20
        Float outFilterScoreMinOverLread = 0.3
        Float outFilterMatchNminOverLread = 0.3
        Int? winAnchorMultimapNmax
        Float? outFilterMismatchNoverReadLmax
        Int? outFilterScoreMin
        String? soloBarcodeMate # 2 for SHARE
        String? clip5pNbases  # 0 34 for SHARE
        # Runtime parameters
        Int? align_cpus
        Float? align_disk_factor
        Float? align_memory_factor
        String? align_docker_image
        
        # Generate h5-specific inputs
        String? gene_naming
        # Runtime parameters
        Float? generate_h5_disk_factor
        Float? generate_h5_memory_factor
        String? generate_h5_docker_image
    
        # QC-specific inputs
        Int? umi_cutoff
        Int? gene_cutoff
        # Runtime parameters
        Int? qc_cpus
        Float? qc_disk_factor
        Float? qc_memory_factor
        String? qc_docker_image

        # Seurat-specific inputs
        Int? seurat_min_features
        Float? seurat_percent_mt
        Int? seurat_min_cells
        Int? seurat_umap_dim
        Float? seurat_umap_resolution
        # Runtime parameters
        Float? seurat_disk_factor
        Float? seurat_memory_factor
        String? seurat_docker_image
    }

    if ( chemistry == "shareseq" && correct_barcodes ) {
        scatter (read_pair in zip(read1, read2)) {
            call share_task_correct_fastq.share_correct_fastq as correct {
                input:
                    fastq_R1 = read_pair.left,
                    fastq_R2 = read_pair.right,
                    whitelist = whitelist,
                    sample_type = "RNA",
                    pkr = subpool,
                    prefix = prefix,
                    paired_rna = paired_rna,
                    cpus = correct_cpus,
                    disk_factor = correct_disk_factor,
                    memory_factor = correct_memory_factor,
                    docker_image = correct_docker_image
            }
        }
    }

    if (  "~{pipeline_modality}" != "no_align" ) {
        call task_starsolo.rna_align as align {
            input:
                chemistry = chemistry,
                fastq_R1 = select_first([correct.corrected_fastq_R1, read1]),
                fastq_R2 = select_first([correct.corrected_fastq_R2, read2]),
                whitelist = whitelist,
                soloMultiMappers = soloMultiMappers,
                soloUMIdedup = soloUMIdedup,
                outFilterMultimapNmax = outFilterMultimapNmax,
                outFilterScoreMinOverLread = outFilterScoreMinOverLread,
                outFilterMatchNminOverLread = outFilterMatchNminOverLread,
                winAnchorMultimapNmax = winAnchorMultimapNmax,
                outFilterMismatchNoverReadLmax = outFilterMismatchNoverReadLmax,
                outFilterScoreMin = outFilterScoreMin,
                soloBarcodeMate = soloBarcodeMate,
                clip5pNbases = clip5pNbases,
                genome_name = genome_name,
                genome_index_tar = idx_tar,
                prefix = prefix,
                cpus = align_cpus,
                disk_factor = align_disk_factor,
                memory_factor = align_memory_factor,
                docker_image = align_docker_image
        }

        call task_generate_h5.generate_h5 as generate_h5 {
            input:
                tar = align.raw_tar,
                genome_name = genome_name,
                prefix = prefix,
                pkr = subpool,
                multimappers = multimappers,
                gene_naming = gene_naming,
                disk_factor = generate_h5_disk_factor,
                memory_factor = generate_h5_memory_factor,
                docker_image = generate_h5_docker_image
        }

        call task_qc_rna.qc_rna as qc_rna {
            input:
                bam = align.output_bam,
                mtx_tar = align.raw_tar,
                umi_cutoff = umi_cutoff,
                gene_cutoff = gene_cutoff,
                subpool = subpool,
                barcode_tag = barcode_tag,
                genome_name = genome_name,
                prefix = prefix,
                cpus = qc_cpus,
                disk_factor = qc_disk_factor,
                memory_factor = qc_memory_factor,
                docker_image = qc_docker_image
        }

        call task_log_rna.log_rna as log_rna {
        input:
            alignment_log = align.log_final_out,
            dups_log = qc_rna.rna_duplicates_log
        }

        if ( "~{pipeline_modality}" == "full") {
            call task_seurat.seurat as seurat {
                input:
                    rna_matrix = generate_h5.h5_matrix,
                    genome_name = genome_name,
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

    output {
        # Correction/trimming
        Array[File]? rna_read1_processed = correct.corrected_fastq_R1
        Array[File]? rna_read2_processed = correct.corrected_fastq_R2
        
        File? task_starsolo_output_bam = align.output_bam
        File? rna_alignment_log = align.log_final_out
        File? task_starsolo_log_out = align.log_out
        File? task_starsolo_log_progress_out = align.log_progress_out
        File? task_starsolo_output_sj = align.output_sj
        File? task_starsolo_barcodes_stats = align.barcodes_stats
        File? task_starsolo_features_stats = align.features_stats
        File? task_starsolo_summary_csv = align.summary_csv
        File? task_starsolo_umi_per_cell = align.umi_per_cell
        File? task_starsolo_raw_tar = align.raw_tar

        File? rna_h5 = generate_h5.h5_matrix

        File? rna_barcode_metadata  = qc_rna.rna_barcode_metadata
        File? rna_duplicates_log = qc_rna.rna_duplicates_log
        File? rna_umi_barcode_rank_plot = qc_rna.rna_umi_barcode_rank_plot
        File? rna_gene_barcode_rank_plot = qc_rna.rna_gene_barcode_rank_plot
        File? rna_gene_umi_scatter_plot = qc_rna.rna_gene_umi_scatter_plot

        File? rna_seurat_notebook_output = seurat.notebook_output
        File? rna_seurat_notebook_log = seurat.notebook_log
        File? rna_seurat_raw_violin_plot = seurat.seurat_raw_violin_plot
        File? rna_seurat_filtered_violin_plot = seurat.seurat_filtered_violin_plot
        File? rna_seurat_raw_qc_scatter_plot = seurat.seurat_raw_qc_scatter_plot
        File? rna_seurat_filtered_qc_scatter_plot = seurat.seurat_filtered_qc_scatter_plot
        File? rna_seurat_variable_genes_plot = seurat.seurat_variable_genes_plot
        File? rna_seurat_PCA_dim_loadings_plot = seurat.seurat_PCA_dim_loadings_plot
        File? rna_seurat_PCA_plot = seurat.seurat_PCA_plot
        File? rna_seurat_heatmap_plot = seurat.seurat_heatmap_plot
        File? rna_seurat_jackstraw_plot = seurat.seurat_jackstraw_plot
        File? rna_seurat_elbow_plot = seurat.seurat_elbow_plot
        File? rna_seurat_umap_cluster_plot = seurat.seurat_umap_cluster_plot
        File? rna_seurat_umap_rna_count_plot = seurat.seurat_umap_rna_count_plot
        File? rna_seurat_umap_gene_count_plot = seurat.seurat_umap_gene_count_plot
        File? rna_seurat_umap_mito_plot = seurat.seurat_umap_mito_plot
        File? rna_seurat_obj = seurat.seurat_filtered_obj
        File? rna_plots_zip = seurat.plots_zip

        File? rna_qc_metrics = log_rna.rna_qc_metrics
    }
}
