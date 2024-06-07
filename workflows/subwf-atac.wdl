version 1.0

# Import the tasks called by the pipeline
import "../tasks/task_chromap_read_format.wdl" as task_chromap_read_format
import "../tasks/task_chromap.wdl" as task_align_chromap
import "../tasks/task_qc_atac.wdl" as task_qc_atac
import "../tasks/task_make_track.wdl" as task_make_track
import "../tasks/task_log_atac.wdl" as task_log_atac
import "../tasks/task_archr.wdl" as task_archr


workflow wf_atac {
    meta {
        version: 'rc-v2.0.0'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) @ Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard combinomics pipeline: Sub-workflow to process the sc-ATAC libraries.'
    }

    input {
        File chrom_sizes
        File genome_index_tar
        File tss_bed
        Int? mapq_threshold = 30
        String? barcode_tag = "CB"
        String? barcode_tag_fragments
        String chemistry
        File? gtf
        String? prefix = "sample"
        String? subpool
        String genome_name
        Int? cutoff
        String pipeline_modality = "full"
        File? barcode_conversion_dict # For 10X multiome

        # Align-specific inputs
        Array[File] read1
        Array[File] read2
        Array[File] fastq_barcode
        Int? align_multimappers
        File reference_fasta
        File whitelist
        Boolean? remove_pcr_duplicates = true
        Boolean? remove_pcr_duplicates_at_cell_level = true
        Boolean? Tn5_shift = true
        Boolean? low_mem = true
        Boolean? bed_output = true
        Boolean? trim_adapters = true
        Int? max_insert_size = 2000
        Int? quality_filter = 0
        Int? bc_error_threshold = 2
        Float? bc_probability_threshold = 0.9
        String? read_format
        # Runtime parameters
        Int? align_cpus
        Float? align_disk_factor = 8.0
        Float? align_memory_factor = 0.15
        String? align_docker_image

        Int? qc_fragment_min_cutoff
        Int? qc_hist_max_fragment = 5000
        Int? qc_hist_min_fragment = 100
        # Runtime parameters
        Int? qc_cpus = 16
        Float? qc_disk_factor = 8.0
        Float? qc_memory_factor = 0.15
        String? qc_docker_image

        # Make track inputs
        # Runtime parameters
        Int make_track_cpus = 8
        Float? make_track_disk_factor = 4
        Float? make_track_memory_factor = 0.3
        String? make_track_docker_image
                
        # ArchR-specific inputs
        File peak_set
        # Runtime parameters
        Float? archr_disk_factor
        Float? archr_memory_factor 
        String? archr_docker_image

    }

    String barcode_tag_fragments_ = if chemistry=="shareseq" then select_first([barcode_tag_fragments, "XC"]) else select_first([barcode_tag_fragments, barcode_tag])

    if ( "~{chemistry}" == "shareseq" && !defined(read_format)) {
        call task_chromap_read_format.get_chromap_read_format as get_chromap_read_format {
            input:
                fastq_path = fastq_barcode[0]
        }
    }

    call task_align_chromap.atac_align_chromap as align {
        input:
            fastq_R1 = read1,
            fastq_R2 = read2,
            fastq_barcode = fastq_barcode,
            reference_fasta = reference_fasta,
            chrom_sizes = chrom_sizes,
            trim_adapters = trim_adapters,
            genome_name = genome_name,
            subpool = subpool,
            multimappers = align_multimappers,
            barcode_inclusion_list = whitelist,
            barcode_conversion_dict = barcode_conversion_dict,
            prefix = prefix,
            disk_factor = align_disk_factor,
            memory_factor = align_memory_factor,
            cpus = align_cpus,
            docker_image = align_docker_image,
            remove_pcr_duplicates = remove_pcr_duplicates,
            remove_pcr_duplicates_at_cell_level = remove_pcr_duplicates_at_cell_level,
            Tn5_shift = Tn5_shift,
            low_mem = low_mem,
            bed_output = bed_output,
            max_insert_size = max_insert_size,
            quality_filter = quality_filter,
            bc_error_threshold = bc_error_threshold,
            bc_probability_threshold = bc_probability_threshold,
            read_format = select_first([get_chromap_read_format.read_format, read_format])
    }

    call task_qc_atac.qc_atac as qc_atac{
        input:
            fragments = align.atac_fragments,
            fragments_index = align.atac_fragments_index,
            barcode_summary = align.atac_align_barcode_statistics,
            tss = tss_bed,
            gtf = gtf,
            subpool = subpool,
            barcode_conversion_dict = barcode_conversion_dict,
            fragment_min_cutoff = qc_fragment_min_cutoff,
            chrom_sizes = chrom_sizes,
            hist_max_fragment = qc_hist_max_fragment,
            hist_min_fragment = qc_hist_min_fragment,
            genome_name = genome_name,
            prefix = prefix,
            cpus = qc_cpus,
            disk_factor = qc_disk_factor,
            docker_image = qc_docker_image,
            memory_factor = qc_memory_factor
        }

    call task_make_track.make_track as track {
        input:
            fragments = align.atac_fragments,
            chrom_sizes = chrom_sizes,
            genome_name = genome_name,
            prefix = prefix,
            cpus = make_track_cpus,
            disk_factor = make_track_disk_factor,
            docker_image = make_track_docker_image,
            memory_factor = make_track_memory_factor
    }

    call task_log_atac.log_atac as log_atac {
        input:
            alignment_log = align.atac_alignment_log,
            barcode_log = align.atac_align_barcode_statistics,
            prefix = prefix
    }

    if (  "~{pipeline_modality}" == "full" ) {
        call task_archr.archr as archr{
            input:
                atac_frag = align.atac_fragments,
                genome = genome_name,
                peak_set = peak_set,
                prefix = prefix,
                memory_factor = archr_memory_factor,
                disk_factor = archr_disk_factor,
                docker_image = archr_docker_image
        }
    }

    output {
        # Align
        File? atac_alignment_log = align.atac_alignment_log
        File? atac_fragments = align.atac_fragments
        File? atac_fragments_index = align.atac_fragments_index

        # QC
        File? atac_qc_chromap_barcode_metadata = qc_atac.atac_qc_chromap_barcode_metadata
        File? atac_qc_snapatac2_barcode_metadata = qc_atac.atac_qc_snapatac2_barcode_metadata
        File? atac_qc_hist_txt = qc_atac.atac_qc_final_hist
        File? atac_qc_tss_enrichment = qc_atac.atac_qc_tss_enrichment_plot
        File? atac_qc_barcode_rank_plot = qc_atac.atac_qc_barcode_rank_plot
        File? atac_qc_insertion_size_histogram = qc_atac.atac_qc_final_hist_png
        File? atac_qc_tsse_fragments_plot = qc_atac.atac_qc_tsse_fragments_plot
        File? atac_qc_fragment_histogram = qc_atac.atac_qc_fragments_histogram
        
        # Track
        File? atac_track_bigwig = track.atac_track_bigwig
        File? atac_track_bigwig_no_nucleosome = track.atac_track_bigwig_no_nucleosome
        File? atac_track_bigwig_mono_nucleosome = track.atac_track_bigwig_mono_nucleosome
        File? atac_track_bigwig_multi_nucleosome = track.atac_track_bigwig_multi_nucleosome

        # Log
        # Int? atac_total_reads = log_atac.atac_total_reads
        # Int? atac_aligned_uniquely = log_atac.atac_aligned_uniquely
        # Int? atac_unaligned = log_atac.atac_unaligned
        # Int? atac_feature_reads = log_atac.atac_feature_reads
        # Int? atac_duplicate_reads = log_atac.atac_duplicate_reads
        # Float? atac_percent_duplicates = log_atac.atac_pct_dup
        File? atac_qc_metrics_csv = log_atac.atac_statistics_csv

        # ArchR
        File? atac_archr_notebook_output = archr.notebook_output
        File? atac_archr_notebook_log = archr.notebook_log
        File? atac_archr_barcode_metadata = archr.archr_barcode_metadata
        File? atac_archr_raw_tss_enrichment = archr.archr_raw_tss_by_uniq_frags_plot
        File? atac_archr_filtered_tss_enrichment = archr.archr_filtered_tss_by_uniq_frags_plot
        File? atac_archr_raw_fragment_size_plot = archr.archr_raw_frag_size_dist_plot
        File? atac_archr_filtered_fragment_size_plot = archr.archr_filtered_frag_size_dist_plot

        File? atac_archr_umap_doublets = archr.archr_umap_doublets
        File? atac_archr_umap_cluster_plot = archr.archr_umap_cluster_plot
        File? atac_archr_umap_num_frags_plot = archr.archr_umap_num_frags_plot
        File? atac_archr_umap_tss_score_plot = archr.archr_umap_tss_score_plot
        File? atac_archr_umap_frip_plot = archr.archr_umap_frip_plot

        File? atac_archr_gene_heatmap_plot = archr.archr_heatmap_plot
        File? atac_archr_arrow = archr.archr_arrow
        File? atac_archr_obj = archr.archr_raw_obj
        File? atac_archr_plots_zip = archr.plots_zip
    }
}
