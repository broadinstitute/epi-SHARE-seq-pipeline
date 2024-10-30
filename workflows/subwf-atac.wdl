version 1.0

# Import the tasks called by the pipeline
import "../tasks/task_chromap_read_format.wdl" as task_chromap_read_format
import "../tasks/task_chromap.wdl" as task_align_chromap
import "../tasks/task_chromap_bam.wdl" as task_align_chromap_bam
import "../tasks/task_qc_atac.wdl" as task_qc_atac
import "../tasks/task_make_track.wdl" as task_make_track

workflow wf_atac {
    meta {
        version: 'rc-v2.0.0'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) @ Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard combinomics pipeline: Sub-workflow to process the sc-ATAC libraries.'
    }

    input {
        File chrom_sizes
        File reference_index_tar_gz
        File? tss_bed
        Int? mapq_threshold = 30
        String chemistry
        File gtf
        String prefix = "combinomics"
        String? subpool
        String genome_name
        File? barcode_conversion_dict # For 10X multiome

        # Align-specific inputs
        Array[File] read1
        Array[File] read2
        Array[File] fastq_barcode
        Int? align_multimappers
        File reference_fasta
        File whitelist
        Boolean? remove_pcr_duplicates = true
        Boolean? remove_pcr_duplicates_at_cell_level = false
        Boolean? remove_pcr_duplicates_at_bulk_level = true
        Boolean? Tn5_shift = false
        Boolean? low_mem = true
        Boolean? bed_output = true
        Boolean? trim_adapters = true
        Int? max_insert_size = 2000
        Int? quality_filter = 0
        Int? bc_error_threshold = 1
        Float? bc_probability_threshold = 0.9
        String? read_format
        # Runtime parameters
        Int? align_bam_cpus
        Float? align_bam_disk_factor = 8.0
        Float? align_bam_memory_factor = 0.15

        Int? align_cpus
        Float? align_disk_factor = 8.0
        Float? align_memory_factor = 0.15
        String? align_docker_image

        Int qc_fragment_min_cutoff = 10
        # Runtime parameters
        Int qc_cpus = 16
        Float qc_disk_factor = 8.0
        Float qc_memory_factor = 0.15
        String qc_docker_image

        # Make track inputs
        # Runtime parameters
        Int make_track_cpus = 8
        Float make_track_disk_factor = 4
        Float make_track_memory_factor = 0.3
        String make_track_docker_image

        Boolean generate_tracks = false
        Boolean generate_bam_alignment = false
    }
    
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
            reference_index_tar_gz = reference_index_tar_gz,
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
            remove_pcr_duplicates_at_bulk_level = remove_pcr_duplicates_at_bulk_level,
            Tn5_shift = Tn5_shift,
            low_mem = low_mem,
            bed_output = bed_output,
            max_insert_size = max_insert_size,
            quality_filter = quality_filter,
            bc_error_threshold = bc_error_threshold,
            bc_probability_threshold = bc_probability_threshold,
            read_format = select_first([read_format, get_chromap_read_format.read_format])
    }

    if ( generate_bam_alignment ) {
        call task_align_chromap_bam.atac_align_chromap as generate_bam {
                input:
                    fastq_R1 = read1,
                    fastq_R2 = read2,
                    fastq_barcode = fastq_barcode,
                    reference_fasta = reference_fasta,
                    reference_index_tar_gz = reference_index_tar_gz,
                    trim_adapters = trim_adapters,
                    genome_name = genome_name,
                    multimappers = align_multimappers,
                    barcode_inclusion_list = whitelist,
                    barcode_conversion_dict = barcode_conversion_dict,
                    prefix = prefix,
                    disk_factor = align_bam_disk_factor,
                    memory_factor = align_bam_memory_factor,
                    cpus = align_bam_cpus,
                    docker_image = align_docker_image,
                    remove_pcr_duplicates = remove_pcr_duplicates,
                    remove_pcr_duplicates_at_cell_level = remove_pcr_duplicates_at_cell_level,
                    remove_pcr_duplicates_at_bulk_level = remove_pcr_duplicates_at_bulk_level,
                    Tn5_shift = Tn5_shift,
                    low_mem = low_mem,
                    max_insert_size = max_insert_size,
                    mapq_threshold = mapq_threshold,
                    bc_error_threshold = bc_error_threshold,
                    bc_probability_threshold = bc_probability_threshold,
                    read_format = select_first([read_format, get_chromap_read_format.read_format])
        }
    }

    call task_qc_atac.qc_atac as qc_atac{
        input:
            fragment_file = align.atac_fragment_file,
            fragment_file_index = align.atac_fragment_file_index,
            chrom_sizes = chrom_sizes,
            gtf = gtf,
            fragment_min_cutoff = qc_fragment_min_cutoff,
            prefix = prefix,
            cpus = qc_cpus,
            disk_factor = qc_disk_factor,
            docker_image = qc_docker_image,
            memory_factor = qc_memory_factor
        }

    if ( generate_tracks ) {
        call task_make_track.make_track as track {
            input:
                fragments = align.atac_fragment_file,
                chrom_sizes = chrom_sizes,
                genome_name = genome_name,
                prefix = prefix,
                cpus = make_track_cpus,
                disk_factor = make_track_disk_factor,
                docker_image = make_track_docker_image,
                memory_factor = make_track_memory_factor
        }
    }

    output {
        # Bam
        File? atac_bam = generate_bam.atac_bam
        File? atac_bam_index = generate_bam.atac_bam_index
        File? atac_bam_alignment_stats = generate_bam.atac_alignment_log

        # Align
        File atac_fragment_file = align.atac_fragment_file
        File atac_fragment_file_index = align.atac_fragment_file_index
        File atac_fragment_file_sorted_by_barcode = align.atac_fragment_file_sorted_by_barcode
        File atac_align_barcode_statistics = align.atac_align_barcode_statistics
        File atac_align_log = align.atac_alignment_log
        Float atac_pcr_duplicates_percentage = align.atac_pcr_duplicates_percentage
        Int atac_reads_count = align.atac_reads_count
        Int atac_mapped_reads = align.atac_mapped_reads
        Int atac_unique_reads = align.atac_unique_reads
        Int atac_multi_mapping_reads = align.atac_multi_mapping_reads
        Int atac_corrected_barcodes = align.atac_corrected_barcodes
        Int atac_unique_mappings_fragments = align.atac_unique_mappings_fragments
        Int atac_multi_mappings_fragments = align.atac_multi_mappings_fragments
        Int atac_final_number_of_fragments = align.atac_final_number_of_fragments
        Int atac_unique_barcodes_unfiltered = align.atac_unique_barcodes_unfiltered
        String atac_alingment_tool_verion = align.atac_chromap_verion

        # QC
        File atac_qc_fragment_size_distribution_plot = qc_atac.atac_fragment_size_distribution_plot
        File atac_qc_tss_enrichment_library_plot = qc_atac.atac_tss_enrichment_library_plot
        File atac_qc_fraction_of_duplicates_distribution_plot = qc_atac.atac_fraction_of_duplicates_distribution_plot
        File atac_qc_fraction_of_mito_distribution_plot = qc_atac.atac_fraction_of_mito_distribution_plot
        File atac_qc_atac_knee_plot= qc_atac.atac_knee_plot
        File atac_qc_n_fragment_vs_tss_enrichment_plot = qc_atac.atac_n_fragment_vs_tss_enrichment_plot
        File atac_qc_atac_n_fragment_vs_tss_enrichment_filtered_plot = qc_atac.atac_n_fragment_vs_tss_enrichment_filtered_plot
        File atac_qc_umap_leiden_plot = qc_atac.atac_umap_leiden_plot
        File atac_qc_snapatac2_h5ad= qc_atac.atac_snapatac2_h5ad
        File atac_qc_barcode_metrics = qc_atac.atac_barcode_metrics
        Float atac_qc_library_tss_overlap = qc_atac.atac_library_tss_overlap
        Float atac_qc_library_tsse = qc_atac.atac_library_tsse
        
        # Track
        File? atac_track_bigwig = track.atac_track_bigwig
        File? atac_track_bigwig_no_nucleosome = track.atac_track_bigwig_no_nucleosome
        File? atac_track_bigwig_mono_nucleosome = track.atac_track_bigwig_mono_nucleosome
        File? atac_track_bigwig_multi_nucleosome = track.atac_track_bigwig_multi_nucleosome
    }
}
