version 1.0

# Import the tasks called by the pipeline
import "../tasks/share_task_bowtie2.wdl" as share_task_align
import "../tasks/share_task_bam2bed.wdl" as share_task_bam2bed
import "../tasks/share_task_count_atac.wdl" as share_task_count
import "../tasks/share_task_qc_atac.wdl" as share_task_qc_atac
import "../tasks/share_task_qc_library.wdl" as share_task_qc_library


workflow wf_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) and Sai Ma @ Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Sub-workflow to process the ATAC portion of SHARE-seq libraries.'
    }

    input {
        # ATAC Sub-worflow inputs
        File read1
        File read2
        File chrom_sizes
        File idx_tar
        File tss_bed
        String prefix = "shareseq-project"
        String genome_name
        Int cutoff
        Int? cpus = 4
        String docker = "polumechanos/share-seq"
    }

    call share_task_align.share_atac_align as align {
        input:
            fastq_R1 = read1,
            fastq_R2 = read2,
            genome_name = genome_name,
            genome_index = idx_tar,
            prefix = prefix
    }

    call share_task_bam2bed.share_atac_bam2bed as bam2bed {
        input:
            bam = align.atac_alignment_alignment,
            bam_index = align.atac_alignment_index,
            genome_name = genome_name,
            chrom_sizes = chrom_sizes,
            prefix = prefix
    }

    call share_task_count.count_reads_atac as count {
        input:
            cutoff = cutoff,
            bedpe = bam2bed.atac_fragments_raw,
            genome_name = genome_name,
            prefix = prefix,
            cpus = cpus
    }

    call share_task_qc_library.qc_library as qc_library{
        input:
            raw_counts = count.atac_counts_unfiltered,
            filtered_counts = count.atac_counts_filtered,
            cutoff = cutoff,
            genome_name = genome_name,
            prefix = prefix,
            assay = "ATAC"
    }

    call share_task_qc_atac.qc_atac as qc_atac{
        input:
            raw_bam = align.atac_alignment,
            raw_bam_index = align.atac_alignment_index,
            filtered_bam = bam2bed.atac_alignment_filtered,
            filtered_bam_index = bam2bed.atac_alignment_filtered_index,
            tss = tss_bed,
            genome_name = genome_name,
            prefix = prefix,
            cpus = cpus
    }

    output {
        File share_atac_alignment_raw = align.atac_alignment
        File share_atac_alignment_raw_index = align.atac_alignment_index
        File share_atac_alignment_log = align.atac_alignment_log

        File share_atac_alignment_filtered = bam2bed.atac_alignment_filtered
        File share_atac_alignment_filtered_index = bam2bed.atac_alignment_filtered_index
        File share_atac_fragments_raw = bam2bed.atac_fragments_raw


        File share_atac_barcodes = count.atac_barcodes
        File share_atac_fragments_filtered = count.atac_fragments_filtered
        File share_atac_counts_raw = count.atac_counts_unfiltered
        File share_atac_counts_filtered = count.atac_counts_filtered

        File share_atac_qc_library_counts = qc_library.lib_size_counts
        File share_atac_qc_library_duplicates = qc_library.lib_size_log
        Array[File] share_atac_qc_library_plots = qc_library.plots

        File share_atac_qc_final = qc_atac.atac_final_stats
        File share_atac_qc_hist_plot = qc_atac.atac_final_hist_pdf
        File share_atac_qc_hist_txt = qc_atac.atac_final_hist
        File share_atac_qc_tss_enrichment = qc_atac.atac_tss_pileup_png

    }
}
