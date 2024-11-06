version 1.0

struct RNA_outputs {
    File rna_qc_metrics
    Int rna_input_reads
    Int rna_aligned_reads
    Int rna_aligned_uniquely
    Int rna_aligned_multimap
    Int rna_unaligned_reads
    Int rna_homopolymer_umis
    Int rna_nonmatch_barcodes
    Int rna_exact_match_barcodes
    Int rna_mismatch_barcodes
    Float rna_frac_valid_barcodes
    Float rna_sequencing_saturation
    Float rna_frac_q30_bases_in_cb_umi
    Float rna_frac_q30_bases_in_read
    Float rna_starsolo_frig
    Int rna_estimated_cells
    Float rna_frac_unique_reads_in_cells
    Int rna_median_reads_per_cell
    Int rna_median_umis_per_cell
    Int rna_genes
    Int rna_unique_reads_mapped_to_genes
    Float rna_qc_rna_frig
    Int rna_duplicate_reads
    Float rna_percent_duplicates
    Float rna_percent_mitochondrial
}