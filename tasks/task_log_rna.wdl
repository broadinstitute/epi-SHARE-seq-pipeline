version 1.0

# TASK
# SHARE-rna-log
# Gather information from log files 


task log_rna {
    meta {
        version: 'v0.1'
        author: 'Neva C. Durand (neva@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: log rna task'
    }

    input {
        # This function takes as input the necessary log files and extracts
        # the quality metrics
        File alignment_log
        File barcode_statistics
        File summary_csv
        File qc_rna_statistics
        String? prefix = "sample"
    }

    command <<<
        # get STARsolo alignment log statistics
        input_reads=$(awk -F "|" '$1~/input reads/{print $2}' ~{alignment_log} | tr -d "\t")
        aligned_uniquely=$(awk -F "|" '$1~/Uniquely mapped reads number/{print $2}' ~{alignment_log} | tr -d "\t")
        aligned_multimap=$(awk -F "|" '$1~/Number of reads mapped to multiple loci/{print $2}' ~{alignment_log} | tr -d "\t")
        aligned_reads=$(($aligned_uniquely + $aligned_multimap))
        unaligned_reads=$(($input_reads - $aligned_reads))

        echo "$input_reads" > input_reads.txt
        echo "$aligned_uniquely" > aligned_uniquely.txt
        echo "$aligned_multimap" > aligned_multimap.txt
        echo "$aligned_reads" > aligned_reads.txt
        echo "$unaligned_reads" > unaligned_reads.txt

        # get STARsolo barcode statistics
        awk -F " " '$1~/noUMIhomopolymer/{print $2}' ~{barcode_statistics} > homopolymer_umis.txt
        awk -F " " '$1~/noNoWLmatch/{print $2}' ~{barcode_statistics} > nonmatch_barcodes.txt
        awk -F " " '$1~/yesWLmatchExact/{print $2}' ~{barcode_statistics} > exact_match_barcodes.txt
        awk -F " " '$1~/yesOneWLmatchWithMM/{print $2}' ~{barcode_statistics} > mismatch_barcodes.txt

        # get STARsolo summary statistics
        awk -F "," '$1~/Reads With Valid Barcodes/{printf "%.2f\n", $2}' ~{summary_csv} > frac_valid_barcodes.txt
        awk -F "," '$1~/Sequencing Saturation/{printf "%.2f\n", $2}' ~{summary_csv} > sequencing_saturation.txt
        awk -F "," '$1~/Q30 Bases in CB\+UMI/{printf "%.2f\n", $2}' ~{summary_csv} > frac_q30_bases_in_cb_umi.txt
        awk -F "," '$1~/Q30 Bases in RNA read/{printf "%.2f\n", $2}' ~{summary_csv} > frac_q30_bases_in_read.txt
        awk -F "," '$1~/Reads Mapped to GeneFull: Unique\+Multiple GeneFull/{printf "%.2f\n", $2}' ~{summary_csv} > starsolo_frig.txt
        awk -F "," '$1~/Estimated Number of Cells/{print $2}' ~{summary_csv} > estimated_cells.txt
        awk -F "," '$1~/Fraction of Unique Reads in Cells/{printf "%.2f\n", $2}' ~{summary_csv} > frac_unique_reads_in_cells.txt
        awk -F "," '$1~/Median Reads per Cell/{print $2}' ~{summary_csv} > median_reads_per_cell.txt
        awk -F "," '$1~/Median UMI per Cell/{print $2}' ~{summary_csv} > median_umis_per_cell.txt
        awk -F "," '$1~/Median GeneFull per Cell/{print $2}' ~{summary_csv} > median_genes_per_cell.txt
        awk -F "," '$1~/Total GeneFull Detected/{print $2}' ~{summary_csv} > genes.txt

        # get qc_rna statistics
        awk -F "," '$1~/RNA_unique_reads_mapped_to_genes/{print $2}' ~{qc_rna_statistics} > unique_reads_mapped_to_genes.txt
        awk -F "," '$1~/RNA_FRIG/{print $2}' ~{qc_rna_statistics} > qc_rna_frig.txt
        awk -F "," '$1~/RNA_duplicate_reads/{print $2}' ~{qc_rna_statistics} > duplicate_reads.txt
        awk -F "," '$1~/RNA_percent_duplicates/{print $2}' ~{qc_rna_statistics} > percent_duplicates.txt
        awk -F "," '$1~/RNA_percent_mitochondrial/{print $2}' ~{qc_rna_statistics} > percent_mitochondrial.txt
        
        # write statistics into metrics file to be used in HTML report
        echo "RNA_input_reads,$input_reads" > ~{prefix}_rna_qc_metrics.csv
        echo "RNA_aligned_reads,$aligned_reads" >> ~{prefix}_rna_qc_metrics.csv
        echo "RNA_uniquely_aligned_reads,$aligned_uniquely" >> ~{prefix}_rna_qc_metrics.csv
        echo "RNA_multimapped_reads,$aligned_multimap" >> ~{prefix}_rna_qc_metrics.csv
        echo "RNA_unaligned_reads,$unaligned_reads" >> ~{prefix}_rna_qc_metrics.csv

        cat ~{qc_rna_statistics} >> ~{prefix}_rna_qc_metrics.csv
    >>>

    output {
        File rna_qc_metrics = "${prefix}_rna_qc_metrics.csv"

        # STARsolo alignment log statistics
        Int rna_input_reads = read_int("input_reads.txt")
        Int rna_aligned_reads = read_int("aligned_reads.txt")
        Int rna_aligned_uniquely = read_int("aligned_uniquely.txt")
        Int rna_aligned_multimap = read_int("aligned_multimap.txt")
        Int rna_unaligned_reads = read_int("unaligned_reads.txt")

        # STARsolo barcode statistics
        Int rna_homopolymer_umis = read_int("homopolymer_umis.txt")
        Int rna_nonmatch_barcodes = read_int("nonmatch_barcodes.txt")
        Int rna_exact_match_barcodes = read_int("exact_match_barcodes.txt")
        Int rna_mismatch_barcodes = read_int("mismatch_barcodes.txt")

        # STARsolo summary statistics
        Float rna_frac_valid_barcodes = read_float("frac_valid_barcodes.txt")
        Float rna_sequencing_saturation = read_float("sequencing_saturation.txt")
        Float rna_frac_q30_bases_in_cb_umi = read_float("frac_q30_bases_in_cb_umi.txt")
        Float rna_frac_q30_bases_in_read = read_float("frac_q30_bases_in_read.txt")
        Float rna_starsolo_frig = read_float("starsolo_frig.txt")
        Int rna_estimated_cells = read_int("estimated_cells.txt")
        Float rna_frac_unique_reads_in_cells = read_float("frac_unique_reads_in_cells.txt")
        Int rna_median_reads_per_cell = read_int("median_reads_per_cell.txt")
        Int rna_median_umis_per_cell = read_int("median_umis_per_cell.txt")
        Int rna_median_genes_per_cell = read_int("median_genes_per_cell.txt")
        Int rna_genes = read_int("genes.txt")

        # qc_rna statistics
        Int rna_unique_reads_mapped_to_genes = read_int("unique_reads_mapped_to_genes.txt")
        Float rna_qc_rna_frig = read_float("qc_rna_frig.txt")
        Int rna_duplicate_reads = read_int("duplicate_reads.txt")
        Float rna_percent_duplicates = read_float("percent_duplicates.txt")
        Float rna_percent_mitochondrial = read_float("percent_mitochondrial.txt")
    }

    runtime {
        docker: 'ubuntu:latest'
    }

    parameter_meta {
        alignment_log: {
            description: 'RNA alignment log file',
	        help: 'STARsolo log file from RNA alignment step.',
            example: 'SS-PKR-30-96-ENTIRE-PLATE.rna.align.hg38.Log.out'
        }

        barcode_statistics: {
            description: 'RNA alignment barcode statistics file',
            help: 'STARsolo barcode statistics file from RNA alignment step.',
            example: 'SS-PKR-30-96-ENTIRE-PLATE.Barcodes.stats'
        }

        summary_csv: {
            description: 'RNA alignment summary csv file',
            help: 'STARsolo summary csv file from RNA alignment step.',
            example: 'SS-PKR-30-96-ENTIRE-PLATE.Summary.csv'
        }

        qc_rna_statistics: {
            description: 'RNA QC statistics file',
            help: 'Statistics file outputted by qc_rna task',
            example: 'SS-PKR-30-96-ENTIRE-PLATE.qc.rna.hg38.qc.statistics.txt'
        }
    }
}