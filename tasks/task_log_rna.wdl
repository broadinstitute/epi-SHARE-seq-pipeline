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
        File dups_log
        String? prefix = "sample"
    }

    command <<<
        total_reads=$(awk -F"|" '$1~/input reads/{print $2}' ~{alignment_log})
        aligned_uniquely=$(awk -F"|" '$1~/Uniquely mapped reads number/{print $2}' ~{alignment_log})
        aligned_multimap=$(awk -F"|" '$1~/Number of reads mapped to multiple loci/{print $2}' ~{alignment_log})
        aligned=$(($aligned_uniquely + $aligned_multimap))
        unaligned=$(($total_reads - $aligned))
        
        echo "RNA_input_reads,$total_reads" > ~{prefix}_rna_qc_metrics.csv
        echo "RNA_aligned_reads,$aligned" >> ~{prefix}_rna_qc_metrics.csv
        echo "RNA_uniquely_aligned_reads,$aligned_uniquely" >> ~{prefix}_rna_qc_metrics.csv
        echo "RNA_multimapped_reads,$aligned_multimap" >> ~{prefix}_rna_qc_metrics.csv
        echo "RNA_unaligned_reads,$unaligned" >> ~{prefix}_rna_qc_metrics.csv

        cat ~{dups_log} >> ~{prefix}_rna_qc_metrics.csv
    >>>

    output {
        File rna_qc_metrics = "~{prefix}_rna_qc_metrics.csv"
    }

    runtime {
        docker: 'ubuntu:latest'
    }
    parameter_meta {
        alignment_log: {
            description: 'RNA alignment log file',
	        help: 'Log file from RNA alignment step.',
            example: 'SS-PKR-30-96-ENTIRE-PLATE.rna.align.hg38.Log.out'
        }

        dups_log: {
            description: 'Group UMI dups log file',
            help: 'Log file from group UMI task',
            example: 'SS-PKR-30-96-ENTIRE-PLATE.rm_dup_barcode.log.txt'
        }
    }
}
