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
    }

    command <<<
        total_reads=$(awk -F"|" '$1~/input reads/{print $2}' ~{alignment_log})
        echo $total_reads > total_reads.txt
        aligned_uniquely=$(awk -F"|" '$1~/Uniquely mapped reads number/{print $2}' ~{alignment_log})
        echo $aligned_uniquely > aligned_uniquely.txt
        aligned_multimap=$(awk -F"|" '$1~/Number of reads mapped to multiple loci/{print $2}' ~{alignment_log})
        echo $aligned_multimap > aligned_multimap.txt
        echo $(($total_reads - $aligned_uniquely - $aligned_multimap)) > unaligned.txt
        awk -F":" '$1~/total reads/{print $2}' ~{dups_log} > feature_reads.txt
        awk -F":" '$1~/duplicate reads/{print $2}' ~{dups_log} > duplicate_reads.txt
    >>>
    output {
        Int rna_total_reads = read_int("total_reads.txt")
        Int rna_aligned_uniquely = read_int("aligned_uniquely.txt")
        Int rna_aligned_multimap = read_int("aligned_multimap.txt")
        Int rna_unaligned = read_int("unaligned.txt")	    
        Int rna_feature_reads = read_int("feature_reads.txt")
        Int rna_duplicate_reads = read_int("duplicate_reads.txt")
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
