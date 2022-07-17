version 1.0

# TASK
# SHARE-atac-log
# Gather information from log files 


task log_atac {
    meta {
        version: 'v0.1'
        author: 'Neva C. Durand (neva@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: log atac task'
    }

    input {
        # This function takes as input the necessary log files and extracts
        # the quality metrics
        File alignment_log
        File dups_log
    }

    command <<<
        total_reads=$(awk 'NR==1{print $1}' ~{alignment_log})
        echo $total_reads > total_reads.txt
        aligned_uniquely=$(awk 'NR==4{print $1}' ~{alignment_log})
        echo $aligned_uniquely > aligned_uniquely.txt
        echo $(($total_reads - $aligned_uniquely)) > unaligned.txt
        awk 'NR==2{print $3}' ~{dups_log} > feature_reads.txt
        awk 'NR==2{print $6}' ~{dups_log} > duplicate_reads.txt
    >>>
    output {
        Int atac_total_reads = read_int("total_reads.txt")
        Int atac_aligned_uniquely = read_int("aligned_uniquely.txt")
        Int atac_unaligned = read_int("unaligned.txt")
        Int atac_feature_reads = read_int("feature_reads.txt")
        Int atac_duplicate_reads = read_int("duplicate_reads.txt")
    }

    runtime {
        docker: 'ubuntu:latest'
    }
    parameter_meta {
        alignment_log: {
            description: 'ATAC alignment log file',
	    help: 'Log file from ATAC alignment step.',
            example: 'SS-PKR-30-96-ENTIRE-PLATE.atac.align.hg38.Log.out'
        }
        dups_log: {
            description: 'ATAC dups log file',
            help: 'Log file from ATAC rmdups step.',
            example: 'SS-PKR-12.atac.counts.mm10.filtered.cs.log'
        }
    }
}
