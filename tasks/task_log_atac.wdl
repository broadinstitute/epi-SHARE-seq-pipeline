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
        File barcode_log
    }
    
    #awk 'NR>1{sum += $2}END{print sum/2}' ~{dups_log} > feature_reads.txt
    #awk 'NR>1{sum += $3}END{print sum/2}' ~{dups_log} > duplicate_reads.txt
    #awk 'NR>1{unique+= $2; dups+=$3}END{printf "%5.1f%", 100*dups/(unique+dups)}' ~{dups_log} > pct_duplicate_reads.txt

    command <<<
        awk -v FS="," 'NR>1{unique+= $2; dups+=$3}END{printf "%5.1f%", 100*dups/(unique+dups)}' ~{barcode_log} > pct_duplicate_reads.txt
    >>>
    output {
        Float? atac_pct_dup = read_float("pct_duplicate_reads.txt")
    }

    runtime {
        docker: 'ubuntu:latest'
    }
    parameter_meta {
        barcode_log: {
            description: 'ATAC alignment log file',
        help: 'Log file from ATAC alignment step.',
            example: 'SS-PKR-30-96-ENTIRE-PLATE.atac.align.hg38.Log.out'
        }
    }
}
