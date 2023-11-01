version 1.0

# TASK
# Create a file with QC metrics
# Gather information from log files


task log_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard IGVF pipeline: log atac task'
    }

    input {
        # This function takes as input the necessary log files and extracts
        # the quality metrics
        File? alignment_log
        File? barcode_log
        String? prefix = "sample"
    }

    command <<<
        # Formatting the output of chromap and extracting statistics
        grep "Number of" ~{alignment_log} | grep -v threads| tr -d '.' | sed 's/ /_/g' | sed 's/:_/,/g'> ~{prefix}_qc_metrics.csv
        grep "#" ~{alignment_log}  | sed 's/, /\n/g' | tr -d '# ' | sed 's/:/,/g' | tr -d '.' >> ~{prefix}_qc_metrics.csv
        # Compute the percentage of duplicates from the barcode log file.
        awk -v FS="," 'NR>1{unique+= $2; dups+=$3}END{printf "percentage_duplicates,%5.1f", 100*dups/(unique+dups)}' ~{barcode_log} >> ~{prefix}_qc_metrics.csv
    >>>
    output {
        File atac_statistics_csv = "~{prefix}_qc_metrics.csv"
    }

    runtime {
        docker: 'ubuntu:latest'
    }
    parameter_meta {
        alignment_log: {
            description: 'ATAC alignment log file',
            help: 'Log file from ATAC alignment step.'
        }
        barcode_log: {
            description: 'ATAC dups log file',
            help: 'Barcode log file from ATAC alignment step.'
        }
    }
}