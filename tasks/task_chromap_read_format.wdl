version 1.0

# TASK
# SHARE-chromap-read-format

task get_chromap_read_format {
    meta {
        version: 'v0.1'
        author: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: get Chromap read format'
    }

    input {
        String fastq_path
    }

    command <<<
        # SHARE R2 FASTQ format is read 2, 15bp linker, 8bp round 1 barcode, 30bp linker,
        # 8bp round 2 barcode, 30bp linker, 8bp round 3 barcode
        # Chromap indexing is 0-based and inclusive

        # Get read2 end position
        r2_end=$(gsutil cp ~{fastq_path} - | gzip -dc | head -n 2 | tail -n 1 | awk '{print length($0) - 100}')
        
        # Calculate barcode positions
        bc1_start=$(( $r2_end + 16 ))
        bc1_end=$(( $bc1_start + 7 ))

        bc2_start=$(( $bc1_end + 31 ))
        bc2_end=$(( $bc2_start + 7 ))

        bc3_start=$(( $bc2_end + 31 ))
        bc3_end=$(( $bc3_start + 7 ))

        echo "r1:0:-1,r2:0:${r2_end},bc:${bc1_start}:${bc1_end},bc:${bc2_start}:${bc2_end},bc:${bc3_start}:${bc3_end}" | tee -a read_format.txt
    >>>

    output {
        String? read_format = read_string("read_format.txt")
    }

    runtime {
        docker: 'google/cloud-sdk:latest'
    }

    parameter_meta {
        fastq_path: {
            description: 'Path to FASTQ file containing cell barcode sequence',
	        help: 'Path to FASTQ file containing cell barcode sequence',
            example: 'SS-PKR-100_R2.fastq.gz'
        }
    }
}
