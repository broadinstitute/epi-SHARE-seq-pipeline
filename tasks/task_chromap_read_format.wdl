version 1.0

# TASK
# SHARE-chromap-read-format

task get_read_format {
    meta {
        version: 'v0.1'
        author: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: get Chromap read format'
    }

    input {
        File fastq_R2 
    }

    command <<<
        # SHARE R2 FASTQ format is read 2, 15bp linker, 8bp round 1 barcode, 30bp linker,
        # 8bp round 2 barcode, 30bp linker, 8bp round 3 barcode
        # Chromap indexing is 0-based and inclusive

        r2_end=$(gzip -dc ~{fastq_R2} | awk 'NR==2 {print length($0)-100}')
        
        bc1_start=$(( $r2_end + 16 ))
        bc1_end=$(( $bc1_start + 7 ))

        bc2_start=$(( $bc1_end + 31 ))
        bc2_end=$(( $bc2_start + 7 ))

        bc3_start=$(( $bc2_end + 31 ))
        bc3_end=$(( $bc3_start + 7 ))

        echo "r1:0:-1,r2:0:${r2_end},bc:${bc1_start}:${bc1_end},bc:${bc2_start}:${bc2_end},bc:${bc3_start}:${bc3_end}" > read_format.txt

        cat read_format.txt
    >>>

    output {
        String? read_format = read_string("read_format.txt")
    }

    runtime {
        docker: 'ubuntu:latest'
    }

    parameter_meta {
        fastq_R2: {
            description: 'R2 FASTQ file',
	        help: 'FASTQ file containing read 2 and cell barcode sequence',
            example: 'SS-PKR-100_R2.fastq.gz'
        }
    }
}
