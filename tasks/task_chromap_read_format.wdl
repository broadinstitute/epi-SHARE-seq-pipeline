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
        File fastq
        
        Float? disk_factor = 1
        Float? memory_factor = 0.15
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fastq, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 1.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(1.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    command <<<
        # SHARE R2 FASTQ format is read 2, 15bp linker, 8bp round 1 barcode, 30bp linker,
        # 8bp round 2 barcode, 30bp linker, 8bp round 3 barcode
        # Chromap indexing is 0-based and inclusive

        r2_end=$(gzip -dc ~{fastq} | awk 'NR==2 {print length($0)-100}')
        
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
        disks: "local-disk ${disk_gb} ${disk_type}"
        memory: "${mem_gb} GB"
        docker: 'ubuntu:latest'
    }

    parameter_meta {
        fastq: {
            description: 'FASTQ file containing cell barcode sequence',
	        help: 'FASTQ file containing cell barcode sequence',
            example: 'SS-PKR-100_R2.fastq.gz'
        }
    }
}
