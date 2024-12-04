version 1.0

# TASK
# SHARE-atac-chromap

task atac_align_chromap {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: align ATAC task using chromap'
    }

    input {
        # This task takes in input the preprocessed ATAC fastqs and align them to the genome.
        Array[File] fastq_R1
        Array[File] fastq_R2
        Array[File] fastq_barcode
        File reference_index_tar_gz
        File reference_fasta
        File? barcode_inclusion_list
        File? barcode_conversion_dict

        Boolean? trim_adapters = true
        Boolean? remove_pcr_duplicates = true
        Boolean? remove_pcr_duplicates_at_cell_level = true
        Boolean? remove_pcr_duplicates_at_bulk_level = false
        Boolean? Tn5_shift = false
        Boolean? low_mem = true
        Boolean? bed_output = true
        Int? max_insert_size = 2000
        Int? quality_filter = 0
        

        Int? multimappers = 4 # As per ENCODE pipeline
        Int? bc_error_threshold = 0
        Float? bc_probability_threshold = 0.9
        #TODO: This should come from a previous task parsing the seqspec.
        String? read_format

        String? subpool = "none"
        String genome_name # GRCh38, mm10
        String prefix = "test-sample"

        Int? cpus = 8
        Float? disk_factor = 1
        #TODO: With this setting it usually caps at 75%.
        Float? memory_factor = 0.15
        #TODO:We need to setup a docker registry.
        String? docker_image = "us.gcr.io/buenrostro-share-seq/task_chromap:dev"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fastq_R1, "G") + size(fastq_R2, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 24.0 + size(reference_fasta, "G") + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    # Define the output names
    String fragment_file = '${prefix}.atac.fragments.${genome_name}.tsv'
    String fragment_file_sorted_by_barcode = '${prefix}.atac.fragments.${genome_name}.sorted.by.barcode.tsv'
    String barcode_log = "${prefix}.atac.align.k${multimappers}.${genome_name}.barcode.summary.csv"
    String alignment_log = "${prefix}.atac.align.k${multimappers}.${genome_name}.log.txt"

    command <<<
        set -e

        bash $(which monitor_script.sh) 1>&2 &

        # Extracting index
        echo '------ Extracting indexing ------' 1>&2
        time tar -xzf ~{reference_index_tar_gz}

        if [[ '~{barcode_inclusion_list}' == *.gz ]]; then
            echo '------ Decompressing the barcode inclusion list ------' 1>&2
            gunzip -c ~{barcode_inclusion_list} > barcode_inclusion_list.txt
        else
            echo '------ No decompression needed for the barcode inclusion list ------' 1>&2
            cat ~{barcode_inclusion_list} > barcode_inclusion_list.txt
        fi

        chromap --version > chromap_version.txt 2>&1
        
        # [r1|r2|bc]:start:end:strand
        # --read-format bc:0:15,r1:16:-1
        # The start and end are inclusive and -1 means the end of the read. User may use multiple fields to specify non-consecutive segments, e.g. bc:0:15,bc:32:-1.
        # The strand is presented by '+' and '-' symbol, if '-' the barcode will be reverse-complemented after extraction
        echo '------ align chromap ------' 1>&2
        chromap -x chromap_index/index \
            ~{true='--trim-adapters ' false='' trim_adapters} \
            ~{true='--remove-pcr-duplicates ' false='' remove_pcr_duplicates} \
            ~{true='--remove-pcr-duplicates-at-cell-level ' false='' remove_pcr_duplicates_at_cell_level} \
            ~{true='--remove-pcr-duplicates-at-bulk-level ' false='' remove_pcr_duplicates_at_bulk_level} \
            ~{true='--Tn5-shift ' false='' Tn5_shift} \
            ~{true='--low-mem ' false='' low_mem} \
            ~{true='--BED ' false='' bed_output} \
            ~{if max_insert_size > 0 then "-l " + "~{max_insert_size}" else "" } \
            ~{"--bc-error-threshold " + bc_error_threshold} \
            ~{"--bc-probability-threshold " + bc_probability_threshold} \
            ~{"--read-format " + read_format} \
            ~{"--allocate-multi-mappings " + multimappers} \
            -r ~{reference_fasta} \
            ~{"-q " + quality_filter} \
            -t ~{cpus} \
            -1 ~{sep="," fastq_R1} \
            -2 ~{sep="," fastq_R2} \
            -b ~{sep="," fastq_barcode} \
            --barcode-whitelist barcode_inclusion_list.txt \
            ~{"--barcode-translate " + barcode_conversion_dict} \
            -o out.fragments.tmp.tsv \
            --summary ~{barcode_log} > ~{alignment_log} 2>&1
        
        if [[ '~{subpool}' != "none" ]]; then
            echo '------  Add subpool to barcode name ------' 1>&2
            awk -v OFS="\t" -v subpool=~{subpool} '{$4=$4"_"subpool; print $0}' out.fragments.tmp.tsv > temp
            mv temp out.fragments.tmp.tsv
            awk -v FS="," -v OFS="," -v subpool=~{subpool} 'NR==1{print $0;next}{$1=$1"_"subpool; print $0}' ~{barcode_log} > temp
            mv temp ~{barcode_log}
        fi
        echo '------ Compute percentage of duplicates ------' 1>&2
        # Compute percentage of duplicates
        awk '{total+=$5}END{printf "%.1f\n", (total-NR)/total*100}' out.fragments.tmp.tsv > duplicates_percentage.txt
        cut -f4 out.fragments.tmp.tsv | sort -u | wc -l > unique_barcodes_unfiltered.txt

        echo '------ Sort fragment file by barcode ------' 1>&2
        # Sort fragments by name
        sort --parallel=~{cpus} -k4,4 out.fragments.tmp.tsv > ~{fragment_file_sorted_by_barcode}

        echo '------ Compress and index fragment file ------' 1>&2
        bgzip -c out.fragments.tmp.tsv > ~{fragment_file}.gz
        tabix --zero-based --preset bed ~{fragment_file}.gz

        grep "Number of reads:" ~{alignment_log} | tr -d '.' | awk '{print $NF}' > reads_count.txt
        grep "Number of mapped reads" ~{alignment_log} | tr -d '.' |awk '{print $NF}' > mapped_reads.txt
        grep "Number of uniquely mapped reads" ~{alignment_log} | tr -d '.' | awk '{print $NF}' > unique_reads.txt
        grep "Number of reads have multi-mappings:" ~{alignment_log} | tr -d '.' | awk '{print $NF}' > multi_mapping_reads.txt
        grep "Number of corrected barcodes" ~{alignment_log} | tr -d '.' | awk '{print $NF}' > corrected_barcodes.txt
        grep "# uni-mappings" ~{alignment_log} | tr -d '.' | tr -d ',' | awk '{print $3}' > uni_mappings_fragments.txt
        grep "# multi-mappings" ~{alignment_log} | tr -d '.' | tr -d ',' | awk '{print $6}' > multi_mappings_fragments.txt
        grep "Number of output mappings" ~{alignment_log} | awk '{print $7}' > final_number_of_fragments.txt
  

    >>>

    output {
        File atac_fragment_file = "~{fragment_file}.gz"
        File atac_fragment_file_index = "~{fragment_file}.gz.tbi"
        File atac_fragment_file_sorted_by_barcode = fragment_file_sorted_by_barcode
        File atac_align_barcode_statistics = barcode_log
        File atac_alignment_log = alignment_log
        Float atac_pcr_duplicates_percentage = read_float("duplicates_percentage.txt")
        Float atac_reads_count = read_float("reads_count.txt")
        Int atac_mapped_reads = read_int("mapped_reads.txt")
        Int atac_unique_reads = read_int("unique_reads.txt")
        Int atac_multi_mapping_reads = read_int("multi_mapping_reads.txt")
        Int atac_corrected_barcodes = read_int("corrected_barcodes.txt")
        Int atac_unique_mappings_fragments = read_int("uni_mappings_fragments.txt")
        Int atac_multi_mappings_fragments = read_int("multi_mappings_fragments.txt")
        Int atac_final_number_of_fragments = read_int("final_number_of_fragments.txt")
        Int atac_unique_barcodes_unfiltered = read_int("unique_barcodes_unfiltered.txt")
        String atac_chromap_version = read_string("chromap_version.txt")

    }


    runtime {
        cpu: cpus
        docker: "${docker_image}"
        singularity: "docker://${docker_image}"
        disks: "local-disk ${disk_gb} ${disk_type}"
        memory: "${mem_gb} GB"
    }

    parameter_meta {
        fastq_R1: {
                description: 'Read1 fastq.',
                help: 'Processed fastq for read1.',
                example: 'input.atac.R1.fq.gz',
            }
        fastq_R2: {
                description: 'Read2 fastq.',
                help: 'Processed fastq for read2.',
                example: 'input.atac.R2.fq.gz',
            }
        fastq_barcode: {
                description: 'Barcode fastq.',
                help: 'Processed fastq for barcode.',
                example: 'input.atac.barcode.fq.gz',
            }
        reference_fasta: {
                description: 'Reference fasta.',
                help: 'Reference fasta file.',
                example: 'reference.fasta',
            }
        barcode_inclusion_list: {
                description: 'Barcode inclusion list.',
                help: 'Barcode inclusion list.',
                example: 'barcode_inclusion_list.txt',
            }
        barcode_conversion_dict: {
                description: 'Barcode conversion dict.',
                help: 'Barcode conversion dict.',
                example: 'barcode_conversion_dict.txt',
            }
        trim_adapters: {
                description: 'Trim adapters.',
                help: 'Trim adapters.',
                example: 'true',
            }
        remove_pcr_duplicates: {
                description: 'Remove PCR duplicates.',
                help: 'Remove PCR duplicates.',
                example: 'true',
            }
        remove_pcr_duplicates_at_cell_level: {
                description: 'Remove PCR duplicates at cell level.',
                help: 'Remove PCR duplicates at cell level.',
                example: 'false',
            }
        remove_pcr_duplicates_at_bulk_level: {
                description: 'Remove PCR duplicates at bulk level.',
                help: 'Remove PCR duplicates at bulk level.',
                example: 'true',
            }
        Tn5_shift: {
                description: 'Tn5 shift.',
                help: 'Tn5 shift.',
                example: 'false',
            }
        low_mem: {
                description: 'Low mem.',
                help: 'Low mem.',
                example: 'true',
            }
        bed_output: {
                description: 'Bed output.',
                help: 'Bed output.',
                example: 'true',
            }
        max_insert_size: {
                description: 'Max insert size.',
                help: 'Max insert size.',
                example: '2000',
            }
        quality_filter: {
                description: 'Quality filter.',
                help: 'Quality filter.',
                example: '0',
            }
        multimappers: {
                description: 'Multimappers.',
                help: 'Multimappers.',
                example: '4',
            }
        bc_error_threshold: {
                description: 'BC error threshold.',
                help: 'BC error threshold.',
                example: '1',
            }
        bc_probability_threshold: {
                description: 'BC probability threshold.',
                help: 'BC probability threshold.',
                example: '0.9',
            }
        read_format: {
                description: 'Read format.',
                help: 'Read format.',
                example: 'bc:0:15,r1:16:-1',
            }
        subpool: {
                description: 'Subpool.',
                help: 'Subpool.',
                example: 'none',
            }
        genome_name: {
                description: 'Genome name.',
                help: 'Genome name.',
                example: 'GRCh38',
            }
        prefix: {
                description: 'Prefix.',
                help: 'Prefix.',
                example: 'test-sample',
            }
        cpus: {
                description: 'CPUs.',
                help: 'CPUs.',
                example: '8',
            }
        disk_factor: {
                description: 'Disk factor.',
                help: 'Disk factor.',
                example: '1',
            }
        memory_factor: {
                description: 'Memory factor.',
                help: 'Memory factor.',
                example: '0.15',
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image.',
                example: 'us.gcr.io/buenrostro-share-seq/task_chromap:dev',
            }
    }
}
