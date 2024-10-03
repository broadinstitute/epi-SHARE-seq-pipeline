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
        File chrom_sizes
        File? barcode_inclusion_list
        File? barcode_conversion_dict

        Boolean? trim_adapters = true
        Boolean? remove_pcr_duplicates = true
        Boolean? remove_pcr_duplicates_at_cell_level = false
        Boolean? remove_pcr_duplicates_at_bulk_level = true
        Boolean? Tn5_shift = false
        Boolean? low_mem = true
        Boolean? bed_output = true
        Int? max_insert_size = 2000
        Int? quality_filter = 0
        

        Int? multimappers = 4 # As per ENCODE pipeline
        Int? bc_error_threshold = 1
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
    String fragments = '${prefix}.atac.fragments.${genome_name}.tsv'
    String barcode_log = "${prefix}.atac.align.k${multimappers}.${genome_name}.barcode.summary.csv"
    String alignment_log = "${prefix}.atac.align.k${multimappers}.${genome_name}.log.txt"

    String monitor_log = "atac_align_monitor.log.txt"

    command <<<
        set -e

        bash $(which monitor_script.sh) 1>&2 &

        # Create index
        mkdir chromap_index
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
        
        # [r1|r2|bc]:start:end:strand
        # --read-format bc:0:15,r1:16:-1
        # The start and end are inclusive and -1 means the end of the read. User may use multiple fields to specify non-consecutive segments, e.g. bc:0:15,bc:32:-1.
        # The strand is presented by '+' and '-' symbol, if '-' the barcode will be reverse-complemented after extraction
        echo '------ align chromap ------' 1>&2
        time chromap -x chromap_index/index \
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
                ~{"--drop-repetitive-reads " + multimappers} \
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

        # Writing non corrected fragments
        awk -v OFS="\t" -v maxinsert=~{max_insert_size} '$3-$2 <= maxinsert' out.fragments.tmp.tsv | bgzip -c > ~{fragments}.gz
        tabix --zero-based --preset bed ~{fragments}.gz

    >>>

    output {
        File atac_fragments = "~{fragments}.gz"
        File atac_fragments_index = "~{fragments}.gz.tbi"
        File atac_align_barcode_statistics = barcode_log
        File atac_alignment_log = alignment_log
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
        chrom_sizes: {
                description: 'Chrom sizes.',
                help: 'Chrom sizes file.',
                example: 'chrom.sizes',
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
