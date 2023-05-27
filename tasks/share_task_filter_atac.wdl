version 1.0

# TASK
# SHARE-atac-filter




task share_atac_filter {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: ATAC filter task'
    }

    input {
        # This task takes in input the aligned bam file and rmeove the low quality reads, the extra chromosomes, marks
        # the duplicats, and convert to a bedpe file.
        Int? shift_plus = 4
        Int? shift_minus = -4
        Int? mapq_threshold = 30
        Int? multimappers = 1
        Int? minimum_fragments_cutoff = 10
        File? bam
        File? bam_index
        File? barcode_conversion_dict
        String? barcode_tag = "CB"
        String? barcode_tag_fragments = "CB"
        String genome_name
        String? prefix = "sample"
        ## Runtime
        Int? cpus = 16
        Float? disk_factor = 8.0
        Float? memory_factor = 0.15
        String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_filter_atac"
    }


    # Determine the size of the input
    Float input_file_size_gb = size(bam, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 32.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    # Determining memory for samtools.
    Int samtools_memory_gb = floor(0.9 * mem_gb) # Samtools has overheads so reducing the memory to 80% of the total.

    # Number of threads to beable to use 4GB of memory per thread seems to be the fastest way
    Int samtools_threads_ = floor(samtools_memory_gb / 4)
    Int samtools_threads =  if samtools_threads_ == 0 then 1 else samtools_threads_

    Int sambamba_threads = floor(cpus/2)

    # Now that we know how many threads we can use to assure 4GB of memory per thread
    # we assign any remaining memory to the threads.
    Int samtools_memory_per_thread_ = floor(samtools_memory_gb * 1024 / samtools_threads) # Computing the memory per thread for samtools in MB.
    Int samtools_memory_per_thread = if samtools_memory_per_thread_ < 768 then 768 else samtools_memory_per_thread_

    # Memory for picard
    Float picard_java_heap_factor = 0.8
    Int picard_java_memory = round(mem_gb * picard_java_heap_factor)

    #String filtering_params = if multimappers == 0 then "-q ${mapq_threshold} -F 1804" else "-F 524"


    String non_mito_bam = "${prefix}.atac.align.k${multimappers}.${genome_name}.nonmito.sorted.bam"

    String tmp_filtered_bam = '${prefix}.filtered.tmp.bam'
    String tmp_fixmate_bam = '${prefix}.filtered.fixmate.tmp.bam'


    String picard_mark_duplicates_metrics = '${prefix}.picard.marksduplicates.metrics.txt'
    String picard_mark_duplicates_log = '${prefix}.picard.marksduplicates.log'

    String final_bam_wdup = '${prefix}.atac.filter.fixmate.wdup.k${multimappers}.${genome_name}.bam'
    String final_bam_wdup_tmp = '${prefix}.atac.filter.fixmate.wdup.k${multimappers}.${genome_name}.tmp.bam'
    String final_bam_wdup_index = '${prefix}.atac.filter.fixmate.wdup.k${multimappers}.${genome_name}.bam.bai'

    String queryname_final_bam = '${prefix}.atac.filter.fixmate.dedup.k${multimappers}.${genome_name}.queryname.final.bam'

    String final_bam = '${prefix}.atac.filter.fixmate.dedup.k${multimappers}.${genome_name}.bam'
    String final_bam_index = '${prefix}.atac.filter.fixmate.dedup.k${multimappers}.${genome_name}.bam.bai'

    String fragments = '${prefix}.atac.filter.fragments.${genome_name}.tsv'

    String monitor_log = "atac_filter_monitor.log"

    command<<<
        set -e

        # I am not writing to a file anymore because Google keeps track of it automatically.
        bash $(which monitor_script.sh) 1>&2 &


        # I need to do this because the bam and bai need to be in the same folder but WDL doesn't allow you to
        # co-localize them in the same path.
        ln -s ~{bam} in.bam
        ln -s ~{bam_index} in.bam.bai

        # "{prefix}.mito.bulk-metrics.tsv"
        # "{prefix}.mito.bc-metrics.tsv"

        echo '------ START: Filtering out mito reads ------' 1>&2
        # The script removes the mithocondrial reads and creates two log file with bulk and barcode statistics.
        time python3 $(which filter_mito_reads.py) -o ~{non_mito_bam} -p ~{cpus} --cutoff ~{minimum_fragments_cutoff} --prefix ~{prefix} --bc_tag ~{barcode_tag_fragments} in.bam

        sambamba index -t ~{cpus} ~{non_mito_bam}

        # Keep only assembled chromosomes
        chrs=$(samtools view -H ~{non_mito_bam}| \
            grep chr | \
            cut -f2 | \
            sed 's/SN://g' | \
            awk '{if(length($0)<6)print}')

        # =============================
        # Remove  unmapped, mate unmapped
        # not primary alignment, reads failing platform
        # Only keep properly paired reads
        # Obtain name sorted BAM file
        # =============================
        echo '------ START: Filter bam -F 524 -f 2 and sort ------' 1>&2
        time sambamba view -h -t ~{sambamba_threads} --num-filter 2/524 -f bam ~{non_mito_bam} $(echo ${chrs}) | \
        sambamba sort -t ~{cpus} -m ~{samtools_memory_gb}G -n -o ~{tmp_filtered_bam} /dev/stdin

        # Assign multimappers if necessary
        if [ ~{multimappers} -le 1 ]; then
            echo '------ START: Fixmate step ------' 1>&2
            time sambamba view -t ~{cpus} -h -f sam ~{tmp_filtered_bam}  | samtools fixmate -@ ~{cpus} -r /dev/stdin ~{tmp_fixmate_bam}
        else
            echo '------ START: Assinging multimappers ------' 1>&2
            time sambamba view -t ~{cpus} -h -f sam ~{tmp_filtered_bam} | \
            python3 $(which assign_multimappers.py) -k ~{multimappers} --paired-end | sambamba view -t ~{cpus} -f sam /dev/stdin | samtools fixmate -@ ~{cpus} -r /dev/stdin ~{tmp_fixmate_bam}
        fi

        # Cleaning up bams we don't need anymore
        rm ~{tmp_filtered_bam}
        # rm ~{non_mito_bam}

        # =============
        # Mark duplicates
        # =============
        echo '------ START: Mark Duplicates ------' 1>&2
        time java -Dsamjdk.compression_level=5 -Xms~{picard_java_memory}G -Xmx~{picard_java_memory}G -jar $(which picard.jar) MarkDuplicates \
        --INPUT ~{tmp_fixmate_bam} --OUTPUT ~{final_bam_wdup_tmp} \
        --METRICS_FILE ~{picard_mark_duplicates_metrics} \
        --VALIDATION_STRINGENCY LENIENT \
        --ASSUME_SORT_ORDER queryname \
        --REMOVE_DUPLICATES false \
        --TAG_DUPLICATE_SET_MEMBERS true \
        --READ_NAME_REGEX NULL \
        --MAX_OPTICAL_DUPLICATE_SET_SIZE -1 \
        --BARCODE_TAG ~{barcode_tag} 2> ~{picard_mark_duplicates_log}

        # Create the final bam removing the duplicates
        echo '------ START: Remove duplicates ------' 1>&2
        time sambamba view -t ~{cpus} -h --num-filter 2/1804 -f bam -o ~{queryname_final_bam} ~{final_bam_wdup_tmp}

        echo '------ START: Sort bam with duplicates by coordinates ------' 1>&2
        time sambamba sort -t ~{cpus} -m ~{samtools_memory_gb}G  -o ~{final_bam_wdup} ~{final_bam_wdup_tmp}
        sambamba index -t ~{cpus} ~{final_bam_wdup}

        echo '------ START: Sort bam without duplicates by coordinates ------' 1>&2
        time sambamba sort -t ~{cpus} -m ~{samtools_memory_gb}G -o ~{final_bam} ~{queryname_final_bam}
        sambamba index -t ~{cpus} ~{final_bam}

        rm ~{tmp_fixmate_bam}

        echo '------ START: Create fragment file ------' 1>&2
        time python3 $(which bam_to_fragments.py) --shift_plus ~{shift_plus} --shift_minus ~{shift_minus} --bc_tag ~{barcode_tag_fragments} -o ~{fragments} ~{final_bam}

        # Change the ATAC 10x barcodes so they can match RNA
        if [ ~{if defined(barcode_conversion_dict) then "true" else "false"} == "true" ];then
            cp ~{fragments} tmp_fragments
            echo '------ START: Convert barcode names for 10x ------' 1>&2
            time awk -F ",|\t" -v OFS="\t" 'FNR==NR{map[$1]=$2; next}{print $1,$2,$3,map[$4],$5}' ~{barcode_conversion_dict} tmp_fragments > ~{fragments}
        fi

        echo '------ START: Compress fragments ------' 1>&2
        time bgzip -c ~{fragments} > ~{fragments}.gz

        # ~{prefix}.fragments.tsv.gz.tbi
        tabix --zero-based --preset bed ~{fragments}.gz
    >>>

    output {
        File? atac_filter_alignment_dedup = final_bam
        File? atac_filter_alignment_dedup_index = final_bam_index
        File? atac_filter_alignment_wdup = final_bam_wdup
        File? atac_filter_alignment_wdup_index = final_bam_wdup_index
        File? atac_filter_alignment_dedup_queryname = queryname_final_bam
        File atac_filter_fragments = "~{fragments}.gz"
        File? atac_filter_fragments_index = "~{fragments}.gz.tbi"
        File? atac_filter_monitor_log = monitor_log
        File? atac_filter_picard_duplicates_metrics = picard_mark_duplicates_metrics
        File? atac_filter_picard_duplicates_log = picard_mark_duplicates_log
        File? atac_filter_mito_metrics_bulk = "~{prefix}.mito.bulk-metrics.tsv"
        File? atac_filter_mito_metrics_barcode = "~{prefix}.mito.bc-metrics.tsv"
    }

    runtime {
        cpu: cpus
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker: "${docker_image}"
        maxRetries: 1
        memory: "${mem_gb} GB"
        memory_retry_multiplier: 2
    }

    parameter_meta {
        bam: {
                description: 'bam file',
                help: 'Aligned reads in bam format',
                example: 'aligned.hg38.bam'
            }
        bam_index: {
            description: 'bai file',
            help: 'Index for the aligned reads in bam format',
            example: 'aligned.hg38.bam.bai'
            }
        barcode_tag: {
            description: 'tag containing the barcode',
            help: 'Which tag inside the bma file contains the cell barcode.',
            examples: ['CB','XC']
            }
        shift_plus: {
                description: 'Integer indicating by how any basepair shift the star position on the positive strand.',
                help: 'Shift the start position of the fragment by the amount of basepairs indicate.',
                example: '4'
            }
        shift_minus: {
                description: 'Integer indicating by how any basepair shift the star position on the negative strand.',
                help: 'Shift the end position of the fragment by the amount of basepairs indicate.',
                example: '-4'
            }
        mapq_threshold: {
            description: 'MAPQ value',
            help: 'Minimum value for MAPQ allowed',
            example: '30'
            }
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus.',
                examples: '4'
            }
        multimappers: {
                    description: 'Specifiy the numbers of multimappers allowed.',
                    help: 'Number of multimppares that have been passed to bowtie2 during alignment',
                    example: [5]
            }
        disk_factor: {
                description: 'Multiplication factor to determine disk required for task filter.',
                help: 'This factor will be multiplied to the size of bams to determine required disk of instance (GCP/AWS) or job (HPCs).',
                default: 8.0
            }
        memory_factor: {
                description: 'Multiplication factor to determine memory required for task filter.',
                help: 'This factor will be multiplied to the size of bams to determine required memory of instance (GCP/AWS) or job (HPCs).',
                default: 0.15
            }
        genome_name: {
                description: 'Reference name.',
                help: 'The name of the reference genome used. This is appended to the output file name.',
                examples: ['GRCh38', 'mm10']
            }
        prefix: {
                description: 'Prefix for output files.',
                help: 'Prefix that will be used to name the output files',
                examples: 'my-experiment'
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for the filtering step.',
                example: ["us.gcr.io/buenrostro-share-seq/share_task_filter"]
            }
    }


}
