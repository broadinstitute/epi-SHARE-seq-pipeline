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
        Int? cpus = 16
        Int? mapq_threshold = 30
        File? bam
        File? bam_index
        Float? disk_factor = 8.0
        Float? memory_factor = 0.15
        Int? multimappers=0
        String? barcode_tag = "CB"
        String docker_image = "polumechanos/share_atac_filter"
        String genome_name
        String? prefix = "sample"

    }


    # Determine the size of the input
    Float input_file_size_gb = size(bam, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 6.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    # Determining memory for samtools.
    Float samtools_memory_gb = 0.8 * mem_gb # Samtools has overheads so reducing the memory to 80% of the total.

    # Number of threads to beable to use 4GB of memory per thread seems to be the fastest way
    Int samtools_threads_ = floor(samtools_memory_gb / 4)
    Int samtools_threads =  if samtools_threads_ == 0 then 1 else samtools_threads_

    # Now that we know how many threads we can use to assure 4GB of memory per thread
    # we assign any remaining memory to the threads.
    Int samtools_memory_per_thread_ = floor(samtools_memory_gb * 1024 / samtools_threads) # Computing the memory per thread for samtools in MB.
    Int samtools_memory_per_thread = if samtools_memory_per_thread_ < 768 then 768 else samtools_memory_per_thread_

    # Memory for picard
    Float picard_java_heap_factor = 0.9
    Int picard_java_memory = round(mem_gb * picard_java_heap_factor)

    String filtering_params = if multimappers == 0 then "-q ${mapq_threshold} -F 1804" else "-F 524"


    String tmp_filtered_bam = '${prefix}.filtered.tmp.bam'
    String tmp_fixmate_bam = '${prefix}.filtered.fixmate.tmp.bam'
    String filtered_bam = '${prefix}.filtered.fixmate.bam'
    String tmp_dup_bam_sorted = '${prefix}.filtered.fixmate.dupmarked.tmp.sorted.bam'
    String picard_mark_duplicates_metrics = '${prefix}.picard.marksduplicates.metrics.txt'
    String picard_mark_duplicates_log = '${prefix}.picard.marksduplicates.log'

    String final_bam_wdup = '${prefix}.atac.filter.fixmate.wdup.k${multimappers}.${genome_name}.bam'
    String final_bam_wdup_index = '${prefix}.atac.filter.fixmate.wdup.k${multimappers}.${genome_name}.bam.bai'
    String final_bam = '${prefix}.atac.filter.fixmate.dedup.k${multimappers}.${genome_name}.bam'
    String final_bam_index = '${prefix}.atac.filter.fixmate.dedup.k${multimappers}.${genome_name}.bam.bai'
    String fragments = '${prefix}.atac.filter.fragments.${genome_name}.tsv.gz'

    String monitor_log = "atac_filter_monitor.log"

    command<<<
        set -e

        bash $(which monitor_script.sh) > ~{monitor_log} 2>&1 &

        # TODO: Fraction of mito reads

        # I need to do this because the bam and bai need to be in the same folder but WDL doesn't allow you to
        # co-localize them in the same path.
        ln -s ~{bam} in.bam
        ln -s ~{bam_index} in.bam.bai

        # Keep only assembled chromosomes
        chrs=$(samtools view -H in.bam | \
            grep chr | \
            cut -f2 | \
            sed 's/SN://g' | \
            grep -v chrM | \
            awk '{if(length($0)<6)print}')

        # =============================
        # Remove  unmapped, mate unmapped
        # not primary alignment, reads failing platform
        # Only keep properly paired reads
        # Obtain name sorted BAM file
        # =============================
        samtools view -@ ~{samtools_threads} ~{filtering_params} -f 2 -u in.bam $(echo ${chrs}) | \
        samtools sort -@ ~{samtools_threads} -m ~{samtools_memory_per_thread}M -n /dev/stdin -o ~{tmp_filtered_bam}

        # Assign multimappers if necessary
        if [ ~{multimappers} -eq 0 ]; then
            samtools fixmate -@ ~{samtools_threads} -r ~{tmp_filtered_bam} ~{tmp_fixmate_bam}
        else
            samtools view -h ~{tmp_filtered_bam} | \
            python3 $(which assign_multimappers.py) -k ~{multimappers} --paired-end | samtools fixmate -r /dev/stdin ~{tmp_fixmate_bam}
        fi

        # Remove orphan reads (pair was removed)
        # and read pairs mapping to different chromosomes
        # Obtain position sorted BAM
        samtools view -F 1804 -f 2 -u ~{tmp_fixmate_bam} | \
        samtools sort -@ ~{samtools_threads} -m ~{samtools_memory_per_thread}M -n /dev/stdin -o ~{filtered_bam}

        # Cleaning up bams we don't need anymore
        rm ~{tmp_fixmate_bam}
        rm ~{tmp_filtered_bam}

        # =============
        # Mark duplicates
        # =============

        java -Xmx~{picard_java_memory}G -jar $(which picard.jar) MarkDuplicates \
        --INPUT ~{filtered_bam} --OUTPUT ~{final_bam_wdup} \
        --METRICS_FILE ~{picard_mark_duplicates_metrics} \
        --VALIDATION_STRINGENCY LENIENT \
        --ASSUME_SORT_ORDER queryname \
        --REMOVE_DUPLICATES false \
        --BARCODE_TAG ~{barcode_tag} 2> ~{picard_mark_duplicates_log}

        # Create the final bam removing the duplicates
        samtools view -F 1804 -f 2 -b ~{final_bam_wdup} | \
        samtools sort -@ ~{samtools_threads} -m ~{samtools_memory_per_thread}M /dev/stdin -o ~{final_bam}

        samtools index ~{final_bam}

        rm ~{filtered_bam}

        # tmp_dup_bam and dinal_bam are the two files necessary for computing qc statistics.

        # TODO: Add library complexity

        # Convert bam to bed.gz and mark duplicates
        # Removing reads that starts and ends at the same position (duplicates)
        bedtools bamtobed -bedpe -i ~{final_bam} | \
            sed 's/_/\t/g' | \
            awk -v OFS="\t" '{if($10=="+"){print $1,$2+4,$6-5,$8}else if($10=="-"){print $1,$2-5,$6+4,$8}}' | \
            sort --parallel=~{cpus} -S 2G  -k4,4 -k1,1 -k2,2n -k3,3n | \
            uniq -c | \
            awk -v OFS="\t" '{print $2, $3, $4, $5, $1}' | \
            sort -k1,1 -k2,2n - | \
            bgzip -c > ~{fragments}

    >>>

    output {
        File? atac_filter_alignment_dedup = final_bam
        File? atac_filter_alignment_dedup_index = final_bam_index
        File? atac_filter_alignment_wdup = final_bam_wdup
        File? atac_filter_fragments = fragments
        File? atac_filter_monitor_log = monitor_log
    }

    runtime {
        cpu : cpus
        memory : "${mem_gb} GB"
        memory_retry_multiplier: 2
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker : "${docker_image}"
        maxRetries:1
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
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz',
                example: ['put link to gcr or dockerhub']
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
