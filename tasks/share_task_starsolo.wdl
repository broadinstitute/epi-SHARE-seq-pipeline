version 1.0

# TASK
# SHARE-atac-STAR
task share_rna_align {
    meta {
        version: 'v0.1'
        author: 'Neva Durand (neva@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: align RNA task'
    }

    input{
        # This function takes in input the pre-processed fastqs 
        # and aligns it to the genome using STAR solo.
        # The CB+UMI is expected to be in read2.
        Array[File] fastq_R1
        Array[File] fastq_R2
        File genome_index_tar
        String genome_name
        String? prefix
        String docker_image = "docker.io/nchernia/share_task_star:1"
        Int cpus = 16
        Float? disk_factor = 50.0
        Float? memory_factor = 2.0
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fastq_R1, "G") + size(fastq_R2, "G")

    # Determining memory size based on the size of the input files.
    Float mem_gb = 5.0 + size(genome_index_tar, "G") + memory_factor * input_file_size_gb

    # Determining disk size based on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    String monitor_log = "monitor.log"

    command <<<
        set -e
     
        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &

        # Untar the genome
        tar xvzf ~{genome_index_tar} --no-same-owner -C ./
        for fq in ~{sep=' ' fastq_R2}
        do
          gunzip -c "${fq}" | awk 'NR%4==2{dict[substr($1,1,24)]}END{for (i in dict){print i}}' >> whitelist.txt
        done
        $(which STAR) \
        --readFilesIn ~{sep=',' fastq_R1} ~{sep=',' fastq_R2}  \
        --soloType CB_UMI_Simple \
        --soloCBstart 1 \
        --soloCBlen 24 \
        --soloUMIstart 25 \
        --soloUMIlen 10 \
        --soloCBmatchWLtype Exact \
        --soloStrand Forward \
        --soloUMIdedup 1MM_All \
        --soloCBwhitelist whitelist.txt \
        --soloFeatures GeneFull  \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes CR UR CY UY CB UB NH HI AS nM MD GX GN \
        --runThreadN ~{cpus} \
        --chimOutType WithinBAM \
        --genomeDir ./ \
        --outFileNamePrefix result/~{prefix}. \
        --outFilterMultimapNmax 20 \
        --outFilterScoreMinOverLread 0.3 \
        --outFilterMatchNminOverLread 0.3 \
        --limitOutSJcollapsed 2000000 \
        --outReadsUnmapped Fastx \
        --readFilesCommand zcat

        cd $(pwd)/result/Solo.out/GeneFull/raw/
        gzip *
        tar -cvf raw.tar *.gz
        gzip raw.tar

    >>>

    output{
        File output_bam = "result/~{prefix}.Aligned.sortedByCoord.out.bam"
        File log_final_out = "result/~{prefix}.Log.final.out"
        File log_out = "result/~{prefix}.Log.out"
        File log_progress_out = "result/~{prefix}.Log.progress.out"
        File output_sj = "result/~{prefix}.SJ.out.tab"
        File barcodes_stats = "result/Solo.out/~{prefix}.Barcodes.stats"
        File features_stats = "result/Solo.out/GeneFull/~{prefix}.Features.stats"
        File summary_csv = "result/Solo.out/GeneFull/~{prefix}.Summary.csv"
        File umi_per_cell = "result/Solo.out/GeneFull/~{prefix}.UMIperCellSorted.txt"
        File raw_tar = "result/Solo.out/GeneFull/raw/~{prefix}.raw.tar.gz"
    }
    runtime{
        cpu : cpus
        memory : "${mem_gb} GB"
        memory_retry_multiplier: 2
        disks: "local-disk ${disk_gb} ${disk_type}"
        maxRetries: 1
        docker: docker_image
    }
}

