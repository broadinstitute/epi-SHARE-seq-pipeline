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
        File whitelist
        Int cpus = 16
    }

    Int mem_gb = 64
    Int disk_gb = 300
    command <<<
        set -e
        # Untar the genome
        tar xvzf ~{genome_index_tar} --no-same-owner -C ./
        for fq in ~{sep=' ' fastq_R2}
        do
          gunzip -c "${fq}" | awk 'NR%4==2{dict[substr($1,1,24)]}END{for (i in dict){print i}}' >> whitelist.txt
        done
        $(which STAR) \
        --readFilesIn ${sep=',' fastq_R1} ${sep=',' fastq_R2}  \
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
        --runThreadN ${cpus} \
        --chimOutType WithinBAM \
        --genomeDir ./ \
        --outFileNamePrefix result/ \
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
        File output_bam = "result/Aligned.sortedByCoord.out.bam"
        File log_final_out = "result/Log.final.out"
        File log_out = "result/Log.out"
        File log_progress_out = "result/Log.progress.out"
        File output_sj = "result/SJ.out.tab"
        File barcodes_stats = "result/Solo.out/Barcodes.stats"
        File features_stats = "result/Solo.out/GeneFull/Features.stats"
        File summary_csv = "result/Solo.out/GeneFull/Summary.csv"
        File umi_per_cell = "result/Solo.out/GeneFull/UMIperCellSorted.txt"
        File raw_tar = "result/Solo.out/GeneFull/raw/raw.tar.gz"
    }
    runtime{
        cpu : cpus
        memory : mem_gb+'G'
        disks : 'local-disk ${disk_gb} SSD'
        maxRetries: 0
        docker: docker_image
    }
}

