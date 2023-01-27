version 1.0

# TASK
# 10X-RNA-STARsolo

task 10x_rna_align {
    meta {
        version: 'v0.1'
        source:  'https://portal.firecloud.org/?return=terra#methods/cumulus/star_solo/7/wdl'
        modified_by: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard Multiome pipeline: align RNA task' # no longer share-seq pipeline...
    }

    input {
        String? prefix
        Array[File] fastq_R1
        Array[File] fastq_R2
        # Chemistry, choosing from tenX_v3 (for 10X V3 chemistry), tenX_v2 (for 10X V2 chemistry), DropSeq, SeqWell and custom
        String chemistry
        # Cell barcode start position (1-based coordinate)
        Int? CBstart
        # Cell barcode length
        Int? CBlen
        # UMI start position (1-based coordinate)
        Int? UMIstart
        # UMI length
        Int? UMIlen
        # Cell barcode white list -- do we want this optional argument ??
        # File? CBwhitelist
        Int cpus = 16
        File genome_index_tar
        String genome
        Int? CBstart
        Int? CBlen
        Int? UMIstart
        Int? UMIlen
        File? CBwhitelist
        File? whitelist
        String docker_image = 'docker.io/nchernia/share_task_star:1'
        String star_version
        Int disk_space
        Int memory
    }


    command <<<
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &

        # Untar the genome
        tar xvzf ~{genome_index_tar} --no-same-owner -C ./

        # Make whitelist
        for fastq in ~{sep=' ' fastq_R2}
        do
            gunzip -c '${fastq}' | awk 'NR%4==2{dict[substr($1,1,24)]}END{for (i in dict){print i}}' >> whitelist.txt
        done

        if ~{chemistry} == 'tenX_v3':
            $(which STAR) \ 
            --readFilesIn ~{sep=',' fastq_R1} ~{sep=',' fastq_R2}  \
            --soloType CB_UMI_Simple \
            --soloCBstart 1
            --soloCBlen 16
            --soloUMIstart 17
            --soloUMIlen 12
            --soloCBwhitelist whitelist.txt
            --genomeDir ./ \
            --runThreadN ~{cpus} \
            --outSAMtype BAM Unsorted \ # do we want sorted by coordinate?
            --outSAMheaderHD \\@HD VN:1.4 SO:unsorted # do we want this?
            --outSAMattributes CR UR CY UY CB UB NH HI AS nM MD GX GN \ # added this
            --outFileNamePrefix result/~{default='10X' prefix}. \

        elif ~{chemistry} == 'tenX_v2':
            $(which STAR) \ 
            --readFilesIn ~{sep=',' fastq_R1} ~{sep=',' fastq_R2}  \
            --soloType CB_UMI_Simple \
            --soloCBstart 1
            --soloCBlen 16
            --soloUMIstart 17
            --soloUMIlen 10
            --soloCBwhitelist whitelist.txt
            --genomeDir ./ \
            --runThreadN ~{cpus} \
            --outSAMtype BAM Unsorted \ # do we want sorted by coordinate?
            --outSAMheaderHD \\@HD VN:1.4 SO:unsorted # do we want this?
            --outSAMattributes CR UR CY UY CB UB NH HI AS nM MD GX GN \ # added this
            --outFileNamePrefix result/~{default='10X' prefix}. \

        elif ~{chemistry} in ['SeqWell', 'DropSeq']:
            $(which STAR) \ 
            --readFilesIn ~{sep=',' fastq_R1} ~{sep=',' fastq_R2}  \
            --soloType CB_UMI_Simple \
            --soloCBstart 1
            --soloCBlen 12
            --soloUMIstart 13
            --soloUMIlen 8
            --soloCBwhitelist None
            --genomeDir ./ \
            --runThreadN ~{cpus} \
            --outSAMtype BAM Unsorted \ # do we want sorted by coordinate?
            --outSAMheaderHD \\@HD VN:1.4 SO:unsorted # do we want this?
            --outSAMattributes CR UR CY UY CB UB NH HI AS nM MD GX GN \ # added this
            --outFileNamePrefix result/~{default='10X' prefix}. \

        elif ~{chemistry} == 'custom':
            $(which STAR) \ 
            --readFilesIn ~{sep=',' fastq_R1} ~{sep=',' fastq_R2}  \
            --soloType CB_UMI_Simple \
            --soloCBstart ~{CBstart}
            --soloCBlen ~{CBlen}
            --soloUMIstart ~{UMIstart}
            --soloUMIlen ~{UMIlen}
            --soloCBwhitelist whitelist.txt
            --genomeDir ./ \
            --runThreadN ~{cpus} \
            --outSAMtype BAM Unsorted \ # do we want sorted by coordinate?
            --outSAMheaderHD \\@HD VN:1.4 SO:unsorted # do we want this?
            --outSAMattributes CR UR CY UY CB UB NH HI AS nM MD GX GN \ # added this
            --outFileNamePrefix result/~{default='10X' prefix}. \
       
        # rename files to include prefix (--outFileNamePrefix doesn't split folder and prefix names in subdirectories) 
        mv result/~{default='10X' prefix}.Solo.out/ result/Solo.out/
        mv result/Solo.out/Barcodes.stats result/Solo.out/~{default='10X' prefix}.Barcodes.stats
        mv result/Solo.out/GeneFull/Features.stats result/Solo.out/GeneFull/~{default='10X' prefix}.Features.stats
        mv result/Solo.out/GeneFull/Summary.csv result/Solo.out/GeneFull/~{default='10X' prefix}.Summary.csv
        mv result/Solo.out/GeneFull/UMIperCellSorted.txt result/Solo.out/GeneFull/~{default='10X' prefix}.UMIperCellSorted.txt

    >>>

    output {
        File output_bam = 'result/~{default='10X' prefix}.Aligned.sortedByCoord.out.bam'
        File log_final_out = 'result/~{default='10X' prefix}.Log.final.out'
        File log_out = 'result/~{default='10X' prefix}.Log.out'
        File log_progress_out = 'result/~{default='10X' prefix}.Log.progress.out'
        File output_sj = 'result/~{default='10X' prefix}.SJ.out.tab'
        File barcodes_stats = 'result/Solo.out/~{default='10X' prefix}.Barcodes.stats'
        File features_stats = 'result/Solo.out/GeneFull/~{default='10X' prefix}.Features.stats'
        File summary_csv = 'result/Solo.out/GeneFull/~{default='10X' prefix}.Summary.csv'
        File umi_per_cell = 'result/Solo.out/GeneFull/~{default='10X' prefix}.UMIperCellSorted.txt'
        File raw_tar = 'result/Solo.out/GeneFull/raw/~{default='10X' prefix}.raw.tar.gz'
    }

    runtime {
        cpu: cpus
        memory : '${mem_gb} GB'
        memory_retry_multiplier: 2
        disks: 'local-disk ${disk_gb} ${disk_type}'
        maxRetries: 1
        docker: docker_image
    }
}

