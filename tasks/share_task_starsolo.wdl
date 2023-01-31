version 1.0

# TASK
# SHARE-rna-STARsolo

task share_rna_align {
    meta {
        version: 'v0.1'
        author: 'Neva Durand (neva@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: align RNA task'
        attribution: '10X multiome STARsolo run parameters and whitelist from Wold Lab at Caltech https://github.com/detrout/woldlab-rna-seq/'
    }

    input{
        # This function takes in input the pre-processed fastqs 
        # and aligns it to the genome using STARsolo.
        # The CB+UMI is expected to be in read2.

        String method # options: 'share-seq', 'tenX' (for 10X multiome)
        Array[File] fastq_R1
        Array[File] fastq_R2
        File? tenX_whitelist = if method == 'tenX' then 'gs://broad-buenrostro-pipeline-genome-annotations/whitelists/737K-arc-v1-GEX.txt.gz' else '' # is there a nicer way to do this?
        File genome_index_tar
        String genome_name
        String prefix
        String docker_image = 'docker.io/nchernia/share_task_star:1'
        Int cpus = 16
        Float? disk_factor = 50.0
        Float? memory_factor = 2.0
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fastq_R1, 'G') + size(fastq_R2, 'G')

    # Determining memory size based on the size of the input files.
    Float mem_gb = 5.0 + size(genome_index_tar, 'G') + memory_factor * input_file_size_gb

    # Determining disk size based on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then 'SSD' else 'LOCAL'

    String monitor_log = 'monitor.log'

    command <<<
        set -e
     
        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &

        # Untar the genome
        tar xvzf ~{genome_index_tar} --no-same-owner -C ./

        # SHARE-seq
        if ~{method} == 'share-seq' 
        then
            # Generate whitelist
            for fq in ~{sep=' ' fastq_R2}
              do
              gunzip -c "${fq}" | awk 'NR%4==2{dict[substr($1,1,24)]}END{for (i in dict){print i}}' >> share_whitelist.txt
            done

            $(which STAR) \
            --readFilesIn ~{sep=',' fastq_R1} ~{sep=',' fastq_R2}  \
            --readFilesCommand zcat
            --runThreadN ~{cpus} \
            --genomeDir ./ \
            --soloType CB_UMI_Simple \
            --soloFeatures GeneFull  \
            --soloStrand Forward \
            --soloCBwhitelist share_whitelist.txt \
            --soloCBmatchWLtype Exact \
            --soloCBstart 1 \
            --soloCBlen 24 \
            --soloUMIstart 25 \
            --soloUMIlen 10 \
            --soloUMIdedup 1MM_All \
            --chimOutType WithinBAM \
            --limitOutSJcollapsed 2000000 \
            --outFilterMultimapNmax 20 \
            --outFilterScoreMinOverLread 0.3 \
            --outFilterMatchNminOverLread 0.3 \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes CR UR CY UY CB UB NH HI AS nM MD GX GN \
            --outReadsUnmapped Fastx \
            --outFileNamePrefix result/~{prefix}. \            

        # 10X multiome
        else
            gunzip -c ~{tenX_whitelist} > tenX_whitelist.txt

            $(which STAR) \
            --readFilesIn ~{sep=',' fastq_R1} ~{sep=',' fastq_R2}  \
            --readFilesCommand zcat \
            --runThreadN ~{cpus} \
            --genomeDir ./ \
            --genomeLoad NoSharedMemory \
            --soloType CB_UMI_Simple \
            --soloFeatures Gene SJ \
            --soloStrand Forward \
            --soloCellFilter EmptyDrops_CR \
            --soloBarcodeReadLength 0 \
            --soloMultiMappers Unique EM \
            --soloCBwhitelist tenX_whitelist.txt \
            --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
            --soloCBlen 16 \
            --soloUMIlen 12 \
            --soloUMIdedup 1MM_CR \
            --soloUMIfiltering MultiGeneUMI_CR \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --sjdbScore 1 \
            --clipAdapterType CellRanger4 \
            --outFilterType BySJout \
            --outFilterMultimapNmax 20 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS NM MD CB CR CY UB UR UY GX GN \
            --outSAMheaderCommentFile COfile.txt \
            --outSAMheaderHD @HD VN:1.4 SO:coordinate \
            --outSAMunmapped Within \
            --outSAMstrandField intronMotif \
            --outFilterScoreMin 30 \
            --outFileNamePrefix result/~{prefix}.
        fi

        # rename files to include prefix (--outFileNamePrefix doesn't split folder and prefix names in subdirectories) 
        mv result/~{prefix}.Solo.out/ result/Solo.out/
        mv result/Solo.out/Barcodes.stats result/Solo.out/~{prefix}.Barcodes.stats
        mv result/Solo.out/GeneFull/Features.stats result/Solo.out/GeneFull/~{prefix}.Features.stats
        mv result/Solo.out/GeneFull/Summary.csv result/Solo.out/GeneFull/~{prefix}.Summary.csv
        mv result/Solo.out/GeneFull/UMIperCellSorted.txt result/Solo.out/GeneFull/~{prefix}.UMIperCellSorted.txt

        cd result/Solo.out/GeneFull/raw/
        gzip *
        tar -cvzf ~{prefix}.raw.tar.gz *.gz

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

