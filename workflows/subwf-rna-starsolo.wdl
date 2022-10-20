version 1.0

# Import the tasks called by the pipeline
workflow wf_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) and Sai Ma @ Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Sub-workflow to process the RNA portion of SHARE-seq libraries.'
    }

    input {
        # RNA Sub-worflow inputs

        # Align
        Array[File] read1
        Array[File] read2
        File idx_tar
        String prefix = "shareseq-project"
        String genome_name
        Int? cpus = 16
        String? docker
        # Update RGID
        Boolean multimappers = false
        # Assign features
        Boolean include_multimappers = false
        Boolean include_introns = false
        File gtf
        String gene_naming = "gene_name"
        # Group UMI
        Boolean remove_single_umi = true
        String mode = "fast"
        Int cutoff = 100
        File whitelist
    }

    call share_task_align_starsolo as align {
        input:
            fastq_R1 = read1,
            fastq_R2 = read2,
            genome_name = genome_name,
            genome_index_tar = idx_tar,
            prefix = prefix,
            whitelist = whitelist,
            cpus = cpus
    }


    output {
        File share_task_starsolo_output_bam = align.output_bam
        File share_task_starsolo_log_final_out = align.log_final_out
        File share_task_starsolo_log_out = align.log_out
        File share_task_starsolo_log_progress_out = align.log_progress_out
        File share_task_starsolo_output_sj = align.output_sj
        File share_task_starsolo_barcodes_stats = align.barcodes_stats
        File share_task_starsolo_features_stats = align.features_stats
        File share_task_starsolo_summary_csv = align.summary_csv
        File share_task_starsolo_umi_per_cell = align.umi_per_cell
        File share_task_starsolo_barcodes_raw = align.barcodes_raw
        File share_task_starsolo_features_raw = align.features_raw
        File share_task_starsolo_matrix_raw = align.matrix_raw
    }
}

task share_task_align_starsolo {
    input{
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
    command{
        set -e
        # Untar the genome
        tar xvzf ${genome_index_tar} --no-same-owner -C ./

        $(which STAR) \
        --readFilesIn ${sep=',' fastq_R1} ${sep=',' fastq_R2}  \
        --soloType CB_UMI_Simple \
        --soloCBstart 1 \
        --soloCBlen 24 \
        --soloUMIstart 25 \
        --soloUMIlen 10 \
        --soloCBmatchWLtype Exact \
	--soloStrand Forward
        --soloUMIfiltering MultiGeneUMI_CR \
        --soloUMIdedup 1MM_CR \
        --soloCBwhitelist ${whitelist} \
        --soloFeatures Gene  \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes CR UR CY UY CB UB NH HI AS nM MD \
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

        gzip $(pwd)/result/Solo.out/Gene/raw/*
    }
    output{
        File output_bam = "result/Aligned.sortedByCoord.out.bam"
        File log_final_out = "result/Log.final.out"
        File log_out = "result/Log.out"
        File log_progress_out = "result/Log.progress.out"
        File output_sj = "result/SJ.out.tab"
        File barcodes_stats = "result/Solo.out/Barcodes.stats"
        File features_stats = "result/Solo.out/Gene/Features.stats"
        File summary_csv = "result/Solo.out/Gene/Summary.csv"
        File umi_per_cell = "result/Solo.out/Gene/UMIperCellSorted.txt"
        File matrix_raw = "result/Solo.out/Gene/raw/matrix.mtx.gz"
        File barcodes_raw = "result/Solo.out/Gene/raw/barcodes.tsv.gz"
        File features_raw = "result/Solo.out/Gene/raw/features.tsv.gz"

    }
    runtime{
        cpu : cpus
        memory : mem_gb+'G'
        disks : 'local-disk ${disk_gb} SSD'
        maxRetries: 0
        docker: docker_image
    }
}
