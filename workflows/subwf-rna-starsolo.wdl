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
    }

    call share_task_align_starsolo as align {
        input:
            fastq_R1 = read1,
            fastq_R2 = read2,
            genome_name = genome_name,
            genome_index_tar = idx_tar,
            prefix = prefix,
            cpus = cpus
    }


    output {
        File share_task_starsolo_output_bam = output_bam
        File share_task_starsolo_log_final_out = log_final_out
        File share_task_starsolo_log_out = log_out
        File share_task_starsolo_log_progress_out = progress_out
        File share_task_starsolo_output_sj = output_sj
        File share_task_starsolo_barcodes_stats = barcodes_stats
        File share_task_starsolo_features_stats = features_stats
        File share_task_starsolo_summary.csv = summary.csv
        File share_task_starsolo_umi_per_cell = umi_per_cell
        File share_task_starsolo_gene_h5_filtered = gene_h5_filtered
        File share_task_starsolo_barcodes_filtered = barcodes_filtered
        File share_task_starsolo_features_fitlered = features_fitlered
        File share_task_starsolo_matrix_filtered = matrix_filtered
        File share_task_starsolo_gene_h5_raw = gene_h5_raw
        File share_task_starsolo_barcodes_raw = barcodes_raw
        File share_task_starsolo_features_raw = features_raw
        File share_task_starsolo_matrix_raw = matrix_raw
    }
}

task share_task_align_starsolo {
    input{
        Array[File] fastq_R1
        Array[File] fastq_R2
        File genome_index_tar
        String genome_name
        String? prefix
        String docker_image = "cumulusprod/starsolo:2.7.10a"
        Int cpus = 16
    }

    Int mem_gb = 64
    Int disk_gb = 250

    command{
         set -e
        # Untar the genome
        tar xvzf ${genome_index_tar} --no-same-owner -C ./

        mkdir out

        $(which STAR) \
        --readFilesIn ${sep=',' fastq_R1} ${sep=',' fastq_R2}  \
        --soloType CB_UMI_Simple \
        --soloCBstart 1 \
        --soloCBlen 24 \
        --soloUMIstart 25 \
        --soloUMIlen 10 \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIfiltering MultiGeneUMI_CR \
        --soloUMIdedup 1MM_CR \
        --soloCBwhitelist None \
        --clipAdapterType CellRanger4 \
        --outFilterScoreMin 30 \
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
    }
    output{
        File output_bam = "result/Aligned.sortedByCoord.out.bam"
        File log_final_out = "result/Log.final.out"
        File log_out = "result/Log.out"
        File log_progress_out = "result/Log.progress.out"
        File output_sj = "result/SJ.out.tab"
        File barcodes_stats = "result/Solo.out/Barcodes.stats"
        File features_stats = "result/Solo.out/Gene/Features.stats"
        File summary.csv = "result/Solo.out/Gene/Summary.csv"
        File umi_per_cell = "result/Solo.out/Gene/UMIperCellSorted.txt"
        File gene_h5_filtered = "result/Solo.out/Gene/filtered/Gene.h5"
        File barcodes_filtered = "result/Solo.out/Gene/filtered/barcodes.tsv"
        File features_fitlered = "result/Solo.out/Gene/filtered/features.tsv"
        File matrix_filtered = "result/Solo.out/Gene/filtered/matrix.mtx"
        File gene_h5_raw = "result/Solo.out/Gene/raw/Gene.h5"
        File barcodes_raw = "result/Solo.out/Gene/raw/barcodes.tsv"
        File features_raw = "result/Solo.out/Gene/raw/features.tsv"
        File matrix_raw = "result/Solo.out/Gene/raw/matrix.mtx"
    }
    runtime{
        cpu : cpus
        memory : mem_gb+'G'
        disks : 'local-disk ${disk_gb} SSD'
        maxRetries: 0
        docker: docker_image
    }
}
