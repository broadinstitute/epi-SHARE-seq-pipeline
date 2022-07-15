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
        Array[File] share_rna_starsolo_ouputs = align.total_outputs
        Array[File]? share_rna_starsolo_logs = align.total_logs
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
            --soloType CB_UMI_Simple \
        --soloCBstart 1 \
        --soloCBlen 24 \
        --soloUMIstart 25 \
        --soloUMIlen 10 \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIfiltering MultiGeneUMI_CR \
        --soloUMIdedup 1MM_CR \
        --clipAdapterType CellRanger4 \
        --outFilterScoreMin 30 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes CR UR CY UY CB UB NH HI AS nM MD
        --runThreadN ${cpus} \
        --chimOutType WithinBAM \
        --genomeDir ./ \
        --readFilesIn ${sep=',' fastq_R1} ${sep=',' fastq_R2}  \
        --outFileNamePrefix out/${default="share-seq" prefix}.rna.align.${genome_name}. \
        #--outFilterMultimapNmax 20 \
        #--outFilterScoreMinOverLread 0.3 \
        #--outFilterMatchNminOverLread 0.3 \
        #--limitOutSJcollapsed 2000000 \
        #--limitIObufferSize 400000000 \
        --outReadsUnmapped Fastx \
        --readFilesCommand zcat
    }
    output{
        Array[File] total_outputs = glob("out/*")
        Array[File]? total_logs = glob("*log*")
    }
    runtime{
        cpu : cpus
        memory : mem_gb+'G'
        disks : 'local-disk ${disk_gb} SSD'
        maxRetries: 0
        docker: docker_image
    }
}
