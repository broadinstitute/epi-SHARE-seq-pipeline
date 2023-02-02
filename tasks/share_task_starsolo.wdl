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

    input {
        # This function takes in input the pre-processed fastqs 
        # and aligns it to the genome using STARsolo.
        # The CB+UMI is expected to be in read2.
        Array[File] fastq_R1
        Array[File] fastq_R2
        String chemistry
        File whitelists_tsv = 'gs://broad-buenrostro-pipeline-genome-annotations/whitelists/whitelists.tsv'
        File? whitelist
        File genome_index_tar
        String genome_name
        String prefix
        String docker_image = 'us.gcr.io/buenrostro-share-seq/share_task_star'
        Int cpus = 16
        Float? disk_factor = 50.0
        Float? memory_factor = 2.0
    }

    Map[String, File] whitelists = read_map(whitelists_tsv)
    File whitelist_ = if (chemistry != 'shareseq') then select_first([whitelist, whitelists[chemistry]]) else '' 

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
        if [ '~{chemistry}' == 'share-seq' ]; then
            # Generate whitelist
            for fq in ~{sep=' ' fastq_R2}
              do
              gunzip -c "${fq}" | awk 'NR%4==2{dict[substr($1,1,24)]}END{for (i in dict){print i}}' >> shareseq_whitelist.txt
            done

            $(which STAR) \
            --readFilesIn ~{sep=',' fastq_R1} ~{sep=',' fastq_R2}  \
            --readFilesCommand zcat
            --runThreadN ~{cpus} \
            --genomeDir ./ \
            --soloType CB_UMI_Simple \
            --soloFeatures GeneFull  \
            --soloStrand Forward \
            --soloCBwhitelist shareseq_whitelist.txt \
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
            --outFileNamePrefix result/ \

            feature_type='GeneFull'

        # 10X v2
        elif [ '~{chemistry}' == '10x_v2' ]; then
            gunzip -c ~{whitelist_} > 10x_v2_whitelist.txt

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
            --soloCBwhitelist 10x_v2_whitelist.txt \
            --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
            --soloCBlen 16 \
            --soloUMIlen 10 \
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
            --outFileNamePrefix result/

            feature_type='Gene'       

        # 10X v3 (multiome)
        elif [ '~{chemistry}' == '10x_v3' ]; then
            gunzip -c ~{whitelist_} > 10x_v3_whitelist.txt

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
            --soloCBwhitelist 10x_v3_whitelist.txt \
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
            --outFileNamePrefix result/

            feature_type='Gene'

        fi

        ls -R result/

        # tar and gzip barcodes, features, and matrix files
        gzip result/Solo.out/$feature_type/raw/*
        tar -cvzf result/raw.tar.gz result/Solo.out/$feature_type/raw/*.gz

        # Move files and rename
        find result -type f -exec mv {} result \;
        ls result/
        for file in $(ls result)
        do 
            mv $file ~{prefix}.$file
        done

    >>>

    output {
        File output_bam = "result/~{prefix}.Aligned.sortedByCoord.out.bam"
        File log_final_out = "result/~{prefix}.Log.final.out"
        File log_out = "result/~{prefix}.Log.out"
        File log_progress_out = "result/~{prefix}.Log.progress.out"
        File output_sj = "result/~{prefix}.SJ.out.tab"
        File barcodes_stats = "result/~{prefix}.Barcodes.stats"
        File features_stats = "result/~{prefix}.Features.stats"
        File summary_csv = "result/~{prefix}.Summary.csv"
        File umi_per_cell = "result/~{prefix}.UMIperCellSorted.txt"
        File raw_tar = "result/~{prefix}.raw.tar.gz"
    }

    runtime {
        cpu : cpus
        memory : "${mem_gb} GB"
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker: docker_image
    }

    parameter_meta {
        fastq_R1: {
            description: 'Read 1 RNA fastq file',
            help: 'Preprocessed RNA fastq for read 1.',
            example: 'rna.R1.fastq.gz'
        }
        fastq_R2: {
            description: 'Read 2 RNA fastq file',
            help: 'Preprocessed RNA fastq for read 2, containing cell barcode and UMI.',
            example: 'rna.R2.fastq.gz'
        }
        chemistry: {
            description: 'Experiment chemistry',
            help: 'Chemistry/method used in the experiment.',
            examples: ['shareseq', '10x_v2', '10x_v3', 'splitseq']
        }
        whitelists_tsv: {
            description: 'TSV with file paths to whitelists',
            help: 'TSV where each row has two columns: chemistry, and file path to corresponding whitelist',
            example: 'gs://broad-buenrostro-pipeline-genome-annotations/whitelists/whitelists.tsv'
        }
        whitelist: {
            description: 'Barcode whitelist',
            help: 'TXT file containing list of known possible barcodes',
            example: 'gs://broad-buenrostro-pipeline-genome-annotations/whitelists/737K-arc-v1-GEX.txt.gz' 
        }
        genome_index_tar: {
            description: 'Genome index files for STARsolo',
            help: 'TAR containing genome index files for STARsolo to use during alignment',
            examples: ['gs://broad-buenrostro-pipeline-genome-annotations/IGVF_human/star/GRCh38_STAR_2.7.10a.tar.gz', 'gs://broad-buenrostro-pipeline-genome-annotations/mm10/star/Mus_musculus.GRCm38.star.2.7.10.tar.gz']
        }
        genome_name: {
            description: 'Reference name',
            help: 'The name of the reference genome used by the aligner.',
            examples: ['hg38', 'mm10']
        }
        prefix: {
            description: 'Prefix for output files',
            help: 'Prefix that will be used to name the output files',
            example: 'SS-PKR-1'
        }
        docker_image: {
            description: 'Docker image',
            help: 'Docker image for alignment step',
            example: ['us.gcr.io/buenrostro-share-seq/share_task_star']
        }
        cpus: {
            description: 'Number of CPUs',
            help: 'Set the number of CPUs to be used by STARsolo',
            default: 16
        }
        disk_factor: {
            description: 'Disk space scaling factor',
            help: 'Scaling factor used to request disk space based on input size',
            example: 10.0
        }
        memory_factor: {
            description: 'Memory scaling factor',
            help: 'Scaling factor used to request memory based on input size',
            example: 2.0
        }
    }
}

