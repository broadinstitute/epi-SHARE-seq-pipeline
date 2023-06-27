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
        Array[File] fastq_R1
        Array[File] fastq_R2
        File whitelist
        File genome_index_tar
        String genome_name
        String prefix
        String chemistry
        
        File? placeholder

        # Runtime parameters
        Int cpus = 16
        Float? disk_factor = 50.0
        Float? memory_factor = 0.5
        String? docker_image = 'us.gcr.io/buenrostro-share-seq/share_task_star'
        
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fastq_R1, 'G') + size(fastq_R2, 'G')

    # Determining memory size based on the size of the input files.
    Float mem_gb = 20.0 + memory_factor * size(genome_index_tar, 'G')

    # Determining disk size based on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then 'SSD' else 'LOCAL'

    String monitor_log = 'monitor.log'

    command <<<
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &

        # Check which fastq contains CB + UMI
        r1_length=$(zcat ~{fastq_R1[0]} | head -2 | tail -1 | awk '{print length($1)}')
        r2_length=$(zcat ~{fastq_R2[0]} | head -2 | tail -1 | awk '{print length($1)}')
        if [ $r1_length -gt $r2_length ]; then
            read_files='~{sep=',' fastq_R1} ~{sep=',' fastq_R2}'
            cb_umi_length=$r2_length
        else
            read_files='~{sep=',' fastq_R2} ~{sep=',' fastq_R1}'
            cb_umi_length=$r1_length
        fi
        echo $r1_length
        echo $r2_length
        echo $read_files

        # Untar the genome
        tar xvzf ~{genome_index_tar} --no-same-owner -C ./

        # SHARE-seq
        if [ '~{chemistry}' == 'shareseq' ]; then
            # Check that CB + UMI length is correct
            if [ $cb_umi_length -ne 34 ]; then
                echo 'CB + UMI length is $cb_umi_length; expected 34'
                exit 1
            fi

            # Generate whitelist
            for fq in ~{sep=' ' fastq_R2}
              do
              gunzip -c "${fq}" | awk 'NR%4==2{dict[substr($1,1,24)]}END{for (i in dict){print i}}' >> shareseq_whitelist.txt
            done

            $(which STAR) \
            --readFilesIn $read_files  \
            --readFilesCommand zcat \
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
            --limitBAMsortRAM 31232551044 \
            --outSAMattributes CR UR CY UY CB UB NH HI AS nM MD GX GN \
            --outReadsUnmapped Fastx \
            --outFileNamePrefix result/ \

            feature_type='GeneFull'

        # 10X v2
        elif [ '~{chemistry}' == '10x_v2' ]; then
            # Check that CB + UMI length is correct
            if [ $cb_umi_length -ne 26 ]; then
                echo 'CB + UMI length is $cb_umi_length; expected 26'
                exit 1
            fi

            if [[ '~{whitelist}' == *.gz ]]; then
                gunzip -c ~{whitelist} > 10x_v2_whitelist.txt
            else
                cat ~{whitelist} > 10x_v2_whitelist.txt
            fi

            $(which STAR) \
            --readFilesIn $read_files \
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
        elif [ '~{chemistry}' == '10x_multiome' ]; then
            # Check that CB + UMI length is correct
            if [ $cb_umi_length -ne 28 ]; then
                echo 'CB + UMI length is $cb_umi_length; expected 28'
                exit 1
            fi

            if [[ '~{whitelist}' == *.gz ]]; then
                gunzip -c ~{whitelist} > 10x_multiome_whitelist.txt
            else
                cat ~{whitelist} > 10x_multiome_whitelist.txt
            fi

            $(which STAR) \
            --readFilesIn $read_files \
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
            --soloCBwhitelist 10x_multiome_whitelist.txt \
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
            --outFileNamePrefix result/ \
            --clipAdapterType CellRanger4 \

            feature_type='Gene'
            # TODO: add the final case in which none of the above is passed.
        fi

        # tar and gzip barcodes, features, and matrix files
        cd result/Solo.out/$feature_type/raw/
        gzip *
        tar -cvzf raw.tar.gz *.gz

        # Move files and rename
        cd ../../../../
        find result -type f -exec mv {} result \;
        cd result
        for file in $(ls)
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
        docker: "${docker_image}"
    }
    parameter_meta {
        fastq_R1: {
            description: 'Read 1 RNA fastq file',
            help: 'Preprocessed RNA fastq for read 1',
            example: 'rna.R1.fastq.gz'
        }
        fastq_R2: {
            description: 'Read 2 RNA fastq file',
            help: 'Preprocessed RNA fastq for read 2',
            example: 'rna.R2.fastq.gz'
        }
        chemistry: {
            description: 'Experiment chemistry',
            help: 'Chemistry/method used in the experiment',
            examples: ['shareseq', '10x_v2', '10x_v3']
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
