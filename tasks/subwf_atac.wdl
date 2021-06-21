version 1.0

workflow wf_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) and Sai Ma @ Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Sub-workflow to process the ATAC portion of SHARE-seq libraries.'
    }
    
    input {
        # ATAC Sub-worflow inputs
        File read1
        File read2
        File chrom_sizes
        File idx_tar
        File tss_bed
        String prefix = "shareseq-project"
        String genome_name
        Int cutoff        
        Int? cpus = 4
        String docker = "polumechanos/share-seq"
    }
    
    call align_atac {
        input:
            fastq_R1 = read1,
            fastq_R2 = read2,
            genome_name = genome_name,
            genome_index = idx_tar,
            prefix = prefix,
            cpus = cpus,
            docker_image = docker
    }
  
    call bam_to_bed_atac {
        input:
            bam = align_atac.atac_bowtie2_align,
            bam_index = align_atac.atac_bowtie2_align_index,
            genome_name = genome_name,
            chrom_sizes = chrom_sizes,
            prefix = prefix,
            cpus = cpus,
            docker_image = docker
    }
    
    call count_reads_atac {
        input:
            cutoff = cutoff,
            bedpe = bam_to_bed_atac.bedpe_cleaned,
            genome_name = genome_name,
            prefix = prefix,
            cpus = cpus,
            docker_image = docker
    }
    
    call qc_lib_size {
        input:
            raw_counts = count_reads_atac.atac_unfiltered_counts,
            filtered_counts = count_reads_atac.atac_filtered_counts,
            cutoff = cutoff,
            genome_name = genome_name,
            prefix = prefix,
            docker_image = docker
    }
    
    call qc_stats_atac {
        input:
            raw_bam = align_atac.atac_bowtie2_align,
            raw_bam_index = align_atac.atac_bowtie2_align_index,
            filtered_bam = bam_to_bed_atac.bam_filtered,
            filtered_bam_index = bam_to_bed_atac.bam_filtered_index,
            tss = tss_bed,
            genome_name = genome_name,
            prefix = prefix,
            cpus = cpus,
            docker_image = docker
    }
    
    output {
        File atac_aligned_raw_bam = align_atac.atac_bowtie2_align 
        File atac_aligned_raw_bai = align_atac.atac_bowtie2_align_index
        File atac_bowtie2_log = align_atac.log
        File atac_aligned_filtered_bam = bam_to_bed_atac.bam_filtered
        File atac_aligned_filtered_bai = bam_to_bed_atac.bam_filtered_index
        File atac_fragment_file_raw = bam_to_bed_atac.bedpe_cleaned
        File atac_fragment_file_filtered = count_reads_atac.bedpe_cleaned_filtered
        File atac_counts_filtered = count_reads_atac.atac_filtered_counts
        File atac_counts_unfiltered = count_reads_atac.atac_unfiltered_counts
        File atac_lib_size_count = qc_lib_size.lib_size_counts
        File atac_duplicates_log = qc_lib_size.lib_size_log
        Array[File] atac_lib_size_plots = qc_lib_size.plots
        File atac_qc_final_stats = qc_stats_atac.final_stats
        File atac_qc_hist_plot = qc_stats_atac.final_hist_stats_pdf
        File atac_qc_hist_stats = qc_stats_atac.final_hist_stats
        File atac_qc_tss_pileup = qc_stats_atac.tss_pileup
    }
}


# TASKS

task align_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: align ATAC task'
    }
    
    input {
        # This task takes in input the preprocessed ATAC fastqs and align them to the genome.
        
        File fastq_R1
        File fastq_R2
        File genome_index       # This is a tar.gz folder with all the index files.
        String genome_name
        String? prefix
        String docker_image
        Int cpus = 16
        
    }
    
    Float input_file_size_gb = size(fastq_R1, "G")
    # This is almost fixed for either mouse or human genome
    Int mem_gb = 2
    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    Int disk_gb = 50

    command {
        set -e
        
        tar zxvf ${genome_index} --no-same-owner -C ./
        genome_prefix=$(basename $(find . -type f -name "*.rev.1.bt2") .rev.1.bt2)
        
        
        bowtie2 -X2000 \
            -p ${cpus} \
            --rg-id ${prefix + "."}atac \
            -x $genome_prefix \
            -1 ${fastq_R1} \
            -2 ${fastq_R2} |\
            samtools view \
                -bS \
                -@ ${cpus} \
                - \
                -o ${prefix + "."}atac.align.${genome_name}.bam \
            2> atac.bowtie2.align.${genome_name}.log
        
        samtools sort \
            -@ ${cpus} \
            -m ${mem_gb}G \
            ${prefix + "."}atac.align.${genome_name}.bam > ${prefix + "."}atac.align.${genome_name}.sorted.bam
        samtools index -@ ${cpus} ${prefix + "."}atac.align.${genome_name}.sorted.bam
        
    }
    
    output {
        File atac_bowtie2_align = glob('*.sorted.bam')[0]
        File atac_bowtie2_align_index = glob('*.sorted.bam.bai')[0]
        File log = 'atac.bowtie2.align.${genome_name}.log'
    }

    runtime {
        #cpu : '${cpus}'
        #memory : '${mem_gb} GB'
        #disks : 'local-disk ${disk_gb} SSD'
        #preemptible : 0 
        docker : docker_image
    }
    
    parameter_meta {
        fastq_R1: {
                description: 'Read1 fastq',
                help: 'Processed fastq for read1.',
                example: 'processed.atac.R1.fq.gz',
            }
        fastq_R2: {
                description: 'Read2 fastq',
                help: 'Processed fastq for read2.',
                example: 'processed.atac.R2.fq.gz'
            }
        genome_index: {
                description: 'Bowtie2 indexes',
                help: 'Index files for bowtie2 to use during alignment.',
                examples: ['hg19.tar.gz']
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.',
                examples: ['hg38', 'mm10', 'both']
            }
        prefix: {
                description: 'Prefix for output files',
                help: 'Prefix that will be used to name the output files',
                examples: 'MyExperiment'
            }
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by bowtie2',
                default: 16
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz',
                example: ['put link to gcr or dockerhub']
            }
    }


}

task bam_to_bed_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: ATAC bam to bed task'
    }
    
    input {
        # This task takes in input the aligned bam file and rmeove the low quality reads, the extra chromosomes, marks 
        # the duplicats, and convert to a bedpe file.
        
        File bam
        File bam_index
        String genome_name
        File chrom_sizes
        String docker_image
        String? prefix
        Int cpus = 4
        
    }
    
    Float input_file_size_gb = size(bam, "G")
    Int mem_gb = 2
    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    Int disk_gb = 50
    
    String filtered_chr_bam = '${prefix + '.'}filtered_chr.bam'
    String bedpe = 'tmp.bedpe'
    String final_bam = '${prefix + '.'}atac.cleaned.${genome_name}.bam'
    String final_bam_index = '${prefix + '.'}atac.cleaned.${genome_name}.bam.bai'
    String final_bedpe = '${prefix + '.'}atac.cleaned.${genome_name}.bedpe.gz'

    command<<<
        set -e
        
        # I need to do this because the bam and bai need to be in the same folder but WDL doesn't allow you to 
        # co-localize them in the same path.
        mv ~{bam} in.bam
        mv ~{bam_index} in.bai
        
        # Remove unwanted chromosomes
        chrs=$(samtools view -H in.bam | \
            grep chr | \
            cut -f2 | \
            sed 's/SN://g' | \
            grep -v chrM | \
            grep -v Y | \
            awk '{if(length($0)<6)print}')
        
        # Sort file by name, remove low quality reads, namesort the input bam
        samtools view -b -q 30 -f 0x2 in.bam $(echo $chrs) | \
        samtools sort -@ ~{cpus} -m ~{mem_gb}G -n -o ~{filtered_chr_bam} -
        
        # Convert bam to bed.gz and mark duplicates
        # Removing reads that starts and ends at the same position (duplicates)
        bedtools bamtobed -bedpe -i ~{filtered_chr_bam} | \
            sed 's/_/\t/g' | \
            awk -v OFS="\t" '{if($10=="+"){print $1,$2+4,$6+4,$8}else if($10=="-"){print $1,$2-5,$6-5,$8}}' | \
            sort --parallel=~{cpus} -S ~{mem_gb}G  -k4,4 -k1,1 -k2,2 -k3,3 | \
            uniq -c | \
            awk -v OFS="\t" '{print $2, $3, $4, $5, $1}' > ~{bedpe}
        
        # Convert the bedpe file to a bam file for QC
        bedToBam -i ~{bedpe} -g ~{chrom_sizes} | \
        samtools sort -@ ~{cpus} -o ~{final_bam} -
        
        # and index the bam
        samtools index -@ ~{cpus} ~{final_bam}
        
        # Compress the bedpe file
        pigz --fast -c -p ~{cpus} ~{bedpe} > ~{final_bedpe}
                
    >>>
    
    output {
        File bam_filtered = final_bam
        File bam_filtered_index = final_bam_index
        File bedpe_cleaned = final_bedpe
    }

    runtime {
      #  cpu : '${cpus}'
      #  memory : '${mem_gb} GB'
      #  disks : 'local-disk ${disk_gb} SSD'
      #  preemptible: 0 
        maxRetries : 0  
        docker: docker_image
    }
    
    parameter_meta {
        bam: {
                description: 'bam file',
                help: 'Aligned reads in bam format',
                example: 'aligned.hg38.bam'
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.',
                examples: ['hg38', 'mm10', 'both']
            }
        chrom_sizes: {
                description: 'Chromosomes size file',
                help: 'File with the length of each chromosome',
                examples: 'hg38.chrom.sizes'
            }
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by bowtie2',
                examples: '4'
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz',
                example: ['put link to gcr or dockerhub']
            }
    }


}

task count_reads_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: ATAC count reads task'
    }
    
    input {
        # This task takes in input the filtered and cleaned align reads in bedpe format and counts the reads per barcode.
        
        Int cpus = 4
        Int cutoff = 300
        File bedpe
        String genome_name
        String? prefix
        String docker_image
        
    }
    
    
    Float input_file_size_gb = size(bedpe, "G")
    Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    Float mem_gb = 8.0
    
    String read_groups_freq = 'read_groups_freq.bed'
    String unfiltered_counts = '${prefix + '.'}atac.${genome_name}.unfiltered.counts.csv'
    String read_groups_freq_rmdup = 'read_groups_freq_rmdup.bed'
    String filtered_counts = '${prefix + '.'}atac.${genome_name}.filtered.counts.csv'
    String filtered_bedpe = '${prefix + '.'}atac.${genome_name}.cleaned.filtered.bedpe.gz'

    command <<<
        set -e
        
        # Count unfiltered reads
        zcat ~{bedpe} | awk -v OFS='\t' '{a[$4] += $5} END{for (i in a) print a[i], i}' | awk -v CUT=~{cutoff} -v OFS='\t' '{if($1 >= $CUT ) print }'> ~{read_groups_freq}
        
        Rscript $(which sum_reads.R) ~{read_groups_freq} ~{unfiltered_counts} --save
    
        # Count filtered reads
        zcat ~{bedpe} | cut -f4 | uniq -c | awk -v CUT=~{cutoff} -v OFS='\t' '{if($1 >= $CUT) print }' > ~{read_groups_freq_rmdup}
        
        Rscript $(which sum_reads.R) ~{read_groups_freq_rmdup} ~{filtered_counts} --save
        
        # Remove barcode with low counts from the fragment file for ATAC
        sed -e 's/,/\t/g' ~{filtered_counts} | awk -v CUT=~{cutoff} -v OFS=',' 'NR>=2 {if($5 >= $CUT) print $1,$2,$3,$4} ' > barcodes.txt
        
        grep -wFf barcodes.txt <(zcat ~{bedpe}) | pigz --fast -p ~{cpus} > ~{filtered_bedpe}        
    >>>
    
    output {
        File bedpe_cleaned_filtered = filtered_bedpe
        File atac_filtered_counts = filtered_counts
        File atac_unfiltered_counts = unfiltered_counts
    }

    runtime {
      #  cpu : '${cpus}'
      #  memory : '${mem_gb} GB'
      #  disks : 'local-disk ${disk_gb} SSD'
      #  preemptible: 0 
        maxRetries : 0  
        docker: docker_image
    }
    
    parameter_meta {
        bedpe: {
                description: 'bedpe file',
                help: 'Aligned reads in bedpe format.',
                example: 'aligned.hg38.bedpe'
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.',
                examples: ['hg38', 'mm10', 'both']
            }
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by bowtie2',
                examples: '4'
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz',
                example: ['put link to gcr or dockerhub']
            }
    }


}

task qc_lib_size {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: ATAC library size task'
    }
    
    input {
        # This task computs the the library size for the library.
        
        File raw_counts
        File filtered_counts
        Int cutoff
        String genome_name
        String? prefix
        String docker_image
        
        
    }
    
    Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    Float input_file_size_gb = size(filtered_counts, "G")
    Float mem_gb = 8.0
    

    command {
        set -e
        
        # TODO remove the hard coded file paths from R scripts
        # TODO create only one R script that uses the parameters to discriminate
        # Estimate lib size
        # both
        #Rscript $(which lib_size_sc_V5_species_mixing.R)./ '${prefix + '.'}atac.${genome_name}' ${cutoff} atac --save
        # hg38/mm10
        Rscript $(which lib_size_sc_V5_single_species.R) ${raw_counts} ${filtered_counts} ${cutoff} ${genome_name} ATAC --save

    }
    
    output {
        File lib_size_counts = glob('*.libsize.counts.csv')[0]
        File lib_size_log = glob('*.dups.log')[0]
        Array[File] plots = glob('*.pdf')
    }

    runtime {
      #  cpu : '${cpus}'
      #  memory : '${mem_gb} GB'
      #  disks : 'local-disk ${disk_gb} SSD'
      #  preemptible: 0 
        maxRetries : 0  
        docker: docker_image
    }
    
    parameter_meta {
        raw_counts: {
                description: 'Barcode count csv',
                help: 'Barcode counts from raw bam in csv format.',
                example: 'raw.counts.csv'
            }
        filtered_counts: {
            description: 'Barcode count csv',
            help: 'Barcode counts from filtered bam in csv format.',
            example: 'filtered.counts.csv'
        }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.',
                examples: ['hg38', 'mm10', 'both']
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz',
                example: ['put link to gcr or dockerhub']
            }
    }


}


task qc_stats_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: ATAC qc statistics task'
    }
    
    input {
        # This function takes in input the raw and filtered bams
        # and compute some alignment metrics along with the TSS
        # enrichment plot.
        
        Int cpus= 4
        File raw_bam
        File raw_bam_index
        File filtered_bam
        File filtered_bam_index
        File tss
        String genome_name
        String? prefix
        String docker_image
        
        
    }
    
    Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    Float input_file_size_gb = size(raw_bam, "G")
    Float mem_gb = 8.0
    
    String stats_log= '${prefix + '.'}atac.stats.${genome_name}.log.txt'
    String hist_log= '${prefix + '.'}atac.hist.${genome_name}.log.txt'
    String hist_log_pdf= '${prefix + '.'}atac.hist.${genome_name}.log.pdf'
    String tss_pileup= '${prefix + '.'}atac.tss.pileup.${genome_name}.log.pdf'
    

    command {
        set -e
        
        mv ${raw_bam} in.raw.bam
        mv ${raw_bam_index} in.raw.bai
        mv ${filtered_bam} in.filtered.bam
        mv ${filtered_bam_index} in.filtered.bai
        
        
        echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Raw" > ${stats_log}
        samtools idxstats in.raw.bam >> ${stats_log}
        
        echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Filtered" >> ${stats_log}
        samtools idxstats in.filtered.bam >> ${stats_log}
        
        echo '' > ${hist_log}
        java -jar $(which picard.jar) CollectInsertSizeMetrics \
            VALIDATION_STRINGENCY=SILENT \
            I=in.raw.bam \
            O=${hist_log} \
            H=${hist_log_pdf} \
            W=1000  2>> picard_run.log
        
        # make TSS pileup fig # original code has a 'set +e' why?
        # the pyMakeVplot is missing
        python $(which make-tss-pileup-jbd.py) \
            -a in.filtered.bam \
            -b ${tss} \
            -e 2000 \
            -p ends \
            -v \
            -u \
            -o ${tss_pileup}

    }
    
    output {
        File final_stats = stats_log
        File final_hist_stats_pdf = hist_log_pdf
        File final_hist_stats = stats_log
        File tss_pileup = tss_pileup
    }

    runtime {
      #  cpu : '${cpus}'
      #  memory : '${mem_gb} GB'
      #  disks : 'local-disk ${disk_gb} SSD'
      #  preemptible: 0 
        maxRetries : 0  
        docker: docker_image
    }
    
    parameter_meta {
        raw_bam: {
                description: 'Unfiltered bam',
                help: 'Not filtered alignment bam file.',
                example: 'aligned.hg38.bam'
            }
        raw_bam: {
                description: 'Filtered bam',
                help: 'Filtered alignment bam file. Typically, no duplicates and quality filtered.',
                example: 'aligned.hg38.rmdup.filtered.bam'
            }
        tss: {
                description: 'TSS bed file',
                help: 'List of TSS in bed format used for the enrichment plot.',
                example: 'refseq.tss.bed'
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.',
                examples: ['hg38', 'mm10', 'both']
            }
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by bowtie2',
                examples: '4'
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz',
                example: ['put link to gcr or dockerhub']
            }
    }


}
