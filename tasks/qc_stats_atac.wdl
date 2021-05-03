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
        File filtered_bam
        File tss
        String genome_name
        String type= "atac"
        String docker_image
        
        
    }
    
    Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    Float input_file_size_gb = size(bam, "G")
    Float mem_gb = 8.0
    
    String stats_log= '${type}.stats.${genome_name}.log.txt'
    String hist_log= '${type}.hist.${genome_name}.log.txt'
    String hist_log_pdf= '${type}.hist.${genome_name}.log.pdf'
    String tss_pileup= '${type}.tss.pileup.${genome_name}.log.pdf'
    

    command {
        set -e
        
        
        echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Raw" > ${stats_log}
        samtools idxstats ${raw_bam} >> ${stats_log}
        
        echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Filtered" >> ${stats_log}
        samtools idxstats ${filtered_bam} >> ${stats_log}
        
        echo '' > ${hist_log}
        java -jar /opt/picard CollectInsertSizeMetrics \
            VALIDATION_STRINGENCY=SILENT \
            I=${raw_bam} \
            O=${hist_log} \
            H=${hist_log_pdf} \
            W=1000  2>> picard_run.log
        
        # make TSS pileup fig # original code has a 'set +e' why?
        # the pyMakeVplot is missing
        /opt/pyMakeVplot.py \
            -a ${filtered_bam} \
            -b ${tss} \
            -e 2000 \
            -p ends \
            -v \
            -u \
            -o ${tss_pileup}

    }
    
    output {
        File final_stats= stats_log
        File final_hist_stats_pdf= stats_log_pdf
        File final_hist_stats= stats_log
        File tss_pileup= tss_pileup
    }

    runtime {
        cpu : ${cpus}
        memory : '${mem_gb} GB'
        disks : 'local-disk ${disk_gb} SSD'
        preemptible: 0 
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
                help: 'Filtered alignment bam file. Typically, no duplicates and quality filtered.'
                example: 'aligned.hg38.rmdup.filtered.bam'
            }
        tss: {
                description: 'TSS bed file',
                help: 'List of TSS in bed format used for the enrichment plot.'
                example: 'refseq.tss.bed'
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.'
                examples: ['hg38', 'mm10', 'both']
            },
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by bowtie2'
                examples: '4'
            },
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz'
                example: ['put link to gcr or dockerhub']
            }
    }


}
