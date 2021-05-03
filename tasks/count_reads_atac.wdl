task count_reads_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: ATAC count reads task'
    }
    
    input {
        # This function takes in input the raw fastqs from Novaseq or
        # Nextseq and perform trimming, adapter removal, appending
        # cell barcode to the read name, and splitting ATAC and RNA.
        
        Int cpus= 4
        Int cutoff= 300
        File bedpe
        String genome_name
        String type= "atac"
        String docker_image
        
        
    }
    
    Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    Float input_file_size_gb = size(bam, "G")
    Float mem_gb = 8.0
    
    String read_groups_freq= 'read_groups_freq.bed'
    String unfiltered_counts= '${type}.unfiltered.counts.${genome_name}.csv'
    String read_groups_freq_rmdup= 'read_groups_freq_rmdup.bed'
    String filtered_counts= '${type}.filtered.counts.${genome_name}.csv'
    String filtered_bedpe= '${type}.cleaned.filtered.${genome_name}.bedpe.gz'

    command {
        set -e
        
        # Count unfiltered reads
        zcat ${bedpe} | \
        awk -v OFS='\t' '{a[$4] += $5} END{for (i in a) print a[i], i}'| \
        awk -v OFS='\t' '{if($1 >= '${cutoff}') print }'> ${read_groups_freq}
        
        Rscript /opt/sum_reads.R ./ ${read_groups_freq.bed} --save
    
        mv '${read_groups_freq.bed}.csv' ${unfiltered_counts}

        # Count filtered reads
        zcat ${bedpe} | \
        cut -f4 | \
        uniq -c | \
        awk -v OFS='\t' '{if($1 >= '${cutoff}') print }' > ${read_groups_freq_rmdup}
        
        Rscript /opt/sum_reads.R ./ ${read_groups_freq_rmdup} --save
        
        mv '${read_groups_freq_rmdup}.csv' ${filtered_counts}
        
		# Remove barcode with low counts from the fragment file for ATAC
        sed -e 's/,/\t/g' ${filtered_counts} | \
        awk -v OFS=',' 'NR>=2 {if($5 >= '${cutoff}') print $1,$2,$3,$4} ' > barcodes.txt
        
        grep -wFf ${barcodes.txt} <(zcat ${bedpe}) | \
        pigz --fast -p ${cpus}  > ${filtered_bedpe}        
    }
    
    output {
        File bedpe_cleaned_filtered= filtered_bedpe
        File stats_filtered= filtered_counts
        File stats_unfiltered= unfiltered_counts
    }

    runtime {
        cpu : ${cpus}
        memory : '${mem_gb} GB'
        disks : 'local-disk ${disk_gb} SSD'
        preemptible: 0 
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
