version 1.0

# TASK
# SHARE-atac-counts

task count_reads_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: ATAC count reads task'
    }

    input {
        # This task takes in input the fragment file and counts the reads per barcode.
        Int? cpus = 2
        Int? memory_gb = 32
        Int? cutoff = 100
        File fragments_raw
        String genome_name
        #String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_count_atac"
        String docker_image = "polumechanos/share_task_count_atac:monitor"
        String? prefix
    }


    Float input_file_size_gb = size(fragments_raw, "G")
    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    Int disk_gb = 100
    Int mem_gb = memory_gb
    Int mem_sort = 8


    String read_groups_freq = 'read_groups_freq.bed'
    String read_groups_freq_rmdup = 'read_groups_freq_rmdup.bed'
    String filtered_counts = '${default="share-seq" prefix}.atac.counts.${genome_name}.filtered.csv'
    String unfiltered_counts = '${default="share-seq" prefix}.atac.counts.${genome_name}.unfiltered.csv'
    String filtered_fragments = '${default="share-seq" prefix}.atac.counts.fragments.filtered.${genome_name}.tsv.gz'

    command <<<
        set -e

        bash $(which monitor_script.sh) > monitoring.log &

        zcat ~{fragments_raw} | cut -f4 | sort --parallel=~{cpus} -S ~{mem_sort}G -u > observed_barcodes_combinations

        # Count unfiltered reads
        zcat ~{fragments_raw} | awk -v OFS='\t' '{a[$4] += $5} END{for (i in a) print a[i], i}' | awk -v CUT=~{cutoff} -v OFS='\t' '{if($1 >= CUT ) print }'> ~{read_groups_freq}

        Rscript $(which sum_reads.R) ~{read_groups_freq} ~{unfiltered_counts} observed_barcodes_combinations --save

        # Count filtered reads
        zcat ~{fragments_raw} | cut -f4 | sort --parallel=~{cpus} -S ~{mem_sort}G | uniq -c | awk -v CUT=~{cutoff} -v OFS='\t' '{if($1 >= CUT) print }' > ~{read_groups_freq_rmdup}

        Rscript $(which sum_reads.R) ~{read_groups_freq_rmdup} ~{filtered_counts} observed_barcodes_combinations --save

        # Remove barcode with low counts from the fragment file for ATAC
        sed -e 's/,/\t/g' ~{filtered_counts} | awk -v CUT=~{cutoff} -v OFS=',' 'NR>=2 {if($NF >= CUT) print $1,$2,$3} ' > barcodes.txt

        grep -wFf barcodes.txt <(zcat ~{fragments_raw}) | sort --parallel=~{cpus} -S ~{mem_sort}G -k1,1 -k2,2n | bgzip -c > ~{filtered_fragments}
    >>>

    output {
        File atac_barcodes = "barcodes.txt"
        File atac_counts_monitor = "monitoring.log"
        File atac_counts_filtered = filtered_counts
        File atac_counts_unfiltered = unfiltered_counts
        File atac_fragments_filtered = filtered_fragments

    }

    runtime {
        cpu : cpus
        memory : mem_gb+'G'
        disks : 'local-disk ${disk_gb} SSD'
        maxRetries : 0
        docker: docker_image
    }

    parameter_meta {
        fragments_raw: {
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
