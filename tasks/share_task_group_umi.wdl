version 1.0

# TASK
# SHARE-rna-group-umi


task group_umi_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: group umi rna task'
    }

    input {
        # This function takes in input the bam file produced by the STAR
        # aligner and group UMI reads whilst filtering duplicates.

        Boolean remove_single_umi = true
        File bam
        Int cpus = 4
        Int cutoff
        String genome_name
        String mode
        String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_group_umi"
        String? prefix
        Int? memory_gb = 16


    }

    #Float input_file_size_gb = size(input[0], "G")
    Int mem_gb = memory_gb
    Int disk_gb = 50
    Int mem_sort = 16
    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)

    String umi_groups_table = "${default="share-seq" prefix}.rna.umi.groups.wdup.${genome_name}.tsv"
    String umi_groups_bed_unfiltered = "${default="share-seq" prefix}.rna.umi.groups.unfiltered.wdup.${genome_name}.bed.gz"
    String umi_groups_bed_filtered = "${default="share-seq" prefix}.rna.umi.groups.filtered.wdup.${genome_name}.bed.gz"
    String umi_counts_unfiltered = "${default="share-seq" prefix}.rna.umi.counts.unfiltered.wdup.${genome_name}.csv"
    String umi_counts_filtered = "${default="share-seq" prefix}.rna.umi.counts.filtered.wdup.${genome_name}.csv"
    String umi_barcodes = "${default="share-seq" prefix}.rna.umi.barcodes.filtered.wdup.${genome_name}.txt"


    command <<<
        set -e

        if [[ '~{mode}' == 'regular' ]]; then
            # Seems to get more slant and fewer UMIs, but get accurate lib size estimation. Slow.
            umi_tools group \
                --extract-umi-method=read_id \
                --per-gene \
                --gene-tag=XT \
                --per-cell \
                -I ~{bam} \
                --output-bam -S ~{prefix + "."}rna.~{genome_name}.grouped.bam \
                --group-out= ~{umi_groups_table} \
                --skip-tags-regex=Unassigned >>./Run.log
        else
            # Custom UMI dedup by matching bc-umi-align position
            samtools view -@ ~{cpus} ~{bam} | \
                grep XT:Z: | \
                sed 's/Unassigned_Ambiguity/discard/g' | \
                sed 's/Unassigned_MappingQuality/discard/g' | \
                sed 's/\_/\t/g' | \
                awk -v OFS='\t' '{if($NF ~/discard/){$NF=$(NF-1)} print $5, $6, $2, $3, $NF}' | \
                sed 's/XT:Z://g' > ~{prefix + "."}rna.~{genome_name}.wdup.bed
            # remove dup reads
            python3 $(which rm_dup_barcode_UMI_v3.py) \
                    -i ~{prefix + "."}rna.~{genome_name}.wdup.bed \
                    -o ~{umi_groups_table} \
                    --m 1
        fi

        cut -f1 ~{umi_groups_table} | sort -u > observed_barcodes_combinations

        # convert groupped UMI to bed file
        if [[ '~{mode}' == 'regular' ]]; then

            ## 3rd column is umi, 4th column is read
            ## note: difference between these two versions of groups.tsv
            ## 1) umitools output keep all the barcode-UMI and don't collapse them
            ## 2) my script already collpsed them at alignment position level
            less ~{umi_groups_table} | \
                sed 's/_/\t/g' | \
                awk -v thr=~{if remove_single_umi then 1 else 0} 'FNR > 1 {if($10 > thr){if($9 != "GGGGGGGGGG"){print}}}' | \
                awk 'NR==1{ id="N"} {if(id != $11 ) {id = $11; print $2, $6, $10}}' | \
                awk -v OFS="\t" 'NR==1{ t1=$1;t2=$2;readsum=0; umisum=0} {if(t1==$1 && t2==$2) {readsum+=$3; umisum+=1} \
                    else {print t1, t2, umisum, readsum; t1=$1;t2=$2;umisum=1;readsum=$3}}' | \
                pigz --fast -p ~{cpus} > ~{umi_groups_bed_unfiltered}
        else

            less ~{umi_groups_table} | \
                sort --parallel=~{cpus} -S ~{mem_sort}G -k1,1 -k2,2 | \
                awk -v thr=~{if remove_single_umi then 1 else 0} -v OFS='\t' '{if($3 > thr){print}}' | \
                awk -v OFS="\t" 'NR==1 { t1=$1;t2=$2;readsum=0; umisum=0} {if(t1==$1 && t2==$2) {readsum+=$3; umisum+=1} else {print t1, t2, umisum, readsum; t1=$1;t2=$2;umisum=1;readsum=$3}}' | \
                pigz --fast -p ~{cpus} > ~{umi_groups_bed_unfiltered}

        fi

        # Count unfiltered reads
        zcat ~{umi_groups_bed_unfiltered} | \
            awk -v OFS='\t' '{a[$1] += $4} END{for (i in a) print a[i], i}' | \
            awk -v thr=~{cutoff} -v OFS='\t' '{if($1 >= thr) print }'> ~{prefix + "."}rna.~{genome_name}.wdup.RG.freq.bed

        Rscript $(which sum_reads.R) ~{prefix + "."}rna.~{genome_name}.wdup.RG.freq.bed ~{umi_counts_unfiltered} observed_barcodes_combinations --save

        # Count filtered reads
        zcat ~{umi_groups_bed_unfiltered} | \
            awk -v OFS='\t' '{a[$1] += $3} END{for (i in a) print a[i], i}' | \
            awk -v thr=~{cutoff} -v OFS='\t' '{if($1 >= thr) print }' > ~{prefix + "."}rna.~{genome_name}.rmdup.RG.freq.bed

        Rscript $(which sum_reads.R) ~{prefix + "."}rna.~{genome_name}.rmdup.RG.freq.bed ~{umi_counts_filtered} observed_barcodes_combinations --save

        # Remove barcode combinations with less then N reads
        awk -v thr=~{cutoff} -v FS=',' -v OFS=',' 'NR>1 && $NF>=thr {NF--; print } ' ~{umi_counts_filtered} > ~{umi_barcodes}
        grep -wFf ~{umi_barcodes} <(zcat ~{umi_groups_bed_unfiltered}) | \
        pigz --fast -p ~{cpus} > ~{umi_groups_bed_filtered}

    >>>

    output {
        File rna_umi_barcodes_filtered = "${umi_barcodes}"
        File rna_umi_bed_filtered = "${umi_groups_bed_filtered}"
        File rna_umi_bed_unfiltered = "${umi_groups_bed_unfiltered}"
        File rna_umi_counts_filtered = "${umi_counts_filtered}"
        File rna_umi_counts_unfiltered = "${umi_counts_unfiltered}"
    }

    runtime {
        cpu : cpus
        memory : mem_gb+'G'
        disks : 'local-disk ${disk_gb} SSD'
        maxRetries : 0
        docker: docker_image
    }

    parameter_meta {
        bam: {
                description: 'Alignment bam file',
                help: 'Aligned reads in bam format.',
                example: 'hg38.aligned.bam'
            }
        remove_single_umi: {
                description: 'Single UMI removal flag',
                help: 'Flag to set if you want to remove the UMI with only one read.',
                default: false,
                example: [true, false]
            }
        mode: {
                description: 'Running mode',
                help: 'Decide if you want ot use umi_toolsFlag to set if you want to include reads overlapping introns.',
                default: false,
                example: [true, false]
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name genome reference used to align.',
                example: ['hg38', 'mm10', 'hg19', 'mm9']
            }
        cutoff: {
                description: 'Read number cutoff',
                help: 'Remove barcodes combination that have less read than the set value (integer).',
                default: 300
            }
        prefix: {
                description: 'Prefix for output files',
                help: 'Prefix that will be used to name the output files',
                examples: 'MyExperiment'
            }
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by the task',
                examples: '4'
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: samtools',
                example: ['put link to gcr or dockerhub']
            }
    }
}
