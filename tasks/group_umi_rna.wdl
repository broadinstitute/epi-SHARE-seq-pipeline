task group_umi_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: group umi rna task'
    }
    
    input {
        # This function takes in input the bam file produced by the STAR
        # aligner and group UMI reads whilst filtering duplicates.
        
        Boolean remove_single_umi= true
        File bam
        String genome_name
        String mode= "fast"
        String? prefix= "rna"
        String docker_image
        Int cpus= 4
        Int cutoff= 300
        
    }
    
    Float input_file_size_gb = size(input[0], "G")
    Float mem_gb = 40.0
    Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    
    command {
        set -e
        
        if [[ '${mode}' == 'regular' ]]; then
        
            # Seems to get more slant and fewer UMIs, but get accurate lib size estimation. Slow.
            umi_tools group \
                --extract-umi-method=read_id \
                --per-gene \
                --gene-tag=XT \
                --per-cell \
                -I ${bam} \
                --output-bam -S ./${prefix}.${genome_name}.grouped.bam \
                --group-out=./${prefix}.${genome_name}.groups.tsv \
                --skip-tags-regex=Unassigned >>./Run.log
        else
        
            # Custom UMI dedup by matching bc-umi-align position
            samtools view -@ ${cpus} ${bam} | \
                grep XT:Z: | \
                sed 's/Unassigned_Ambiguity/discard/g' | \
                sed 's/Unassigned_MappingQuality/discard/g' | \
                sed 's/\_/\t/g' | \
                awk -v OFS='\t' '{if($NF ~/discard/){$NF=$(NF-1)} print $5, $6, $2, $3, $NF}' | \
                sed 's/XT:Z://g' > ${prefix}.${genome_name}.wdup.bed
            # remove dup reads
            python3 /opt/rm_dup_barcode_UMI_v3.py \
                    -i ${prefix}.${genome_name}.wdup.bed \
                    -o ${prefix}.${genome_name}.groups.tsv \
                    --m 1
        fi
            
        # convert groupped UMI to bed file
        if [[ '${mode}' == 'regular' ]]; then
            
            ## 3rd column is umi, 4th column is read
            ## note: difference between these two versions of groups.tsv
            ## 1) umitools output keep all the barcode-UMI and don't collapse them
            ## 2) my script already collpsed them at alignment position level
            less ${prefix}.${genome_name}.groups.tsv | \
                sed 's/_/\t/g' | \
                awk -v thr=${if remove_single_umi then 1 else 0} 'FNR > 1 {if($10 > thr){if($9 != "GGGGGGGGGG"){print}}}' | \
                awk 'NR==1{ id="N"} {if(id != $11 ) {id = $11; print $2, $6, $10}}' | \
                awk -v OFS="\t" 'NR==1{ t1=$1;t2=$2;readsum=0; umisum=0} {if(t1==$1 && t2==$2) {readsum+=$3; umisum+=1} \
                    else {print t1, t2, umisum, readsum; t1=$1;t2=$2;umisum=1;readsum=$3}}' | \
                pigz --fast -p $cores > ${prefix}.${genome_name}.bed.gz
        else
        
            less ${prefix}.${genome_name}.groups.tsv | \
                awk -v thr=${if remove_single_umi then 1 else 0} -v OFS='\t' '{if($3 > thr){print}}' | \
                awk -v OFS="\t" 'NR==1{ t1=$1;t2=$2; t3=$3} {if(t1==$1 && t2==$2) {readsum+=$3} else {print t1, t2, t3, readsum; t1=$1;t2=$2;t3=$3;readsum=$3}}' | \
                pigz --fast -p $cores > ${prefix}.${genome_name}.bed.gz
        fi

        # Count unfiltered reads
        zcat ${prefix}.${genome_name}.bed.gz | \
            awk -v OFS='\t' '{a[$1] += $4} END{for (i in a) print a[i], i}' | \
            awk -v OFS='\t' '{if($1 >= '${cutoff}') print }'> ${prefix}.${genome_name}.wdup.RG.freq.bed
        Rscript /opt/sum_reads.R ./ ${prefix}.${genome_name}.wdup.RG.freq.bed --save\
        mv ${prefix}.${genome_name}.wdup.RG.freq.bed.csv ${prefix}.${genome_name}.unfiltered.counts.csv

        # Count filtered reads
        zcat ${prefix}.${genome_name}.bed.gz | \
            awk -v OFS='\t' '{a[$1] += $3} END{for (i in a) print a[i], i}' | \
            awk -v OFS='\t' '{if($1 >= '${cutoff}') print }' > ${prefix}.${genome_name}.rmdup.RG.freq.bed
        Rscript $myPATH/sum_reads.R $dir/fastqs/ ${prefix}.${genome_name}.rmdup.RG.freq.bed --save
        mv ${prefix}.${genome_name}.rmdup.RG.freq.bed.csv ${prefix}.${genome_name}.filtered.counts.csv
        
        # Remove barcode combinations with less then N reads
        sed -e 's/,/\t/g' ${prefix}.${genome_name}.filtered.counts.csv | \
            awk -v OFS=',' 'NR>=2 {if($5 >= '${cutoff}') print $1,$2,$3,$4} ' > ${prefix}.${genome_name}.barcodes.txt
        grep -wFf ${prefix}.${genome_name}.barcodes.txt <(zcat ${prefix}.${genome_name}.bed.gz) | pigz --fast -p $cores > ${prefix}.${genome_name}.cutoff.bed.gz
        
    }
    
    output {
        File filtered_counts_rna= ${prefix}.${genome_name}.filtered.counts.csv
        File unfiltered_counts_rna= ${prefix}.${genome_name}.unfiltered.counts.csv
        File groupped_umi_rna= ${prefix}.${genome_name}.bed.gz
        File groupped_umi_filtered_rna= ${prefix}.${genome_name}.cutoff.bed.gz
    }

    runtime {
        cpu : ${cpus}
        memory : '${mem_gb} GB'
        disks : 'local-disk ${disk_gb} SSD'
        preemptible: 0 
        docker: docker_image
    }
    
    parameter_meta {
        bam: {
                description: 'Alignment bam file',
                help: 'Aligned reads in bam format.',
                example: 'hg38.aligned.bam'
            },
        remove_single_umi: {
                description: 'Single UMI removal flag',
                help: 'Flag to set if you want to remove the UMI with only one read.'
                default: false
                examples: [true, false]
            },
        mode: {
                description: 'Running mode',
                help: 'Decide if you want ot use umi_toolsFlag to set if you want to include reads overlapping introns.'
                default: false
                examples: [true, false]
            },
        genome_name: {
                description: 'Reference name',
                help: 'The name genome reference used to align.'
                examples: ['hg38', 'mm10', 'hg19', 'mm9']
            },
        cutoff: {
                description: 'Read number cutoff',
                help: 'Remove barcodes combination that have less read than the set value (integer).'
                default: 300
            },
        prefix: {
                description: 'Prefix for output files',
                help: 'Prefix that will be used to name the output files'
                examples: 'MyExperiment'
            },
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by bowtie2'
                examples: '4'
            },
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: samtools'
                example: ['put link to gcr or dockerhub']
            }
    }


}
