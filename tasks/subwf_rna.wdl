version 1.0

struct Annotation{
    String genome_name
    File gene_gtf
    File gene_bed
}

workflow wf_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) and Sai Ma @ Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Sub-workflow to process the RNA portion of SHARE-seq libraries.'
    }
    
    input {
        # RNA Sub-worflow inputs

        # Align
        File read1
        File idx_tar
        String prefix = "shareseq-project"
        String genome_name
        Int? cpus = 4
        String docker = "polumechanos/share-seq"
        # Update RGID
        Boolean multimappers = false
        # Split mixed genome
        String? genome_name2
        # Assign features
        Boolean include_multimappers = false
        Boolean include_introns = false
        File gtf
        File? gtf2
        String gene_naming = "gene_name"
        # Group UMI
        Boolean remove_single_umi = true
        String mode = "fast"
        Int cutoff = 300
        # Lib_size QC
        Boolean qc = false
        File genes_annotation_bed
        File? genes_annotation_bed2
    }

    Annotation genome1 = object{
                                genome_name : genome_name,
                                gene_gtf : gtf,
                                gene_bed : genes_annotation_bed
                                }

    Annotation? genome2 = object{
                                genome_name : genome_name2,
                                gene_gtf : gtf2,
                                gene_bed : genes_annotation_bed2
                                 }

    Array[Annotation] annotations = select_all([genome1, genome2])

    call align_rna {
        input:
            fastq_R1 = read1,
            genome_name = genome_name,
            genome_index_tar = idx_tar,
            prefix = prefix,
            cpus = cpus,
            docker_image = docker
    }

    call update_rgid_rna{
        input:
            bam = align_rna.rna_align_bam,
            multimapper = multimappers,
            genome_name = genome_name,
            prefix = prefix,
            docker_image = docker
    }

    if ( defined(genome_name2) ){
        call split_mixed_alignments_rna{
            input:
                bam = update_rgid_rna.rna_rgid_updated_bam,
                bai = update_rgid_rna.rna_rgid_updated_bai,
                genome1_name = genome_name,
                genome2_name = genome_name2,
                prefix = prefix,
                docker_image = docker
        }
        Array[Int]? range = [1,2]
    }

    Array[Int] indexes = select_first([range,[1]])

    Array[File] assign_feature_inputs = select_first([split_mixed_alignments_rna.rna_splitted_bam,
                                                     [update_rgid_rna.rna_rgid_updated_bam]])
    
    scatter( index in indexes ){
        call assign_features_rna{
            input:
                multimapper = include_multimappers,
                intron = include_introns,
                bam = assign_feature_inputs[index],
                gtf = annotations[index]["gene_gtf"],
                gene_naming = gene_naming,
                genome_name = annotations[index]["genome_name"],
                prefix = prefix,
                docker_image = docker
        }

        call group_umi_rna{
            input:
                bam = assign_features_rna.assigned_features_rna_bam,
                mode = mode,
                cutoff = cutoff,
                remove_single_umi = remove_single_umi,
                genome_name = annotations[index]["genome_name"],
                prefix = prefix,
                docker_image = docker
        }

        # TODO: the genes annotation needs to be passed just like the gtf
        # Check which one is the bam in input
        call qc_libsize_rna{
            input:
                qc = qc,
                bam = assign_features_rna.assigned_features_rna_bam,
                genes_annotations_bed = annotations[index]["gene_bed"],
                genome_name = annotations[index]["genome_name"],
                docker_image = docker,
                prefix = prefix,
        }
    }

    output {
        File rna_aligned_raw_bam = align_rna.rna_align_bam
        File rna_aligned_raw_bai = align_rna.rna_align_bai
        File rna_align_log = align_rna.rna_align_log
        Array[File] rna_rgid_updated_bam = assign_feature_inputs
        Array[File] rna_assigned_features_bam = assign_features_rna.assigned_features_rna_bam
        Array[File] rna_assigned_features_bai = assign_features_rna.assigned_features_rna_bai
        Array[File] rna_filtered_counts = group_umi_rna.filtered_counts_rna
        Array[File] rna_unfiltered_counts = group_umi_rna.unfiltered_counts_rna
        Array[File] rna_groupped_umi_unfiltered = group_umi_rna.groupped_umi_rna
        Array[File] rna_groupped_umi_filtered = group_umi_rna.groupped_umi_filtered_rna
        Array[File] rna_read_distribution_extra = qc_libsize_rna.read_distribution2_rna
        Array[File] rna_read_distribution_plot = qc_libsize_rna.read_distribution_plot_rna
    }
}

# TASKS

task align_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: align RNA task'
    }
    
    input {
        # This function takes in input the pre-processed fastq and align it to the genome
        # using STAR.

        File fastq_R1
        File genome_index_tar
        String genome_name
        String? prefix
        String docker_image = "polumechanos/share-seq"
        Int cpus = 4
    }
    #Float input_file_size_gb = size(input[0], "G")
    Int samtools_mem_gb = 1
    #Float mem_gb = 40.0
    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)

    command {
        set -e
        # Untar the genome
        tar xvzf ${genome_index_tar} --no-same-owner -C ./

        mkdir out

        $(which STAR) \
            --runThreadN ${cpus} \
            --chimOutType WithinBAM \
            --genomeDir ./ \
            --readFilesIn ${fastq_R1}  \
            --outFileNamePrefix out/${prefix + "."}rna.${genome_name}. \
            --outFilterMultimapNmax 20 \
            --outFilterMismatchNoverLmax 0.06 \
            --limitOutSJcollapsed 2000000 \
            --outSAMtype BAM Unsorted \
            --limitIObufferSize 400000000 \
            --outReadsUnmapped Fastx \
            --readFilesCommand zcat

        $(which samtools) sort \
            -@ ${cpus} \
            -m ${samtools_mem_gb}G \
            -o out/${prefix + "."}rna.${genome_name}.Aligned.out.sorted.bam \
            out/${prefix + "."}rna.${genome_name}.Aligned.out.bam

        $(which samtools) index \
            -@ ${cpus} \
            out/${prefix + "."}rna.${genome_name}.Aligned.out.sorted.bam
    }
    
    output {
        File rna_align_bam= glob('out/*.sorted.bam')[0]
        File rna_align_bai= glob('out/*.sorted.bam.bai')[0]
        File rna_align_log= glob('out/*.Log.final.out')[0]
    }

    runtime {
#        cpu : ${cpus}
#        memory : '${mem_gb} GB'
#        disks : 'local-disk ${disk_gb} SSD'
#        preemptible: 0
        maxRetries: 0
        docker: docker_image
    }
    
    parameter_meta {
        fastq_R1: {
                description: 'Read1 fastq',
                help: 'Processed fastq for read1.',
                example: 'processed.atac.R1.fq.gz'
            }
        genome_index_tar: {
                description: 'STAR indexes',
                help: 'Index files for STAR to use during alignment in tar.gz.',
                example: ['']
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.',
                example: ['hg38', 'mm10', 'both']
            }
        prefix: {
                description: 'Prefix for output files',
                help: 'Prefix that will be used to name the output files',
                example: 'MyExperiment'
            }
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by bowtie2',
                example: '4'
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: STAR',
                example: ['put link to gcr or dockerhub']
            }
    }
}

task update_rgid_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: update RGID task'
    }
    
    input {
        # This function takes in input the aligned bam and update the read group names
        Boolean multimapper = false
        File bam
        String genome_name
        String? prefix
        String docker_image
        Int cpus = 4
    }
    
    #Float input_file_size_gb = size(input[0], "G")
    #Float mem_gb = 40.0
    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)

    String updated_bam = "${prefix + "."}${genome_name}.rigid.reheader.${if multimapper then "multi" else "unique"}.st.bam"
    String updated_bam_index = "${prefix + "."}${genome_name}.rigid.reheader.${if multimapper then "multi" else "unique"}.st.bam.bai"

    command <<<
        set -e
        # Update RGID, remove low quality reads and unwanted chrs
        # If keepig multimappers, keep primary aligned reads only,
        # otherwise filter by quality score (-q 30)
        
        $(which samtools) view -h -@ ~{cpus} ~{bam} | \
            sed 's/chrMT/chrM/g' | \
            awk -v OFS='\t' '{$1=substr($1,1,length($1)-34)""substr($1,length($1)-22,23)""substr($1,length($1)-34,11); print $0}' | \
            $(which samtools) view -@ ~{cpus} -bS ~{if multimapper then "-F 256" else "-q 30"} > ~{updated_bam}
        $(which samtools) index -@ ~{cpus} ~{updated_bam}
    >>>

    output {
        File rna_rgid_updated_bam = updated_bam
        File rna_rgid_updated_bai = updated_bam_index
    }

    runtime {
        #cpu : ${cpus}
        #memory : '${mem_gb} GB'
        #disks : 'local-disk ${disk_gb} SSD'
        #preemptible: 0
        maxRetries: 0
        docker: docker_image
    }
    
    parameter_meta {
        bam: {
                description: 'Alignment bam file',
                help: 'Aligned reads in bam format.',
                example: 'hg38.aligned.bam'
            }
        multimapper: {
                description: 'Multimappers flag',
                help: 'Flag to set if you want to keep the multimapping reads.',
                default: false,
                example: [true, false]
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.',
                example: ['hg38', 'mm10', 'both']
            }
        prefix: {
                description: 'Prefix for output files',
                help: 'Prefix that will be used to name the output files',
                example: 'MyExperiment'
            }
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by bowtie2',
                example: '4'
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: samtools',
                example: ['put link to gcr or dockerhub']
            }
    }
}

task split_mixed_alignments_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: split mixed alignments task'
    }
    
    input {
        # This function takes in input the bam file produced by the STAR
        # aligner run on a mixed index (e.g. mouse + human) and split
        # the reads into the two genomes
        
        File bam
        File bai
        String genome1_name
        String? genome2_name
        String? prefix
        String docker_image
        Int cpus= 4
        
    }
    
    #Float input_file_size_gb = size(bam, "G")
    #Float mem_gb = 40.0
    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    
    String genome1_bam= "${prefix + "."}rna.mixed.${genome1_name}.rigid.reheader.unique.st.bam"
    String genome1_index= "${prefix + "."}.rna.mixed.${genome1_name}.rigid.reheader.unique.st.bam.bai"
    String genome2_bam= "${prefix + "."}rna.mixed.${genome2_name}.rigid.reheader.unique.st.bam"
    String genome2_index= "${prefix + "."}rna.mixed.${genome2_name}.rigid.reheader.unique.st.bam.bai"

    command <<<
        set -e
        
        # Split reads aligned onto a mixed species index into two files
        # one for ech of the indexes
        
        chrs1=`samtools view -H ~{bam}| grep ~{genome1_name} | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<8)print}'`
        chrs2=`samtools view -H ~{bam}| grep ~{genome2_name} | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<8)print}'`
        
        samtools view -@ ~{cpus} -b ~{bam} -o temp1.bam `echo ${chrs1[@]}`
        samtools view -@ ~{cpus} -b ~{bam} -o temp2.bam `echo ${chrs2[@]}`
        samtools view -@ ~{cpus} -h temp1.bam | sed 's/~{genome1_name}_/chr/g' | sed 's/chrMT/chrM/g'| samtools view -@ ~{cpus} -b -o ~{genome1_bam}
        samtools view -@ ~{cpus} -h temp2.bam | sed 's/~{genome2_name}_/chr/g' | sed 's/chrMT/chrM/g'| samtools view -@ ~{cpus} -b -o ~{genome2_bam}

        #samtools index -@ ${cpus} ~{genome1_bam}
        #samtools index -@ ${cpus} ~{genome2_bam}
    >>>
    
    output {
        Array[File] rna_splitted_bam = [genome1_bam, genome2_bam]
    }

    runtime {
        #cpu : ${cpus}
        #memory : '${mem_gb} GB'
        #disks : 'local-disk ${disk_gb} SSD'
        #preemptible : 0
        maxRetries : 0
        docker: docker_image
    }
    
    parameter_meta {
        bam: {
                description: 'Alignment bam file',
                help: 'Aligned reads in bam format.',
                example: 'hg38.aligned.bam'
            }
        genome1_name: {
                description: 'Reference name',
                help: 'The name of the first genome reference used to create the mixed genome index.',
                example: ['hg38', 'mm10', 'hg19', 'mm9']
            }
        genome2_name: {
                description: 'Reference name',
                help: 'The name of the second genome reference used to create the mixed genome index.',
                example: ['hg38', 'mm10', 'hg19', 'mm9']
            }
        prefix: {
                description: 'Prefix for output files',
                help: 'Prefix that will be used to name the output files',
                example: 'MyExperiment'
            }
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by bowtie2',
                example: '4'
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: samtools',
                example: ['put link to gcr or dockerhub']
            }
    }
}

task assign_features_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: assign features rna task'
    }
    
    input {
        # This function takes in input the bam file produced by the STAR
        # aligner run on a mixed index (e.g. mouse + human) and split
        # the reads into the two genomes
        
        Boolean multimapper = false
        Boolean intron = false
        File bam
        File gtf
        String gene_naming = "gene_name"
        String genome_name
        String? prefix
        String docker_image
        Int cpus= 4
    }
    
    #Float input_file_size_gb = size(input[0], "G")
    #Float mem_gb = 40.0
    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    String out_bam = "${prefix + "."}rna.featureCounts.${if multimapper then "multi." else "unique."}${if intron then "intron." else "exon."}${genome_name}.wdup.bam"
    String out_bai = "${prefix + "."}rna.featureCounts.${if multimapper then "multi." else "unique."}${if intron then "intron." else "exon."}${genome_name}.wdup.bam.bai"

    
    command {
        set -e

        mv ${bam} temp_input.bam
        
        # Count reads in exons
        # If multimappers are selected use '-Q 0 -M' options.
        # For unique mappers use '-Q 30'
        featureCounts -T ${cpus} \
            -Q ${if multimapper then "0 -M " else "30"} \
            -a ${gtf} \
            -t exon \
            -g ${gene_naming} \
            -o ${prefix + "."}rna_featureCount.${genome_name}.feature.count.txt \
            -R BAM \
            temp_input.bam >> featureCount.log
        
        temp_filename=temp_input.bam.featureCounts.bam
        
        # Extract reads that assigned to genes
        if [[ ${intron} == "true" ]]; then
            featureCounts -T ${cpus} \
            -Q ${if multimapper then "0 -M " else "30"} \
            -a ${gtf} \
            -t gene \
            -g ${gene_naming} \
            -o ${prefix + "."}rna.featureCount.${genome_name}.feature.count.txt \
            -R BAM \
            $temp_filename >> featureCount.log
            
            temp_filename=$temp_filename.featureCounts.bam
        fi
        
        samtools sort -@ ${cpus} -m 2G -o ${out_bam} $temp_filename
        samtools index -@ ${cpus} ${out_bam}
        
    }
    
    output {
        File assigned_features_rna_bam = out_bam
        File assigned_features_rna_bai = out_bai
        File featureCounts_log = "featureCount.log"
    }

    runtime {
        #cpu : ${cpus}
        #memory : '${mem_gb} GB'
        #disks : 'local-disk ${disk_gb} SSD'
        #preemptible: 0
        maxRetries : 0
        docker: docker_image
    }
    
    parameter_meta {
        bam: {
                description: 'Alignment bam file',
                help: 'Aligned reads in bam format.',
                example: 'hg38.aligned.bam'
            }
        gtf: {
                description: 'GTF file',
                help: 'Genes definitions in GTF format.',
                example: 'hg38.refseq.gtf'
            }
        multimapper: {
                description: 'Multimappers flag',
                help: 'Flag to set if you want to keep the multimapping reads.',
                default: false,
                example: [true, false]
            }
        intron: {
                description: 'Introns flag',
                help: 'Flag to set if you want to include reads overlapping introns.',
                default: false,
                example: [true, false]
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name genome reference used to align.',
                examples: ['hg38', 'mm10', 'hg19', 'mm9'],
            }
        gene_naming: {
                description: 'Gene nomenclature',
                help: 'Choose if you want to use the official gene symbols (gene_name) or ensemble gene names (gene_id).',
                default: 'gene_name',
                examples: ['gene_name', 'gene_id']
            }
        prefix: {
                description: 'Prefix for output files',
                help: 'Prefix that will be used to name the output files',
                example: 'MyExperiment'
            }
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by bowtie2',
                example: '4'
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: samtools',
                example: ['put link to gcr or dockerhub']
            }
    }
}

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
        String genome_name
        String mode = "fast"
        String? prefix
        String docker_image
        Int cpus = 4
        Int cutoff = 300
        
    }
    
    #Float input_file_size_gb = size(input[0], "G")
    #Float mem_gb = 40.0
    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    
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
                --group-out= ~{prefix + "."}rna.~{genome_name}.groups.tsv \
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
                    -o ~{prefix + "."}rna.~{genome_name}.groups.tsv \
                    --m 1
        fi
            
        # convert groupped UMI to bed file
        if [[ '~{mode}' == 'regular' ]]; then
            
            ## 3rd column is umi, 4th column is read
            ## note: difference between these two versions of groups.tsv
            ## 1) umitools output keep all the barcode-UMI and don't collapse them
            ## 2) my script already collpsed them at alignment position level
            less ~{prefix + "."}rna.~{genome_name}.groups.tsv | \
                sed 's/_/\t/g' | \
                awk -v thr=~{if remove_single_umi then 1 else 0} 'FNR > 1 {if($10 > thr){if($9 != "GGGGGGGGGG"){print}}}' | \
                awk 'NR==1{ id="N"} {if(id != $11 ) {id = $11; print $2, $6, $10}}' | \
                awk -v OFS="\t" 'NR==1{ t1=$1;t2=$2;readsum=0; umisum=0} {if(t1==$1 && t2==$2) {readsum+=$3; umisum+=1} \
                    else {print t1, t2, umisum, readsum; t1=$1;t2=$2;umisum=1;readsum=$3}}' | \
                pigz --fast -p ~{cpus} > ~{prefix + "."}rna.~{genome_name}.bed.gz
        else
        
            less ~{prefix + "."}rna.~{genome_name}.groups.tsv | \
                awk -v thr=~{if remove_single_umi then 1 else 0} -v OFS='\t' '{if($3 > thr){print}}' | \
                awk -v OFS="\t" 'NR==1{ t1=$1;t2=$2; t3=$3} {if(t1==$1 && t2==$2) {readsum+=$3} else {print t1, t2, t3, readsum; t1=$1;t2=$2;t3=$3;readsum=$3}}' | \
                pigz --fast -p ~{cpus} > ~{prefix + "."}rna.~{genome_name}.bed.gz
        fi

        # Count unfiltered reads
        zcat ~{prefix + "."}rna.~{genome_name}.bed.gz | \
            awk -v OFS='\t' '{a[$1] += $4} END{for (i in a) print a[i], i}' | \
            awk -v thr=~{cutoff} -v OFS='\t' '{if($1 >= thr) print }'> ~{prefix + "."}rna.~{genome_name}.wdup.RG.freq.bed

        Rscript $(which sum_reads.R) ~{prefix + "."}rna.~{genome_name}.wdup.RG.freq.bed ~{prefix + "."}rna.~{genome_name}.unfiltered.counts.csv --save

        # Count filtered reads
        zcat ~{prefix + "."}rna.~{genome_name}.bed.gz | \
            awk -v OFS='\t' '{a[$1] += $3} END{for (i in a) print a[i], i}' | \
            awk -v thr=~{cutoff} -v OFS='\t' '{if($1 >= thr) print }' > ~{prefix + "."}rna.~{genome_name}.rmdup.RG.freq.bed

        Rscript $(which sum_reads.R) ~{prefix + "."}rna.~{genome_name}.rmdup.RG.freq.bed ~{prefix + "."}rna.~{genome_name}.filtered.counts.csv --save

        # Remove barcode combinations with less then N reads
        sed -e 's/,/\t/g' ~{prefix + "."}rna.~{genome_name}.filtered.counts.csv | \
            awk -v thr=~{cutoff} -v OFS=',' 'NR>=2 {if($5 >= thr) print $1,$2,$3,$4} ' > ~{prefix + "."}rna.~{genome_name}.barcodes.txt
        grep -wFf ~{prefix + "."}rna.~{genome_name}.barcodes.txt <(zcat ~{prefix + "."}rna.~{genome_name}.bed.gz) | \
        pigz --fast -p ~{cpus} > ~{prefix + "."}rna.~{genome_name}.cutoff.bed.gz
        
    >>>
    
    output {
        File filtered_counts_rna = "${prefix + "."}rna.${genome_name}.filtered.counts.csv"
        File unfiltered_counts_rna = "${prefix + "."}rna.${genome_name}.unfiltered.counts.csv"
        File groupped_umi_rna = "${prefix + "."}rna.${genome_name}.bed.gz"
        File groupped_umi_filtered_rna = "${prefix + "."}rna.${genome_name}.cutoff.bed.gz"
    }

    runtime {
     #   cpu : ${cpus}
     #   memory : '${mem_gb} GB'
     #   disks : 'local-disk ${disk_gb} SSD'
     #   preemptible: 0
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

task qc_libsize_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: group umi rna task'
    }
    
    input {
        # This function takes in input the bam file produced by the STAR
        # aligner and group UMI reads whilst filtering duplicates.
        Boolean qc = false
        File bam
        File genes_annotations_bed
        String genome_name
        String? prefix
        String docker_image
    }

    #Float input_file_size_gb = size(input[0], "G")
    #Float mem_gb = 40.0
    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    
    command {
        set -e
        # Calculate gene body coverage and reads distribution
        INPUT=${bam}

        if [[ '${qc}' == 'true' ]]; then
            $(which samtools) view -s 0.01 -o temp.1pct.bam ${bam}
            INPUT="temp.1pct.bam"
        fi
        
        ## TODO: Where is this python script?
        # Calculate read distribution
        # python3 $(which read_distribution.py) -i $INPUT -r ${genes_annotations_bed} > ${prefix + "."}rna.${genome_name}.reads_distribution.txt 2>>./Run.log
        
        # The two files created here are necessary for the
        # Read_distribution.R script.
        tail -n +5 ${prefix + "."}rna.${genome_name}.reads_distribution.txt | head -n -1 > temp1.txt
        head -n 3  ${prefix + "."}rna.${genome_name}.reads_distribution.txt | grep Assigned | sed 's/Total Assigned Tags//g' | sed 's/ //g' > temp2.txt
        
        # Plot reads distribution
        Rscript $(which Read_distribution.R) . ${prefix + "."}rna.${genome_name} --save

    }
    
    output {
        #File read_distribution_rna= "${prefix + "."}rna.${genome_name}.reads_distribution.txt"
        File read_distribution2_rna= "${prefix + "."}rna.${genome_name}.reads_distribution2.txt"
        File read_distribution_plot_rna= "${prefix + "."}rna.${genome_name}.reads_distribution.pdf"
    }

    runtime {
     #   cpu : ${cpus}
     #   memory : '${mem_gb} GB'
     #   disks : 'local-disk ${disk_gb} SSD'
     #   preemptible: 0
        maxRetries : 0
        docker: docker_image
    }
    
    parameter_meta {
        bam: {
                description: 'Alignment bam file',
                help: 'Aligned reads in bam format.',
                example: 'hg38.aligned.bam'
            }
        qc: {
                description: 'QC run flag',
                help: 'Flag to set if you are doing a qc run.',
                default: false,
                example: [true, false]
            }
        genes_annotations_bed: {
                description: 'Genes annotations in bed format',
                help: 'The genes annotations to use when calculating the reads distributions.',
                example: 'hg38.UCSC_RefSeq.bed'
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name genome reference used to align.',
                example: ['hg38', 'mm10', 'hg19', 'mm9']
            }
        prefix: {
                description: 'Prefix for output files',
                help: 'Prefix that will be used to name the output files',
                example: 'MyExperiment'
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: samtools',
                example: ['put link to gcr or dockerhub']
            }
    }
}
