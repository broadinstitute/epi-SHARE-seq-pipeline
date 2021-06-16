version 1.0

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

    }
    
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
                genome_name_1 = genome_name,
                genome_name_2 = genome_name2,
                prefix = prefix,
                docker_image = docker
        }
    }

    Array[Array[File]] assign_feature_input = select_first([split_mixed_alignments_rna.rna_splitted_genomes_bam,
                                                           [[update_rgid_rna.rna_rgid_updated_bam, update_rgid_rna.rna_rgid_updated_bai]]])

    output {
        File rna_aligned_raw_bam = align_rna.rna_align_bam
        File rna_aligned_raw_bai = align_rna.rna_align_bai
        File rna_align_log = align_rna.rna_align_log
        Array[Array[File]] rna_rgid_updated_bam = assign_feature_input
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
        String genome_name_1
        String? genome_name_2
        String? prefix
        String docker_image
        Int cpus= 4
        
    }
    
    #Float input_file_size_gb = size(bam, "G")
    #Float mem_gb = 40.0
    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    
    String genome1_bam= "${prefix + "."}rna.mixed.${genome_name_1}.rigid.reheader.unique.st.bam"
    String genome1_index= "${prefix + "."}.rna.mixed.${genome_name_1}.rigid.reheader.unique.st.bam.bai"
    String genome2_bam= "${prefix + "."}rna.mixed.${genome_name_2}.rigid.reheader.unique.st.bam"
    String genome2_index= "${prefix + "."}rna.mixed.${genome_name_2}.rigid.reheader.unique.st.bam.bai"

    command <<<
        set -e
        
        # Split reads aligned onto a mixed species index into two files
        # one for ech of the indexes
        
        chrs1=`samtools view -H ~{bam}| grep ~{genome_name_1} | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<8)print}'`
        chrs2=`samtools view -H ~{bam}| grep ~{genome_name_2} | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<8)print}'`
        
        samtools view -@ ~{cpus} -b ~{bam} -o temp1.bam `echo ${chrs1[@]}`
        samtools view -@ ~{cpus} -b ~{bam} -o temp2.bam `echo ${chrs2[@]}`
        samtools view -@ ~{cpus} -h temp1.bam | sed 's/~{genome_name_1}_/chr/g' | sed 's/chrMT/chrM/g'| samtools view -@ ~{cpus} -b -o ~{genome1_bam}
        samtools view -@ ~{cpus} -h temp2.bam | sed 's/~{genome_name_2}_/chr/g' | sed 's/chrMT/chrM/g'| samtools view -@ ~{cpus} -b -o ~{genome2_bam}

        samtools index -@ ${cpus} ~{genome1_bam}
        samtools index -@ ${cpus} ~{genome2_bam}
    >>>
    
    output {
        Array[Array[File]] rna_splitted_genomes_bam = [[genome1_bam, genome1_index], [genome2_bam, genome2_index]]
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
        genome_name_1: {
                description: 'Reference name',
                help: 'The name of the first genome reference used to create the mixed genome index.',
                example: ['hg38', 'mm10', 'hg19', 'mm9']
            }
        genome_name_2: {
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
