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
        String genome_name_1
        String genome_name_2
        String? prefix= "rna"
        String docker_image
        Int cpus= 4
        
    }
    
    Float input_file_size_gb = size(input[0], "G")
    Float mem_gb = 40.0
    Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    
    String genome1_bam= "${prefix}.${genome_name_1}.rigid.reheader.unique.st.bam"
    String genome1_index= "${prefix}.${genome_name_1}.rigid.reheader.unique.st.bam.bai"
    String genome2_bam= "${prefix}.${genome_name_2}.rigid.reheader.unique.st.bam"
    String genome2_index= "${prefix}.${genome_name_2}.rigid.reheader.unique.st.bam.bai"

    command {
        set -e
        
        # Split reads aligned onto a mixed species index into two files
        # one for ech of the indexes
        
        chrs1=`samtools view -H $Name.$Species.rigid.reheader.unique.st.bam | grep ${genome_name_1} | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<8)print}'`
        chrs2=`samtools view -H $Name.$Species.rigid.reheader.unique.st.bam | grep ${genome_name_2} | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<8)print}'`
        
        samtools view -@ ${cpus} -b ${bam} -o temp1.bam `echo ${chrs1[@]}`
        samtools view -@ ${cpus} -b ${bam} -o temp2.bam `echo ${chrs2[@]}`
        samtools view -@ ${cpus} -h temp1.bam | sed 's/${genome_name_1}_/chr/g' | sed 's/chrMT/chrM/g'| samtools view -@ ${cpus} -b -o ${genome1_bam}
        samtools view -@ ${cpus} -h temp2.bam | sed 's/${genome_name_2}_/chr/g' | sed 's/chrMT/chrM/g'| samtools view -@ ${cpus} -b -o ${genome2_bam}

        samtools index -@ ${cpus} ${genome1_bam}
        samtools index -@ ${cpus} ${genome2_bam}
        
    }
    
    output {
        Array[File] rna_splitted_genomes_bam= [${genome1_bam}, ${genome1_bam_index}, ${genome2_bam}, ${genome2_bam_index}]
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
        genome_name_1: {
                description: 'Reference name',
                help: 'The name of the first genome reference used to create the mixed genome index.'
                examples: ['hg38', 'mm10', 'hg19', 'mm9']
            },
        genome_name_2: {
                description: 'Reference name',
                help: 'The name of the second genome reference used to create the mixed genome index.'
                examples: ['hg38', 'mm10', 'hg19', 'mm9']
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
