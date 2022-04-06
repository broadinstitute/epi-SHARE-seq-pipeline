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
    Int mem_gb = 8
    Int disk_gb = 50
    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)

    String genome1_bam= "${prefix + "."}rna.mixed.${genome1_name}.rigid.reheader.unique.st.bam"
    String genome1_index= "${prefix + "."}.rna.mixed.${genome1_name}.rigid.reheader.unique.st.bam.bai"
    String genome2_bam= "${prefix + "."}rna.mixed.${genome2_name}.rigid.reheader.unique.st.bam"
    String genome2_index= "${prefix + "."}rna.mixed.${genome2_name}.rigid.reheader.unique.st.bam.bai"

    command <<<
        set -e

        # Split reads aligned onto a mixed species index into two files
        # one for ech of the indexes

        chrs1=`samtools view -H ~{bam}| grep ~{genome1_name} | cut -f2 | grep chr | sed 's/SN://g' | awk '{if(length($0)<12)print}'`
        chrs2=`samtools view -H ~{bam}| grep ~{genome2_name} | cut -f2 | grep chr | sed 's/SN://g' | awk '{if(length($0)<12)print}'`

        samtools view -@ ~{cpus} -b ~{bam} -o temp1.bam `echo ${chrs1[@]}`
        samtools view -@ ~{cpus} -b ~{bam} -o temp2.bam `echo ${chrs2[@]}`
        samtools view -@ ~{cpus} -h temp1.bam | sed 's/~{genome1_name}_//g' | sed 's/chrMT/chrM/g'| samtools view -@ ~{cpus} -b -o ~{genome1_bam}
        samtools view -@ ~{cpus} -h temp2.bam | sed 's/~{genome2_name}_//g' | sed 's/chrMT/chrM/g'| samtools view -@ ~{cpus} -b -o ~{genome2_bam}

        #samtools index -@ ${cpus} ~{genome1_bam}
        #samtools index -@ ${cpus} ~{genome2_bam}
    >>>

    output {
        Array[File] rna_splitted_bam = [genome1_bam, genome2_bam]
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
