task group_umi_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: group umi rna task'
    }
    
    input {
        # This function takes in input the bam file produced by the STAR
        # aligner and group UMI reads whilst filtering duplicates.
        
        Boolean qc= true
        File bam
        File genes_annotations_bed
        String genome_name
        String? prefix= "rna"
        String docker_image
        
    }

    Float input_file_size_gb = size(input[0], "G")
    Float mem_gb = 40.0
    Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    
    command {
        set -e
        
        # Calculate gene body coverage and reads distribution
        
        INPUT=${bam}
        
        if [[ '${qc}' == 'false' ]]; then
            samtools view -s 0.01 -o ./temp.1pct.bam ${bam}
            INPUT="./temp.1pct.bam"
        fi
        
        # Calculate read distribution
        python3 /opt/read_distribution.py -i $INPUT -r ${genes_annotations_bed} > ${prefix}.${genome_name}.reads_distribution.txt 2>>./Run.log
        
        # The two files created here are necessary for the
        # Read_distribution.R script.
        tail -n +5 ${prefix}.${genome_name}.reads_distribution.txt | head -n -1 > temp1.txt
        head -n 3  ${prefix}.${genome_name}.reads_distribution.txt | grep Assigned | sed 's/Total Assigned Tags//g' | sed 's/ //g' > temp2.txt
        
        # Plot reads distribution
        Rscript/opt/Read_distribution.R ./ ${prefix}.${genome_name} --save

    }
    
    output {
        File read_distribution_rna= "${prefix}.${genome_name}.reads_distribution.txt"
        File read_distribution2_rna= "${prefix}.${genome_name}.reads_distribution2.txt"
        File read_distribution_plot_rna= "${prefix}.${genome_name}.reads_distribution.pdf"
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
        qc: {
                description: 'QC run flag',
                help: 'Flag to set if you are doing a qc run.'
                default: false
                examples: [true, false]
            },
        genes_annotations_bed: {
                description: 'Genes annotations in bed format',
                help: 'The genes annotations to use when calculating the reads distributions.'
                examples: 'hg38.UCSC_RefSeq.bed'
            },
        genome_name: {
                description: 'Reference name',
                help: 'The name genome reference used to align.'
                examples: ['hg38', 'mm10', 'hg19', 'mm9']
            },
        prefix: {
                description: 'Prefix for output files',
                help: 'Prefix that will be used to name the output files'
                examples: 'MyExperiment'
            },
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: samtools'
                example: ['put link to gcr or dockerhub']
            }
    }
}
