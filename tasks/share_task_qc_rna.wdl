version 1.0

# TASK
# SHARE-rna-qc


task qc_rna {
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
        String docker_image = "polumechanos/share_task_qc_rna"
    }

    #Float input_file_size_gb = size(input[0], "G")
    Int mem_gb = 8
    Int disk_gb = 50
    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)

    String reads_distribution = "${default="share-seq" prefix}.rna.qc.${genome_name}.reads_distribution.txt"
    # Generated automatically inside the R scripts
    String reads_distribution2 = "${default="share-seq" prefix}.rna.qc.${genome_name}.reads_distribution.txt"
    String reads_distribution_plot = "${default="share-seq" prefix}.rna.qc.${genome_name}.reads_distribution.pdf"

    command {
        set -e
        # Calculate gene body coverage and reads distribution
        ln -s ${bam} in.bam
        INPUT=in.bam

        ln -s ${genes_annotations_bed} gene.bed.gz

        gunzip gene.bed.gz

        samtools index in.bam

        if [[ '${qc}' == 'true' ]]; then
            $(which samtools) view -s 0.01 -o temp.1pct.bam ${bam}
            INPUT="temp.1pct.bam"
        fi

        ## TODO: Where is this python script?
        # Calculate read distribution
        python3 $(which read_distribution.py) -i $INPUT -r gene.bed > ${reads_distribution} 2>>./Run.log

        # The two files created here are necessary for the
        # Read_distribution.R script.
        tail -n +5 ${reads_distribution} | head -n -1 > temp1.txt
        head -n 3  ${reads_distribution} | grep Assigned | sed 's/Total Assigned Tags//g' | sed 's/ //g' > temp2.txt

        # Plot reads distribution
        Rscript $(which Read_distribution.R) . ${default="share-seq" prefix}.rna.qc.${genome_name} --save

    }

    output {
        File rna_qc_reads_distribution = reads_distribution
        File rna_qc_reads_distribution2 = reads_distribution2
        File rna_qc_reads_distribution_plot = reads_distribution_plot
    }

    runtime {
        #cpu : cpus
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
