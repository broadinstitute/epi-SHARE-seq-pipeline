version 1.0

# TASK
# SHARE-rna-update-rgid


task share_rna_update_rgid {
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
        String docker_image = "polumechanos/share_task_update_rgid"
        Int cpus = 4
    }

    #Float input_file_size_gb = size(input[0], "G")
    Int mem_gb = 8
    Int disk_gb = 50
    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)

    String updated_bam = "${default="share-seq" prefix}.rna.reheaded.alignment.${if multimapper then "multi" else "unique"}.${genome_name}.bam"
    String updated_bam_index = "${default="share-seq" prefix}.rna.reheaded.alignment.${if multimapper then "multi" else "unique"}.${genome_name}.bam.bai"

    command <<<
        set -e
        # Update RGID, remove low quality reads and unwanted chrs
        # If keeping multimappers, keep primary aligned reads only,
        # otherwise filter by quality score (-q 30)

        $(which samtools) view -H ~{bam} | sed 's/chrMT/chrM/g' > header.sam

        cat header.sam <($(which samtools) view -@ ~{cpus} ~{bam} | sed 's/chrMT/chrM/g') | \
            $(which samtools) view -@ ~{cpus} -bS ~{if multimapper then "-F 256" else "-q 30"} > ~{updated_bam}

        $(which samtools) index -@ ~{cpus} ~{updated_bam}
    >>>

    output {
        File rna_reheaded_alignment = updated_bam
        File rna_reheaded_alignment_index = updated_bam_index
    }

    runtime {
        cpu : cpus
        memory : mem_gb+'G'
        disks : 'local-disk ${disk_gb} SSD'
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
