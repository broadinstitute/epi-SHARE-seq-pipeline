task qc_lib_size {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: ATAC library size task'
    }
    
    input {
        # This function takes in input the raw and filtered bams
        # and compute some alignment metrics along with the TSS
        # enrichment plot.
        
        File raw_counts
        File filtered_counts
        Int cutoff
        String genome_name
        String type
        String docker_image
        
        
    }
    
    Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    Float input_file_size_gb = size(bam, "G")
    Float mem_gb = 8.0
    

    command {
        set -e
        
        # TODO remove the hard coded file paths from R scripts
        # TODO create only one R script that uses the parameters to discriminate
        # Estimate lib size
        # both
        Rscript /opt/lib_size_sc_V5_species_mixing.R ./ ${this is the prefix to the counts} ${cutoff} ${type} --save
        # hg38/mm10
        Rscript /opt/lib_size_sc_V5_single_species.R ./ ${this is the prefix to the counts} ${cutoff} ${genome_name} ${type} --save

    }
    
    output {
        Array[File] lib_size_unfiltered= glob('*.unfiltered.counts.csv')
        Array[File] lib_size_filtered= glob('*.filtered.counts.csv')
    }

    runtime {
        cpu : ${cpus}
        memory : '${mem_gb} GB'
        disks : 'local-disk ${disk_gb} SSD'
        preemptible: 0 
        docker: docker_image
    }
    
    parameter_meta {
        raw_bam: {
                description: 'Unfiltered bam',
                help: 'Not filtered alignment bam file.',
                example: 'aligned.hg38.bam'
            }
        raw_bam: {
                description: 'Filtered bam',
                help: 'Filtered alignment bam file. Typically, no duplicates and quality filtered.'
                example: 'aligned.hg38.rmdup.filtered.bam'
            }
        tss: {
                description: 'TSS bed file',
                help: 'List of TSS in bed format used for the enrichment plot.'
                example: 'refseq.tss.bed'
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.'
                examples: ['hg38', 'mm10', 'both']
            },
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by bowtie2'
                examples: '4'
            },
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz'
                example: ['put link to gcr or dockerhub']
            }
    }


}
