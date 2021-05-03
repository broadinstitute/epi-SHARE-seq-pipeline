task pre_process {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: preprocess task'
    }
    
    input {
        # This function takes in input the raw fastqs from Novaseq or
        # Nextseq and perform trimming, adapter removal, appending
        # cell barcode to the read name, and splitting ATAC and RNA.
        
        Array[File] input
        File config
        String? prefix
        Boolean? qc= false
        String docker_image
        
    }
    
    Float input_file_size_gb = size(input[0], "G")
    Float mem_gb = 8.0
    Int disk_gb = round(20.0 + 4 * input_file_size_gb)

    command {
        set -e
        
        python3 /opt/fastq.process.py3.v0.6.py \
            -a ${input[0]} \
            -b ${input[1]} \
            ${"--c=" + input[2]} \
            ${"--d=" + input[3]} \
            --out ./processed-fastqs/${default="processed" prefix} \
            -y '${config}'
        
        # Compressing the fastqs
        pigz --fast -p 4 ./processed-fastqs/*.fq
        
    }
    
    output {
        # Useful for ATAC and RNA so I don't have to copy twice
        File atac_processed_fastq_R1= glob('./processed-fastqs/*ATAC.R1.fq.gz')
        File atac_processed_fastq_R2= glob('./processed-fastqs/*ATAC.R2.fq.gz')
        File rna_processed_fastq_R1= glob('./processed-fastqs/*RNA.R1.fq.gz')
    }

    runtime {
        cpu : 4
        memory : '${mem_gb} GB'
        disks : 'local-disk ${disk_gb} SSD'
        preemptible: 0 
        docker: docker_image
    }
    
    parameter_meta {
        input: {
                description: 'Read1, Read2 FASTQs. Optional: Index1, Index2',
                help: '[Read1, Read2, Index1, Index2] or [Read1, Read2] depending if your sequencer is Nextseq or Novaseq.',
                example: [
                    'R1',
                    'R2',
                    'I1',
                    'I2'
                ]
            }
        config: {
                description: 'Config YAML file.',
                help: 'Configuration file that contains the adapter information'
                examples: ['https://raw.githubusercontent.com/masai1116/SHARE-seq-alignment/main/config_example.ymal']
            },
        qc: {
                description: 'Boolean used to switch to QC run.'
            },
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz'
                example: ['put link to gcr or dockerhub']
            }
    }


}
