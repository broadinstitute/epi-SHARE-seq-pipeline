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
        
        # TODO: Create a log file inside the python script that can be
        # reported to the user
        
        File R1
        File R2
        File? I1
        File? I2
        Boolean? qc= false
        Int cpus= 4
        String? prefix
        String rna_primers
        String atac_primers
        String docker_image
        
    }
    
    Float input_file_size_gb = size(input[0], "G")
    Float mem_gb = 8.0
    Int disk_gb = round(20.0 + 4 * input_file_size_gb)

    command {
        set -e
        
        python3 /opt/fastq.process.py3.v0.6.py \
            #-a ${input[0]} \
            #-b ${input[1]} \
            #${if length(input) == 4 then "--c=" + input[2] else ""} \
            #${if length(input) == 4 then "--d=" + input[3] else ""} \
            -a ${R1} \
            -b ${R2} \
            ${"--c" + I1} \
            ${"--d" + I2} \
            ${if qc then "--qc"} \
            --rna_primers ${rna_primers} \
            --atac_primers ${atac_primers} \
            --out ./processed-fastqs/"${default="shareseq-project." prefix+"."}preprocessed
        
        # Compressing the fastqs
        pigz --fast -p ${cpus} ./processed-fastqs/*.fq
        
    }
    
    output {
        # Useful for ATAC and RNA so I don't have to copy twice
        File atac_processed_fastq_R1= glob('./processed-fastqs/*atac.R1.fq.gz')[0]
        File atac_processed_fastq_R2= glob('./processed-fastqs/*atac.R2.fq.gz')[0]
        File rna_processed_fastq_R1= glob('./processed-fastqs/*rna.R1.fq.gz')[0]
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
            },
        rna_primers: {
                description: 'Primers used for RNA.',
                help: 'Comma separated list of primers used for the RNA part of SHARE-seq'
                examples: ["P1.52,P1.53,P1.54,P1.55,P1.56"]
            },
        atac_primers: {
            description: 'Primers used for ATAC.',
            help: 'Comma separated list of primers used for the ATAC part of SHARE-seq'
            examples: ["P1.04,P1.05,P1.06,P1.07,P1.08"]
        },
        qc: {
                description: 'Boolean used to switch to QC run.'
                default: false
                example: [true, false]
            },
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz'
                example: ['put link to gcr or dockerhub']
            }
    }


}
