version 1.0

workflow wf_preprocess {
  input {
    # Preprocess inputs
    File read1
    File read2
    File? index1
    File? index2
    Boolean? qc
    String atac_primers
    String rna_primers
    String prefix = "shareseq-project"
    Int? cpus = 4
  }

  call preprocess {
    input:
      R1 = read1,
      R2 = read2,
      I1 = index1,
      I2 = index2,
      qc = qc,
      atac_primers = atac_primers,
      rna_primers = rna_primers,
      prefix = prefix      
  }

  output {
    File atac_processed_fastq_R1 = preprocess.atac_processed_fastq_R1
    File atac_processed_fastq_R2 = preprocess.atac_processed_fastq_R2
    File rna_processed_fastq_R1 = preprocess.rna_processed_fastq_R1
  }
}

task preprocess {
  meta {
    version: 'v0.1'
    author: 'Eugenio Mattei (emattei@broadinstitute.org) and Sai Ma @ Broad Institute of MIT and Harvard'
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
    Boolean? qc = false
    Int cpus= 4
    String? prefix
    String rna_primers
    String atac_primers
    String docker_image = "polumechanos/share-seq"
  }
  
  Float input_file_size_gb = size(R1, "G")
  Float mem_gb = 8.0
  Int disk_gb = round(20.0 + 4 * input_file_size_gb)

  command {
    set -e
    
    python3 $(which fastq.process.py3.v0.6.py) \
      -a ${R1} \
      -b ${R2} \
      ${if defined(I1) then "--c " + I1 else ""} \
      ${if defined(I2) then "--d " + I2 else ""} \
      ${if defined(qc) then "--qc" else ""} \
      --rna_primers ${rna_primers} \
      --atac_primers ${atac_primers} \
      --out ${default="shareseq-project." prefix+"."}preprocessed
        
    # Compressing the fastqs
    pigz --fast -p ${cpus} *.fq
    
  }
    
  output {
    # Useful for ATAC and RNA so I don't have to copy twice
    File atac_processed_fastq_R1= glob('*atac.R1.fq.gz')[0]
    File atac_processed_fastq_R2= glob('*atac.R2.fq.gz')[0]
    File rna_processed_fastq_R1= glob('*rna.R1.fq.gz')[0]
  }

  runtime {
    cpu : 4
    memory : '${mem_gb} GB'
    disks : 'local-disk ${disk_gb} SSD'
    preemptible: 0 
    docker: docker_image
  }
  
  parameter_meta {
    R1: {
      description: 'Read1 fastq',
      help: 'Array containing read1 for different lanes.',
      example: ['R1.L1', 'R1.L2', 'R1.L3', 'R1.L4']
    }
    R2: {
      description: 'Read2 fastq',
      help: 'Array containing read2 for different lanes.',
      example: ['R2.L1', 'R2.L2', 'R2.L3', 'R2.L4']
    }
    I1: {
      description: 'Index1 FASTQ',
      help: 'Optional array containing first index file (check the correct name for this) for different lanes.',
      example: ['I1.L1', 'I1.L2', 'I1.L3', 'I1.L4']
    }
    I2: {
      description: 'Index2 FASTQ',
      help: 'Optional array containing first index file (check the correct name for this) for different lanes.',
      example: ['I2.L1', 'I2.L2', 'I2.L3', 'I2.L4']
    }
    prefix: {
      description: 'Prefix for output.',
      help: 'String to use as prefix for the outputs',
      examples: "my-experiment"
    }
    rna_primers: {
      description: 'Primers used for RNA.',
      help: 'Comma separated list of primers used for the RNA part of SHARE-seq',
      examples: ["P1.52,P1.53,P1.54,P1.55,P1.56"]
    }
    atac_primers: {
      description: 'Primers used for ATAC.',
      help: 'Comma separated list of primers used for the ATAC part of SHARE-seq',
      examples: ["P1.04,P1.05,P1.06,P1.07,P1.08"]
    }
    qc: {
      description: 'Boolean used to switch to QC run.',
      default: false,
      example: [true, false]
    }
    docker_image: {
      description: 'Docker image.',
      help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein Bio; apt install pigz.',
      example: ['put link to gcr or dockerhub']
    }
  }
}
