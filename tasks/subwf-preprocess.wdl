version 1.0

workflow wf_preprocess {
    
  input {
    # Preprocess inputs
    File read1
    File read2
    File? index1
    File? index2
    File yaml
    Boolean? qc
    String prefix = "shareseq-project"
    String output_table_name
    String terra_project
    String workspace_name
    Int? cpus = 4
  }

  call preprocess {
    input:
      R1 = read1,
      R2 = read2,
      I1 = index1,
      I2 = index2,
      yaml = yaml,
      qc = qc,
      prefix = prefix      
  }
  
  call gather_outputs{
    input:
      yaml = yaml,
      read_paths = preprocess.demux_fastqs,
      new_table_name = output_table_name
  }

  call entities_batch_upsert{
    input:
      tsv_file = gather_outputs.load_tsv,
      terra_project = terra_project,
      workspace_name = workspace_name
  }

  output {
    File response = entities_batch_upsert.upsert_entities_response
    File discarded_R1 = preprocess.discarded_fastq_R1
    File discarded_R2 = preprocess.discarded_fastq_R2
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
    File yaml
    Boolean? qc
    Int cpus= 4
    String? prefix
    String docker_image = "polumechanos/share-seq-undetermined"
  }
  
  Float input_file_size_gb = size(R1, "G")
  Float mem_gb = 8.0
  Int disk_gb = round(20.0 + 6 * input_file_size_gb)

  command {
    set -e
    
    mkdir out
    mkdir discard

    python3 $(which fastq.process.py3.v0.7.py) \
      -a ${R1} \
      -b ${R2} \
      ${if defined(I1) then "--c " + I1 else ""} \
      ${if defined(I2) then "--d " + I2 else ""} \
      ${if defined(qc) then "--qc" else ""} \
      -y ${yaml} \
      --out ${default="shareseq-project." prefix+"."}preprocessed
    
  }
    
  output {
    Array[File] demux_fastqs = glob('out/*.fq.gz')
    File discarded_fastq_R1 = glob('discard/*R1.fq.gz')[0]
    File discarded_fastq_R2 = glob('discard/*R2.fq.gz')[0]
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
    yaml: {
      description: 'YAML file',
      help: 'YAML-formatted file containing the metadata describing the experiment.'
    }
    prefix: {
      description: 'Prefix for output.',
      help: 'String to use as prefix for the outputs',
      examples: "my-experiment"
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
  


task gather_outputs {
    input {
        Array[String] read_paths
        File yaml
        String new_table_name
        
        
    }

    command <<<
        python3 /software/gather_outputs.py -n ~{new_table_name} -y ~{yaml} -p ~{sep="," read_paths}
    >>>

    runtime {
        docker: "polumechanos/share-seq-undetermined"

    }

    output {
        File load_tsv = "results_load_table.tsv"
    }
}

task entities_batch_upsert {
  input {
    File      tsv_file
    String    terra_project
    String    workspace_name
  }
    String    docker="polumechanos/update_terra_table"

  command {
    set -e

    python3 /software/batch_upsert_entities_standard.py \
        -t "~{tsv_file}" \
        -p "~{terra_project}" \
        -w "~{workspace_name}"
  }

  runtime {
    docker : docker
    memory: "2 GB"
    cpu: 1
  }

  output {
    String upsert_entities_response = stdout()
  }
}
