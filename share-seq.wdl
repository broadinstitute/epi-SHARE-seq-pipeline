version 1.0

# Import the sub-workflow for preprocessing the fastqs.
import "tasks/subwf_preprocess.wdl" as step1

# WDLworkflow for SHARE-seq

workflow ShareSeq {
    input {
        # Preprocess inputs
        Array[File] R1
        Array[File] R2
        Array[File]? I1
        Array[File]? I2
        Boolean qc = false
        String atac_primers
        String rna_primers
        String prefix = "shareseq-project"
        Int preprocess_cpus = 4
    }
    
    # TODO: find a better way to do the next 3 checks
    if ( length(I1) != length(I2) ) {
        call raise_exception as error_input_data { 
            input:
            msg = 'Length of inputs for Index_1 and Index_2 are different. Please check your inputs.'
        }
    }
    
    if ( length(R1) != length(R2) ) {
        call raise_exception as error_input_data { 
            input:
            msg = 'Length of inputs for Read_1 and Read_2 are different. Please check your inputs.'
        }
    }
    
    if ( length(I1) > 0 ) {
        if ( length(R1) != length(I1) ) {
            call raise_exception as error_input_data { 
                input:
                msg = 'Different numbers of reads and indexes. Please check your inputs.'
            }
        }
    }
    
    # Preprocess the input fastqs lane by lane.
    scatter( i in range(length(R1)) ){
      # All of this workaround is to make the optional I1,I2 work without
      # asking the user to input a complicated structure like
      # Array[Array[File]].
      
      # I cannot access or use 'length' on Array[File]?
      # As a workaround I am first converting Array[File]? into Array[File].
      Array[File] I1_not_optional= if defined(I1) then select_first([I1]) else []
      Array[File] I2_not_optional= if defined(I2) then select_first([I2]) else []
      
      # Now that I have a type on which use length, I can check if the
      # file exists and assign it 
      if (length(I1_not_optional) >= i ) {
          File? I1_tmp = I1_not_optional[i]
          File? I2_tmp = I2_not_optional[i]
      }
      
      # Call the preprocess step for each lane
      call step1.wf_preprocess {
        input:
          R1 = R1[i],
          R2 = R2[i],
          I1 = I1_tmp,
          I2 = I2_tmp,
          qc = qc,
          atac_primers = atac_primers,
          rna_primers = rna_primers,
          prefix = prefix,
          cpus = preprocess_cpus
      }
    }
    
    

    output {
        Array[File] output1 = wf_preprocess.atac_processed_fastq_R1
        Array[File] output2 = wf_preprocess.atac_processed_fastq_R2
        Array[File] output3 = wf_preprocess.rna_processed_fastq_R1
    }
    
}


# Task to report errors to user.
# From https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/chip.wdl
task raise_exception {
  input {
      String msg
      Array[String]? vals
  }
  command {
      echo -e "\n* Error: ${msg}\n" >&2
      echo -e "* Vals: ${sep=',' vals}\n" >&2
      exit 2
  }
  output {
      String error_msg = '${msg}'
  }
  runtime {
      maxRetries : 0
      cpu : 1
      memory : '2 GB'
      time : 1
      disks : 'local-disk 10 SSD'
  }
}
