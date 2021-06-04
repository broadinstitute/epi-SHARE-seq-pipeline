version 1.0

# Import the sub-workflow for preprocessing the fastqs.
import "tasks/subwf_preprocess.wdl" as step1

# WDLworkflow for SHARE-seq

workflow ShareSeq {
    input {
        # Preprocess inputs
        Array[File] R1
        Array[File] R2
        # The index files are not requested but putting them as optional creates a lot of problem in the logic of WDL.
        Array[File] I1 = []
        Array[File] I2 = []
        Boolean? qc = false
        String atac_primers
        String rna_primers
        String? prefix
        Int preprocess_cpus = 4
        # WDL does not allow to assign None|null|nill if I want to explicitly pass an empty value so I am creating it
        # myself and use this as a substitute. Look inside the wf_preprocess call to see how I am using it.
        File? empty
    }
    
    # Check if the number of reads in input is correct
    if ( length(R1) != length(R2) ) {
        call raise_exception as error_input_data { 
            input:
                msg = 'Length of inputs for Read_1 and Read_2 are different. Please check your inputs.'
        }
    }
    
    
    # Preprocess the input fastqs lane by lane.
    scatter( i in range(length(R1)) ){
        # Call the preprocess step for each lane
        call step1.wf_preprocess {
            input:
                read1 = R1[i],
                read2 = R2[i],
                # If the user passed the index files than get the file from the array. If not the task expect an 
                # optional value. WDL does not allow to use None, so I created an optional values that I use instead
                # of that. TODO: Is there a bettere way to do?
                index1 = if length(I1) > 0 then I1[i] else empty,
                index2 = if length(I2) > 0 then I2[i] else empty,
                qc = qc,
                atac_primers = atac_primers,
                rna_primers = rna_primers,
                prefix = prefix,
                cpus = preprocess_cpus
        }
    }
    
    # TODO: If I just have one lane is there a way to skip this?
    call merge_fastqs{
        input:
            atac_read1 = wf_preprocess.atac_processed_fastq_R1,
            atac_read2 = wf_preprocess.atac_processed_fastq_R2,
            rna_read1 = wf_preprocess.rna_processed_fastq_R1,
            prefix = prefix
    }
    

    output {
        File output1 = merge_fastqs.merged_atac_fastq_R1
        File output2 = merge_fastqs.merged_atac_fastq_R2
        File output3 = merge_fastqs.merged_rna_fastq_R1
    }
    
}


task merge_fastqs{
    input{
        Array[File] atac_read1
        Array[File] atac_read2
        Array[File] rna_read1
        String? prefix
    }
    command{
        cat ${sep=' ' atac_read1} > ${prefix + "."}merged.atac.R1.fq.gz
        cat ${sep=' ' atac_read2} > ${prefix + "."}merged.atac.R2.fq.gz
        cat ${sep=' ' rna_read1} > ${prefix + "."}merged.rna.R1.fq.gz
    }
    output{
        File merged_atac_fastq_R1 = glob('*.merged.atac.R1.fq.gz')[0]
        File merged_atac_fastq_R2 = glob('*.merged.atac.R2.fq.gz')[0]
        File merged_rna_fastq_R1 = glob('*.merged.rna.R1.fq.gz')[0]
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
