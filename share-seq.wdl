version 1.0

import tasks/preprocess.wdl
import tasks/raise_exception.wdl

# WDLworkflow for SHARE-seq

workflow ShareSeq {
    input {
        # Preprocess inputs
        Array[File] R1
        Array[File] R2
        Array[File]? I1
        Array[File]? I2        
        Boolean qc= false
        String atac_primers
        String rna_primers
        String prefix= "shareseq-project"
        Int preprocess_cpus= 4
    }
    
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
    
    call preprocess {
        input:
            R1= R1,
            R2= R2,
            I1= I1,
            I2= I2,
            qc= qc,
            atac_primers= atac_primers,
            rna_primers= rna_primers,
            cpus= preprocess_cpus
    }

    output {
        File output1= preprocess.atac_processed_fastq_R1,
        File output2= preprocess.atac_processed_fastq_R2,
        File output3= rna_processed_fastq_R1
    }
    
}
