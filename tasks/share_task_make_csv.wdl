version 1.0

# Task 
# Write a CSV of the numbers, files, and logs from a run 

task write_csv {
    meta {
        version: 'v0.1'
        author: 'Michael Shriver (shriverm@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Take in arrays that have the relevant information from a run of the pipeline, and generate a csv file with the name of some output from the pipeline, and     the output itself' 
    }
    
    input {
        # takes in four arrays, one for names and the rest for the data 
        # associated with that name
        String? prefix
        Array[String] names_data
        Array[Int] numeric_data
        Array[String] image_data 
        Array[String] log_data
    }

    String output_file = "${default="share-seq-run-values" prefix}.csv"
    
    command <<< 
        touch names.txt
        touch numeric.txt 
        touch image.txt 
        touch log.txt
        echo ${~{names_data}[@]} >> names.txt
        echo $~{numeric_data}[@] >> numeric.txt
        echo $~{image_data}[@] >> image.txt
        echo $~{log_data}[@] >> log.txt
        python3 /software/write_csv.py ~{output_file} names.txt numeric.txt image.txt log.txt
        >>>

    output {
        File out_file = "~{output_file}"
    }

}