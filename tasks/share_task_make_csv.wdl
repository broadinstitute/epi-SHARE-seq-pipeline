version 1.0

# Task 
# Write a CSV of the numbers, files, and logs from a run 

task write_csv {
    meta {
        version: 'v0.1'
        author: 'Michael Shriver (shriverm@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Take in arrays that have the relevant information from a run of the pipeline,
                      and generate a csv file with the name of some output from the pipeline, and 
                      the output itself' 
    }
    
    input {
        # takes in four arrays, one for names and the rest for the data 
        # associated with that name
        Array[String] names_data
        Array[Int] numeric_data
        Array[String] inage_data 
        Array[String] log_data
    }

    command <<< 
        touch run_data.csv
        >>>

}