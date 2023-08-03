version 1.0

# TASK
# SHARE-html-report
# Gather information from log files

#updating this to also generate the file with the data

task html_report {
    meta {
        version: 'v0.1'
        author: 'Neva C. Durand (neva@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: create html report task'
    }

    input {
        # This function takes as input the files to append to the report
        # and the metrics and writes out an html file

        String? prefix

        # Stats for ATAC and RNA, will go at top of html
        Int? atac_total_reads
        Int? atac_aligned_uniquely
        Int? atac_unaligned
        Int? atac_feature_reads
        Int? atac_duplicate_reads
        Float? atac_nrf
        Float? atac_pbc1
        Float? atac_pbc2
        Float? atac_percent_duplicates
        Int? rna_total_reads
        Int? rna_aligned_uniquely
        Int? rna_aligned_multimap
        Int? rna_unaligned
        Int? rna_feature_reads
        Int? rna_duplicate_reads
       
        ## JPEG files to be encoded and appended to html
        Array[File?] image_files
        Array[String?] image_file_paths
        ## Raw text logs to append to end of html
        Array[String?] log_files

        File? joint_barcode_stats

    }

    String output_file = "${default="share-seq" prefix}.html"
    # need to select from valid files since some are optional
    Array[File] valid_image_files = select_all(image_files)
    Array[String] valid_log_files = select_all(log_files)
    # new stuff im adding for the dictionary functionality 
    String other_output_file = "${default="share-seq-data" prefix}.csv"
    
    # made these strings local variables so the files that are generated from 
    # each of them can be pipeline outputs for debugging purposes
    String names_for_csv = "${default="share-seq-data" prefix}names.txt"
    String images_for_csv = "${default="share-seq-data" prefix}images.txt"
    String logs_for_csv = "${default="share-seq-data" prefix}logs.txt"
    Array[String] names_numeric_fields = ['atac_total_reads', 'atac_aligned_uniquely','atac_unaligned', 'atac_feature_reads','atac_duplicate_reads', 'atac_percent_duplicates', 'atac_nrf','atac_pbc1','atac_pbc2', 'rna_total_reads','rna_aligned_uniquely','rna_aligned_multimap','rna_unaligned','rna_feature_reads','rna_duplicate_reads', 'rna_percent_duplicates']
    String? joint_stats = joint_barcode_stats
    command <<<

        echo "~{sep="\n" valid_image_files}" > ~{images_for_csv}
        echo "~{sep="\n" valid_log_files}" > ~{logs_for_csv}
        echo "~{sep="\n" names_numeric_fields}" > ~{names_for_csv}
        echo "~{sep="\n" valid_image_files}" >> ~{names_for_csv}
        echo "~{sep="\n" valid_log_files}" >> ~{names_for_csv}
        
        echo ~{atac_total_reads}
        echo ~{atac_aligned_uniquely} 
        echo ~{atac_unaligned} 
        echo ~{atac_feature_reads} 
        echo ~{atac_duplicate_reads}
        echo ~{atac_percent_duplicates} 
        echo ~{atac_nrf} 
        echo ~{atac_pbc1}
        echo ~{atac_pbc2}
        echo ~{rna_total_reads} 
        echo ~{rna_aligned_uniquely}
        echo ~{rna_aligned_multimap} 
        echo ~{rna_unaligned}
        echo ~{rna_feature_reads}
        echo ~{rna_duplicate_reads}
        echo "end of outputs"

        
        echo ~{atac_total_reads} >> output.txt
        echo ~{atac_aligned_uniquely} >> output.txt
        echo ~{atac_unaligned} >> output.txt
        echo ~{atac_feature_reads} >> output.txt
        echo ~{atac_duplicate_reads} >> output.txt
        echo ~{atac_percent_duplicates} >> output.txt
        echo ~{atac_nrf} >> output.txt
        echo ~{atac_pbc1} >> output.txt
        echo ~{atac_pbc2} >> output.txt
        echo ~{rna_total_reads} >> output.txt
        echo ~{rna_aligned_uniquely} >> output.txt
        echo ~{rna_aligned_multimap} >> output.txt
        echo ~{rna_unaligned} >> output.txt
        echo ~{rna_feature_reads} >> output.txt
        echo ~{rna_duplicate_reads} >> output.txt

        percent=$(( ~{default=0 rna_duplicate_reads}*100/~{default=1 rna_feature_reads} ))
        echo $percent >> output.txt
        PYTHONIOENCODING=utf-8 python3 /software/write_html.py ~{output_file} ~{images_for_csv} output.txt ~{logs_for_csv} ~{joint_barcode_stats}
        
        python3 /software/write_csv.py ~{other_output_file} ~{names_for_csv} output.txt ~{images_for_csv} ~{logs_for_csv}
    >>>
    
    output {
        File html_report_file = "~{output_file}"
        File csv_report_file = "~{other_output_file}"

    }

    runtime {
        #docker: 'nchernia/share_task_html_report:14'
        docker: 'mshriver01/share_task_html_report:latest'
    }
}
