version 1.0

# TASK
# SHARE-html-report
# Gather information from log files


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

        # Stats for ATAC and RNA, will go at top of html from terra table
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

        ## Raw text logs to append to end of html
        Array[String?] log_files

        ##values from tasks
        File? joint_qc_vals
        File? archr_vals

    }

    String output_file = "${default="share-seq" prefix}.html"
    # need to select from valid files since some are optional
    #comment
    Array[File] valid_image_files = select_all(image_files)
    Array[String] valid_log_files = select_all(log_files)
    String output_csv_file = "${default="share-seq" prefix}.csv"

    command <<<

        echo "~{sep="\n" valid_image_files}" > image_list.txt
        echo "~{sep="\n" valid_log_files}" > log_list.txt

        echo "<h3>Summary Statistics</h3><p><table><tr><td colspan=2>ATAC</td></tr><tr><td>Total reads</td><td> id ='ATAC Total reads' " ~{atac_total_reads} "</td></tr>" > output.txt
        echo "<tr><td>Aligned uniquely</td><td> id = 'Aligned uniquely'" ~{atac_aligned_uniquely} "</td></tr>" >> output.txt
        echo "<tr><td>Unaligned</td><td> id = 'ATAC Unaligned'" ~{atac_unaligned} "</td></tr>" >> output.txt
        echo "<tr><td>Unique Reads</td><td> id = 'Unique Reads'" ~{atac_feature_reads} "</td></tr>" >> output.txt
        echo "<tr><td>Duplicate Reads</td><td> id = 'ATAC Duplicate Reads'" ~{atac_duplicate_reads} "</td></tr>" >> output.txt
        echo "<tr><td>Percent Duplicates</td><td> id = 'ATAC Percent Duplicates'" ~{atac_percent_duplicates} "</td></tr>" >> output.txt
        echo "<tr><td>NRF=Distinct/Total</td><td> id = 'NRF=Distinct/Total'" ~{atac_nrf} "</td></tr>" >> output.txt
        echo "<tr><td>PBC1=OnePair/Distinct</td><td> id = 'PBC1=OnePair/Distinct'" ~{atac_pbc1} "</td></tr>" >> output.txt
        echo "<tr><td>PBC2=OnePair/TwoPair</td><td> id = 'PBC2=OnePair/TwoPair'" ~{atac_pbc2} "</td></tr>" >> output.txt
        echo "<td colspan=2>RNA</td></tr><tr><td>Total reads</td><td> id = 'RNA Total reads'" ~{rna_total_reads} "</td></tr>" >> output.txt
        echo "<tr><td>Aligned uniquely</td><td> id = 'RNA Aligned uniquely'" ~{rna_aligned_uniquely} "</td></tr>" >> output.txt
        echo "<tr><td>Aligned multimap</td><td> id = 'Aligned multimap'" ~{rna_aligned_multimap} "</td></tr>" >> output.txt
        echo "<tr><td>Unaligned</td><td> id = 'RNA Unaligned'" ~{rna_unaligned} "</td></tr>" >> output.txt
        echo "<tr><td>Filtered (feature) Reads</td><td> id = 'Filtered (feature) Reads'" ~{rna_feature_reads} "</td></tr>" >> output.txt
        echo "<tr><td>Duplicate Reads</td><td> id = 'RNA Duplicate Reads'" ~{rna_duplicate_reads} "</td></tr>" >> output.txt
        percent=$(( ~{default=0 rna_duplicate_reads}*100/~{default=1 rna_feature_reads} ))
        echo "<tr><td>Percent Duplicates</td><td> id = 'RNA Percent Duplicates'" $percent "</td></tr></table>" >> output.txt
        
        echo "atac_total_reads, " ~{atac_total_reads} > csv_in.txt
        echo "atac_aligned_uniquely, " ~{atac_aligned_uniquely} >> csv_in.txt
        echo "atac_unaligned, " ~{atac_unaligned} >> csv_in.txt
        echo "atac_feature_reads, "~{atac_feature_reads} >> csv_in.txt
        echo "atac_duplicate_reads, " ~{atac_duplicate_reads} >> csv_in.txt
        echo "atac_percent_duplicates, " ~{atac_percent_duplicates} >> csv_in.txt
        echo "atac_nrf, " ~{atac_nrf} >> csv_in.txt
        echo "atac_pbc1, " ~{atac_pbc1} >> csv_in.txt
        echo "atac_pbc2, " ~{atac_pbc2} >> csv_in.txt
        echo "rna_total_reads, "~{rna_total_reads} >> csv_in.txt
        echo "rna_aligned_uniquely, " ~{rna_aligned_uniquely} >> csv_in.txt
        echo "rna_aligned_multimap, " ~{rna_aligned_multimap} >> csv_in.txt
        echo "rna_unaligned, " ~{rna_unaligned} >> csv_in.txt
        echo "rna_feature_reads, " ~{rna_feature_reads} >> csv_in.txt
        echo "rna_duplicate_reads, " ~{rna_duplicate_reads} >> csv_in.txt
        

        PYTHONIOENCODING=utf-8 python3 /software/write_html.py ~{output_file} image_list.txt log_list.txt ~{output_csv_file} csv_in.txt ~{joint_qc_vals} ~{archr_vals} --input_file_name output.txt
    
    

    
    >>>
    output {
        File html_report_file = "~{output_file}"
        File csv_summary_file = "~{output_csv_file}"
    }

    runtime {
        #docker: 'nchernia/share_task_html_report:14'
        docker: 'mshriver01/share_task_html_report:latest'
    }
}
