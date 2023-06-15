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

        ## Raw text logs to append to end of html
        Array[String?] log_files

    }

    String output_file = "${default="share-seq" prefix}.html"
    # need to select from valid files since some are optional
    Array[File] valid_image_files = select_all(image_files)
    Array[String] valid_log_files = select_all(log_files)

    # TODO: command needs fixing in two ways. Firstly, all of the echoes are printing formatted summary stats to an output file, that the write_html script then prints the 
    # rest of the report into. The stats aren't all that well formated. Two thoughts for this. Firstly, the names and numbers could be written into a text file that is then 
    # passed to write_html as an argument, so write_html can format the stats well. Secondly, could see what happens if the wdl list is passed to the command line, but that 
    # seems unlikely and would require some weird stuff with lists given the need to move strings, ints, and floats. Second thing that needs to happen is that, in the last line, 
    # the call needs to be modified to work with the amount of parameters. Also inputs will need to be updated to reflect the three lists of image files moving around
    
    command <<<

        echo "~{sep="\n" valid_image_files}" > image_list.txt
        echo "~{sep="\n" valid_log_files}" > log_list.txt

        echo "<h3>Summary Statistics</h3><p><table><tr><td colspan=2>ATAC</td></tr><tr><td>Total reads</td><td>" ~{atac_total_reads} "</td></tr>" > output.txt
        echo "<tr><td>Aligned uniquely</td><td>" ~{atac_aligned_uniquely} "</td></tr>" >> output.txt
        echo "<tr><td>Unaligned</td><td>" ~{atac_unaligned} "</td></tr>" >> output.txt
        echo "<tr><td>Unique Reads</td><td>" ~{atac_feature_reads} "</td></tr>" >> output.txt
        echo "<tr><td>Duplicate Reads</td><td>" ~{atac_duplicate_reads} "</td></tr>" >> output.txt
        echo "<tr><td>Percent Duplicates</td><td>" ~{atac_percent_duplicates} "</td></tr>" >> output.txt
        echo "<tr><td>NRF=Distinct/Total</td><td>" ~{atac_nrf} "</td></tr>" >> output.txt
        echo "<tr><td>PBC1=OnePair/Distinct</td><td>" ~{atac_pbc1} "</td></tr>" >> output.txt
        echo "<tr><td>PBC2=OnePair/TwoPair</td><td>" ~{atac_pbc2} "</td></tr>" >> output.txt
        echo "<td colspan=2>RNA</td></tr><tr><td>Total reads</td><td>" ~{rna_total_reads} "</td></tr>" >> output.txt
        echo "<tr><td>Aligned uniquely</td><td>" ~{rna_aligned_uniquely} "</td></tr>" >> output.txt
        echo "<tr><td>Aligned multimap</td><td>" ~{rna_aligned_multimap} "</td></tr>" >> output.txt
        echo "<tr><td>Unaligned</td><td>" ~{rna_unaligned} "</td></tr>" >> output.txt
        echo "<tr><td>Filtered (feature) Reads</td><td>" ~{rna_feature_reads} "</td></tr>" >> output.txt
        echo "<tr><td>Duplicate Reads</td><td>" ~{rna_duplicate_reads} "</td></tr>" >> output.txt
        percent=$(( ~{default=0 rna_duplicate_reads}*100/~{default=1 rna_feature_reads} ))
        echo "<tr><td>Percent Duplicates</td><td>" $percent "</td></tr></table>" >> output.txt
        PYTHONIOENCODING=utf-8 python3 /software/write_html.py ~{output_file} image_list.txt log_list.txt --input_file_name output.txt
    >>>
    output {
        File html_report_file = "~{output_file}"
    }

    runtime {
        docker: 'nchernia/share_task_html_report:14'
    }
}
