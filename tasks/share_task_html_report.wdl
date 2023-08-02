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

        # Stats for ATAC and RNA, will go at top of html
        Int? atac_total_reads
        Int? atac_aligned_uniquely
        Int? atac_unaligned
        Int? atac_feature_reads
        Int? atac_duplicate_reads
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
    # need to select from valid files since some are optional
    Array[File] valid_image_files = select_all(image_files)
    Array[String] valid_log_files = select_all(log_files)

    command <<<

        echo "~{sep="\n" valid_image_files}" > image_list.txt
        echo "~{sep="\n" valid_log_files}" > log_list.txt

        echo "<h3>Summary Statistics</h3><p><table><tr><td colspan=2>ATAC</td></tr><tr><td>Total reads</td><td>" ~{atac_total_reads} "</td></tr>" > output.txt
        echo "<tr><td>Aligned uniquely</td><td>" ~{atac_aligned_uniquely} "</td></tr>" >> output.txt
        echo "<tr><td>Unaligned</td><td>" ~{atac_unaligned} "</td></tr>" >> output.txt
        echo "<tr><td>Filtered Reads</td><td>" ~{atac_feature_reads} "</td></tr>" >> output.txt
        echo "<tr><td>Duplicate Reads</td><td>" ~{atac_duplicate_reads} "</td></tr>" >> output.txt
        percent=$(( ~{default=0 atac_duplicate_reads}*100/~{default=1 atac_feature_reads} ))
        echo "<tr><td>Percent Duplicates</td><td>" $percent "</td></tr>" >> output.txt
        echo "<td colspan=2>RNA</td></tr><tr><td>Total reads</td><td>" ~{rna_total_reads} "</td></tr>" >> output.txt
        echo "<tr><td>Aligned uniquely</td><td>" ~{rna_aligned_uniquely} "</td></tr>" >> output.txt
        echo "<tr><td>Aligned multimap</td><td>" ~{rna_aligned_multimap} "</td></tr>" >> output.txt
        echo "<tr><td>Unaligned</td><td>" ~{rna_unaligned} "</td></tr>" >> output.txt
        echo "<tr><td>Filtered (feature) Reads</td><td>" ~{rna_feature_reads} "</td></tr>" >> output.txt
        echo "<tr><td>Duplicate Reads</td><td>" ~{rna_duplicate_reads} "</td></tr>" >> output.txt
        percent=$(( ~{default=0 rna_duplicate_reads}*100/~{default=1 rna_feature_reads} ))
        echo "<tr><td>Percent Duplicates</td><td>" $percent "</td></tr></table>" >> output.txt
        PYTHONIOENCODING=utf-8 python3 /software/write_html.py output.html image_list.txt log_list.txt --input_file_name output.txt
    >>>
    output {
        File html_report_file = glob('*.html')[0]
    }

    runtime {
        docker: 'us.gcr.io/buenrostro-share-seq/share_task_html_report:v0.0.1'
    }
}
