version 1.0

# TASK
# SHARE-html-report
# Gather information from log files


task html_report {
    meta {
        version: 'v1.0'
        author: 'Neva C. Durand (neva@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: create html report task'
    }

    input {
        # This function takes as input the files to append to the report
        # and the metrics and writes out an html file

        String? prefix

        File? atac_metrics
        File? rna_metrics

        
        Int? rna_total_reads
        Int? rna_aligned_uniquely
        Int? rna_aligned_multimap
        Int? rna_unaligned
        Int? rna_feature_reads
        Int? rna_duplicate_reads
        Float? rna_frig

        ## JPEG files to be encoded and appended to html
        Array[File?] image_files

        ## Raw text logs to append to end of html
        Array[String?] log_files

        String docker_image = 'us.gcr.io/buenrostro-share-seq/task_html_report:dev'

    }

    String output_file = "${default="combinomics" prefix}.html"
    # need to select from valid files since some are optional
    Array[File]? valid_image_files = select_all(image_files)
    Array[String]? valid_log_files = select_all(log_files)
    String output_csv_file = "${default="share-seq" prefix}.csv"

    command <<<

        echo "~{sep="\n" valid_image_files}" > image_list.txt
        echo "~{sep="\n" valid_log_files}" > log_list.txt
        
        echo "rna_total_reads," ~{if defined(rna_total_reads) then rna_total_reads else 0} > rna_metrics.csv
        echo "rna_aligned_uniquely," ~{if defined(rna_aligned_uniquely) then rna_aligned_uniquely else 0} >> rna_metrics.csv
        echo "rna_aligned_multimap," ~{if defined(rna_aligned_multimap) then rna_aligned_multimap else 0} >> rna_metrics.csv
        echo "rna_unaligned," ~{if defined(rna_unaligned) then rna_unaligned else 0} >> rna_metrics.csv
        echo "rna_feature_reads," ~{if defined(rna_feature_reads) then rna_feature_reads else 0} >> rna_metrics.csv
        echo "rna_duplicate_reads," ~{if defined(rna_duplicate_reads) then rna_duplicate_reads else 0} >> rna_metrics.csv
        echo "rna_frig," ~{if defined(rna_frig) then rna_frig else 0} >> rna_metrics.csv
        
        # This can be uncommented once we get the rna_metrics file sorted out
        #PYTHONIOENCODING=utf-8 python3 /software/write_html.py ~{output_file} image_list.txt log_list.txt ~{output_csv_file} ~{atac_metrics} ~{rna_metrics}
        PYTHONIOENCODING=utf-8 python3 /software/write_html.py \
        ~{output_file} \
        image_list.txt \
        log_list.txt \
        ~{output_csv_file} \
        --rna_metrics rna_metrics.csv \
        ~{if defined(atac_metrics) then '--atac_metrics ~{atac_metrics}' else ''}

        # TODO: this needs to be fix.
        echo 'PKR,~{prefix}' | cat - ~{output_csv_file} > temp && mv temp ~{output_csv_file} 
    >>>
    output {
        File html_report_file = "~{output_file}"
        File csv_summary_file = "~{output_csv_file}"
    }

    runtime {
        docker: "${docker_image}"
    }
}
