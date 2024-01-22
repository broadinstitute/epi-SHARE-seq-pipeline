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
        
        # This can be uncommented once we get the rna_metrics file sorted out
        #PYTHONIOENCODING=utf-8 python3 /software/write_html.py ~{output_file} image_list.txt log_list.txt ~{output_csv_file} ~{atac_metrics} ~{rna_metrics}
        PYTHONIOENCODING=utf-8 python3 /software/write_html.py \
        ~{output_file} \
        image_list.txt \
        log_list.txt \
        ~{output_csv_file} \
        ~{if defined(rna_metrics) then '--rna_metrics ~{rna_metrics}' else ''} \
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
