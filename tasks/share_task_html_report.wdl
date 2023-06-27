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

        ## Raw text logs to append to end of html
        Array[String?] log_files


    }

    String output_file = "${default="share-seq" prefix}.html"
    # need to select from valid files since some are optional
    Array[File] valid_image_files = select_all(image_files)
    Array[String] valid_log_files = select_all(log_files)
    String stuff = "push through changes"

    # new stuff im adding for the dictionary functionality

    # try to implicity coerce them all to strings here 
    String other_output_file = "${default="share-seq-data" prefix}.csv"
    #issue with this line at the moment, commenting it out to see if hard coded numbers are ok
    Array[String] names_of_data = ['atac_total_reads', 'atac_aligned_uniquely','atac_unaligned', 'atac_feature_reads','atac_duplicate_reads','atac_nrf','atac_pbc1','atac_pbc2','atac_percent_duplicates','rna_total_reads','rna_aligned_uniquely','rna_aligned_multimap','rna_unaligned','rna_feature_reads','rna_duplicate_reads', 'atac_total_reads', 'atac_aligned_uniquely', 'atac_unaligned', 'atac_feature_reads', 'atac_duplicate_reads', 'atac_nrf', 'atac_pbc1', 'atac_pbc2', 'atac_percent_duplicates', 'rna_total_reads', 'rna_aligned_uniquely', 'rna_aligned_multimap', 'rna_unaligned', 'rna_feature_reads', 'rna_duplicate_reads', 'joint_qc_plot', 'joint_density_plot', 'share_rna_umi_barcode_rank_plot', 'share_rna_gene_barcode_rank_plot', 'share_rna_gene_umi_scatter_plot', 'share_rna_seurat_raw_violin_plot', 'share_rna_seurat_raw_qc_scatter_plot', 'share_rna_seurat_filtered_violin_plot', 'share_rna_seurat_filtered_qc_scatter_plot', 'share_rna_seurat_variable_genes_plot', 'share_rna_seurat_PCA_dim_loadings_plot', 'share_rna_seurat_PCA_plot', 'share_rna_seurat_heatmap_plot', 'share_rna_seurat_jackstraw_plot', 'share_rna_seurat_elbow_plot', 'share_rna_seurat_umap_cluster_plot', 'share_rna_seurat_umap_rna_count_plot', 'share_rna_seurat_umap_gene_count_plot', 'share_rna_seurat_umap_mito_plot', 'share_atac_qc_barcode_rank_plot', 'share_atac_qc_hist_plot', 'share_atac_qc_tss_enrichment', 'share_atac_archr_raw_tss_enrichment', 'share_atac_archr_filtered_tss_enrichment', 'share_atac_archr_raw_fragment_size_plot', 'share_atac_archr_filtered_fragment_size_plot', 'share_atac_archr_umap_cluster', 'share_atac_archr_umap_num_frags_plot', 'share_atac_archr_umap_tss_score_plot', 'share_atac_archr_umap_frip_plot', 'share_atac_archr_raw_tss_by_unique_frags', 'share_atac_archr_filtered_tss_by_unique_frags', 'share_atac_archr_unfiltered_tss_by_unique_frags', 'share_atac_archr_filtered_tss_by_unique_frags','share_atac_archr_umap_cluster_plot', 'share_atac_archr_strict_umap_num_frags_plot', 'share_atac_archr_strict_umap_tss_score_plot', 'share_atac_archr_strict_umap_frip_plot', 'share_rna_alignment_log',  'share_task_starsolo_barcodes_stats', 'share_task_starsolo_features_stats', 'share_task_starsolo_summary_csv', 'share_task_starsolo_umi_per_cell', 'share_task_starsolo_raw_tar','share_rna_seurat_notebook_log', 'share_atac_alignment_log', 'share_atac_archr_notebook_log', 'dorcs_notebook_log']
    Array[Int] nums_to_save = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    Array[Int?] first_ints = [atac_total_reads, atac_aligned_uniquely, atac_unaligned, atac_feature_reads]
    Array[Float?] first_floats = [atac_nrf, atac_pbc1, atac_pbc2, atac_percent_duplicates]
    Array[Int?] second_ints = [rna_total_reads, rna_aligned_uniquely, rna_aligned_multimap, rna_unaligned, rna_feature_reads, rna_duplicate_reads]
    # added echo command to get the names and the numbers in bash, and added 
    # a call to the make csv python script with the relevant files
    # comment to push changes through
    command <<<

        echo "~{sep="\n" names_of_data}" > names_list.txt
        echo "~{sep="\n" valid_image_files}" > image_list.txt
        echo "~{sep="\n" valid_log_files}" > log_list.txt
        echo "~{sep="\n" nums_to_save}" > numeric_list.txt
        


        
        echo ~{atac_total_reads} '/n' >> output.txt
        echo ~{atac_aligned_uniquely} '/n' >> output.txt
        echo ~{atac_unaligned} '/n' >> output.txt
        echo ~{atac_feature_reads} '/n' >> output.txt
        echo ~{atac_duplicate_reads} '/n' >> output.txt
        echo ~{atac_percent_duplicates} '/n' >> output.txt
        echo ~{atac_nrf} '/n' >> output.txt
        echo ~{atac_pbc1} '/n' >> output.txt
        echo ~{atac_pbc2} '/n' >> output.txt
        echo ~{rna_total_reads} '/n' >> output.txt
        echo ~{rna_aligned_uniquely} '/n' >> output.txt
        echo ~{rna_aligned_multimap} '/n' >> output.txt
        echo ~{rna_unaligned} '/n' >> output.txt
        echo ~{rna_feature_reads} '/n' >> output.txt
        echo ~{rna_duplicate_reads} '/n' >> output.txt
        percent=$(( ~{default=0 rna_duplicate_reads}*100/~{default=1 rna_feature_reads} ))
        echo $percent '/n' >> output.txt
        PYTHONIOENCODING=utf-8 python3 /software/write_html.py ~{output_file} image_list.txt output.txt log_list.txt
        
        python3 /software/write_csv.py ~{other_output_file} names_list.txt numeric_list.txt image_list.txt log_list.txt
    >>>
    
    output {
        File html_report_file = "~{output_file}"
        #added output for the overall csv
        File csv_report_file = "~{other_output_file}"
    }

    runtime {
        #docker: 'nchernia/share_task_html_report:14'
        docker: 'mshriver01/share_task_html_report:latest'
    }
}



#echo "<h3>Summary Statistics</h3><p><table><tr><td colspan=2>ATAC</td></tr><tr><td>Total reads</td><td>" ~{atac_total_reads} "</td></tr>" > output.txt
 #       echo "<tr><td>Aligned uniquely</td><td>" ~{atac_aligned_uniquely} "</td></tr>" >> output.txt
  #      echo "<tr><td>Unaligned</td><td>" ~{atac_unaligned} "</td></tr>" >> output.txt
  #      echo "<tr><td>Unique Reads</td><td>" ~{atac_feature_reads} "</td></tr>" >> output.txt
   #     echo "<tr><td>Duplicate Reads</td><td>" ~{atac_duplicate_reads} "</td></tr>" >> output.txt
    #    echo "<tr><td>Percent Duplicates</td><td>" ~{atac_percent_duplicates} "</td></tr>" >> output.txt
     #   echo "<tr><td>NRF=Distinct/Total</td><td>" ~{atac_nrf} "</td></tr>" >> output.txt
     #   echo "<tr><td>PBC1=OnePair/Distinct</td><td>" ~{atac_pbc1} "</td></tr>" >> output.txt
      #  echo "<tr><td>PBC2=OnePair/TwoPair</td><td>" ~{atac_pbc2} "</td></tr>" >> output.txt
       # echo "<td colspan=2>RNA</td></tr><tr><td>Total reads</td><td>" ~{rna_total_reads} "</td></tr>" >> output.txt
        #echo "<tr><td>Aligned uniquely</td><td>" ~{rna_aligned_uniquely} "</td></tr>" >> output.txt
       # echo "<tr><td>Aligned multimap</td><td>" ~{rna_aligned_multimap} "</td></tr>" >> output.txt
       # echo "<tr><td>Unaligned</td><td>" ~{rna_unaligned} "</td></tr>" >> output.txt
       # echo "<tr><td>Filtered (feature) Reads</td><td>" ~{rna_feature_reads} "</td></tr>" >> output.txt
       # echo "<tr><td>Duplicate Reads</td><td>" ~{rna_duplicate_reads} "</td></tr>" >> output.txt
       # percent=$(( ~{default=0 rna_duplicate_reads}*100/~{default=1 rna_feature_reads} ))
       # echo "<tr><td>Percent Duplicates</td><td>" $percent "</td></tr></table>" >> output.txt