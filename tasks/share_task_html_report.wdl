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
        Int atac_total_reads
        Int atac_aligned_uniquely 
        Int atac_unaligned 
        Int atac_feature_reads
        Int atac_duplicate_reads 
        Int rna_total_reads
        Int rna_aligned_uniquely 
	Int rna_aligned_multimap 
        Int rna_unaligned 
        Int rna_feature_reads 
        Int rna_duplicate_reads 

        ## JPEG files to be encoded and appended to html
        # RNA plots
        File share_rna_qc_library_plot
        Array[File] share_rna_umi_qc_plots
        File? share_rna_seurat_violin_plot 
        File? share_rna_seurat_mitochondria_qc_plot 
        File? share_rna_seurat_features_plot 
        File? share_rna_seurat_PCA_dim_loadings_plot 
        File? share_rna_seurat_PCA_plot 
        File? share_rna_seurat_heatmap_plot 
        File? share_rna_seurat_jackstraw_plot 
        File? share_rna_seurat_elbow_plot 
        File? share_rna_seurat_umap_plot 
        # ATAC plots
        File share_atac_qc_library_plot
        File share_atac_qc_hist_plot
        File share_atac_qc_tss_enrichment
        File? share_atac_archr_gene_heatmap_plot 
        File? share_atac_archr_tss_enrichment_raw 
        File? share_atac_archr_tss_enrichment_filtered 
        File? share_atac_archr_fragment_size_plot 
        File? share_atac_archr_doublet_plot 
        File? share_atac_archr_umap_plot 
        # DORC plots
        File? seurat_violin_plot 
        File? j_plot 

        ## Raw text logs to append to end of html
        # RNA logs
#        File share_rna_alignment_log
#        File share_rna_featurecount_exon_txt
#        File? share_rna_featurecount_intron_txt
#        File share_rna_qc_reads_distribution
#        File share_rna_qc_reads_distribution2
#        File share_rna_umi_rm_dup_log 
#        File share_rna_seurat_notebook_log
        # ATAC logs
#        File share_atac_alignment_log
#        File share_atac_archr_notebook_output 
#        File share_atac_archr_notebook_log 
#        File share_atac_archr_papermill_log 
        # DORCs logs
#        File dorcs_notebook_log 
    }

    command <<<
	fnames=$(ls -m *jpg)
        text=$(echo "<h3>Summary Statistics</h3><p><table><tr><td colspan=2>ATAC</td></tr><tr><td>Total reads</td><td>" ~{atac_total_reads} "</td></tr>")
	text=$(echo $text "<tr><td>Aligned uniquely</td><td>" ~{atac_aligned_uniquely} "</td></tr>")
        text=$(echo $text "<tr><td>Unaligned</td><td>" ~{atac_unaligned} "</td></tr>")
        text=$(echo $text "<tr><td>Filtered Reads</td><td>" ~{atac_feature_reads} "</td></tr>")
        text=$(echo $text "<tr><td>Duplicate Reads</td><td>" ~{atac_duplicate_reads} "</td></tr>")
        percent=$(( ~{atac_duplicate_reads}*100/~{atac_feature_reads} ))
        text=$(echo $text "<tr><td>Percent Duplicates</td><td>" ~{percent} "</td></tr>")
        text=$(echo $text "<td colspan=2>RNA</td></tr><tr><td>Total reads</td><td>" ~{rna_total_reads} "</td></tr>")
        text=$(echo $text "<tr><td>Aligned uniquely</td><td>" ~{rna_aligned_uniquely} "</td></tr>")
        text=$(echo $text "<tr><td>Aligned multimap</td><td>" ~{rna_aligned_multimap} "</td></tr>")
        text=$(echo $text "<tr><td>Unaligned</td><td>" ~{rna_unaligned} "</td></tr>")
        text=$(echo $text "<tr><td>Filtered (feature) Reads</td><td>" ~{rna_feature_reads} "</td></tr>")
        text=$(echo $text "<tr><td>Duplicate Reads</td><td>" ~{rna_duplicate_reads} "</td></tr>")
        percent=$(( ~{rna_duplicate_reads}*100/~{rna_feature_reads} ))
        text=$(echo $text "<tr><td>Percent Duplicates</td><td>" ~{percent} "</td></tr>")       
        python3 $(which write_html.py) $fnames $text	
    >>>
    output {
	File html_report = output.html
    }

    runtime {
        docker: 'ubuntu:latest'
    }
    parameter_meta {
        alignment_log: {
            description: 'ATAC alignment log file',
	    help: 'Log file from ATAC alignment step.',
            example: 'SS-PKR-30-96-ENTIRE-PLATE.atac.align.hg38.Log.out'
        }
        dups_log: {
            description: 'ATAC dups log file',
            help: 'Log file from ATAC rmdups step.',
            example: 'SS-PKR-12.atac.counts.mm10.filtered.cs.log'
        }
    }
}
