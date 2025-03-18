version 1.0

# Import the sub-workflow for preprocessing the fastqs.
import "https://raw.githubusercontent.com/IGVF/single-cell-pipeline/52b4ddfef587aa2af771fe1a1c698d2e97e14a28/single_cell_pipeline.wdl" as preprocessing
import "structs/atac_output_struct.wdl"
import "structs/rna_output_struct.wdl"
import "structs/joint_output_struct.wdl"

# WDL workflow for SHARE-seq
struct Combinomics_output{
    Atac_outputs? atac_struct_output
    RNA_outputs? rna_struct_output
    Joint_outputs? joint_struct_output
}

workflow combinomics {

    input {
        # Commond inputs
        Boolean create_onlist_mapping = false
        String prefix
        String? subpool = "none"
        File genome_tsv

        # ATAC-specific inputs
        Array[File] atac_read1
        Array[File] atac_read2
        Array[File] fastq_barcode
        File atac_barcode_inclusion_list
        File? chromap_genome_index_tar_gz
        File? genome_fasta
        String atac_read_format

        # RNA-specific inputs
        Array[File] rna_read1
        Array[File] rna_read2
        Array[File] fastq_barcode_rna = []
        File rna_barcode_inclusion_list
        String kb_mode = "nac"
        String rna_read_format
        File? kb_genome_index_tar_gz
    }


    if ( process_atac && process_rna ) {
        if ( read1_atac[0] != "" && read1_rna[0] != "" ) {
            call joint_qc.joint_qc_plotting as joint_qc {
                input:
                    atac_barcode_metadata = atac.atac_qc_barcode_metrics,
                    rna_barcode_metadata = rna.rna_barcode_metadata,
                    prefix = prefix,
                    genome_name = genome_name_
            }
        }
    }

    output{
        # RNA outputs
        File? rna_final_bam = rna.task_starsolo_output_bam
        File? rna_bam_index = rna.task_starsolo_output_bam_index
        File? rna_starsolo_raw_tar = rna.task_starsolo_mtx_unique_tar
        File? rna_align_output_folder_tar = rna.task_starsolo_output_folder_tar

        File? rna_h5 = rna.rna_h5
        File? rna_barcode_metadata  = rna.rna_barcode_metadata
        File? rna_seurat_notebook_output = rna.rna_seurat_notebook_output
        File? rna_seurat_obj = rna.rna_seurat_obj

        # ATAC ouputs
        File? atac_final_bam = atac.atac_bam
        File? atac_bam_index = atac.atac_bam_index
        File? atac_bam_log = atac.atac_bam_alignment_stats
        File? atac_fragments = atac.atac_fragment_file
        File? atac_fragments_index = atac.atac_fragment_file_index
        File? atac_barcode_alignment_stats = atac.atac_align_barcode_statistics
        File? atac_barcode_metrics = atac.atac_qc_barcode_metrics
        File? atac_h5ad = atac.atac_qc_snapatac2_h5ad

        # Joint outputs
        File? joint_barcode_metadata = joint_qc.joint_barcode_metadata

        # Report
        File? html_summary = html_report.html_report_file
        File? csv_summary_file = html_report.csv_summary_file

        # Combined outputs
        Combinomics_output combinomics_struct_output = object{
            atac_struct_output: atac.atac_struct_output,
            rna_struct_output: rna.rna_struct_output,
            joint_struct_output: object{
                joint_qc_plot: atac.atac_qc_barcode_metrics,
                joint_density_plot: joint_qc.joint_density_plot,
                joint_barcode_metadata: joint_qc.joint_barcode_metadata
            }
        }
    }

}

