version 1.0

# Import the sub-workflow for preprocessing the fastqs.
import "workflows/subwf-atac-single-organism.wdl" as share_atac
import "workflows/subwf-rna-single-organism.wdl" as share_rna
import "task/dorcs_task_find_dorcs.wdl" as find_dorcs

# WDL workflow for SHARE-seq

workflow ShareSeq {

    input {
        # Common inputs
        String prefix = "shareseq-project"
        String genome_name
        Int? cpus = 16


        # ATAC specific inputs
        File chrom_sizes
        File read1_atac
        File read2_atac
        File idx_tar_atac
        File tss_bed
        Int? cpus_atac
        Int? cutoff_atac

        # RNA specific inputs
        Boolean multimappers = false
        Boolean include_multimappers = false
        Boolean include_introns = false
        File read1_rna
        File genes_annotation_bed
        File gtf
        File idx_tar_rna
        Int? cpus_rna
        String? gene_naming = "gene_name"

        # Group UMI
        Boolean remove_single_umi = true
        String mode = "fast"
        Int cutoff_rna = 100

        # Lib_size QC
        Boolean qc = false


        # DORCs specific inputs
        Int? cpus_dorcs = 4
        String save_plots_to_dir = "TRUE"
        String output_filename = "output.ipynb"

        # Seurat filters
        Int minFeature_RNA = 200
        Int maxFeature_RNA = 2500
        Float percentMT_RNA = 5
        Int minCells_RNA = 3

        # DORCs filter
        Int dorcGeneCutOff = 10
        Float fripCutOff = 0.3
        Float corrPVal = 0.05
        Int topNGene = 20

        # Regulatory region around TSS. Default is +/- 50Kb
        Int windowPadSize = 50000
        Int bootstraps = 100

        String docker_image_dorcs = "polumechanos/dorcs_task_find_dorcs"
        Int mem_gb_dorcs = 16
    }

    call share_rna.wf_rna as rna{
        File read1 = read1_rna
        File idx_tar = idx_tar_rna
        String prefix = prefix
        String genome_name = genome_name
        Int? cpus = cpus_rna
        # Update RGID
        Boolean multimappers = multimappers
        # Assign features
        Boolean include_multimappers = include_multimappers
        Boolean include_introns = include_introns
        File gtf = gtf
        String gene_naming = gene_naming
        # Group UMI
        Boolean remove_single_umi = remove_single_umi
        String mode = mode
        Int cutoff = cutoff_rna
        # Lib_size QC
        Boolean qc = qc
        File genes_annotation_bed = genes_annotation_bed
    }
    call share_atac.wf_atac as atac{
        File read1 = read1_atac
        File read2 = read2_atac
        File chrom_sizes = chrom_sizes
        File idx_tar = idx_tar_atac
        File tss_bed = tss_bed
        String prefix = prefix
        String genome_name = genome_name
        Int cutoff = cutfoff_atac
        Int? cpus = cpus_atac
    }
    call find_dorcs{
        File rna_matrix = rna.share_rna_h5_matrix
        File atac_fragments = atac.share_atac_alignment_filtered
        File peak_file

        String genome = genome_name
        Int n_cores = cpus_dorcs
        String save_plots_to_dir = save_plots_to_dir
        String output_filename = output_filename

        Int minFeature_RNA = minFeature_RNA
        Int maxFeature_RNA = maxFeature_RNA
        Float percentMT_RNA = percentMT_RNA
        Int minCells_RNA = minCells_RNA

        Int dorcGeneCutOff = dorcGeneCutOff
        Float fripCutOff = fripCutOff
        Float corrPVal = corrPVal
        Int topNGene = topNGene

        # Regulatory region around TSS. Default is +/- 50Kb
        Int windowPadSize = windowPadSize
        Int bootstraps = bootstraps

        String docker_image = docker_image_dorcs
        Int mem_gb = mem_gb_dorcs
    }

}


task merge_fastqs{
    input{
        Array[File] atac_read1
        Array[File] atac_read2
        Array[File] rna_read1
        String? prefix
    }
    command{
        cat ${sep=' ' atac_read1} > ${prefix + "."}merged.atac.R1.fq.gz
        cat ${sep=' ' atac_read2} > ${prefix + "."}merged.atac.R2.fq.gz
        cat ${sep=' ' rna_read1} > ${prefix + "."}merged.rna.R1.fq.gz
    }
    output{
        File merged_atac_fastq_R1 = glob('*.merged.atac.R1.fq.gz')[0]
        File merged_atac_fastq_R2 = glob('*.merged.atac.R2.fq.gz')[0]
        File merged_rna_fastq_R1 = glob('*.merged.rna.R1.fq.gz')[0]
    }
}


# Task to report errors to user.
# From https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/chip.wdl
task raise_exception {
  input {
    String msg
    Array[String]? vals
  }
  command {
    echo -e "\n* Error: ${msg}\n" >&2
    echo -e "* Vals: ${sep=',' vals}\n" >&2
    exit 2
  }
  output {
    String error_msg = '${msg}'
  }
  runtime {
    maxRetries : 0
    cpu : 1
    memory : '2 GB'
    time : 1
    disks : 'local-disk 10 SSD'
  }
}
