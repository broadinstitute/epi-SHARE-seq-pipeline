version 1.0

task run_fastqc{
    input{
      Array[File] fastq
      String? sample_prefix = "report"
      Int memory = 8
      Int disk_space = 300
      Int? num_threads = 4
    }

    command<<<
        mkdir output

        for fq in ~{sep=" " fastq}
        do
            name=${fq##*/}
            prefix=${name%%.*}
            gzip -dc ${fq} | fastqc stdin:$prefix -c contaminants.tsv --noextract --threads ~{num_threads} --outdir '.'
        done

        multiqc '.' --force -o ~{sample_prefix}
        
        cp ~{sample_prefix}/~{sample_prefix}_multiqc_report.html .

        tar cvzf ~{sample_prefix}_multiqc_output.tar.gz output
    >>>

  output {
    Array[File] fastqc_html=glob("*_fastqc.html")
    Array[File] fastqc_zip=glob("*_fastqc.zip")
    File multiqc_report="~{sample_prefix}_multiqc_report.html"
    File multiqc_tar="~{sample_prefix}_multiqc_output.tar.gz"
  }

  runtime {
    docker: "ctzouana/cutadapt7:latest"
    memory: "${memory}GB"
    disks: "local-disk ${disk_space} HDD"
    cpu: "${num_threads}"
  }

  meta {
    author: "Eugenio Mattei"
  }

}