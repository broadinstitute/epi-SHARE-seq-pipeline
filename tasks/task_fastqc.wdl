version 1.0

task run_fastqc{
    input{
      Array[File] fastq
      String? sample_prefix = "multiqc_output"
      Int memory = 8
      Int disk_space = 300
      Int? num_threads = 4
    }

    command<<<

        for fq in ~{sep=" " fastq}
        do
            name=${fq##*/}
            prefix=${name%%.*}
            gzip -dc ${fq} | fastqc stdin:$prefix --contaminants /common/share_contaminants.tsv --adapters /common/share_contaminants.tsv --noextract --threads ~{num_threads} --outdir '.'
        done

        multiqc '.' --force -o ~{sample_prefix}
        
        cp ~{sample_prefix}/multiqc_report.html ~{sample_prefix}_multiqc_report.html

        tar cvzf ~{sample_prefix}_multiqc_output.tar.gz ~{sample_prefix}
    >>>

  output {
    File multiqc_report="~{sample_prefix}_multiqc_report.html"
    File multiqc_tar="~{sample_prefix}_multiqc_output.tar.gz"
  }

  runtime {
    docker: "polumechanos/multiqc"
    memory: "${memory}GB"
    disks: "local-disk ${disk_space} HDD"
    cpu: "${num_threads}"
  }

  meta {
    author: "Eugenio Mattei"
  }

}