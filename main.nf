nextflow.enable.dsl=2

reads_ch = Channel.fromPath('data/ggal/*.fq')
        .map{ [ it.name.tokenize("_")[0], it] }
        .groupTuple()
        .view()

process foo {
  input:
    tuple val(sample_id), file(sample_file)
  output:
    tuple val(sample_id), file("${sample_id}.bam")
  script:
  """
    echo your_command_here --reads $sample_file > ${sample_id}.bam
  """
}


workflow{

    bam_ch = foo(reads_ch)
    bam_ch.view()

}


