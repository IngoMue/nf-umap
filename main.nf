#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """\

         ================================================
         NF-μmap - A nextflow mapping workflow for museomics data
         https://github.com/IngoMue/nf-umap     
         Author: Ingo A. Müller
         ================================================
         |refseq       : ${params.refseq}
         |refprefix    : ${params.refprefix}
         |readpairs    : ${params.read_pairs}
         |mergedreads  : ${params.merged_reads}
         |outdir       : ${params.outdir}
         |
         |Include read pairs?....${params.inclRdPrs}
         |Include merged reads?..${params.inclMrgRds}
         |
         |Map with bwa-mem2?.....${params.usebwamem2}
         |
         |Skip merging?..........${params.skipmerge}
         |
         |X11 unix server?.......${params.X11}
         |
         |Print full QC report?..${params.fullreport}
         |Run DamageProfiler?....${params.runDMGprof}
         ================================================
         """
         .stripIndent()

refseq_ch = Channel.fromPath(params.refseq)
        .ifEmpty { error "No reference sequence found from: ${params.refseq}" }

read_prs = Channel.fromPath(params.read_pairs)
        .ifEmpty { error "Can't find read pairs matching pattern: ${params.read_pairs}" }
        .map { [ it.name.tokenize("_")[0..-2].join("_"), it] }
        .groupTuple()

mrg_reads = Channel.fromPath(params.merged_reads)
        .map { [ it.name.tokenize("_")[0..-2].join("_"), it] }
       	.ifEmpty { error "Can't find merged reads matching pattern: ${params.merged_reads}" }
        .groupTuple()


process faidx_ref {

    publishDir "${params.outdir}/01_RefIndex", mode:'copy'
    tag "Indexing $ref_file with samtools faidx"

    input:
      file ref_file

    output:
      file("*.fai")

    script:
    """
      samtools faidx $ref_file
    """
}

process mem2_index_ref {
    label 'High_CPU'

    publishDir "${params.outdir}/01_RefIndex", mode:'copy'
    tag "Indexing $ref_file with bwa-mem2 index"

    input:
      file ref_file

    output:
      file '*'

    script:
    """
      bwa-mem2 index -p $params.refprefix $ref_file
    """
}

process index_ref {
    label 'High_CPU'

    publishDir "${params.outdir}/01_RefIndex", mode:'copy'
    tag "Indexing $ref_file with bwa index"

    input:
      file ref_file

    output:
      file '*'

    script:
    """
      bwa index -p $params.refprefix $ref_file
    """
}


process bwa_mem2_PRS {
    label 'High_CPU'

    tag "Mapping read pairs for $sample_id onto $params.refprefix using bwa-mem2"

    input:
      tuple val(sample_id), file(sample_file), val(ref_ID), file(index), file(faidx)

    output:
      tuple val(sample_id), file("${sample_id}_PRS.sam")

    when:
      params.inclRdPrs == true

    script:
    """
      bwa-mem2 mem $params.refprefix $sample_file -t ${task.cpus} > ${sample_id}_PRS.sam
    """
}

process bwa_mem_PRS {
    label 'High_CPU'

    tag "Mapping read pairs for $sample_id onto $params.refprefix using bwa mem"

    input:
      tuple val(sample_id), file(sample_file), val(ref_ID), file(index), file(faidx)

    output:
      tuple val(sample_id), file("${sample_id}_PRS.sam")

    when:
      params.inclRdPrs == true

    script:
    """
      bwa mem $params.refprefix $sample_file -t ${task.cpus} > ${sample_id}_PRS.sam
    """
}

process quickcheck_PRS_sams {
    publishDir "${params.outdir}/00_quickcheck/", pattern: '*fofn', mode: 'copy'
    tag "Running quickcheck for all .sam files (read pairs)"

    input:
      file(sams)

    output:
      file("bad_PRS_sams.fofn")

    script:
    """
      samtools quickcheck ${sams} -v > bad_PRS_sams.fofn && echo 'All good' || echo 'Some files may be truncated, check bad_PRS_sams.fofn for a list of files'
    """

}

process bwa_mem2_MRG {
    label 'High_CPU'

    tag "Mapping merged reads for $mrg_id onto $params.refprefix using bwa-mem2"

    input:
      tuple val(mrg_id), file(merged_file), val(ref_ID), file(index), file(faidx)

    output:
      tuple val(mrg_id), file("${mrg_id}_MRG.sam")

    when:
      params.inclMrgRds == true

    script:
    """
      bwa-mem2 mem $params.refprefix $merged_file -t ${task.cpus} > ${mrg_id}_MRG.sam
    """
}

process bwa_mem_MRG {
    label 'High_CPU'

    tag "Mapping merged reads for $mrg_id onto $params.refprefix using bwa mem"

    input:
      tuple val(mrg_id), file(merged_file), val(ref_ID), file(index), file(faidx)

    output:
      tuple val(mrg_id), file("${mrg_id}_MRG.sam")

    when:
      params.inclMrgRds == true

    script:
    """
      bwa mem $params.refprefix $merged_file -t ${task.cpus} > ${mrg_id}_MRG.sam
    """
}

process quickcheck_MRG_sams {
    publishDir "${params.outdir}/00_quickcheck/", pattern: '*fofn', mode: 'copy'
    tag "Running quickcheck for all .sam files (merged reads)"

    input:
      file(sams)

    output:
      file("bad_MRG_sams.fofn")

    script:
    """
      samtools quickcheck ${sams} -v > bad_MRG_sams.fofn && echo 'All good' || echo 'Some files may be truncated, check bad_MRG_sams.fofn for a list of files'
    """

}


process samtools_bam_srt_idx_PRS {
    publishDir "${params.outdir}/02_bams/ReadPairs", pattern: '*_sorted{.bam,.bam.bai}', mode: 'copy'
    tag "Bam conversion, sorting and indexing for $sample_id (read pairs)"

    input:
      tuple val(sample_id), file(sam_file)

    output:
      file("${sample_id}_PRS_sorted.bam")
      file("${sample_id}_PRS_sorted.bam.bai")

    when:
      params.inclRdPrs == true

    script:
    """
      samtools view -@ ${task.cpus} -bS $sam_file > ${sample_id}_PRS.bam

      samtools sort -@ ${task.cpus} -o ${sample_id}_PRS_sorted.bam ${sample_id}_PRS.bam

      samtools index -@ ${task.cpus} ${sample_id}_PRS_sorted.bam
    """
}

process samtools_bam_srt_idx_MRG {
    publishDir "${params.outdir}/02_bams/MergedReads", pattern: '*_sorted{.bam,.bam.bai}', mode: 'copy'
    tag "Bam conversion, sorting and indexing for $sample_id (merged reads)"

    input:
      tuple val(sample_id), file(sam_file)

    output:
      file("${sample_id}_MRG_sorted.bam")
      file("${sample_id}_MRG_sorted.bam.bai")

    when:
      params.inclMrgRds == true

    script:
    """
      samtools view -@ ${task.cpus} -bS $sam_file > ${sample_id}_MRG.bam

      samtools sort -@ ${task.cpus} -o ${sample_id}_MRG_sorted.bam ${sample_id}_MRG.bam

      samtools index -@ ${task.cpus} ${sample_id}_MRG_sorted.bam
    """
}


process samtools_merge {
    publishDir "${params.outdir}/03_merged_bams/", pattern: '*{.bam,.bam.bai}', mode: 'copy'
    tag "Merging all .bam files for $sample_id"

    input:
      tuple val(sample_id), file(bam_file)

    output:
      tuple val(sample_id), file("${sample_id}.bam")
      file '*.bai'

    when:
      params.skipmerge == false

    script:
    """
      samtools merge -o ${sample_id}.bam $bam_file

      samtools index -@ ${task.cpus} ${sample_id}.bam
    """
}

process quickcheck_bams {
    publishDir "${params.outdir}/00_quickcheck/", pattern: '*fofn', mode: 'copy'
    tag "Running quickcheck for all .bam files"

    input:
      file(bams)

    output:
      file("bad_bams.fofn")

    script:
    """
      samtools quickcheck ${bams} -v > bad_bams.fofn && echo 'All good' || echo 'Some files may be truncated, check bad_bams.fofn for a list of files'
    """

}

process qualimap_bamqc {
    publishDir "${params.outdir}/04_mapping_stats/${sample_id}", pattern: '*.{txt,pdf}', mode: 'move'
    tag "Generating mapping QC report for $sample_id"
    label 'High_mem'

    input:
      tuple val(sample_id), file(bam_file)

    output:
      file '*'

    script:
        if (params.X11 == true && params.fullreport == false)
            """
              qualimap bamqc -bam ${bam_file} -outdir ${sample_id}/ -c --java-mem-size=${task.memory.toGiga()}G -nt ${task.cpus}
              cp ./${sample_id}/genome_results.txt ./${sample_id}_genome_results.txt
	    """
        else if (params.X11 == false && params.fullreport == false)
            """
              unset DISPLAY
              qualimap bamqc -bam ${bam_file} -outdir ${sample_id}/ -c --java-mem-size=${task.memory.toGiga()}G -nt ${task.cpus}
              cp ./${sample_id}/genome_results.txt ./${sample_id}_genome_results.txt
            """
        else if (params.X11 == true && params.fullreport == true)
            """
              qualimap bamqc -bam ${bam_file} -outdir ${sample_id}/ -outfile ${sample_id}_report.pdf -c --java-mem-size=${task.memory.toGiga()}G -nt ${task.cpus}
              cp ./${sample_id}/* .
            """
        else if (params.X11 == false && params.fullreport == true)
            """
              unset DISPLAY
              qualimap bamqc -bam ${bam_file} -outdir ${sample_id}/ -outfile ${sample_id}_report.pdf -c --java-mem-size=${task.memory.toGiga()}G -nt ${task.cpus}
              cp ./${sample_id}/* .
            """

}

process dmgprof {
    publishDir "${params.outdir}/05_damage/${sample_id}", pattern: '*.pdf', mode: 'move'
    tag "Generating damage report for $sample_id"
    label 'High_mem'

    input:
      tuple val(sample_id), file(bam_file), val(ref_ID), file(index), file(faidx)

    output:
      file '*'

    script:
        """
          damageprofiler -Xmx${task.memory.toGiga()}g -i $bam_file -r ${params.refseq} -o .
        """

}


workflow {
    faidx_ref(refseq_ch)
    if (params.usebwamem2) {
        mem2_index_ref(refseq_ch)
        bwa_index = mem2_index_ref.out
    } else {
        index_ref(refseq_ch)
        bwa_index = index_ref.out
    }
    refindex = refseq_ch.combine(bwa_index).map { it -> tuple(it.simpleName, it) }.combine(faidx_ref.out)

    if (params.inclRdPrs) {
        read_PRS_idx = read_prs.combine(refindex)
        if (params.usebwamem2) {
            bwa_mem2_PRS(read_PRS_idx)
            mapped_PRS = bwa_mem2_PRS.out
        } else {
            bwa_mem_PRS(read_PRS_idx)
            mapped_PRS = bwa_mem_PRS.out
        }
        quickcheck_PRS_sams(mapped_PRS.flatten().filter( ~/.*sam$/ ).toList())
        samtools_bam_srt_idx_PRS(mapped_PRS)
    }

    if (params.inclMrgRds) {
        MRG_reads_idx = mrg_reads.combine(refindex)
        if (params.usebwamem2) {
            bwa_mem2_MRG(MRG_reads_idx)
            mapped_MRG = bwa_mem2_MRG.out
        } else {
            bwa_mem_MRG(MRG_reads_idx)
            mapped_MRG = bwa_mem_MRG.out
        }
        quickcheck_MRG_sams(mapped_MRG.flatten().filter( ~/.*sam$/ ).toList())
        samtools_bam_srt_idx_MRG(mapped_MRG)
    }

    if (params.inclRdPrs == true && params.inclMrgRds == true) {
        samtools_merge(samtools_bam_srt_idx_PRS.out[0].mix(samtools_bam_srt_idx_MRG.out[0]).map { [ it.name.split(/_L0\d+/)[0], it] }.groupTuple())
    } else if (params.inclRdPrs == true && params.inclMrgRds == false) {
        samtools_merge(samtools_bam_srt_idx_PRS.out[0].map { [ it.name.split(/_L0\d+/)[0], it] }.groupTuple())
    } else if (params.inclRdPrs == false && params.inclMrgRds == true) {
        samtools_merge(samtools_bam_srt_idx_MRG.out[0].map { [ it.name.split(/_L0\d+/)[0], it] }.groupTuple())
    }

    quickcheck_bams(samtools_merge.out[0].flatten().filter( ~/.*bam$/ ).toList())

    qualimap_bamqc(samtools_merge.out[0])
    if (params.runDMGprof) {
        dmgprof(samtools_merge.out[0].combine(refindex))
    }

}
