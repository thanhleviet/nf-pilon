#!/usr/bin/env nextflow
nextflow.enable.dsl=2

if (params.sample_sheet) {
    ch_pilon = Channel.fromPath(params.sample_sheet)
                        .splitCsv(header: true)
                        .map {row -> tuple(row.sample_id,[row.sr1,row.sr2],row.contigs)}
} else {
    ch_contigs = Channel.fromPath(params.contigs, checkIfExists: true)
                        .map { it -> tuple(it.baseName.replaceAll(/_assembly|_medaka/,''), file(it)) }

    ch_short_reads = Channel.fromFilePairs(params.short_reads, checkIfExists: true)
                            .map {it -> tuple(it[0].replaceAll(/_S[0-9]{1,3}/,''), it[1])}

    ch_pilon = ch_short_reads.join(ch_contigs)
}


process pilon {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    cpus 8

    input:
    tuple val(sample_id), path(reads), path(contigs)
    output:
    path "${sample_id}_pilon.fasta"
    path "*.changes", optional: true
    
    script:
    """
    pilon.py -t ${task.cpus} -n ${params.iterations} ${reads[0]} ${reads[1]} ${contigs}
    mv final.polished.fasta  ${sample_id}_pilon.fasta
    """
}


workflow {
    main:
    pilon(ch_pilon)
}