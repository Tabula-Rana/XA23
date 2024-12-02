process fasterq_dump {
    tag "fasterq_dump on $sra"

    queue "short"
    memory params.sra_mem
    time params.sra_time
    cpus params.sra_threads

    publishDir "${params.output_root_dir}/${subdir}/${sra}", mode: 'copy'

    input:
        tuple val(subdir), val(sra)
        val fasterq_dump_cl

    output:
        tuple val(sra), path("*.fastq")

    script:
    """
    ${params.fasterq_dump_bin} -e $params.sra_threads $fasterq_dump_cl $sra
    """
}

// process gff2gff {
//     publishDir "${annotation_gtf}", mode: 'copy'

//     input:
//     val annotation_gtf
//     val gffread_bin
//     path annotation_gff3

//     output:
//     path annotation_gtf

//     script:
//     """
//     // gff3 -> gtf
//     $gffread_bin -T ${annotation_gff3} -o ${annotation_gtf}
//     """
// }



process fastqc {
    tag "fastqc on $fastq"
    label "conda"

    queue "short"
    memory params.fastqc_mem
    time params.fastqc_time
    cpus params.fastqc_threads

    publishDir "${params.output_root_dir}", mode: 'copy'

    input:
        path(fastq)

    output:
        path("*")

    script:
    """
    fastqc -t $params.fastqc_threads $fastq
    """
}


process multiqc {
    tag "multiqc"
    label "conda"

    queue "short"
    memory params.multiqc_mem
    time params.multiqc_time
    cpus params.multiqc_threads

    publishDir "${params.output_root_dir}", mode: 'copy'

    input:
        path(fastqcs)

    output:
        path("*")

    script:
    """
    multiqc .
    """

}


process fastp {
    tag "fastp on $fastq"
    label "conda"

    memory params.fastp_mem
    time params.fastp_time
    cpus params.fastp_threads
    queue "short"
    
    publishDir "${params.output_root_dir}/${subdir}", mode: 'move'


    input:
        tuple val(subdir), path(fastq)
        val fastp_cl
    
    output:
        path("*")

    script:
    
    fastp_fastq = fastq.toString() - ~/\.\w+$/ + ".fastp.fastq"  

    """
    fastp \
        -w $params.fastp_threads \
        -i $fastq \
        -o $fastp_fastq \
        -V \
        $fastp_cl
    """

}

process fastq_paired {
    tag "fastq_paired on $subdir"

    memory params.fastq_paired_mem
    time params.fastq_paired_time
    cpus params.fastq_paired_threads
    queue "short"
    
    publishDir "${params.output_root_dir}/${subdir}", mode: 'move'


    input:
        tuple val(subdir), path(fastq1), path(fastq2)
        val fastq_paired_bin
        val fastq_pair_cl
    
    output:
        path("*")

    script:

    """
    $fastq_paired_bin $fastq_pair_cl $fastq1 $fastq2
    """
}


process FilterOutShortReads {
    memory "32 GB"
    queue "short"
    time "60m"
    cpus params.num_threads
    
    publishDir "${params.output_root_dir}/${sra}", mode: 'copy'

    tag "FilterOutShortReads on $sra"

    input:
    val seqtk_bin
    tuple val(sra), path(reads)
    val length
    
    output:
    tuple val(sra), path("*_{1,2}.filtered.fastq")

 
    script:
    
    filtered_ids = reads[1].simpleName + "_gte" + length + ".list"

    """
    $seqtk_bin comp ${reads[1]} | awk '{ if (\$2 >= ${length}) { print} }' | cut --fields 1 > $filtered_ids
    $seqtk_bin subseq ${reads[0]} $filtered_ids > ${reads[0].simpleName}".filtered.fastq"
    $seqtk_bin subseq ${reads[1]} $filtered_ids > ${reads[1].simpleName}".filtered.fastq"
    """

}


process SubsampleReads {
    memory "32 GB"
    queue "short"
    time "60m"
    cpus params.num_threads
    
    publishDir "${params.output_root_dir}/${sra}", mode: 'copy'

    tag "SubsampleReads on $sra"

    input:
    val seqtk_bin
    tuple val(sra), path(reads)
    val sample_size
    
    output:
    tuple val(sra), path("*_{1,2}.subsampled.fastq")

 
    script:

    """
    $seqtk_bin sample -s100 ${reads[0]} $sample_size > ${reads[0].simpleName}".subsampled.fastq"
    $seqtk_bin sample -s100 ${reads[1]} $sample_size > ${reads[1].simpleName}".subsampled.fastq"
    """

}


process RemoveGenome {
    input:
    val STAR_bin
    val genome_dir
    val results

    script:
    """
    ${STAR_bin} \
        --genomeLoad Remove \
        --genomeDir $genome_dir
    """
}

process demultiplex_barcodes_inDrops_v3 {
    publishDir "${params.output_root_dir}/${params.lane}/demultiplexed", mode: 'copy'
    
    queue "short"
    cpus 1
    memory params.seqtk_mem
    time params.seqtk_time
    
    tag "demultiplex_barcodes_inDrops_v3 on ${demultiplexed_lib_barcodes_fq}"

    input:
    val seqtk_bin
    tuple val(sample_name), path(barcode_reads_fq_1), path(barcode_reads_fq_2), path (library_read_fq), path(demultiplexed_lib_barcodes_fq)

    output:
    path "${barcode_reads_fq_1.simpleName}.${demultiplexed_lib_barcodes_fq.simpleName}.fastq", emit: barcode_reads_fq_1
    path "${barcode_reads_fq_2.simpleName}.${demultiplexed_lib_barcodes_fq.simpleName}.fastq", emit: barcode_reads_fq_2
    path "${library_read_fq.simpleName}.${demultiplexed_lib_barcodes_fq.simpleName}.fastq", emit: library_read_fq

    script:
    
    """
    $seqtk_bin comp $demultiplexed_lib_barcodes_fq | cut --fields 1 > demultiplexed_ids
    $seqtk_bin subseq ${barcode_reads_fq_1} demultiplexed_ids > "${barcode_reads_fq_1.simpleName}.${demultiplexed_lib_barcodes_fq.simpleName}.fastq"
    $seqtk_bin subseq ${barcode_reads_fq_2} demultiplexed_ids > "${barcode_reads_fq_2.simpleName}.${demultiplexed_lib_barcodes_fq.simpleName}.fastq"
    $seqtk_bin subseq ${library_read_fq} demultiplexed_ids > "${library_read_fq.simpleName}.${demultiplexed_lib_barcodes_fq.simpleName}.fastq"
    """

}

process demultiplex_libraries_inDrops_v3 {
    publishDir "${params.output_root_dir}/${params.lane}/demultiplexed", mode: 'copy'
    
    queue "short"
    cpus params.cutadapt_threads
    memory params.cutadapt_mem
    time params.cutadapt_time
    
    label "conda"
    tag "demultiplex_libraries_inDrops_v3 on ${lib_barcode_read}"

    
    input:
    tuple val(sample_name), path(lib_barcode_read)
    val barcodes_fa

    output:
    tuple val(sample_name), path("*.fastq.gz")

    script:
    """
    cutadapt -e 1 -j ${params.cutadapt_threads} -g ^file:${barcodes_fa} -o "{name}.fastq.gz" --action=none ${lib_barcode_read}
    """
}
