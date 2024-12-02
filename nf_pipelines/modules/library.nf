process concat_fastq {
    tag "concat_fastq on ${sample_name}.${suffix}"

    queue "short"
    memory params.concat_mem
    time params.concat_time
    cpus 1

    input:
        tuple val(sample_name), val(fastqs)
        val suffix
    
    output:
        tuple val(sample_name), path("${sample_name}.${suffix}.fastq")
    
    script:
    
    def fastqs_str = fastqs instanceof List ? fastqs.join(" ") : fastqs
    
    """
    cat $fastqs_str > ${sample_name}.${suffix}.fastq
    """
}

process FilterOutShortReads_single {
    memory "16 GB"
    queue "short"
    time "1h"
    cpus 1
    publishDir "${input_dir}", mode: 'copy'

    tag "FilterOutShortReads_single on $fastq"

    input:
    val seqtk_bin
    tuple val(input_dir), path(fastq)
    val len
    
    output:
    path("*.filtered.fastq.gz")


    script:
    
    filtered_ids = fastq.simpleName + "_gte${len}.list"
    filtered_fastq = "${fastq.simpleName}.filtered.fastq"


    """
    $seqtk_bin comp ${fastq} | awk '{ if (\$2 >= ${len}) { print} }' | cut --fields 1 > $filtered_ids
    $seqtk_bin subseq ${fastq} $filtered_ids > ${filtered_fastq}
    gzip -f ${filtered_fastq}
    """

}


// process cutadapt {
//     tag "cutadapt on $fastq"
//     label "conda"

//     queue "short"
//     memory params.cutadapt_mem
//     time params.cutadapt_time
//     cpus params.cutadapt_threads

//     publishDir "${params.output_root_dir}/${output_folder}", mode: 'copy'


//     input:
//         tuple val(sample_name), path(fastq)
//         val cutadapt_cl
//         val output_folder
    
//     output:
//         path("*.fastq.gz")
//         path("*.cutadapt.stats.txt")


//     script:
    
//     trimmed_fastq = "${sample_name}.trimmed.fastq.gz"  
//     cutadapt_stats = "${sample_name}.cutadapt.stats.txt"


//     """
//     cutadapt \
//         -j 0 \
//         -o $trimmed_fastq \
//         $cutadapt_cl \
//         --rename='{id} adapter={adapter_name} {match_sequence}' \
//         $fastq \
//         > $cutadapt_stats
//     """

// }


process cutadapt {
    tag "cutadapt on $fastq"
    label "conda"

    queue "short"
    memory params.cutadapt_mem
    time params.cutadapt_time
    cpus params.cutadapt_threads

    publishDir "${params.output_root_dir}/${subdir}", mode: 'move'


    input:
        tuple val(subdir), val(output_name), path(fastq)
        val cutadapt_cl
    
    output:
        path(output_name)
        path("*.cutadapt.stats.txt")


    script:
    
    cutadapt_stats = "${output_name}.cutadapt.stats.txt"

    """
    cutadapt \
        -j 0 \
        -o $output_name \
        $cutadapt_cl \
        $fastq \
        > $cutadapt_stats
    """

}