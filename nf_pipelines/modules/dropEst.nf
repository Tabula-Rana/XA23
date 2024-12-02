#! /usr/bin/env nextflow
nextflow.enable.dsl=2

/* droptag:
    inDrop v1 & v2
    File 1: barcode reads. (R2)
        Structure:
        Cell barcode, part 1
        Spacer
        Cell barcode, part 2
        UMI
    File 2: gene reads (R1)

    inDrop v3
    File 1: cell barcode, part 1 (default length: 8bp) (R2)
    File 2: cell barcode + UMI, part 1 (default length: >= 14bp) (R4)
    File 3: gene read (R1)
    File 4 (optional): library tag (R3)
*/


process dropTag_v2() {
    tag "dropTag_v2 on $subdir"
    label "dropEst"
    
    queue "short"
    memory params.droptag_mem  //"24 GB"
    time params.droptag_time // "60m"
    cpus params.droptag_threads
    
    publishDir "${params.output_root_dir}/${subdir}", mode: 'copy' 

    input:
        tuple val(subdir), path(read1), path(read2)
        path inDrops_v2_xml
    
    output:
        tuple val(subdir), path("*.fastq.gz"), emit: tagged_reads
        tuple val(subdir), path("*.params.gz"), emit: tagged_params
        tuple val(subdir), path("*.tagged.rds"), emit: tagged_rds
        path("*")

    script:
    """
    droptag -c $inDrops_v2_xml -S -s -p $params.droptag_threads $read2 $read1
    """
}

process dropTag_microwell() {
    tag "dropTag_microwell on $sample_name"
    label "dropEst"

    queue "short"
    memory params.droptag_mem  //"24 GB"
    time params.droptag_time // "60m"
    cpus params.droptag_threads
    
    publishDir "${params.output_root_dir}/${sample_name}", mode: 'copy' 


    input:
        tuple val(sample_name), path(reads)
        path microwell_xml
    
    output:
        tuple val(sample_name), path("*.fastq.gz"), emit: tagged_reads
        tuple val(sample_name), path("*.params.gz"), emit: tagged_params
        tuple val(sample_name), path("*.tagged.rds"), emit: tagged_rds
        path("*")

    script:
    """
    droptag -c $microwell_xml -S -s -p ${params.droptag_threads} ${reads[0]} ${reads[1]}
    """
}

process dropTag_v3() {
    tag "dropTag_v3 on ${sample_name}"
    label "dropEst"

    queue "short"
    cpus params.droptag_threads
    memory params.droptag_mem
    time params.droptag_time
    publishDir "${params.output_root_dir}/${sample_name}", mode: 'copy' 

    input:
        tuple val(sample_name), path(r1), path(r2), path(r4)
        path inDrops_v3_xml
    
    output:
        tuple val(sample_name), path("*.fastq.gz"), emit: tagged_reads
        tuple val(sample_name), path("*.params.gz"), emit: tagged_params
        tuple val(sample_name), path("*.rds"), emit: tagged_rds
        path("*")

    script:

    """
    droptag -c $inDrops_v3_xml -S -s -p ${params.droptag_threads} -n $sample_name $r2 $r4 $r1
    """
}

// process dropTag_microwell() {
//     tag "dropTag_microwell on $sample_name"
//     queue "short"
//     publishDir "${params.output_root_dir}/${subdir}/${sample_name}", mode: 'copy' 
//     cpus params.droptag_threads
//     label "dropEst"

//     input:
    
//     tuple val(sample_name), path(r1), path(r2)
//     val Dropest_xml
    
//     output:
//     tuple val(sample_name), val(sample_name), path("*.fastq.gz"), emit: tagged_reads
//     tuple val(sample_name), val(sample_name), path("*.params.gz"), emit: tagged_params
//     tuple val(sample_name), path("*.tagged.rds"), emit: tagged_rds

//     script:
//     """
//     droptag -c $Dropest_xml -S -s -p ${params.droptag_threads} $r2 $r1
//     """
// }

process dropEst {
    tag "dropEst on $subdir"
    publishDir "${params.output_root_dir}/${subdir}", mode: 'move'
    queue "short"
    cpus 1
    time params.dropest_time
    memory params.dropest_mem
    label "dropEst"

    input:
        // tuple val(subdir), path(bam_fs, stageAs: "bam?/*"), path(params_fs, stageAs: "param?/*")
        tuple val(subdir), path(bam_f), path(params_f)
        path annotation_gtf
        path Dropest_xml
        path barcodes
        val gene_parts
    
    output:
        path("*")

    script:

    """
    dropest \
        -r $params_f \
        -w \
        -V \
        -m \
        -L $gene_parts \
        -g $annotation_gtf \
        -c $Dropest_xml \
        $bam_f
    """
}
