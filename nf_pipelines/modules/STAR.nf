process IndexTranscriptome {
    publishDir "${genome_dir}", mode: 'copy'
    cpus params.star_threads
    memory params.star_memory  // "34 GB"
    time params.star_time  // "1h 30m"
    queue "short"
    label "STAR"

    input:
    val genome_fasta
    val genome_dir
    val annotation_gtf
    val STAR_overhang
    val STAR_bin
    val sjdbGTFtagExonParentGene  // geneID
    val sjdbGTFtagExonParentTranscript  // Parent

    output:
    path "*"
    
    script:
    """
    
    ${STAR_bin} \
        --runThreadN $params.star_threads \
        --runMode genomeGenerate \
        --genomeDir "." \
        --genomeFastaFiles $genome_fasta \
        --sjdbGTFtagExonParentTranscript $sjdbGTFtagExonParentTranscript \
        --sjdbGTFtagExonParentGene $sjdbGTFtagExonParentGene \
        --sjdbGTFfile $annotation_gtf \
        --sjdbOverhang $STAR_overhang
    """
}


process StarAlignAndCount_inDrops_v2 {
    publishDir "${params.output_root_dir}/${sra}", mode: 'copy'
    cpus params.n_threads
    maxForks 1
    label "STAR"

    input:
    val STAR_bin
    val genome_dir
    val indrops_cell_barcode1
    val indrops_cell_barcode2
    tuple val(sra), path(reads)
    
    output:
    path "STAR/*"
    
    script:
    """
    mkdir STAR

    ${STAR_bin} \
        --runThreadN $task.cpus \
        --limitBAMsortRAM 30000000000 \
        --genomeLoad LoadAndKeep \
        --genomeDir $genome_dir \
        --readFilesIn ${reads[0]} ${reads[1]} \
        --soloFeatures Gene GeneFull \
        --soloType CB_UMI_Complex \
        --soloCBmatchWLtype 1MM \
        --soloAdapterSequence GAGTGATTGCTTGTGACGCCTT \
        --soloAdapterMismatchesNmax 2 \
        --soloCBwhitelist $indrops_cell_barcode1 $indrops_cell_barcode2 \
        --soloCBposition 0_0_2_-1 3_1_3_8 \
        --soloUMIposition 3_9_3_14 \
        --soloUMIdedup 1MM_CR \
        --soloUMIfiltering MultiGeneUMI_CR \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes NH HI AS nM CR CY UR UY CB UB sM \
        --soloMultiMappers EM \
        --outFileNamePrefix STAR/
    """
}


process StarRemoveGenome {
    input:
    val STAR_bin
    val genome_dir
    val results
    label "STAR"

    script:
    """
    ${STAR_bin} \
        --genomeLoad Remove \
        --genomeDir $genome_dir
    """
}

// process StarAlignAndCount_inDrops_v3 {
//     publishDir "${params.output_root_dir}/${sra}", mode: 'copy'
//     cpus params.n_threads
//     maxForks 1
    // label "STAR"

//     input:
//     val STAR_bin
//     path genome_dir
//     val indrops_cell_barcode1
//     val indrops_cell_barcode2
//     tuple val(sra), path(reads)
    
//     output:
//     path "STAR/*"
    
//     script:
//     """
//     mkdir STAR

//     ${STAR_bin} \
//         --runThreadN $task.cpus \
//         --limitBAMsortRAM 30000000000 \
//         --genomeLoad LoadAndKeep \
//         --genomeDir $genome_dir \
//         --readFilesIn ${reads[0]} ${reads[1]} \
//         --soloFeatures Gene GeneFull \
//         --soloType CB_UMI_Complex \
//         --soloCBmatchWLtype 1MM \
//         --soloAdapterSequence GAGTGATTGCTTGTGACGCCTT \
//         --soloAdapterMismatchesNmax 2 \
//         --soloCBwhitelist $indrops_cell_barcode1 $indrops_cell_barcode2 \
//         --soloCBposition 0_0_2_-1 3_1_3_8 \
//         --soloUMIposition 3_9_3_14 \
//         --soloUMIdedup 1MM_CR \
//         --soloUMIfiltering MultiGeneUMI_CR \
//         --outSAMtype BAM SortedByCoordinate \
//         --outSAMunmapped Within \
//         --outSAMattributes NH HI AS nM CR CY UR UY CB UB sM \
//         --soloMultiMappers EM \
//         --outFileNamePrefix STAR/
//     """
// }

process StarAlign {
    tag "StarAlign on $subdir"
    label "STAR"

    queue "short"
    memory params.star_mem  // "32 GB"
    time params.star_time  // "1h"
    cpus params.star_threads
    
    publishDir "${params.output_root_dir}/${subdir}", mode: 'copy'


    input:
    val STAR_bin
    val genome_dir
    tuple val(subdir), path(fastqs)
    val outFilterMultimapNmax
    val readFilesCommand
    val quantMode
    val outFilterScoreMinOverLread
    val outFilterMatchNminOverLread
    val outFilterType
    val alignSJoverhangMin
    val alignSJDBoverhangMin
    val outFilterMismatchNmax
    val outFilterMismatchNoverReadLmax
    val alignIntronMin
    val alignIntronMax
    val outSJfilterIntronMaxVsReadN 
    val genomeLoad
    val genomeLoaded
    val join_str // for paired reads use space, for several one-sided reads use ","
    

    output:
        path("STAR_out/*")
        tuple val(subdir), path("STAR_out/Aligned.sortedByCoord.out.bam"), emit: aligned_bams
    

    script:
    
    reads = fastqs.join(join_str)

    if(readFilesCommand.length()) {
        readFilesCommand_line = "--readFilesCommand $readFilesCommand"
    }
    else {
        readFilesCommand_line = ""
    }

    // log.info "sample_name: ${sample_name}"
    // log.info "subdir: ${subdir}"
    // module load gcc/6.2.0
    // module load star/2.7.9a

    """

    mkdir -p STAR_out

    $STAR_bin \
        --readFilesIn $reads \
        --runThreadN $params.star_threads \
        --limitBAMsortRAM 30000000000 \
        --genomeLoad $genomeLoad ${readFilesCommand_line} \
        --outFilterMultimapNmax $outFilterMultimapNmax \
        --quantMode $quantMode \
        --genomeDir $genome_dir \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outFileNamePrefix STAR_out/ \
        --outFilterScoreMinOverLread $outFilterScoreMinOverLread \
        --outFilterMatchNminOverLread $outFilterMatchNminOverLread \
        --outFilterType $outFilterType \
        --alignSJoverhangMin $alignSJoverhangMin \
        --alignSJDBoverhangMin $alignSJDBoverhangMin \
        --outFilterMismatchNmax $outFilterMismatchNmax \
        --outFilterMismatchNoverReadLmax $outFilterMismatchNoverReadLmax \
        --alignIntronMin $alignIntronMin \
        --alignIntronMax $alignIntronMax \
        --outSJfilterIntronMaxVsReadN $outSJfilterIntronMaxVsReadN
    """
}


process StarLoadGenome {
    memory params.star_mem  // "32 GB"
    cpus params.star_threads
    maxForks 1
    tag "StarLoadGenome on $genome_dir"
    label "STAR"

    input:
    val STAR_bin
    val genome_dir

    output:
    val 1
    
    script:
    """
    $STAR_bin \
    --genomeLoad LoadAndExit \
    --genomeDir $genome_dir

    """
}