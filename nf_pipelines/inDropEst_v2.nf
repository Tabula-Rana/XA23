// export NXF_CONTAINER_ENTRYPOINT_OVERRIDE=true

/* droptag:
    inDrop v1 & v2
    File 1: barcode reads. (R1)
        Structure:
        Cell barcode, part 1
        Spacer
        Cell barcode, part 2
        UMI
    File 2: gene reads (R2)
*/

// resources
params.star_threads = 20
params.star_mem = "32 GB"
params.star_time = "1h"

params.droptag_threads = 20
params.droptag_mem = "32 GB"
params.droptag_time = "1h"

params.dropest_mem = "32 GB"
params.dropest_time = "1h"

// input
params.read1 = "OUTPUT_DIR/inDropest_v2/fastq/*/SRR*/SRR*_1.cutadapt.fastq.paired.fq" // gene read
params.read2 = "OUTPUT_DIR/inDropest_v2/fastq/*/SRR*/SRR*_2.fastq.paired.fq" // barcode read

// reference
params.reference_root_dir = "ROOT_DIR/xenopus/reference/Xentr10.0"
params.annotation_gtf = "/XENTR_10.0_Xenbase.corrected.gtf"
params.genome_dir = "/STAR_2.7.9a_XENTR_10.0_Xenbase_33"

// programs' params
// droptag
params.indrops_v2_xml = "ROOT_DIR/ksp503/projs/dropEst/configs/indrop_v1_2.xml"
params.indrops_v2_barcodes = "ROOT_DIR/ksp503/projs/dropEst/data/barcodes/indrop_v1_2"
// STAR
params.STAR_bin = "STAR"
params.outFilterMultimapNmax = 10
params.readFilesCommand = "zcat" // zcat
params.quantMode = "-"
params.outFilterScoreMinOverLread = 0.66 // alignment will be output only if its score is higher than or equal to this value, normalized to read length 
params.outFilterMatchNminOverLread = 0.66
params.outFilterType = "BySJout" // reduces the number of spurious junctions
params.alignSJoverhangMin = 8 // minimum overhang for unannotated junctions
params.alignSJDBoverhangMin = 1 // minimum overhang for annotated junctions
params.outFilterMismatchNmax = 999 // maximum number of mismatches per pair, large number switches off this filter
params.outFilterMismatchNoverReadLmax = 0.05 // max number of mismatches per pair relative to read length: for 2x100b, max number of mis-matches is 0.04*200=8 for the paired read
// 0.05 * 120 = 6
params.alignIntronMin = 20 // minimum intron length
params.alignIntronMax = 100000 // maximum intron length 
params.outSJfilterIntronMaxVsReadN = "50000 100000 200000"
params.genomeLoad = "NoSharedMemory" // for slurm/standalone
params.join_str = ","
// dropEst
params.gene_parts = "eiEIBA" // do count UMIs with both exon and inton reads

params.genome_remove = false // cluster

// output
params.output_root_dir = "OUTPUT_DIR/inDropest_v2_cutadapt/results"


include {
    StarAlign;
    } from "modules/STAR.nf"

include {
    dropTag_v2;
    dropEst
    } from "modules/dropEst.nf"


workflow {

    genome_dir = params.reference_root_dir + params.genome_dir
    annotation_gtf = params.reference_root_dir + params.annotation_gtf

    read1 = Channel.fromPath(params.read1)
            .map { it -> tuple("${it.parent.parent.baseName}/${it.parent.baseName}", it) } // stage/srr, file
        
    read2 = Channel.fromPath(params.read2)
        .map { it -> tuple("${it.parent.parent.baseName}/${it.parent.baseName}", it) } // stage/srr, file
    
    data_ch = read1.join(read2)

    // data_ch.view()

    dropTag_v2(
        data_ch,
        params.indrops_v2_xml
    )
    
    StarAlign(
        params.STAR_bin,
        genome_dir,
        dropTag_v2.out.tagged_reads,
        params.outFilterMultimapNmax,
        params.readFilesCommand,
        params.quantMode,
        params.outFilterScoreMinOverLread,
        params.outFilterMatchNminOverLread,
        params.outFilterType,
        params.alignSJoverhangMin,
        params.alignSJDBoverhangMin,
        params.outFilterMismatchNmax,
        params.outFilterMismatchNoverReadLmax,
        params.alignIntronMin,
        params.alignIntronMax,
        params.outSJfilterIntronMaxVsReadN,
        params.genomeLoad,
        false,
        params.join_str
    )

    dropEst_inp = StarAlign.out.aligned_bams
        .join(dropTag_v2.out.tagged_params)

    // log.info "annotation_gtf: ${annotation_gtf}"
    
    dropEst(
        dropEst_inp,
        annotation_gtf,
        params.indrops_v2_xml,
        params.indrops_v2_barcodes,
        params.gene_parts
    )
}