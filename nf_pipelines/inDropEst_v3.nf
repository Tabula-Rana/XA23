// export NXF_CONTAINER_ENTRYPOINT_OVERRIDE=true

/* droptag:
    inDrop v3
    File 1: cell barcode, part 1 (default length: 8bp) (R2)
    File 2: cell barcode + UMI, part 1 (default length: >= 14bp) (R4)
    File 3: gene read (R1)
    File 4 (optional): library tag (R3)
*/

// resources
params.concat_mem = "1 GB"
params.concat_time = "3h"

params.droptag_threads = 20
params.droptag_mem = "32 GB"
params.droptag_time = "1h"

params.star_mem = "32 GB"
params.star_threads = 20
params.star_time = "3h"

params.dropest_time = "8h"
params.dropest_mem = "80 GB"

// input
params.read1_nextseq = "OUTPUT_DIR/inDropest_v3/fastq/*/NextSeq/SRR*_1.cutadapt.fastq.paired.fq"
params.read2_nextseq = "OUTPUT_DIR/inDropest_v3/fastq/*/NextSeq/SRR*_2.fastq.paired.fq"
params.read4_nextseq = "OUTPUT_DIR/inDropest_v3/fastq/*/NextSeq/SRR*_4.fastq.paired.fq"
params.read1_novaseq = "OUTPUT_DIR/inDropest_v3/fastq/*/SL-NV[B|K]/Undetermined*_R1_001.{S,N}*.cutadapt.fastq.paired.fq"
params.read2_novaseq = "OUTPUT_DIR/inDropest_v3/fastq/*/SL-NV[B|K]/Undetermined*_R2_001.{S,N}*.fastq.paired.fq"
params.read4_novaseq = "OUTPUT_DIR/inDropest_v3/fastq/*/SL-NV[B|K]/Undetermined*_R4_001.{S,N}*.fastq.paired.fq"

params.reference_root_dir = "ROOT_DIR/xenopus/reference/Xentr10.0"
params.annotation_gtf = "/XENTR_10.0_Xenbase.corrected.gtf"
params.genome_dir = "/STAR_2.7.9a_XENTR_10.0_Xenbase_94"

// params
// droptag
params.indrops_v3_xml = "ROOT_DIR/ksp503/projs/dropEst/configs/indrop_v3.xml"
params.indrops_v3_barcodes = "ROOT_DIR/ksp503/projs/dropEst/data/barcodes/indrop_v3"
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
params.output_root_dir = "OUTPUT_DIR/inDropest_v3_cutadapt/results"


include {
    concat_fastq as concat_fastq_r1;
    concat_fastq as concat_fastq_r2;
    concat_fastq as concat_fastq_r4;
    } from "modules/library.nf"

include {
    StarAlign;
    } from "modules/STAR.nf"

include {
    dropTag_v3;
    dropEst
    } from "modules/dropEst.nf"


workflow {
    genome_dir = params.reference_root_dir + params.genome_dir
    annotation_gtf = params.reference_root_dir + params.annotation_gtf

    // novaseq reads
    Channel.fromPath(params.read1_novaseq)
        .map { infile -> tuple( infile.parent.parent.baseName, infile ) } // lib_name, fastq
        .set { lib_ch_novaseq }
    
    Channel.fromPath(params.read2_novaseq)
        .map { infile -> tuple( infile.parent.parent.baseName, infile ) } // lib_name, fastq
        .set { r2_ch_novaseq }
    
    Channel.fromPath(params.read4_novaseq)
        .map { infile -> tuple( infile.parent.parent.baseName, infile ) } // lib_name, fastq
        .set { r4_ch_novaseq }

    // nextseq reads
    Channel.fromPath(params.read1_nextseq)
        .map { infile -> tuple( infile.parent.parent.baseName, infile ) } // lib_name, fastq
        .set { lib_ch_nextseq }
    
    Channel.fromPath(params.read2_nextseq)
        .map { infile -> tuple( infile.parent.parent.baseName, infile ) } // lib_name, fastq
        .set { r2_ch_nextseq }
    
    Channel.fromPath(params.read4_nextseq)
        .map { infile -> tuple( infile.parent.parent.baseName, infile ) } // lib_name, fastq
        .set { r4_ch_nextseq }

    // runs combined

    lib_ch_novaseq.concat( lib_ch_nextseq )
        .groupTuple()
        .set { lib_ch }
    // lib_ch.view()
    concat_fastq_r1(lib_ch, "R1")

    r2_ch_novaseq.concat( r2_ch_nextseq )
        .groupTuple()
        .set { r2_ch }

    // r2_ch_novaseq.view()

    concat_fastq_r2(r2_ch, "R2")
    
    r4_ch_novaseq.concat( r4_ch_nextseq )
        .groupTuple()
        .set { r4_ch }

    concat_fastq_r4(r4_ch, "R4")
        
    lib_r2_r4_ch = concat_fastq_r1.out.join(concat_fastq_r2.out).join(concat_fastq_r4.out)

    // lib_r2_r4_ch.view()

    dropTag_v3(
        lib_r2_r4_ch,
        params.indrops_v3_xml
        )
    
    // dropTag_v3.out.tagged_reads.groupTuple(size: 3, by: 0).set { tagged_reads }

    // tagged_reads = tagged_reads.map(it -> tuple(it[0], it[1].flatten()))

    // tagged_reads.view()

    StarAlign(
        params.STAR_bin,
        genome_dir,
        dropTag_v3.out.tagged_reads,
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

    // dropTag_v3.out.tagged_params.groupTuple(size: 3, by: 0).set { tagged_params }

    dropEst_inp = StarAlign.out.aligned_bams
        .join(dropTag_v3.out.tagged_params)
   
    dropEst(
        dropEst_inp,
        annotation_gtf,
        params.indrops_v3_xml,
        params.indrops_v3_barcodes,
        params.gene_parts
        )
}