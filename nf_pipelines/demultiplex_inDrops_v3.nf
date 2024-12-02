/*
* File 1: cell barcode, part 1 (default length: 8bp)
* File 2: cell barcode + UMI, part 1 (default length: >= 14bp)
* File 3: gene read
* File 4 (optional): library tag
*/

// resources
params.cutadapt_threads = 20
params.cutadapt_time = "40m"
params.cutadapt_mem = "12 GB"
params.seqtk_mem = "62 GB"
params.seqtk_time = "6h"

// tools
params.bin_root_dir = "ROOT_DIR/tools"

// input
params.data_root_dir = "ROOT_DIR/XenopusCellAtlas/20180611_XTr_rep2_NovaSeq/raw_data2/seq/illumina/proc/SL-NVK/180716_SL-NVK_0046_BHCYVKDMXX/Data/Intensities/BaseCalls"
params.lane = "L002"
params.barcodes = "CODE_DIR/nf_pipelines/inputs/pool2.fa"

// output
params.output_root_dir = "ROOT_DIR/xenopus/old_atlas/inDropest_v3/SL-NVK"

include {
    demultiplex_barcodes_inDrops_v3;
    demultiplex_libraries_inDrops_v3
    } from "modules/scRNAseq.nf"


workflow {
    seqtk_bin = params.bin_root_dir + "/seqtk/seqtk"

    Channel.fromPath("${params.data_root_dir}/Undetermined*${params.lane}_R3_001.fastq.gz")
        .map { infile -> tuple( params.lane, infile ) }
        .set { demutli_ch }
    
    // demutli_ch.view()
    
    demultiplex_libraries_inDrops_v3(demutli_ch, params.barcodes)
    
    Channel.fromPath("${params.data_root_dir}/Undetermined*${params.lane}_R1_001.fastq.gz")
        .map { infile -> tuple( params.lane, infile ) }
        .groupTuple()
        .set { lib_ch }
    
    Channel.fromPath("${params.data_root_dir}/Undetermined*${params.lane}_R2_001.fastq.gz")
        .map { infile -> tuple( params.lane, infile ) }
        .groupTuple()
        .set { r2_ch }
    
    Channel.fromPath("${params.data_root_dir}/Undetermined*${params.lane}_R4_001.fastq.gz")
        .map { infile -> tuple( params.lane, infile ) }
        .groupTuple()
        .set { r4_ch }

    r4_r2_lib_ch = r4_ch.join(r2_ch).join(lib_ch)

    // demultiplex_libraries_inDrops_v3.out.view()
    inp_ch = r4_r2_lib_ch.combine(
        demultiplex_libraries_inDrops_v3.out.transpose(), by: 0)

    // inp_ch.view()

    demultiplex_barcodes_inDrops_v3(
        seqtk_bin,
        inp_ch
    )
}