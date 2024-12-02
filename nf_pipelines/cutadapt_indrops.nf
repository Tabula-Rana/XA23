// module load java/jdk-17.0.7

/*
    inDrop v1 & v2
    File 1: barcode reads. (R1)
        Structure:
        Cell barcode, part 1
        Spacer
        Cell barcode, part 2
        UMI
    File 2: gene reads (R2)

    inDrop v3
    File 1: cell barcode, part 1 (default length: 8bp) (R2)
    File 2: cell barcode + UMI, part 1 (default length: >= 14bp) (R4)
    File 3: gene read (R1)
    File 4 (optional): library tag (R3)
*/


// input/output
// // v2
// params.input_type = "v2"
// params.data_root_dir = "ROOT_DIR/XenopusCellAtlas/Briggs_v2/"
// params.input = "*/*/*_1.fastq"  // R1 - lib, R2 - barcode
// params.output_root_dir = "OUTPUT_DIR/inDropest_v2/fastq"

//v3 NextSeq
params.input_type = "v3_nextseq"
params.data_root_dir = "ROOT_DIR/XenopusCellAtlas/20180109_XTr_rep2_NextSeq/"
params.input = "*/*/*_1.fastq"
params.output_root_dir = "OUTPUT_DIR/inDropest_v3/fastq"

// //v3 NovaSeq
// params.input_type = "v3_novaseq"
// params.data_root_dir = "ROOT_DIR/xenopus/old_atlas/inDropest_v3/fastq/demultiplexed/*/*/"
// params.input = "Undetermined_S0_L00[1|2]_R1_001.[S|N]*.fastq"
// params.output_root_dir = "OUTPUT_DIR/inDropest_v3/fastq"


// params

// polyA
// polyG + quality triming
// adapter trimming
// min_length threshold

params.cutadapt_cl = "-a file:CODE_DIR/nf_pipelines/inputs/illumina_nextseq_p7.fasta \
    --nextseq-trim=20 \
    --poly-a \
    -m 20"

// resources
params.cutadapt_mem = "8 GB"
params.cutadapt_time = "40m"
params.cutadapt_threads = "20"


include {
    cutadapt;
    } from "modules/library.nf"


workflow {
    

    if (params.input_type == "v2") {
        
        fastqs = Channel.fromPath("${params.data_root_dir}${params.input}")
        .map { it -> tuple(it.parent.toString().replace("${params.data_root_dir}", ""), "${it.baseName}.cutadapt.fastq", it) } // subdir (stage/srr), trimmed_filename, filepath
    }
    else if (params.input_type == "v3_novaseq") {
        fastqs = Channel.fromPath("${params.data_root_dir}${params.input}")
        .map { it -> tuple("${it.baseName.split('\\.')[1]}/${it.parent.parent.baseName}", "${it.baseName}.cutadapt.fastq", it) } // subdir (lib_name/run), trimmed_filename, filepath
    }
    else if (params.input_type == "v3_nextseq") {
        fastqs = Channel.fromPath("${params.data_root_dir}${params.input}")
        .map { it -> tuple("${it.parent.baseName}/NextSeq", "${it.baseName}.cutadapt.fastq", it) } // subdir (lib_name/run), trimmed_filename, filepath
    }
    else {
        fastqs = null
    }

    // fastqs.view()

    cutadapt(
        fastqs,
        params.cutadapt_cl
    )
}
