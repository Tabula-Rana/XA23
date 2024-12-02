// module load java/jdk-17.0.7
/*
* takes 2 lists of fastq files
* combines them by file name
* returns files with paired reads
*/


// input/output
// // v2
// params.input_type = "v2"
// params.fastq1 = "OUTPUT_DIR/inDropest_v2/fastq/*/*/*_1.cutadapt.fastq"
// params.fastq2 = "ROOT_DIR/XenopusCellAtlas/Briggs_v2/*/*/*_2.fastq"
// params.output_root_dir = "OUTPUT_DIR/inDropest_v2/fastq"
// // v3_nextseq
// params.input_type = "v3_nextseq"
// params.fastq1 = "OUTPUT_DIR/inDropest_v3/fastq/*/NextSeq/*_1.cutadapt.fastq"
// params.fastq2 = "ROOT_DIR/XenopusCellAtlas/20180109_XTr_rep2_NextSeq/*/*/*_4.fastq"
// params.output_root_dir = "OUTPUT_DIR/inDropest_v3/fastq"
// v3_novaseq
params.input_type = "v3_novaseq"
params.fastq1 = "OUTPUT_DIR/inDropest_v3/fastq/*/SL-NV[B|K]/Undetermined_S0_L00[1|2]_R1_001.[S|N]*.cutadapt.fastq"
params.fastq2 = "ROOT_DIR/xenopus/old_atlas/inDropest_v3/fastq/demultiplexed/*/*/Undetermined_S0_L00[1|2]_R2_001.[S|N]*.fastq"
params.output_root_dir = "OUTPUT_DIR/inDropest_v3/fastq"

// params
params.fastq_pair_cl = "-t 10000000 -s"
params.fastq_pair_bin = "TOOLS_DIR/bin/fastq_pair"

// resources
params.fastq_paired_mem = "64 GB"
params.fastq_paired_time = "3h"
params.fastq_paired_threads = "1"


include {
    fastq_paired;
    } from "modules/scRNAseq.nf"


workflow {
    
    if (params.input_type == "v2") {
        fastq1s = Channel.fromPath(params.fastq1)
            .map { it -> tuple("${it.parent.parent.baseName}/${it.parent.baseName}", it) } // stage/srr, file
        
        fastq2s = Channel.fromPath(params.fastq2)
            .map { it -> tuple("${it.parent.parent.baseName}/${it.parent.baseName}", it) } // stage/srr, file
    }
    else if (params.input_type == "v3_nextseq") {
        fastq1s = Channel.fromPath(params.fastq1)
            .map { it -> tuple("${it.parent.parent.baseName}/NextSeq", it) } // lib_name/NextSeq, file
        
        fastq2s = Channel.fromPath(params.fastq2)
            .map { it -> tuple("${it.parent.baseName}/NextSeq", it) } // lib_name/NextSeq, file
    }
    else if (params.input_type == "v3_novaseq") {
        fastq1s = Channel.fromPath(params.fastq1)
            .map { it -> tuple("${it.parent.parent.baseName}/${it.parent.baseName}", it) } // lib_name/run, file
        
        fastq2s = Channel.fromPath(params.fastq2)
            .map { it -> tuple("${it.baseName.split('\\.')[1]}/${it.parent.parent.baseName}", it) } // lib_name/run, file
    }
    else {
        fastq1s = null
        fastq2s = null
    }

    fastqs = fastq1s.join(fastq2s)

    // fastqs.view()

    fastq_paired(
        fastqs,
        params.fastq_pair_bin,
        params.fastq_pair_cl
    )
}
