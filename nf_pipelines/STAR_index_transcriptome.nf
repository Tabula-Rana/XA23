// env
params.star_threads = 20
params.star_time = "60m"
params.star_memory = "64 GB"

params.reference_root_dir = "ROOT_DIR/xenopus/reference/Xentr10.0"
params.STAR_overhang = "33"
params.genome_fasta = "/XENTR_10.0_genome.fasta"
params.genome_dir = "/STAR_2.7.9a_XENTR_10.0_Xenbase"
params.annotation_gtf = "/XENTR_10.0_Xenbase.corrected.gtf"
params.STAR_bin = "STAR"
params.sjdbGTFtagExonParentGene = "gene_id"
params.sjdbGTFtagExonParentTranscript = "transcript_id"

// params.reference_root_dir = "ROOT_DIR/xenopus/reference/hg38_HCMV"
// params.STAR_overhang = "74"
// params.genome_fasta = "/hg38_HCMV.fa"
// params.genome_dir = "/STAR_2.7.9a_hg38_HCMV.knownGene"
// params.annotation_gtf = "/hg38_HCMV.knownGene.gtf"
// params.STAR_bin = "STAR"
// params.sjdbGTFtagExonParentGene = "gene_id"
// params.sjdbGTFtagExonParentTranscript = "transcript_id"

include {
    IndexTranscriptome;
    } from "modules/STAR.nf"

workflow {
    annotation_gtf = params.reference_root_dir + params.annotation_gtf
    genome_dir = params.reference_root_dir + params.genome_dir + "_" + params.STAR_overhang
    genome_fasta = params.reference_root_dir + params.genome_fasta
    
    IndexTranscriptome(
        genome_fasta,
        genome_dir,
        annotation_gtf,
        params.STAR_overhang,
        params.STAR_bin,
        params.sjdbGTFtagExonParentGene,
        params.sjdbGTFtagExonParentTranscript
    )
}