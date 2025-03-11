suppressMessages(library(minfi))
# BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
# library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
# library(IlluminaHumanMethylation450kmanifest)
suppressMessages(library(ggplot2))
suppressMessages(library(MethylCIBERSORT))

# MethylCIBERSORT
# https://yuantian1991.github.io/notes/Using-MethylCIBERSORT-for-Cell-Type-Deconvolution

# load signatures
data("V2_Signatures")

get_cibersort_upload_files <- function(path2betamatrix, mixture_matrix_name, base_sig_matrix, base_sig_matrix_filename, sample_name) {
    ###########################################################################
    # MethylCIBERSORT
    print("==================================")
    print("Start MethylCIBERSORT processing")

    # Get beta values for matrix samples
    # mat <- getBeta(mSetSq)
    # mat <- getBeta(gset)

    mat <- read.csv(path2betamatrix)

    # mat <- mat[!duplicated(mat), ]  # this removes dupilcated whole rows
    # Note this results in duplicate values, pretty sure due to different values on opposite
    # DNA strands... but would need to investigate to make sure.
    # mat[duplicated(mat$probe) | duplicated(mat$probe, fromLast = TRUE), ]  # - this checks for duplicated probes

    # mat <- mat[!duplicated(mat$probe), ]  # This removes duplications for probes regardless of value
    rownames(mat) <- mat$probe
    mat$probe <- NULL
    # colnames(mat) = "Test2"
    colnames(mat) <- sample_name
    mat = as.matrix(mat)

    ###########################################################################
    # Export signature to file

    # Move into results dir
    # setwd(path2_results)

    write.table(
        base_sig_matrix,
        file = base_sig_matrix_filename,
        sep = "\t", row.names = FALSE, quote = FALSE
    )

    # Export mixture matrix to file
    # print(head(mat))
    # print(head(base_sig_matrix))
    # print(base_sig_matrix$NAME)
    load("/home/dejlu879/ProjectProtocol/bm_pipeline_dev/resources/mat.RData")
    Prep.CancerType(
        Beta = mat,
        Probes = base_sig_matrix$NAME,
        fname = mixture_matrix_name
    )

    # Log where files are
    #     print(paste0(
    #       "CIBERSORT signature file: ", file.path(
    #         path2_results,
    #         paste0(sigName, ".txt")
    #       )
    #     ))

    #     print(paste0(
    #       "CIBERSORT matrix file: ", file.path(
    #         path2_results,
    #         paste0(mixture_matrix, ".txt")
    #       )
    #     ))
}

# path2targets <- file.path(cwd, paste0("raw_data/TCGA_", tissue))
# path2data <- file.path(cwd, paste0("raw_data/TCGA_", tissue, "/rdata"))
# targets_fp <- paste0(tissue, '_target_sheet.csv')
# path2_results <- file.path(cwd, paste0('results/TCGA_', tissue))
# qc_report_pdf_name = paste0("TCGA_", tissue, "_qcReport.pdf")
# logfile = paste0(date, "-metharray.", tissue, ".log.txt")

# Log parameters

# path2betamatrix = "/home/dejlu879/ProjectProtocol/ext_bm_pipeline_dev/results/test2/mod_calling/test2.methatlas.csv"

# Get arguments from commandline
args <- commandArgs(trailingOnly = TRUE)
input_beta_matrix <- args[1]
cibersort_mix_matrix_upload_fname <- args[2]
cibersort_bladder_ref_upload_fname <- args[3]
sample_name <- args[4]

get_cibersort_upload_files(
    path2betamatrix = input_beta_matrix,
    mixture_matrix_name = cibersort_mix_matrix_upload_fname,
    base_sig_matrix = Signatures$bladder_v2_Signature.txt,
    base_sig_matrix_filename = cibersort_bladder_ref_upload_fname,
    sample_name
)
