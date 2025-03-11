# Set up
library(minfi)
# BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kmanifest)

# set up paths
cwd <- "/home/dejlu879/ProjectProtocol/20240123-ProjectProtocolDev/BLCA_example_idats"
setwd(cwd)
path2samplesheet_BLCA <- "gdc_sample_sheet.2022-06-16.tsv"
# path2samplesheet_SKCM <- "raw_data/TCGA_SKCM/gdc_sample_sheet.2022-06-20.tsv"
path2BLCAdata <- file.path(cwd)
# path2SKCMdata <- file.path(cwd, "raw_data/TCGA_SKCM")

# Process data
# base<-file.path("/Users/lucy/Documents/immune_infiltrate/raw_data/TCGA_BLCA/rdata")
base <- file.path(cwd)
# base<-file.path("/Users/lucy/Documents/immune_infiltrate/raw_data/TCGA_SKCM")

# set data
meth.anno <- IlluminaHumanMethylation450kanno.ilmn12.hg19

###########################################################################
# Minfi process SKCM data

# Load in sample_sheet
# sample_sheet_SKCM = read.table(path2samplesheet_SKCM, header=T, sep='\t')
targets <- read.metharray.sheet(base = path2BLCAdata, pattern = "BLCA_target_sheet.csv")

# create Large RGChannelSet - raw data from array
# failed because no matching in SKCM data
# using python, I moved all SKCM files into
# /Users/lucy/Documents/immune_infiltrate/raw_data/TCGA_SKCM/rdata
# and generated a targets file at the same time
# NOTE: below line worked when targets=NULL
# Make sure targets exists, otherwise downstream will break
if (is.NULL(targets)) {
    print("targets file does not exist, have not made rgSet.")
} else {
    rgSet <- read.metharray.exp(base, targets = targets, force = TRUE, recursive = TRUE)
}

# rgSet rows = probes
# rgSet cols = samples

library(minfi)

# Set the directory containing the IDAT files
# idat_dir <- "path/to/idat/files"

# Get the base name (without _Grn.idat or _Red.idat)
# ext_bm_pipeline_dev/BLCA_example_idats/9565ff11-956a-4c92-83f4-114e7557e1e1_noid_Grn.idat
sample_base <- "9565ff11-956a-4c92-83f4-114e7557e1e1_noid"

# Read IDAT files into an RGChannelSet
rgSet <- read.metharray(basenames = sample_base)

# Check the RGset object

rgSet

# Explore
green_rgSet <- getGreen(rgSet)

# ================================================== #
# Minfi preprocessing steps
# Preprocess (normalisation within array -?)

# ===== Remove failed samples ====== #
# Step 1
detP <- detectionP(rgSet)

# Step 2: look at detP info - want to record this
table(colMeans(detP) < 0.01)
# First try - 0.05
# TRUE
# 475
# Second try - 0.01
# FALSE  TRUE
# 1   474

head(detP)
colSums(detP)
dim(detP)

# Step 3: remove poor quality samples (where too many probes didnt work)
# from both rgSet and targets file
keep <- colMeans(detP) < 0.01 # returns TRUE/FALSE list
rgSet_filt <- rgSet[, keep] # samples to filter are in COLUMNS
# NOTE: did not remove anything during devpt, possibly has already been filtered by TCGA
targets <- targets[keep, ] # samples to filter are in ROWS

# Step 4: Check what targets now looks like
targets[, 1:5]

# Step 5: do preprocessing, outputs a MethylSet
mSetSq <- preprocessNoob(rgSet_filt) # methylcibersort paper used noob

# Step 6: remove cross-reacting probes etc
# ONLY available on EPIC arrays

###################################
# NOTE check this block of code
# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq), rownames(detP)), ]

# remove any probes that have failed in 50% or more samples
keep <- rowSums(detP < 0.5) == ncol(mSetSq)
table(keep) # look at keep result
mSetSqFlt <- mSetSq[keep, ]
mSetSqFlt

# NOTE: sex check, may not be available
keep <- !(featureNames(mSetSq) %in% meth.anno$Name[meth.anno$CHR %in% c("X", "Y")])
table(keep)
mSetSqFlt <- mSetSq[keep, ]

# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSq)
mSetSqFlt

###################################

##### CHECKED HISTORY OF ALL BELOW #####

# ====== Plot P values ====== #
# pal <- brewer.pal(8,"Dark2")
# barplot(colSums(detP), las=2,
#         cex.names=0.8, ylab="Mean detection p-values")
#
# colnames(detP)<-targets$Name
# greycols<-gray.colors(4, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL)
# data<-colMeans(detP)
# labels_vec<-colnames(detP)
# rot_angle<-45
# barplot(colMeans(detP), las=2,
#         cex.names=0.8, ylab="Mean detection p-values")
# rotate_x <- function(data,labels_vec, rot_angle) {
#   plt <- barplot(data, xaxt="n")
#   text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1,1), xpd = TRUE, cex=1)
# }
# rotate_x(data,colnames(detP), 45)
#
# abline(h=0.05,col="red")
# legend("top", legend=levels(factor(targets$Exposure)), fill=greycols,
#        bg="white",horiz = FALSE,bty='n', inset=c(1.50,0), xpd = TRUE)

# ====== QC report====== #
# qc_report_pdf_name = "TCGA_SKCM_qcReport.pdf"
# qcReport(rgSet, sampNames=targets$Name, sampGroups=NULL,
#          pdf=qc_report_pdf_name)

# ======================= #

# Quantile normalisation - don't do for cibersort (recommended by methylcibersort paper)
# This function implements stratified quantile normalization preprocessing for
# Illumina methylation microarrays.
# Probes are stratified by region (CpG island, shore, etc.)

# mSetSq <- preprocessQuantile(rgSet, fixOutliers = TRUE, quantileNormalize = TRUE)

# ================================================== #
# Genomic locations
# load in annotation info
data(IlluminaHumanMethylation450kanno.ilmn12.hg19) # (?)
# IlluminaHumanMethylation450kanno.ilmn12.hg19@defaults  # etc to look at stuff

# get annotations of genomic locations of probes
gset <- mapToGenome(mSetSq)
annotation <- getAnnotation(gset)
# ================================================== #

###########################################################################
# MethylCIBERSORT
# https://yuantian1991.github.io/notes/Using-MethylCIBERSORT-for-Cell-Type-Deconvolution
library(MethylCIBERSORT)

# select reference dataset for methylcibersort, loaded as Stromal_v2
data("StromalMatrix_V2")

# Get beta values for matrix samples
Mat <- getBeta(mSetSq)
write.table(Mat, file = "Mat_raw.csv", sep = ",", quote = FALSE)
length(Mat)
# [1] 485512

###########################################################################
# Do the code I have in the pipeline methylcibersort.r

mat <- read.csv("Mat_raw.csv")
rownames(mat) <- mat$probe
mat$probe <- NULL
# colnames(mat) = "Test2"
colnames(mat) <- sample_name



###########################################################################

# select only the probes we have in our matrix and are available in
# methylcibersort stromal ref data
Int <- intersect(rownames(Mat), rownames(Stromal_v2))
Mat_matched <- Mat[match(Int, rownames(Mat)), ]
Stromal_matched <- Stromal_v2[match(Int, rownames(Stromal_v2)), ]

# Check how many probes were removed
dim(Int)
dim(Stromal_v2)
dim(Stromal_matched)
dim(Mat)
dim(Mat_matched)

# Export Signature matrix file for reference on CIBERSORT
# Note that instead of doing this, we can use the pre-prepared signature files
# of each cancer/disease supplied by authors
sigName <- "ExampleType"
RefData <- Stromal_matched # reference matrix matched to our filtered data
RefPheno <- Stromal_v2.pheno # list of cell types we are looking for
Signature <- FeatureSelect.V4(
    CellLines.matrix = NULL,
    Heatmap = FALSE,
    export = TRUE,
    sigName = sigName,
    Stroma.matrix = RefData,
    deltaBeta = 0.2,
    FDR = 0.01,
    MaxDMRs = 100,
    Phenotype.stroma = RefPheno
)

# Below code will generate the "mixture matrix" that is uploaded to CIBERSORT
fname <- "ExportData" # will make a file called "ExportData.txt"
Prep.CancerType(
    Beta = Mat,
    Probes = rownames(Signature$SignatureMatrix),
    fname = fname
)


# Checks
dim(Int)
dim(Mat) #
dim(Stromal_v2) #
dim(Mat_matched) # [1] 485512    475
dim(Stromal_matched) # [1] 485512     56

# ======================= #
# Note that there are cancer specific signatures available in methylcibersort
# check paper for appropriateness of use

# load signatures
data("V2_Signatures")
# look at available signatures
names(Signatures)
