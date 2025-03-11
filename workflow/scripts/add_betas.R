################################################################################
# LUCY PICARD
# ---------------------------------- #
# add_betas.r
################################################################################

suppressMessages(library(minfi))
suppressMessages(library(data.table))

fp <- snakemake@input[["pre_beta"]]
# fp = "/home/dejlu879/ProjectProtocol/ext_bm_pipeline_dev/results/test2/mod_calling/test2.pre_beta.csv"

df <- read.csv(fp)

# Add unmeth column
df$unmeth <- df$N - df$X

# Convert to matrix directly (one column for one sample)
Meth <- matrix(df$X, ncol = 1)
Unmeth <- matrix(df$unmeth, ncol = 1)

# Add row and column names
rownames(Meth) <- rownames(Unmeth) <- df$position
colnames(Meth) <- colnames(Unmeth) <- "Sample_1"

# Create MethylSet for a single sample
mset <- MethylSet(Meth = Meth, Unmeth = Unmeth)

# Calculate beta values
beta_values <- getBeta(mset)

# add beta column
df$beta = beta_values

# drop unmeth column
df$unmeth = NULL

write.csv(df, snakemake@output[["post_beta"]], quote = FALSE, row.names = FALSE, )
