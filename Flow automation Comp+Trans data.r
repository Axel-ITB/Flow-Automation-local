#Flow automation Comp+Trans data

# Minimal example: load compensated + transformed FCS files and plot FSC vs SSC

# Required libraries
if (!requireNamespace("flowCore", quietly = TRUE)) {
  stop("Package 'flowCore' is required. Please install it from Bioconductor.")
}
if (!requireNamespace("ggcyto", quietly = TRUE)) {
  stop("Package 'ggcyto' is required. Please install it from Bioconductor.")
}

library(flowCore)
library(ggcyto)

# 1. Load all FCS files from the 'comp+trans data' directory
#    Assumes files are already compensated and transformed
fcs_path <- "comp+trans data"
fs_raw <- read.flowSet(path = fcs_path, pattern = "\\.fcs$", alter.names = TRUE)

# 1.1 Check channels
pData(parameters(fs_raw[[1]]))

# 2. Basic filtering to remove debris, doublets and extreme values
fs_filtered <- fsApply(fs_raw, function(ff) {
  expr <- exprs(ff)
  # Remove events with negative FSC/SSC
  keep_pos <- expr[, "FSC.A"] >= 0 & expr[, "SSC.A"] >= 0
  expr <- expr[keep_pos, ]
  # Debris filter
  debris_keep <- expr[, "FSC.A"] > 50000 & expr[, "SSC.A"] > 10000
  # Doublet filter using FSC-A/FSC-H ratio
  ratio <- expr[, "FSC.A"] / (expr[, "FSC.H"] + 1e-6)
  doublet_keep <- ratio > 0.85 & ratio < 1.15
  # Upper bounds for extremely high values
  high_keep <- expr[, "FSC.A"] < 250000 & expr[, "SSC.A"] < 250000
  keep <- debris_keep & doublet_keep & high_keep
  ff[keep, ]
})

# 3. Plot FSC.A vs SSC.A for the first sample
p <- autoplot(fs_filtered[[1]], x = "FITC.A", y = "SSC.A", bins = 128) +
  ggplot2::ggtitle("FSC vs SSC after filtering")
print(p)
