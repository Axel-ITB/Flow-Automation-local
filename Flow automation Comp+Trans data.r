#Flow automation Comp+Trans data

# Minimal example: load compensated  FCS files and plot FSC vs SSC

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
#    These files are compensated but not transformed
fcs_path <- "comp+trans data"
fs_raw <- read.flowSet(path = fcs_path, pattern = "\\.fcs$", alter.names = TRUE, truncate_max_range = FALSE)


# 1.1 Check channels
pData(parameters(fs_raw[[1]]))



# 3. Plot FSC.A vs SSC.A for the first sample
p <- autoplot(fs_raw[[1]], x = "FSC.A", y = "SSC.A", bins = 128) +
  ggplot2::ggtitle("FSC vs SSC after filtering")
print(p)

# 3.1 Histogram of APC.A
p_hist <- autoplot(fs_raw[[1]], "APC.A", bins = 128) +
  ggplot2::ggtitle("Histogram of APC.A after filtering")
print(p_hist)

# 4. Transform fluorescent channels (e.g. FITC.A, APC.A, etc.)
#    Identify fluorescent channels by excluding scatter and time parameters
param_names <- colnames(fs_raw[[1]])
fluor_channels <- grep("\\.A$", param_names, value = TRUE)
fluor_channels <- setdiff(
  fluor_channels,
  c("FSC.A", "FSC.H", "SSC.A", "SSC.H", "Time")
)

# Estimate logicle transformation from the first sample and apply to all
# `estimateLogicle` already returns a transformation list, so we can use it
# directly with `transform()`
logicle <- estimateLogicle(fs_filtered[[1]], fluor_channels)
fs_trans <- transform(fs_filtered, logicle)

# 4. Plot an example fluorescence vs SSC after transformation
example_chan <- fluor_channels[1]
p <- autoplot(fs_trans[[1]], x = example_chan, y = "SSC.A", bins = 128) +
  ggplot2::ggtitle(paste("FSC vs SSC after", example_chan, "transformation"))
print(p)

p_hist <- autoplot(fs_trans[[1]], "FITC.A", bins = 128) +
  ggplot2::ggtitle("Histogram of APC.A after filtering")
print(p_hist)
