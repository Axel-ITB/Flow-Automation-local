# Load the flowCore and flowStats libraries for flow cytometry analysis
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("flowCore", quietly = TRUE)) {
  BiocManager::install("flowCore")
}
if (!requireNamespace("flowStats", quietly = TRUE)) {
  BiocManager::install("flowStats")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}
if (!requireNamespace("ggcyto", quietly = TRUE)) {
  BiocManager::install("ggcyto")
}
if (!requireNamespace("flowDensity", quietly = TRUE)) {
  BiocManager::install("flowDensity")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  BiocManager::install("pheatmap")
}
if (!requireNamespace("MASS", quietly = TRUE)) {
  install.packages("MASS")
}
library(flowCore)
library(flowStats)
library(flowViz)
library(ggplot2)
library(ggcyto)
library(flowDensity)
library(MASS)
library(CATALYST)
library(writexl)

#setwd("C:/Users/AxelBergenstråle/OneDrive - ITB-MED AB/Desktop/R Projects/FACS Automation local")

#Steps
# 1. Load and simple filtering data. Done.
# 2. Load controls. Done.
# Transform flourescence channels for samples. done.
# On sample: Calculate APC.A cutoff for back gating, using density peak detection gating.
# Using APC.A back gate the SSC and FSC channels to remove debris.
# 3. Compute spillover matrix, on nontransformed data.
# 4. Apply compensation matrix to experimental samples
# 5. Transform data
# 6. Convert to CATALYST format
# 7. Run CATALYST clustering and visualization


# STEP 1: Load your experimental data
# -----------------------------------

# Read all FCS files (samples) from the directory 'fcs_files/' into a flowSet object # nolint: line_length_linter. # nolint
fcs_path <- "fcs_files"

fs_raw <- read.flowSet(path = fcs_path, pattern = "\\.fcs$", alter.names = TRUE)
#Filter out events with negative FSC.A or SSC.A in all samples of the flowSet
#This is not a must but these events are garbage, and it will mess with future plots.
fs_positive <- fsApply(fs_raw, function(fs_raw) {
  expr <- exprs(fs_raw)
  keep <- expr[, "FSC.A"] >= 0 & expr[, "SSC.A"] >= 0
  fs_raw[keep, ]
})
#Automated debris and doublet filtering for all samples in your flowSet
fs <- fsApply(fs_positive, function(ff) {
  expr <- exprs(ff)
  # 1. Remove debris (adjust thresholds as needed)
  debris_keep <- expr[, "FSC.A"] > 50000 & expr[, "SSC.A"] > 10000
  # 2. Remove doublets (FSC-A/FSC-H ratio close to 1)
  ratio <- expr[, "FSC.A"] / (expr[, "FSC.H"] + 1e-6)
  doublet_keep <- ratio > 0.85 & ratio < 1.15
  # 3. Remove extreme high values (hard upper threshold)
  high_keep <- expr[, "FSC.A"] < 250000 & expr[, "SSC.A"] < 250000   # NEW: Remove extreme high values
  # Combine filters
  keep <- debris_keep & doublet_keep & high_keep
  ff[keep, ]
})
#Some code to test the filtering. Printed values will most likely go down.
stats <- list(
    sum_FSC_A_raw_neg = sum(exprs(fs_raw[[5]])[, "FSC.A"] < 0),
    sum_FSC_A_fs_neg  = sum(exprs(fs[[5]])[, "FSC.A"] < 0),
    mean_FSC_A_raw    = mean(exprs(fs_raw[[5]])[, "FSC.A"]),
    mean_FSC_A_fs     = mean(exprs(fs[[5]])[, "FSC.A"]),
    sum_SSC_A_raw_neg = sum(exprs(fs_raw[[5]])[, "SSC.A"] < 0),
    sum_SSC_A_fs_neg  = sum(exprs(fs[[5]])[, "SSC.A"] < 0),
    sum_FSC_A_raw     = sum(exprs(fs_raw[[5]])[, "FSC.A"]),
    sum_FSC_A_fs      = sum(exprs(fs[[5]])[, "FSC.A"])
  )
stats
#Use fs for all downstream analysis!
#Remove fs_raw and fs_positive from the global environment
rm(fs_raw, fs_positive, fcs_path)


# STEP 2: Load compensation control files
# ---------------------------------------
# Automatically list and load all FCS control files (unstained + single-stained controls)
control_files <- list.files(
  path       = "fcs_files/controls/",
  pattern    = "\\.fcs$",
  full.names = TRUE
)

# Create a flowSet containing only the control samples
fs_controls <- read.flowSet(files = control_files, alter.names = TRUE)

# STEP 3: Assign meaningful sample names
# ---------------------------------------
# Provide clear, ordered names for your control tubes:
# Check current names
sampleNames(fs_controls)

# Define the mapping (pattern-based)
new_names <- rep(NA, length(fs_controls))
new_names[grep("Unstained", sampleNames(fs_controls), ignore.case = TRUE)] <- "Unstained"
new_names[grep("APC", sampleNames(fs_controls), ignore.case = TRUE)] <- "APC"
new_names[grep("BV786", sampleNames(fs_controls), ignore.case = TRUE)] <- "BV786"
new_names[grep("BV510", sampleNames(fs_controls), ignore.case = TRUE)] <- "BV510"
new_names[grep("BB515", sampleNames(fs_controls), ignore.case = TRUE)] <- "BB515"
new_names[grep("\\bPE", sampleNames(fs_controls), ignore.case = TRUE)] <- "PE"

# Check if assignment was successful
print(new_names)

# Assign back to sampleNames safely
sampleNames(fs_controls) <- new_names

#Need to start gating for FSC and SSC here,

channels <- c("APC.A", "BV786.A", "BV510.A", "BB515.A", "PE.A")

logicle <- logicleTransform(w = 0.5, t = 262143, m = 4.5, a = 0)
trans_list <- transformList(channels, logicle)
fs_noncomp_trans <- transform(fs, trans_list)
fs_controls_noncomp_trans <- transform(fs_controls, trans_list)

#Plot APC.A (CD4)
autoplot(fs_noncomp_trans[[3]], "APC.A") +
  ggtitle("APC-A Histogram") +
  xlab("APC-A") +
  ylab("Count")

# On sample: Calculate APC.A cutoff for back gating, using density peak detection gating.
# Initialize a list to store gated flowFrames
APC_positive_list <- list()
thresholds <- numeric(length(fs_noncomp_trans))

# Loop over each sample in the flowSet
for (i in seq_along(fs_noncomp_trans)) {
  fr <- fs_noncomp_trans[[i]]
  sample_name <- sampleNames(fs_noncomp_trans)[i]

  # Automatic threshold using density gating
  threshold <- deGate(fr, channel = "APC.A", all.cuts = FALSE)
  thresholds[i] <- threshold  # store threshold for reference

  # Subset to APC+ cells
  fr_APC_pos <- fr[exprs(fr)[, "APC.A"] > threshold, ]

  APC_positive_list[[sample_name]] <- fr_APC_pos

  dens <- density(exprs(fr)[,"APC.A"])
  plot(dens$x, dens$y, type = "l", main = paste("APC.A Density:", sample_name))
  abline(v = threshold, col = "red", lty = 2)

}
fs_APC_pos <- flowSet(APC_positive_list)

autoplot(fs_APC_pos[[1]], x = "FSC.A", y = "SSC.A", bins = 128)
autoplot(fs_noncomp_trans[[5]], x = "FSC.A", y = "SSC.A", bins = 128)
# Using APC.A back gate the SSC and FSC channels to remove debris.

# fs_noncomp_trans: original (transformed, non-compensated)
# fs_APC_pos: APC+ gated cells

# List to hold the new back-gated populations
backgated_list <- list()
backgated_list_cells_in_polygon <- list()

for (i in seq_along(fs_APC_pos)) {
  apc_pos_fr <- fs_APC_pos[[i]]
  name       <- sampleNames(fs_APC_pos)[i]
  expr       <- exprs(apc_pos_fr)
  fsc        <- expr[, "FSC.A"]
  ssc        <- expr[, "SSC.A"]

  #--- density-based backgate ----------------------------------------------
  dens <- MASS::kde2d(fsc, ssc, n = 200)
  z_sorted  <- sort(as.vector(dens$z), decreasing = TRUE)
  cum_z     <- cumsum(z_sorted)
  level_idx <- which(cum_z >= 0.90 * sum(z_sorted))[1]   # contour enclosing ~90%
  dens_lvl  <- z_sorted[level_idx]

  contour_list <- contourLines(dens$x, dens$y, dens$z, levels = dens_lvl)
  if (length(contour_list) > 0) {
    cl <- contour_list[[which.max(sapply(contour_list, function(x) length(x$x)))]]
    gate_coords <- cbind(FSC.A = cl$x, SSC.A = cl$y)
  } else {
    # fallback: convex hull
    coords      <- cbind(FSC.A = fsc, SSC.A = ssc)
    gate_coords <- coords[chull(coords), ]
  }
  pg <- polygonGate(filterId = "APC_backgate", .gate = gate_coords)
  backgated_fr <- Subset(apc_pos_fr, pg)
  backgated_list[[name]] <- backgated_fr
  cells_in_polygon <- Subset(fr, pg)
  backgated_list_cells_in_polygon[[name]] <- cells_in_polygon

  png(filename = paste0("backgating_", gsub("[^A-Za-z0-9]", "_", name), ".png"),
      width = 900, height = 900)
  plot(fsc, ssc, col = rgb(0.7, 0.7, 0.7, 0.3), pch = 16, cex = 0.4,
       xlab = "FSC.A", ylab = "SSC.A",
       main = paste("Density Back-Gating:", name))
  points(exprs(backgated_fr)[, "FSC.A"], exprs(backgated_fr)[, "SSC.A"],
         col = rgb(1, 0, 0, 0.5), pch = 16, cex = 0.5)
  polygon(gate_coords[c(seq_len(nrow(gate_coords)), 1), 1],
          gate_coords[c(seq_len(nrow(gate_coords)), 1), 2],
          border = "blue", lwd = 2)
  dev.off()
}

# Optional: convert back to flowSet
fs_backgated_APC_pos <- flowSet(backgated_list)
fs_backgated <- flowSet(backgated_list_cells_in_polygon)

# STEP 4: Compute spillover (compensation) matrix from single-stained controls
# ----------------------------------------------------------------------------
# Calculate the spillover matrix using the specified fluorescence channels.
# The "patt" parameter selects which fluorescent channels to include in the calculation.
spill_list <- spillover(
  fs_controls,
  patt = "APC.A|BV786.A|BV510.A|BB515.A|PE.A",  # regex for selecting fluorescence channels # nolint: line_length_linter.
  unstained = "Unstained",                 # specify your unstained control
  fsc = "FSC.A",                           # forward scatter channel (exclude from compensation) # nolint
  ssc = "SSC.A",                           # side scatter channel (exclude from compensation)
  method = "median"                        # use median intensity for spillover calculation
)


spill_mat <- spill_list

# Verify matrix
print(spill_mat)

#Check spill_mat class. Should be a Matrix / Matrix array
class(spill_mat)


# STEP 6: Inspect loaded control and sample names (optional, for verification)
# ----------------------------------------------------------------------------
sampleNames(fs_controls)  # Check assigned control names
fs                        # Quickly view your experimental samples in the flowSet

# STEP 7: Apply compensation matrix to your experimental samples

#Inverse the transformation to apply the spillover matrix
inv_logicle <- inverseLogicleTransform(logicle)

trans_list_inv <- transformList(channels, inv_logicle)
fs_backgated_untransformed_noncomp <- transform(fs_backgated, trans_list_inv)



# --------------------------------------------------------------
# Compensate all samples in your experimental flowSet (fs) using the spillover matrix


fs_backgated_untransformed_comp <- compensate(fs, spill_mat)

# Your data (fs_comp) is now compensated and ready for transformation and gating.



# STEP 8: Transform
# ------------------------------------------------------------------------
# Transform the data using the `transform` function from the flowCore package
# Define the transformation function
channels <- c("APC.A", "BV786.A", "BV510.A", "BB515.A", "PE.A")


logicle <- logicleTransform(w = 0.5, t = 262143, m = 4.5, a = 0)
trans_list <- transformList(channels, logicle)
fs_backgated_comp_trans <- transform(fs_comp, trans_list)



# Visualize FSC-A vs SSC-A density
autoplot(fs_comp_trans[[5]], x = "FSC.A", y = "SSC.A", bins = 128)
autoplot(fs[[5]], x = "FSC.A", y = "SSC.A", bins = 128)
autoplot(fs_backgated[[5]], x = "FSC.A", y = "SSC.A", bins = 128)
autoplot(fs_backgated_APC_pos[[5]], x = "FSC.A", y = "SSC.A", bins = 128)
autoplot(fs_backgated_comp_trans[[5]], x = "PE.A", y = "SSC.A", bins = 128)




#Clustering and visualization using CATALYST
# Define the marker panel used for clustering
panel <- data.frame(
  channel      = c("APC.A", "BV786.A", "BV510.A", "BB515.A", "PE.A"),
  antigen      = c("CD4", "CD2", "CD3", "CD8", "CD34"),
  marker_class = rep("type", 5),  # or "state"/"none" as appropriate
  stringsAsFactors = FALSE
)

# Metadata describing each sample
md <- data.frame(
  file_name = sampleNames(fs_backgated_comp_trans),
  sample_id = sampleNames(fs_backgated_comp_trans),
  condition = "sample",
  stringsAsFactors = FALSE
)



#------------------------------------- Clump together experimental samples for CATALYST analysis-------
# -----------------Run either this or the section below to run all samples separately. Both are not needed.
# Convert the flowSet to a SingleCellExperiment
sce <- CATALYST::prepData(
  fs_backgated_comp_trans,
  panel = panel,
  md = md,
  panel_cols = list(channel = "channel", antigen = "antigen", class = "marker_class"),
  md_cols = list(file = "file_name", id = "sample_id", factors = "condition")
)

# Perform FlowSOM clustering
sce <- CATALYST::cluster(
  sce,
  features = panel$antigen,
  xdim = 10, ydim = 10, maxK = 7
)

# Explore cluster expression to annotate cell types
plotClusterExprs(sce, k = "meta7")
# Alternatively:
plotMultiHeatmap(sce, k = "meta20")
# Use known marker expression patterns to interpret clusters as specific cell types.
#------------------------------^^^^^^^^Clump together experimental samples for CATALYST analysis-------

#--------------------------------------------Run all samples seperately

# Convert and cluster each sample independently
sce_list <- lapply(seq_along(fs_backgated_comp_trans), function(i) {
  fr <- fs_backgated_comp_trans[[i]]
  sample_name <- sampleNames(fs_backgated_comp_trans)[i]

  fs_single <- flowSet(list(fr))      # create one-sample flowSet
  sampleNames(fs_single) <- sample_name

  md_i <- data.frame(
    file_name = sample_name,
    sample_id = sample_name,
    condition = "sample",
    stringsAsFactors = FALSE
  )

  sce_i <- CATALYST::prepData(
    fs_single,                         # use the one-sample flowSet
    panel = panel,
    md = md_i,
    panel_cols = list(channel = "channel", antigen = "antigen", class = "marker_class"),
    md_cols = list(file = "file_name", id = "sample_id", factors = "condition")
  )

  sce_i <- CATALYST::cluster(
    sce_i,
    features = panel$antigen,
    xdim = 10, ydim = 10, maxK = 10
  )
  sce_i
})

# sce_list now contains a clustered SingleCellExperiment for each sample
names(sce_list) <- sampleNames(fs_backgated_comp_trans)

#--------------------------------------------^^^^Run all samples seperately


plotClusterExprs(sce_list[[1]], k = "meta10")


#-------------------------------------- Classify marker expression per cluster ----
# Generate a summary table of marker levels for each independently clustered
# sample. Levels are assigned based on the median expression of each marker per
# cluster, where values below the 33rd percentile are "Negative", above the 67th
# percentile are "Positive" and the rest are "Mid". In addition, report the
# number and percentage of events belonging to each cluster within its sample.
cluster_levels <- lapply(seq_along(sce_list), function(i) {
  sce_i <- sce_list[[i]]
  sample_name <- names(sce_list)[i]

  cl_ids <- CATALYST::cluster_ids(sce_i, k = "meta10")
  expr_mat <- assay(sce_i, "exprs")[panel$antigen, , drop = FALSE]
  uniq_cl <- sort(unique(cl_ids))

  medians <- sapply(uniq_cl, function(cl) {
    apply(expr_mat[, cl_ids == cl, drop = FALSE], 1, median, na.rm = TRUE)
  })
  medians <- t(medians)
  colnames(medians) <- panel$antigen
  rownames(medians) <- paste0("Cluster", uniq_cl)

  low  <- apply(medians, 2, quantile, probs = 0.33)
  high <- apply(medians, 2, quantile, probs = 0.67)

  classify <- function(value, marker) {
    if (value < low[marker]) {
      "Negative"
    } else if (value > high[marker]) {
      "Positive"
    } else {
      "Mid"
    }
  }

  labels <- t(apply(medians, 1, function(row) {
    sapply(seq_along(row), function(j) classify(row[j], names(row)[j]))
  }))
  colnames(labels) <- panel$antigen

  n_events <- sapply(uniq_cl, function(cl) sum(cl_ids == cl))
  pct_events <- n_events / length(cl_ids) * 100

  data.frame(
    Sample = sample_name,
    Cluster = rownames(medians),
    Events = n_events,
    Percent = pct_events,
    labels,
    stringsAsFactors = FALSE
  )
})

cluster_levels_df <- do.call(rbind, cluster_levels)
print(cluster_levels_df)
writexl::write_xlsx(cluster_levels_df, "cluster_summary.xlsx")

