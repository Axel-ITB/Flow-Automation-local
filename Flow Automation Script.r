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

setwd("C:/Users/AxelBergenstråle/OneDrive - ITB-MED AB/Desktop/R Projects/FACS Automation local")

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
  level_idx <- which(cum_z >= 0.95 * sum(z_sorted))[1]   # contour enclosing ~95%
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
fs_backgated <- flowSet(backgated_list)

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
# --------------------------------------------------------------
# Compensate all samples in your experimental flowSet (fs) using the spillover matrix
fs_comp <- compensate(fs, spill_mat)

# Your data (fs_comp) is now compensated and ready for transformation and gating.



# STEP 8: Transform
# ------------------------------------------------------------------------
# Transform the data using the `transform` function from the flowCore package
# Define the transformation function
channels <- c("APC.A", "BV786.A", "BV510.A", "BB515.A", "PE.A")


logicle <- logicleTransform(w = 0.5, t = 262143, m = 4.5, a = 0)
trans_list <- transformList(channels, logicle)
fs_comp_trans <- transform(fs_comp, trans_list)

# # Estimate transformation parameters from one sample
# trans_auto <- estimateLogicle(fs_comp[[1]], channels)
# # Apply to all samples in flowSet
# fs_comp_trans <- transform(fs_comp, trans_auto)

# STEP 8: Visualize the transformed data
# Load ggcyto for visualization


# Visualize FSC-A vs SSC-A density
autoplot(fs_comp_trans[[5]], x = "FSC.A", y = "SSC.A", bins = 128)
autoplot(fs[[5]], x = "FSC.A", y = "SSC.A", bins = 128)




# ---- CATALYST workflow from here ----
# Assumes fs_comp_trans is your compensated and transformed flowSet

# Load CATALYST and SingleCellExperiment
if (!requireNamespace("CATALYST", quietly = TRUE)) BiocManager::install("CATALYST")
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment")
library(CATALYST)
library(SingleCellExperiment)


#Converting to CATALYST, with md.
#-----------------------------------------------------------------------------------------------

# This entire comment block is for converting to catalyst format using already transformed data. but that fucks with FSC and SSC, so im trying to run it on untransformed data. 
# Convert flowSet to a CATALYST-compatible SingleCellExperiment (SCE) object
# You need to specify the marker/channel names for panel and md
# panel <- data.frame(
#   fcs_colname = colnames(exprs(fs_comp_trans[[1]])),
#   antigen = colnames(exprs(fs_comp_trans[[1]])),
#   stringsAsFactors = FALSE
# )
# features <- panel$fcs_colname


# md <- data.frame(
#   file_name = as.character(fs_comp_trans@phenoData@data$name),
#   sample_id = as.character(fs_comp_trans@phenoData@data$name),
#   condition = rep("sample", length(fs_comp_trans)),
#   stringsAsFactors = FALSE,
#   check.names = FALSE
# )

# #Debugging
# str(md)
# print(head(md))


#Converting to CATALYST, without md.
#-----------------------------------------------------------------------------------------------
# # Running the code with md=md does not work. Don't know why. running it without md works.

panel <- data.frame(
  fcs_colname = colnames(exprs(fs_comp_trans[[1]])),
  antigen = colnames(exprs(fs_comp_trans[[1]])),
  stringsAsFactors = FALSE
)
features <- panel$fcs_colname

sce <- prepData(
  fs_comp_trans,
  panel = panel,
  features = features,
  cofactor = 1,
  transform = FALSE
)
assay_data <- assay(sce, "counts")
summary(assay_data["FSC.A", ])
summary(assay_data["SSC.A", ])
plotScatter(sce, chs = c("FSC.A", "SSC.A"), assay = "counts")

#Gives errors chatgpt thinks its some kind of column naming missmatch. 
# Use compensated, but not transformed, data
# panel2 <- data.frame(
#   fcs_colname = colnames(exprs(fs_comp[[1]])),
#   antigen = colnames(exprs(fs_comp[[1]])),
#   stringsAsFactors = FALSE
# )
# panel2$marker_class <- ifelse(grepl("FSC|SSC", panel2$fcs_colname), "none", "type")


# features <- panel2$fcs_colname
# sce2 <- prepData(
#   fs_comp,
#   panel = panel2,
#   features = features,
#   cofactor = 5 # or your preferred cofactor for arcsinh
# )
# plotScatter(sce2, chs = c("FSC.A", "SSC.A"), assay = "counts")

# #debugging attempt
# print(colnames(exprs(fs_comp[[1]])))
# print(panel2$fcs_colname)

# #Attempt to fix
# # Get the actual column names and their descriptions
# col_names <- colnames(fs_comp[[1]])
# desc_names <- pData(parameters(fs_comp[[1]]))$desc

# # If desc is NA, fall back to col_names
# desc_names[is.na(desc_names)] <- col_names[is.na(desc_names)]

# panel2 <- data.frame(
#   fcs_colname = col_names,
#   antigen = desc_names,
#   marker_class = ifelse(grepl("FSC|SSC", desc_names), "none", "type"),
#   stringsAsFactors = FALSE
# )

# features <- panel2$fcs_colname

# sce2 <- prepData(
#   fs_comp,
#   panel = panel2,
#   features = features,
#   cofactor = 5
# )
# plotScatter(sce2, chs = c("FSC.A", "SSC.A"), assay = "counts")
# plotScatter(sce2, chs = c(panel2$fcs_colname[panel2$antigen == "FSC.A"],
#                           panel2$fcs_colname[panel2$antigen == "SSC.A"]),
#             assay = "counts")


# # Minimal panel with only scatter
# panel_test <- data.frame(
#   fcs_colname = colnames(fs_comp[[1]])[1:4], # e.g., $P1N, $P2N, $P3N, $P4N
#   antigen = pData(parameters(fs_comp[[1]]))$desc[1:4],
#   marker_class = "none",
#   stringsAsFactors = FALSE
# )
# panel_test$antigen[is.na(panel_test$antigen)] <- panel_test$fcs_colname[is.na(panel_test$antigen)]

# sce_test <- prepData(
#   fs_comp,
#   panel = panel_test,
#   features = panel_test$fcs_colname,
#   cofactor = 5
# )
# plotScatter(sce_test, chs = panel_test$antigen[1:2], assay = "counts")




#End of debugging attempt




# Quick QC plot: scatter of FSC-A vs SSC-A

CATALYST::plotScatter(sce, c("FSC.A", "SSC.A"), assay = "counts")
plotScatter(sce, chs = c("FSC.A", "SSC.A"), assay = "counts")



# Example: clustering with FlowSOM (adjust markers as needed)
# Define which channels to use for clustering (e.g., all fluorescence channels)
cluster_channels <- c("APC.A", "BV786.A", "BV510.A", "BB515.A", "PE.A")
# Add this before running CATALYST::cluster
assay(sce, "exprs") <- assay(sce, "counts")
sce <- CATALYST::cluster(sce, features = cluster_channels, xdim = 10, ydim = 10, maxK = 10, seed = 123)
# After running CATALYST::cluster()
table(sce$cluster_id)
# Cluster counts per sample
table(sce$sample_id, sce$cluster_id)


# ------------- View Meta-clusters -------------
CATALYST::plotExprHeatmap(
  sce,
  features = cluster_channels,
  by = "cluster_id",
  scale = "last",
  fun = "median",
  k = "meta10"
)


# ------------- View Meta-clusters for specific sample -------------
# Subset SCE to a specific sample (e.g., sample 5)
sce_sample_to_plot <- sce[, sce$sample_id == unique(sce$sample_id)[5]]

CATALYST::plotExprHeatmap(
  sce_sample_to_plot,
  features = cluster_channels,
  by = "cluster_id",
  scale = "last",
  fun = "median",
  k = "meta10"
)



  # "APC.A"    = "APC.A CD4",
  # "BV786.A"  = "BV786.A CD2",
  # "BV510.A"  = "BV510.A CD3",
  # "BB515.A"  = "BB515.A CD8",
  # "PE.A"     = "PE.A CD34"



# Show the number of cells in each meta-cluster
table(sce$cluster_id)
# Show the number of cells per sample and meta-cluster
clust_to_meta <- metadata(sce)$cluster_codes$meta10
head(clust_to_meta)


# 1. Pull cluster-to-meta mapping (meta10 assumed here)
meta_map <- metadata(sce)$cluster_codes$meta10

# 2. Map it back to each cell using cluster_id
sce$meta10_cell <- meta_map[as.character(sce$cluster_id)]
table(sample_id = sce$sample_id, meta_cluster = sce$meta10_cell)





# ------------- Figure out what the different clusters are -------------
# Calculate median expression of each marker per cluster
cluster_medians <- aggregate(
  t(assay(sce, "exprs")[cluster_channels, ]),
  by = list(cluster = sce$cluster_id),
  FUN = median
)
# Install if needed: install.packages("pheatmap")
library(pheatmap)
# Remove the cluster column for heatmap input
mat <- as.matrix(cluster_medians[,-1])
rownames(mat) <- paste("Cluster", cluster_medians$cluster)
# Plot heatmap
pheatmap(mat, cluster_rows = FALSE, cluster_cols = TRUE, 
         main = "Median marker expression per cluster")











# ------------------------ Visualizing clusters ------------------------
# Visualize clusters on a t-SNE/UMAP map
# Need to change the runDR. Takes forever to run.
sce <- CATALYST::runDR(sce, dr = "UMAP", features = cluster_channels)
CATALYST::plotDR(sce, "UMAP", color_by = "cluster_id")

# Downsample to 5000 cells per sample using filterSCE (works in older CATALYST)
sce_sub <- CATALYST::filterSCE(sce, n = 5000, by = "sample_id", seed = 123)
# Run UMAP on the downsampled SCE
sce_sub <- CATALYST::runDR(sce_sub, dr = "UMAP", features = cluster_channels)
CATALYST::plotDR(sce_sub, "UMAP", color_by = "cluster_id")

# ...continue with CATALYST downstream analysis as needed...
#OBS No gating in FSC or SSC performed! Need to add using openCYTO or flowdensity.