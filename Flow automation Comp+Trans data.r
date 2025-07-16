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
library(flowDensity)
#setwd("C:/Users/AxelBergenstråle/OneDrive - ITB-MED AB/Desktop/R Projects/FACS Automation local")


# 1. Load all FCS files from the 'comp+trans data' directory
#    These files are compensated but not transformed
fcs_path <- "comp+trans data"
fs_raw <- read.flowSet(path = fcs_path, pattern = "\\.fcs$", alter.names = TRUE, truncate_max_range = FALSE)


# 1.1 Check channels
pData(parameters(fs_raw[[1]]))





# 3. Plot FSC.A vs SSC.A for the [[]] sample
p <- autoplot(fs_raw[[1]], x = "FSC.A", y = "SSC.A", bins = 128) +
  ggplot2::ggtitle("FSC vs SSC after filtering")
print(p)

# 3.1 Histogram of APC.A
p_hist <- autoplot(fs_raw[[1]], "FSC.A", bins = 128) +
  ggplot2::ggtitle("Histogram of APC.A after filtering")
print(p_hist)

# 3.2 Automatically gate on FSC.A to remove debris
dir.create("Step 3.2 plots", showWarnings = FALSE)
fsc_filtered_list <- lapply(seq_along(fs_raw), function(i) {
  fr <- fs_raw[[i]]
   fsc_vals <- exprs(fr)[, "FSC.A"]

  ## Density-based threshold: find the valley after the first peak
  dens <- density(fsc_vals)
  y <- dens$y
  x <- dens$x
  dy <- diff(y)
  sign_changes <- diff(sign(dy))
  peaks <- which(sign_changes == -2) + 1
  valleys <- which(sign_changes == 2) + 1

  threshold <- deGate(
    fr,
    channel   = "FSC.A",
    use.upper = TRUE,
    upper     = FALSE,
    all.cuts  = FALSE
  )

  if (length(peaks) > 0) {
    first_peak <- peaks[1]
    valley_after <- valleys[valleys > first_peak][1]
    if (!is.na(valley_after)) {
      threshold <- x[valley_after]
    }
  }

 plot_file <- file.path(
    "Step 3.2 plots",
    paste0(sampleNames(fs_raw)[i], "_FSC_A_density.png")
  )
  png(plot_file)
  plot(
    dens$x,
    dens$y,
    type = "l",
    main = paste("FSC.A Density:", sampleNames(fs_raw)[i])
  )
  abline(v = threshold, col = "red", lty = 2)
  dev.off()

  fr[fsc_vals > threshold, ]
})
fs_filtered <- flowSet(fsc_filtered_list)

# 3.3 Check
p_hist <- autoplot(fs_filtered[[1]], "FSC.A", bins = 128) +
  ggplot2::ggtitle("Histogram of APC.A after filtering")
print(p_hist)


# 3.4 Plotting x = "SSC.W", y = "FSC.H"
p <- autoplot(fs_filtered[[3]], x = "FSC.W", y = "FSC.H", bins = 1000) +
  ggplot2::ggtitle("FSC vs SSC after filtering")
print(p)

p_hist <- autoplot(fs_filtered[[3]], "FSC.W", bins = 1000) +
  ggplot2::ggtitle("Histogram of FSC.A after filtering")
print(p_hist)

# 3.5 Gating fsc.W and FSC.H
# Helper function ---------------------------------------------------------------
# Automatically gate around the dominant density peak of a numeric vector.
# Returns the nearest valley thresholds on both sides. Set `plot = TRUE` to
dir.create("Step 3.5 plots", showWarnings = FALSE)

gate_main_peak <- function(vals, plot = FALSE, main = "gate_main_peak",
                           valley_max_y = Inf, bw = NULL) {
  if (length(vals) < 2 || all(is.na(vals))) {
    rng <- range(vals, na.rm = TRUE)
    return(list(left = rng[1], right = rng[2]))
  }
  dens <- if (!is.null(bw)) density(vals, bw = bw) else density(vals)

  y <- dens$y
  x <- dens$x
  dy <- diff(y)
  sc <- diff(sign(dy))
  peaks <- which(sc == -2) + 1
  valleys <- which(sc == 2) + 1
   if (!is.infinite(valley_max_y)) {
    valleys <- valleys[y[valleys] <= valley_max_y]
  }

  if (length(peaks) == 0) {
    rng <- range(vals, na.rm = TRUE)
    return(list(left = rng[1], right = rng[2]))
  }
  main_peak <- peaks[which.max(y[peaks])]

  # Add distance and density thresholds to avoid valleys inside the main peak
min_valley_distance <- 0.05 * diff(range(x))  # 5% of the range, adjust as needed
valley_density_threshold <- 0.5 * y[main_peak]  # Only valleys less than 50% of peak height

left_idx <- valleys[
  (valleys < main_peak) &
  (abs(x[valleys] - x[main_peak]) > min_valley_distance) &
  (y[valleys] < valley_density_threshold)
]
right_idx <- valleys[
  (valleys > main_peak) &
  (abs(x[valleys] - x[main_peak]) > min_valley_distance) &
  (y[valleys] < valley_density_threshold)
]

left <- if (length(left_idx) > 0) max(left_idx) else integer(0)
right <- if (length(right_idx) > 0) min(right_idx) else integer(0)

  if (plot) {
    plot(dens$x, dens$y, type = "l", main = main)
    abline(v = x[main_peak], col = "blue", lty = 2)
    if (length(left))  abline(v = x[left],  col = "red", lty = 2)
    if (length(right)) abline(v = x[right], col = "red", lty = 2)
  }

 rng <- range(vals, na.rm = TRUE)
  list(
    left  = if (length(left))  x[left]  else rng[1],
    right = if (length(right)) x[right] else rng[2]
  )
}
#------------End of helper function------------------------------------------------

# Use density-based thresholds around the dominant peak for both FSC.H and fsc.W
# to remove outliers. The gate_main_peak helper finds the main peak and the
# nearest valleys on each side.
gated_list <- lapply(seq_along(fs_filtered), function(i) {
  fr <- fs_filtered[[i]]
  fsc_h <- exprs(fr)[, "FSC.H"]
  fsc_w <- exprs(fr)[, "FSC.W"]
   th_fsc <- gate_main_peak(fsc_h, valley_max_y = 2e6)
  th_fsc <- gate_main_peak(fsc_w, valley_max_y = 2e6)
  keep <- (fsc_h >= th_fsc$left & fsc_h <= th_fsc$right) &
          (fsc_w >= th_fsc$left & fsc_w <= th_fsc$right)

  # Optional diagnostic plots
   png(file.path(
    "Step 3.5 plots",
    paste0(sampleNames(fs_filtered)[i], "_FSC_H_density.png")
  ))
  gate_main_peak(
    fsc_h,
    plot = TRUE,
    main = paste("FSC.H Density:", sampleNames(fs_filtered)[i]),
    valley_max_y = 2e6
  )
  dev.off()
  png(file.path(
    "Step 3.5 plots",
    paste0(sampleNames(fs_filtered)[i], "_FSC_W_density.png")
  ))
  gate_main_peak(
    fsc_w,
    plot = TRUE,
    main = paste("FSC.W Density:", sampleNames(fs_filtered)[i]),
    valley_max_y = 2e6
  )
  dev.off()
  fr[keep, ]
})
fs_filtered_2x <- flowSet(gated_list)

# 3.6 plotting fsc.W and FSC.H after gating

p <- autoplot(fs_filtered_2x[[3]], x = "FSC.W", y = "FSC.H", bins = 1000) +
  ggplot2::ggtitle("FSC vs SSC after filtering")
print(p)

# 3.7 ploitting SSC.W and SSC-H Before gating

p <- autoplot(fs_filtered_2x[[3]], x = "SSC.W", y = "SSC.H", bins = 1000) +
  ggplot2::ggtitle("FSC vs SSC after filtering")
print(p)

p_hist <- autoplot(fs_filtered_2x[[3]], "SSC.W", bins = 1000) +
  ggplot2::ggtitle("Histogram of SSC.W before filtering")
print(p_hist)

# 3.8 Gating SSC.W and SSC.H
gated_ssc <- lapply(seq_along(fs_filtered_2x), function(i) {
  fr <- fs_filtered_2x[[i]]
  ssc_h <- exprs(fr)[, "SSC.H"]
  ssc_w <- exprs(fr)[, "SSC.W"]
  th_h <- gate_main_peak(ssc_h, valley_max_y = 2e6)
  th_w <- gate_main_peak(ssc_w, valley_max_y = 2e6)
  keep <- (ssc_h >= th_h$left & ssc_h <= th_h$right) &
          (ssc_w >= th_w$left & ssc_w <= th_w$right)

  gate_main_peak(
    ssc_h,
    plot = TRUE,
    main = paste("SSC.H Density:", sampleNames(fs_filtered_2x)[i]),
    valley_max_y = 2e6
  )
  gate_main_peak(
    ssc_w,
    plot = TRUE,
    main = paste("SSC.W Density:", sampleNames(fs_filtered_2x)[i]),
    valley_max_y = 2e6
  )
  fr[keep, ]
})
fs_filtered_3x <- flowSet(gated_ssc)

# 3.9 plotting SSC.W and SSC-H after gating
p <- autoplot(fs_filtered_3x[[3]], x = "SSC.W", y = "SSC.H", bins = 1000) +
  ggplot2::ggtitle("FSC vs SSC after filtering")
print(p)

p_hist <- autoplot(fs_filtered_3x[[3]], "SSC.W", bins = 1000) +
  ggplot2::ggtitle("Histogram of SSC.W before filtering")
print(p_hist)




# 3.10 Plot CD45 Alexa.Fluor.532.A and SSC-H. Only include >-10000 Alexa.Fluor.532.A values
fr_pos <- fs_filtered_3x[[3]][exprs(fs_filtered_3x[[3]])[, "Alexa.Fluor.532.A"] > -100000, ]
p <- autoplot(fr_pos, x = "Alexa.Fluor.532.A", y = "SSC.A", bins = 1000) +
  ggplot2::ggtitle("FSC vs SSC after filtering")
print(p)



# 4. Transform fluorescent channels (e.g. FITC.A, APC.A, etc.)
#    Identify fluorescent channels and excluding scatter and time parameters
param_names <- colnames(fs_raw[[3]])
fluor_channels <- grep("\\.A$", param_names, value = TRUE)
fluor_channels <- setdiff(
  fluor_channels,
  c("FSC.A", "FSC.H", "SSC.A", "SSC.H", "SSC.B.A")
)
fluor_channels


# Estimate logicle transformation from the first sample and apply to all
# `estimateLogicle` already returns a transformation list, so we can use it
# directly with `transform()`
logicle <- estimateLogicle(fs_filtered_3x[[3]], fluor_channels)
fs_trans <- transform(fs_filtered_3x, logicle)

# 4. Plot an example fluorescence vs SSC after transformation
example_chan <- fluor_channels[3]
p <- autoplot(fs_trans[[3]], x = example_chan, y = "SSC.A", bins = 128) +
  ggplot2::ggtitle(paste("FSC vs SSC after", example_chan, "transformation"))
print(p)

p_hist <- autoplot(fs_trans[[3]], "Alexa.Fluor.532.A", bins = 128) +
  ggplot2::ggtitle("Histogram of Alexa.Fluor.532.A before gating")
print(p_hist)


p <- autoplot(fs_trans[[2]], x = "Alexa.Fluor.532.A", y = "SSC.h", bins = 1000) +
  ggplot2::ggtitle("FSC vs SSC before gating")
print(p)


# 4.1 Gate Alexa.Fluor.532.A to isolate positive population
#    Identify valley between two peaks on transformed data
#    Keep events above the threshold

# Directory for diagnostic plots
dir.create("Step 4.1 Alexa gating", showWarnings = FALSE)

alexa_pos_list <- lapply(seq_along(fs_trans), function(i) {
  fr <- fs_trans[[i]]
  alexa_vals <- exprs(fr)[, "Alexa.Fluor.532.A"]

  dens <- density(alexa_vals)
  y <- dens$y
  x <- dens$x
  dy <- diff(y)
  sign_changes <- diff(sign(dy))
  peaks <- which(sign_changes == -2) + 1
  valleys <- which(sign_changes == 2) + 1

  threshold <- deGate(
    fr,
    channel   = "Alexa.Fluor.532.A",
    use.upper = TRUE,
    upper     = FALSE,
    all.cuts  = FALSE
  )

    if (length(peaks) >= 2) {
    # choose the two most prominent peaks and cut at the lowest valley
    top_two <- order(y[peaks], decreasing = TRUE)[1:2]
    sel_peaks <- sort(peaks[top_two])
    valley_between <- valleys[valleys > sel_peaks[1] & valleys < sel_peaks[2]]
    if (length(valley_between) > 0) {
      valley_idx <- valley_between[which.min(y[valley_between])]
      threshold <- x[valley_idx]
    }
  } else if (length(peaks) == 1) {
     # fall back to valley after the single detected peak
    valley_after <- valleys[valleys > peaks[1]]
    if (length(valley_after) > 0) {
      valley_idx <- valley_after[which.min(y[valley_after])]
      threshold <- x[valley_idx]
    }
  }

  plot_file <- file.path(
    "Step 4.1 Alexa gating",
    paste0(sampleNames(fs_trans)[i], "_AlexaFluor532_density.png")
  )
  png(plot_file)
  plot(
    dens$x,
    dens$y,
    type = "l",
    main = paste("Alexa.Fluor.532.A Density:", sampleNames(fs_trans)[i])
  )
  abline(v = threshold, col = "red", lty = 2)
  dev.off()

  fr[alexa_vals > threshold, ]
})
fs_alexa_pos <- flowSet(alexa_pos_list)

# 4.2 Check the Alexa.Fluor.532.A histogram after gating
p_hist <- autoplot(fs_alexa_pos[[3]], "Alexa.Fluor.532.A", bins = 128) +
  ggplot2::ggtitle("Histogram of Alexa.Fluor.532.A after filtering")
print(p_hist)


p <- autoplot(fs_alexa_pos[[2]], x = "Alexa.Fluor.532.A", y = "SSC.h", bins = 100) +
  ggplot2::ggtitle("FSC vs SSC after filtering")
print(p)


# 5.1 check alexa fluor 700 before gating
p_hist <- autoplot(fs_alexa_pos[[3]], "Alexa.Fluor.700.A", bins = 128) +
  ggplot2::ggtitle("Histogram of Alexa.Fluor.700.A before gating")  
print(p_hist)

p <- autoplot(fs_alexa_pos[[3]], x = "Alexa.Fluor.700.A", y = "SSC.a", bins = 100) +
  ggplot2::ggtitle("FSC vs SSC after filtering")
print(p)

# 5.2 Gate Alexa.Fluor.700.A 
dir.create("Step 5.2 Alexa700 gating", showWarnings = FALSE)

alexa700_low_list <- lapply(seq_along(fs_alexa_pos), function(i) {
  fr <- fs_alexa_pos[[i]]
  alexa_vals <- exprs(fr)[, "Alexa.Fluor.700.A"]

  dens <- density(alexa_vals)
  y <- dens$y
  x <- dens$x
  dy <- diff(y)
  sign_changes <- diff(sign(dy))
  peaks <- which(sign_changes == -2) + 1
  valleys <- which(sign_changes == 2) + 1

  threshold <- deGate(
    fr,
    channel   = "Alexa.Fluor.700.A",
    use.upper = FALSE,
    upper     = TRUE,
    all.cuts  = FALSE
  )

  if (length(peaks) >= 2) {
    top_two <- order(y[peaks], decreasing = TRUE)[1:2]
    sel_peaks <- sort(peaks[top_two])
    valley_between <- valleys[valleys > sel_peaks[1] & valleys < sel_peaks[2]]
    if (length(valley_between) > 0) {
      valley_idx <- valley_between[which.min(y[valley_between])]
      threshold <- x[valley_idx]
    }
  } else if (length(peaks) == 1) {
    valley_before <- valleys[valleys < peaks[1]]
    if (length(valley_before) > 0) {
      valley_idx <- valley_before[which.min(y[valley_before])]
      threshold <- x[valley_idx]
    }
  }

  plot_file <- file.path(
    "Step 5.2 Alexa700 gating",
    paste0(sampleNames(fs_alexa_pos)[i], "_AlexaFluor700_density.png")
  )
  png(plot_file)
  plot(
    dens$x,
    dens$y,
    type = "l",
    main = paste("Alexa.Fluor.700.A Density:", sampleNames(fs_alexa_pos)[i])
  )
  abline(v = threshold, col = "red", lty = 2)
  dev.off()

  fr[alexa_vals < threshold, ]
})
fs_alexa700_low <- flowSet(alexa700_low_list)

# 5.3 check alexa fluor 700 after gating
p_hist <- autoplot(fs_alexa700_low[[3]], "Alexa.Fluor.700.A", bins = 128) +
  ggplot2::ggtitle("Histogram of Alexa.Fluor.700.A after gating")  
print(p_hist)

p <- autoplot(fs_alexa700_low[[3]], x = "Alexa.Fluor.700.A", y = "SSC.a", bins = 100) +
  ggplot2::ggtitle("FSC vs SSC after gating")
print(p)

# 6.1 check AlexaFluor 532 before gating. Figure 6. 
p_hist <- autoplot(fs_alexa700_low[[3]], "Alexa.Fluor.532.A", bins = 128) +
  ggplot2::ggtitle("Histogram of Alexa.Fluor.700.A after gating")  
print(p_hist)

p <- autoplot(fs_alexa700_low[[3]], x = "Alexa.Fluor.532.A", y = "SSC.a", bins = 100) +
  ggplot2::ggtitle("FSC vs SSC after gating")
print(p)



# 6.2 Gate Lymphs from the Lymph Mono population
#    Using CD45 (Alexa.Fluor.532.A) versus SSC.A. We approximate the
#    polygon gate in the reference figure with 1D thresholds on both
#    parameters.
dir.create("Step 6 Lymph gating", showWarnings = FALSE)

lymph_list <- lapply(seq_along(fs_alexa700_low), function(i) {
  fr <- fs_alexa700_low[[i]]
  cd45_vals <- exprs(fr)[, "Alexa.Fluor.532.A"]
  ssc_vals  <- exprs(fr)[, "SSC.A"]

  th_ssc  <- gate_main_peak(ssc_vals)
  th_cd45 <- gate_main_peak(cd45_vals)

  plot_file <- file.path(
    "Step 6 Lymph gating",
    paste0(sampleNames(fs_alexa700_low)[i], "_CD45_vs_SSC_A.png")
  )
  png(plot_file)
  smoothScatter(
    cd45_vals,
    ssc_vals,
    xlab = "Alexa.Fluor.532.A",
    ylab = "SSC.A",
    main = paste("CD45 vs SSC.A:", sampleNames(fs_alexa700_low)[i])
  )
  abline(h = c(th_ssc$left, th_ssc$right), col = "red", lty = 2)
  abline(v = c(th_cd45$left, th_cd45$right), col = "red", lty = 2)
  dev.off()

  keep <- (ssc_vals >= th_ssc$left & ssc_vals <= th_ssc$right) &
          (cd45_vals >= th_cd45$left & cd45_vals <= th_cd45$right)
  fr[keep, ]
})
fs_lymphs <- flowSet(lymph_list)



