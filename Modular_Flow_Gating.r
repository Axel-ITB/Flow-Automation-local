# Modular Flow Gating Script
#setwd("C:/Users/AxelBergenstråle/OneDrive - ITB-MED AB/Desktop/R Projects/FACS Automation local")

# This script provides a set of small, composable gating methods that can be
# chained together. Each method isolates one peak from a density distribution by
# cutting off the surrounding valleys. The methods are numbered so that users can
# specify a sequence like c(1, 3, 5) to apply several gating strategies in
# order.

# Method 1: isolate the dominant peak (works with any number of peaks)
# Method 2: two peaks, keep the left one
# Method 3: two peaks, keep the right one
# Method 4: three peaks, keep the left one
# Method 5: three peaks, keep the middle one
# Method 6: three peaks, keep the right one

# ---- Utility functions ------------------------------------------------------

#' Find peaks and valleys in a numeric vector using a kernel density estimate.
#'
#' @param vals Numeric vector.
#' @param bw   Optional bandwidth passed to `density`.
#' @return A list with elements `x`, `y`, `peaks`, and `valleys` (indices in
#'         the density object).
find_peaks_valleys <- function(vals, bw = NULL) {
  dens <- if (!is.null(bw)) density(vals, bw = bw) else density(vals)
  y <- dens$y
  x <- dens$x
  dy <- diff(y)
  sign_changes <- diff(sign(dy))
  list(
    x = x,
    y = y,
    peaks = which(sign_changes == -2) + 1,
    valleys = which(sign_changes ==  2) + 1
  )
}

#' Gate around a specific peak.
#'
#' @param vals         Numeric vector.
#' @param total_peaks  Expected number of peaks in the distribution.
#' @param peak_number  Which peak to select (1 = leftmost, 2 = middle, 3 = rightmost).
#' @param plot         If TRUE, plot the density and thresholds.
#' @param main         Plot title when `plot = TRUE`.
#' @return A list with `left` and `right` thresholds.
 gate_specific_peak <- function(vals, total_peaks, peak_number, plot = FALSE, main = "") {
  pv <- find_peaks_valleys(vals)
  x <- pv$x
  peaks <- pv$peaks
  valleys <- pv$valleys

  if (length(peaks) < total_peaks) {
    stop(sprintf("Expected %d peaks but found %d", total_peaks, length(peaks)))
  }

  # Order peaks by position and select the requested peak
  peaks_sorted <- peaks[order(x[peaks])]
  peak_idx <- peaks_sorted[peak_number]

  # Find nearest valleys around the peak
  left_candidates  <- valleys[valleys < peak_idx]
  right_candidates <- valleys[valleys > peak_idx]
  left_valley  <- if (length(left_candidates))  max(left_candidates)  else NA
  right_valley <- if (length(right_candidates)) min(right_candidates) else NA

  left_threshold  <- if (!is.na(left_valley))  x[left_valley]  else min(vals, na.rm = TRUE)
  right_threshold <- if (!is.na(right_valley)) x[right_valley] else max(vals, na.rm = TRUE)

  if (plot) {
    plot(x, pv$y, type = "l", main = main)
    abline(v = x[peak_idx], col = "blue", lty = 2)
    abline(v = c(left_threshold, right_threshold), col = "red", lty = 2)
  }

  list(left = left_threshold, right = right_threshold)
}

#' Gate around the dominant peak regardless of how many peaks exist.
#'
#' @param vals Numeric vector.
#' @param plot If TRUE, plot the density.
#' @param main Plot title when `plot = TRUE`.
#' @return A list with `left` and `right` thresholds.
 gate_main_peak <- function(vals, plot = FALSE, main = "") {
  pv <- find_peaks_valleys(vals)
  x <- pv$x
  y <- pv$y
  peaks <- pv$peaks
  valleys <- pv$valleys

  if (!length(peaks)) {
    rng <- range(vals, na.rm = TRUE)
    return(list(left = rng[1], right = rng[2]))
  }

  main_peak <- peaks[which.max(y[peaks])]
  left_candidates  <- valleys[valleys < main_peak]
  right_candidates <- valleys[valleys > main_peak]
  left_valley  <- if (length(left_candidates))  max(left_candidates)  else NA
  right_valley <- if (length(right_candidates)) min(right_candidates) else NA

  left_threshold  <- if (!is.na(left_valley))  x[left_valley]  else min(vals, na.rm = TRUE)
  right_threshold <- if (!is.na(right_valley)) x[right_valley] else max(vals, na.rm = TRUE)

  if (plot) {
    plot(x, y, type = "l", main = main)
    abline(v = x[main_peak], col = "blue", lty = 2)
    abline(v = c(left_threshold, right_threshold), col = "red", lty = 2)
  }

  list(left = left_threshold, right = right_threshold)
}

# ---- Individual gating methods ---------------------------------------------

# Method 1: isolate the dominant peak (works with any number of peaks)
 gate_method1 <- function(vals, plot = FALSE) {
  gate_main_peak(vals, plot = plot, main = "Method 1: dominant peak")
}

# Method 2: two peaks, keep the left one
 gate_method2 <- function(vals, plot = FALSE) {
  gate_specific_peak(vals, total_peaks = 2, peak_number = 1, plot = plot,
                     main = "Method 2: left of two peaks")
}

# Method 3: two peaks, keep the right one
 gate_method3 <- function(vals, plot = FALSE) {
  gate_specific_peak(vals, total_peaks = 2, peak_number = 2, plot = plot,
                     main = "Method 3: right of two peaks")
}

# Method 4: three peaks, keep the left one
 gate_method4 <- function(vals, plot = FALSE) {
  gate_specific_peak(vals, total_peaks = 3, peak_number = 1, plot = plot,
                     main = "Method 4: left of three peaks")
}

# Method 5: three peaks, keep the middle one
 gate_method5 <- function(vals, plot = FALSE) {
  gate_specific_peak(vals, total_peaks = 3, peak_number = 2, plot = plot,
                     main = "Method 5: middle of three peaks")
}

# Method 6: three peaks, keep the right one
 gate_method6 <- function(vals, plot = FALSE) {
  gate_specific_peak(vals, total_peaks = 3, peak_number = 3, plot = plot,
                     main = "Method 6: right of three peaks")
}

# A lookup table for the numbered methods
 gating_methods <- list(
  `1` = gate_method1,
  `2` = gate_method2,
  `3` = gate_method3,
  `4` = gate_method4,
  `5` = gate_method5,
  `6` = gate_method6
)

#' Apply a sequence of gating methods to a numeric vector.
#'
#' @param vals    Numeric vector.
#' @param methods Integer vector with method numbers.
#' @param plot    If TRUE, plotting is enabled for every step.
#' @return Filtered numeric vector after all gating steps.
 run_gating_sequence <- function(vals, methods, plot = FALSE) {
  current_vals <- vals
  for (m in methods) {
    fn <- gating_methods[[as.character(m)]]
    if (is.null(fn)) stop(sprintf("Unknown method %s", m))
    th <- fn(current_vals, plot = plot)
    current_vals <- current_vals[current_vals >= th$left & current_vals <= th$right]
  }
  current_vals
}

# ---- Example ----------------------------------------------------------------
# The example is executed only when this file is run directly, not when sourced.
if (sys.nframe() == 0) {
  set.seed(1)
  # Three-peak synthetic distribution
  vals <- c(rnorm(1000, 1), rnorm(1000, 3), rnorm(1000, 5))
  # Apply Method 1 followed by Method 3 (dominant peak, then right of two peaks)
  gated <- run_gating_sequence(vals, c(1, 3))
  cat("Remaining events:", length(gated), "\n")
}

