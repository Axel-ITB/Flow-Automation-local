#!/usr/bin/env Rscript

# Lightweight CLI to run modular gating on .fcs files in a folder.
# - Sources Modular_Flow_Gating.r for the gating methods (1..6)
# - Reads all FCS files in a directory
# - Applies a sequence of methods to one or more channels (1D density gates)
# - Supports same methods for all channels, or per-channel rules
# - Prints per-file counts; optionally writes gated FCS + per-step density plots

suppressWarnings(suppressMessages({
  # flowCore is required to read/write FCS
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required. Install via Bioconductor: BiocManager::install('flowCore')")
  }
  library(flowCore)
}))

# Source the modular gating functions
source("Modular_Flow_Gating.r")

cat_usage <- function() {
  cat(
    "\nUsage:\n",
    "  # Single channel (methods applied to that channel)\n",
    "  Rscript Run_Modular_Flow_Gating_FCS.r dir=\"comp+trans data\" channel=\"FSC.A\" methods=1,3 [plots=false] [out=\"Gated_Output\"]\n\n",
    "  # Multiple channels (same methods for all)\n",
    "  Rscript Run_Modular_Flow_Gating_FCS.r dir=\"comp+trans data\" channels=\"FSC.A,SSC.A\" methods=1,3 plots=true out=\"Gated_Output\"\n\n",
    "  # Per-channel methods (rules syntax: <channel>:<m1,m2; channel2:m3...>)\n",
    "  Rscript Run_Modular_Flow_Gating_FCS.r dir=\"comp+trans data\" rules=\"FSC.A:1,3;SSC.A:2,6\" plots=true out=\"Gated_Output\"\n\n",
    "  # Ordered pipeline across channels (steps run in sequence)\n",
    "  Rscript Run_Modular_Flow_Gating_FCS.r dir=\"comp+trans data\" pipeline=\"FSC.A:1,3; SSC.A:2; FSC.A:5\" plots=true out=\"Gated_Output\"\n\n",
    "  # From config file (one step per line)\n",
    "  Rscript Run_Modular_Flow_Gating_FCS.r config=\"gating_plan.txt\"\n\n",
    "Arguments:\n",
    "  dir:       Folder containing .fcs files\n",
    "  channel:   Single channel name to gate (e.g., 'FSC.A')\n",
    "  channels:  Comma-separated channels (e.g., 'FSC.A,SSC.A')\n",
    "  methods:   Comma-separated method numbers (applied to channel(s))\n",
    "  rules:     Per-channel methods, e.g. 'FSC.A:1,3;SSC.A:2' (overrides channels/methods)\n",
    "  pipeline:  Ordered steps across channels, e.g. 'FSC.A:1,3;SSC.A:2' (aliases: sequence, steps, plan)\n",
    "  config:    Path to a text file with ordered steps and options\n",
    "  plots:     true/false to save step plots (default: false)\n",
    "  out:       Optional output folder to write gated FCS + summary + plots\n",
    "\nConfig file format (order matters):\n",
    "  # Comments and blank lines are ignored\n",
    "  dir = comp+trans data          # optional; can also be passed via CLI\n",
    "  out = Gated_Output             # optional\n",
    "  plots = true                   # optional\n",
    "  \n",
    "  FSC.A: 1,3                     # channel: methods (comma-separated)\n",
    "  SSC.A: 2                       # next step(s)\n",
    "  FSC.A: 5                       # next step(s)\n",
    sep = ""
  )
}

parse_args <- function(args) {
  # Simple key=value parsing; booleans recognized for plot
  kv <- strsplit(args, "=", fixed = TRUE)
  lst <- list()
  for (p in kv) {
    if (length(p) == 2) {
      key <- trimws(p[[1]])
      val <- trimws(p[[2]])
      # strip optional quotes
      val <- sub('^"', '', sub('"$', '', val))
      val <- sub("^'", '', sub("'$", '', val))
      lst[[key]] <- val
    }
  }
  # Do not require dir here; a config file may provide it. We'll validate later in main.
  # Normalize flags
  lst$plots <- if (!is.null(lst$plots)) tolower(lst$plots) %in% c("true", "1", "t", "yes", "y") else FALSE
  # Normalize channels/methods/rules
  if (!is.null(lst$methods)) lst$methods <- as.integer(strsplit(lst$methods, ",")[[1]])
  if (!is.null(lst$channels)) lst$channels <- trimws(unlist(strsplit(lst$channels, ",")))
  # Allow pipeline aliases
  if (!is.null(lst$sequence) && is.null(lst$pipeline)) lst$pipeline <- lst$sequence
  if (!is.null(lst$steps)    && is.null(lst$pipeline)) lst$pipeline <- lst$steps
  if (!is.null(lst$plan)     && is.null(lst$pipeline)) lst$pipeline <- lst$plan
  lst
}

plot_density_thresholds <- function(vals, thresholds, title = "", out_png) {
  # Basic density plot with thresholds
  png(out_png, width = 1200, height = 800)
  dens <- density(vals)
  plot(dens$x, dens$y, type = "l", main = title, xlab = "Value", ylab = "Density")
  abline(v = c(thresholds$left, thresholds$right), col = "red", lty = 2)
  dev.off()
}

# Given original vector vals and a sequence of methods, return logical mask of kept indices
gating_mask <- function(vals, methods, plot_dir = NULL, sample_id = NULL, channel = NULL, keep_init = NULL, step_start = 1) {
  n <- length(vals)
  keep <- if (is.null(keep_init)) rep(TRUE, n) else keep_init
  idx <- which(keep)
  current_vals <- vals[idx]
  step <- 1
  for (m in methods) {
    fn <- gating_methods[[as.character(m)]]
    if (is.null(fn)) stop(sprintf("Unknown method %s", m))
    th <- fn(current_vals, plot = FALSE)

    # Optional plotting to file
    if (!is.null(plot_dir)) {
      dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
      global_step <- step_start + step - 1
      title <- sprintf("%s | %s | step %02d | method %s", sample_id %||% "sample", channel %||% "channel", global_step, m)
      out_png <- file.path(plot_dir, sprintf("%s__%s__step%02d__method%s.png", sample_id %||% "sample", channel %||% "channel", global_step, m))
      plot_density_thresholds(current_vals, th, title, out_png)
    }

    pass <- current_vals >= th$left & current_vals <= th$right
    keep[idx] <- pass
    idx <- which(keep)
    current_vals <- vals[idx]
    step <- step + 1
  }
  keep
}

# null-coalescing helper
`%||%` <- function(a, b) if (!is.null(a)) a else b

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    cat_usage()
    quit(status = 1)
  }
  cfg <- parse_args(args)

  dir <- cfg$dir
  out_dir <- cfg$out
  plots <- cfg$plots
  plots_dir <- NULL
  if (!is.null(out_dir) && isTRUE(plots)) plots_dir <- file.path(out_dir, "plots") else if (isTRUE(plots)) plots_dir <- "Gating_Plots"

  # If a config file is provided, merge basic options (dir/out/plots) early
  if (!is.null(cfg$config)) {
    if (!file.exists(cfg$config)) stop(sprintf("Config file not found: %s", cfg$config))
    lines_early <- readLines(cfg$config, warn = FALSE)
    for (ln in lines_early) {
      s <- trimws(ln)
      if (nchar(s) == 0 || substr(s,1,1) %in% c("#",";")) next
      if (!grepl("=", s, fixed = TRUE)) next
      kv <- strsplit(s, "=", fixed = TRUE)[[1]]
      if (length(kv) < 2) next
      key <- tolower(trimws(kv[1]))
      val <- trimws(paste(kv[-1], collapse = "="))
      val <- sub('^"', '', sub('"$', '', val))
      val <- sub("^'", '', sub("'$", '', val))
      if (key == "dir" && is.null(dir)) dir <- val
      if (key == "out" && is.null(out_dir)) out_dir <- val
      if (key == "plots" && !isTRUE(plots)) plots <- tolower(val) %in% c("true","1","t","yes","y")
    }
    if (!is.null(out_dir) && isTRUE(plots)) plots_dir <- file.path(out_dir, "plots") else if (isTRUE(plots) && is.null(out_dir)) plots_dir <- "Gating_Plots"
  }

  # Early merge from config (dir/out/plots), if provided
  if (!is.null(cfg$config)) {
    lines <- readLines(cfg$config, warn = FALSE)
    for (ln in lines) {
      s <- trimws(ln)
      if (nchar(s) == 0 || substr(s,1,1) %in% c("#",";")) next
      if (grepl("=", s, fixed = TRUE)) {
        kv <- strsplit(s, "=", fixed = TRUE)[[1]]
        if (length(kv) >= 2) {
          key <- tolower(trimws(kv[1]))
          val <- trimws(paste(kv[-1], collapse = "="))
          val <- sub('^"', '', sub('"$', '', val))
          val <- sub("^'", '', sub("'$", '', val))
          if (key == "dir" && is.null(dir)) dir <- val
          if (key == "out" && is.null(out_dir)) out_dir <- val
          if (key == "plots" && !isTRUE(plots)) plots <- tolower(val) %in% c("true","1","t","yes","y")
        }
      }
    }
    if (!is.null(out_dir) && isTRUE(plots)) plots_dir <- file.path(out_dir, "plots") else if (isTRUE(plots) && is.null(out_dir)) plots_dir <- "Gating_Plots"
  }

  if (is.null(dir)) { cat_usage(); stop("Missing required arg: dir (in CLI or config)") }

  if (is.null(dir)) { cat_usage(); stop("Missing required arg: dir (in CLI or config)") }
  if (!dir.exists(dir)) stop(sprintf("Directory not found: %s", dir))

  fs <- read.flowSet(path = dir, pattern = "\\.fcs$", alter.names = TRUE, truncate_max_range = FALSE)
  if (length(fs) == 0) stop(sprintf("No .fcs files found under: %s", dir))

  # Build channel->methods mapping OR ordered pipeline
  param_names <- colnames(exprs(fs[[1]]))
  rules <- list()
  pipeline_steps <- NULL

  parse_pipeline <- function(txt) {
    # Accept separators ';' or '|'
    txt <- gsub("\\|", ";", txt)
    parts <- unlist(strsplit(txt, ";"))
    steps <- list()
    for (part in parts) {
      part <- trimws(part)
      if (nchar(part) == 0) next
      kv <- strsplit(part, ":", fixed = TRUE)[[1]]
      if (length(kv) != 2) stop(sprintf("Invalid pipeline segment: '%s'", part))
      ch <- trimws(kv[1])
      meth <- as.integer(strsplit(trimws(kv[2]), ",")[[1]])
      steps[[length(steps) + 1]] <- list(channel = ch, methods = meth)
    }
    steps
  }

  parse_config_file <- function(path) {
    if (!file.exists(path)) stop(sprintf("Config file not found: %s", path))
    lines <- readLines(path, warn = FALSE)
    steps <- list()
    opts <- list()
    for (ln in lines) {
      s <- trimws(ln)
      if (nchar(s) == 0) next
      if (substr(s, 1, 1) %in% c("#", ";")) next
      if (grepl("=", s, fixed = TRUE)) {
        kv <- strsplit(s, "=", fixed = TRUE)[[1]]
        if (length(kv) >= 2) {
          key <- trimws(kv[1])
          val <- trimws(paste(kv[-1], collapse = "="))
          val <- sub('^"', '', sub('"$', '', val))
          val <- sub("^'", '', sub("'$", '', val))
          if (tolower(key) == "plots") opts$plots <- tolower(val) %in% c("true", "1", "t", "yes", "y")
          else if (tolower(key) == "out") opts$out <- val
          else if (tolower(key) == "dir") opts$dir <- val
        }
        next
      }
      if (grepl(":", s, fixed = TRUE)) {
        kv <- strsplit(s, ":", fixed = TRUE)[[1]]
        if (length(kv) == 2) {
          ch <- trimws(kv[1])
          meth <- as.integer(strsplit(trimws(kv[2]), ",")[[1]])
          steps[[length(steps) + 1]] <- list(channel = ch, methods = meth)
          next
        }
      }
      stop(sprintf("Unrecognized config line: %s", s))
    }
    list(steps = steps, options = opts)
  }

  # If a config file is provided, parse steps into pipeline (options already merged earlier)
  if (!is.null(cfg$config)) {
    conf <- parse_config_file(cfg$config)
    if (length(conf$steps) > 0) pipeline_steps <- conf$steps
  }

  # Build selection from CLI flags if present
  if (!is.null(cfg$pipeline)) {
    pipeline_steps <- parse_pipeline(cfg$pipeline)
  }
  if (!is.null(cfg$rules)) {
    # Parse rules like "FSC.A:1,3;SSC.A:2,6"
    parts <- unlist(strsplit(cfg$rules, ";"))
    for (part in parts) {
      part <- trimws(part)
      if (nchar(part) == 0) next
      kv <- strsplit(part, ":", fixed = TRUE)[[1]]
      if (length(kv) != 2) stop(sprintf("Invalid rule segment: '%s'", part))
      ch <- trimws(kv[1])
      meth <- as.integer(strsplit(trimws(kv[2]), ",")[[1]])
      rules[[ch]] <- meth
    }
  }
  # Fallback to channels+methods only if no pipeline or rules provided
  if (is.null(pipeline_steps) && length(rules) == 0) {
    chs <- NULL
    if (!is.null(cfg$channels)) chs <- cfg$channels
    if (!is.null(cfg$channel)) chs <- unique(c(chs, cfg$channel))
    if (is.null(chs) || is.null(cfg$methods)) {
      cat_usage(); stop("Provide either config=... OR pipeline=... OR rules=... OR channel(s)+methods")
    }
    for (ch in chs) rules[[ch]] <- cfg$methods
  }

  # Validate channels exist
  if (!is.null(pipeline_steps)) {
    chs <- vapply(pipeline_steps, function(s) s$channel, character(1))
    for (ch in unique(chs)) {
      if (!(ch %in% param_names)) {
        cat("Available channels (first file):\n"); print(param_names)
        stop(sprintf("Channel '%s' not found. See above for available names.", ch))
      }
    }
  } else {
    for (ch in names(rules)) {
      if (!(ch %in% param_names)) {
        cat("Available channels (first file):\n"); print(param_names)
        stop(sprintf("Channel '%s' not found. See above for available names.", ch))
      }
    }
  }

  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  }

  summary <- data.frame(
    file = character(0),
    total = integer(0),
    kept = integer(0),
    kept_pct = numeric(0),
    stringsAsFactors = FALSE
  )

  res_list <- vector("list", length(fs))
  for (i in seq_along(fs)) {
    fr <- fs[[i]]
    smpl <- sampleNames(fs)[i]

    keep_all <- rep(TRUE, nrow(fr))

    if (!is.null(pipeline_steps)) {
      # Apply ordered pipeline across channels
      global_step <- 1
      for (st in pipeline_steps) {
        ch <- st$channel
        vals <- exprs(fr)[, ch]
        ch_plot_dir <- if (!is.null(plots_dir) && isTRUE(plots)) file.path(plots_dir, smpl, ch) else NULL
        keep_all <- gating_mask(vals, methods = st$methods, plot_dir = ch_plot_dir, sample_id = smpl, channel = ch, keep_init = keep_all, step_start = global_step)
        global_step <- global_step + length(st$methods)
      }
    } else {
      # Combine masks across channels (AND). Apply channels in the order provided.
      for (ch in names(rules)) {
        vals <- exprs(fr)[, ch]
        ch_plot_dir <- if (!is.null(plots_dir) && isTRUE(plots)) file.path(plots_dir, smpl, ch) else NULL
        mask_ch <- gating_mask(vals, methods = rules[[ch]], plot_dir = ch_plot_dir, sample_id = smpl, channel = ch, keep_init = keep_all)
        keep_all <- keep_all & mask_ch
      }
    }

    kept_fr <- fr[keep_all, ]
    res_list[[i]] <- kept_fr

    tot <- nrow(fr)
    kept <- sum(keep_all)
    summary[nrow(summary) + 1, ] <- list(smpl, tot, kept, round(100 * kept/tot, 2))

    if (!is.null(out_dir)) {
      out_file <- file.path(out_dir, paste0(smpl, "_gated.fcs"))
      suppressWarnings(write.FCS(kept_fr, out_file))
    }
  }

  # Print summary
  cat("\nGating summary\n", sep="")
  print(summary, row.names = FALSE)

  # Write summary CSV if out_dir is set
  if (!is.null(out_dir)) {
    write.csv(summary, file = file.path(out_dir, "gating_summary.csv"), row.names = FALSE)
  }

  invisible(res_list)
}

if (sys.nframe() == 0) main()
