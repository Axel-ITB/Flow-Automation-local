#!/usr/bin/env Rscript

# Extract text from a PDF into .txt files using pdftools
# Usage examples:
#   Rscript Extract_PDF_Text.r                       # uses default file
#   Rscript Extract_PDF_Text.r file="path/to/file.pdf"
#   Rscript Extract_PDF_Text.r file="Gating_Strategy_Pages_17_to_33.pdf" out="Gating_Strategy_Text"

args <- commandArgs(trailingOnly = TRUE)

kv <- strsplit(args, "=", fixed = TRUE)
cfg <- list(file = "Gating_Strategy_Pages_17_to_33.pdf", out = NULL)
for (p in kv) {
  if (length(p) == 2) {
    key <- trimws(p[[1]]); val <- trimws(p[[2]])
    val <- sub('^"', '', sub('"$', '', val))
    val <- sub("^'", '', sub("'$", '', val))
    cfg[[key]] <- val
  }
}

if (!requireNamespace("pdftools", quietly = TRUE)) {
  stop("Package 'pdftools' is required. Install with install.packages('pdftools') and rerun.")
}

pdf_file <- cfg$file
if (!file.exists(pdf_file)) stop(sprintf("PDF not found: %s", pdf_file))

message("Reading PDF: ", pdf_file)
txt_pages <- pdftools::pdf_text(pdf_file)

base_out <- cfg$out
if (is.null(base_out) || nchar(base_out) == 0) {
  # default base name without extension
  base_out <- sub("\\.pdf$", "", basename(pdf_file), ignore.case = TRUE)
}

dir.create(base_out, showWarnings = FALSE)

# Write combined text
combined_path <- file.path(base_out, paste0(basename(base_out), "_text.txt"))
writeLines(txt_pages, con = combined_path)

# Write per-page files
for (i in seq_along(txt_pages)) {
  page_path <- file.path(base_out, sprintf("%s_page_%02d.txt", basename(base_out), i))
  writeLines(txt_pages[[i]], con = page_path)
}

message("Wrote: ", combined_path)
message("Per-page files under: ", base_out)

