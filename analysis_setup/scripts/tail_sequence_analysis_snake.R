# ---- CLI args: in_ends out_tsv [rfuns_dir] [pure_A_threshold] [min_run] ----
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2)

in_ends  <- args[[1]]
out_tsv  <- args[[2]]
rfuns_dir <- if (length(args) >= 3 && nzchar(args[[3]])) args[[3]] else "scripts/R_functions"
pure_A    <- if (length(args) >= 4 && nzchar(args[[4]])) as.numeric(args[[4]]) else 0.98
min_run   <- if (length(args) >= 5 && nzchar(args[[5]])) as.integer(args[[5]]) else 5L

# ensure output dir exists
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(stringi)
})

# 1) load tail end sequences
tail_end <- read_tsv(
  in_ends,
  col_select = c("read_id","hit_mode","extracted"),
  col_types  = cols(
    read_id  = col_character(),
    hit_mode = col_character(),
    extracted = col_character()
  )
)

################################################################################
# analyse 3'-end sequences
################################################################################

# 1) strip adapter at the 3' end in the original orientation
tail_end <- tail_end %>%
  mutate(
    extracted = toupper(extracted),
    # remove any trailing stuff after CTGTAG if present
    extracted = str_remove(extracted, "CTGTAG.*$"),
    # for deletions use 6N+CACCATCA, otherwise 7N+CACCATCA
    extracted = if_else(
      hit_mode == "del",
      str_remove(extracted, "[ACGTN]{6}CACCATCA.*$"),
      str_remove(extracted, "[ACGTN]{7}CACCATCA.*$")
    )
  )

# 2) reverse so the original 3' end is at the left
tail_end <- tail_end %>%
  mutate(extracted_rev = stri_reverse(extracted))

# 3) run main function analysing soft-clipped tail sequence composition
source("scripts/R_functions/run_AT_tail_analysis.R")
tail_end <- run_AT_tail_analysis(df = tail_end)

# 4) categorise tail types
source("scripts/R_functions/classify_from_counts.R")
tail_end <- classify_from_counts(tail_end, min_run = 5L)

# 5) extra patch for classification
#    - define N-insensitive length for purity checks
pure_A_threshold <- 0.9
tail_end <- tail_end %>%
  mutate(
    prop_A        = if_else(tail_len > 0, A_count / tail_len, NA_real_),
    class = case_when(
      # A-only (≥ threshold of A over non-Ns)
      class == "A" & prop_A > pure_A_threshold ~ "A_pure",
      # T/G/C-only (strict equality over non-Ns)
      class == "T" & T_count == tail_len   ~ "T_pure",
      class == "G" & G_count == tail_len   ~ "G_pure",
      class == "C" & C_count == tail_len   ~ "C_pure",
      # legacy remap
      class == "A_ten_A"                   ~ "A_pure",
      TRUE                                 ~ class
    )
  )

###############################################################################
# write the file (to current working directory)
###############################################################################
write_tsv(tail_end, out_tsv)

message(sprintf("Wrote: %s (rows: %s) in %s", out_tsv, nrow(tail_end), getwd()))
