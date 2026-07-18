# ---- inputs/outputs from CLI ----
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 4)
in_file        <- args[[1]]   # TX_SUMMARY
in_names       <- args[[2]]   # READ_IDS_TSV
classification <- args[[3]]   # TAIL_CLASS_TSV
output         <- args[[4]]   # OUTDIR/<prefix>.final_table.tsv
logdir         <- if (length(args) >= 5 && nzchar(args[[5]])) args[[5]] else file.path(dirname(output), "R_logs")

# make sure output dir exists
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)
################################################################################
# libraries
################################################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

################################################################################
# read data
################################################################################
desired_cols <- c(
  "read_id", "alignment_genome", "alias", "sequence_length_template",
  "poly_tail_length", "alignment_mapping_quality", "alignment_genome_start",
  "alignment_genome_end", "alignment_strand_start", "alignment_strand_end",
  "alignment_direction"
)

df <- read_tsv(
  in_file,
  col_select = all_of(desired_cols),
  col_types = cols(
    read_id                  = col_character(),
    alignment_genome         = col_character(),
    alias		  = col_character(),
    sequence_length_template = col_integer(),
    poly_tail_length         = col_integer(),
    alignment_mapping_quality  = col_integer(),
    alignment_genome_start   = col_integer(),
    alignment_genome_end     = col_integer(),
    alignment_strand_start   = col_integer(),
    alignment_strand_end     = col_integer(),
    alignment_direction      = col_character()
  ),
  show_col_types = FALSE
) %>% mutate(alignment_length = alignment_genome_end - alignment_genome_start) %>%
	rename(barcode = alias, alignment_mapq = alignment_mapping_quality)

#renaming and alignment_length for compatibility, column names of dorado summary were changed
#somewhere between version 0.9.X and 2.X 

names_filtered <- read_tsv(
  in_names,
  col_types = cols(read_id = col_character()),
  show_col_types = FALSE
) %>% distinct(read_id, .keep_all = TRUE)

################################################################################
# stats BEFORE filtering
################################################################################
data_bf <- df %>%
  summarise(
    counts = n(),
    mapped = sum(alignment_genome != "*", na.rm = TRUE),
    tailed = sum(poly_tail_length > 0, na.rm = TRUE)
  )

source("scripts/R_functions/histogram_length.R")
hist_before_filter <- histogram_length(
  df, length_col = "sequence_length_template",
  max_len = 4000, bin = 50
)

write_tsv(data_bf, file = file.path(logdir, "data_bf.tsv"))
write_tsv(hist_before_filter, file = file.path(logdir, "hist_before_filter.tsv"))

################################################################################
# FILTERING (keep full rows; do NOT call count() here)
################################################################################
df_filt <- df %>%
  semi_join(names_filtered, by = "read_id") %>%
  filter(!is.na(alignment_genome), alignment_genome != "*")

################################################################################
# stats AFTER filtering
################################################################################
data_af <- df_filt %>%
  summarise(
    counts = n(),
    mapped = sum(alignment_genome != "*", na.rm = TRUE),  # will equal counts
    tailed = sum(poly_tail_length > 0, na.rm = TRUE)
  )

hist_after_filter <- histogram_length(
  df_filt, length_col = "sequence_length_template",
  max_len = 4000, bin = 50
)

write_tsv(data_af, file = file.path(logdir, "data_af.tsv"))
write_tsv(hist_after_filter, file = file.path(logdir, "hist_after_filter.tsv"))

################################################################################
# add 3'-end annotations (left join keeps all reads from df_filt)
################################################################################
end_seq <- read_tsv(
  classification,
  col_select = c("read_id", "class", "tail_len", "T_count", "G_count", "C_count", "prop_A"),
  col_types = cols(
    read_id  = col_character(),
    class    = col_character(),
    tail_len = col_integer(),
    T_count  = col_integer(),
    G_count  = col_integer(),
    C_count  = col_integer(),
    prop_A   = col_number()
  ),
  show_col_types = FALSE
)

end_seq <- left_join(df_filt, end_seq, by = "read_id")

################################################################################
# save the annotated and filtered table
################################################################################
write_tsv(end_seq, output)


