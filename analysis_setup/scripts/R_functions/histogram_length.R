# This function works on transcriptome summary file and prepares data frame
# with number of reads in different length bins to check overall real length
# distribution in my Nanopore data
histogram_length <- function(data,
                             length_col = "sequence_length_template",
                             max_len = 6000,
                             bin = 50) {
  # prepare
  breaks <- seq(0, max_len, by = bin)
  len <- dplyr::pull(data, !!ensym(length_col))
  
  # make a bins table to join on (so we keep empty bins in order)
  bins_df <- tibble(
    bin = cut(breaks[-length(breaks)],
              breaks = breaks, right = FALSE, include.lowest = TRUE),
    bin_start = head(breaks, -1),
    bin_end   = tail(breaks, -1)
  )
  
  # count in-range bins
  in_range_counts <- tibble(
    bin = cut(len, breaks = breaks, right = FALSE, include.lowest = TRUE)
  ) %>%
    count(bin, name = "n")
  
  # overflow (> max_len)
  overflow_n <- sum(len > max_len, na.rm = TRUE)
  
  # assemble histogram table
  hist_tbl <- bins_df %>%
    left_join(in_range_counts, by = "bin") %>%
    mutate(n = tidyr::replace_na(n, 0L)) %>%
    select(bin_start, bin_end, n)
  
  # append overflow row (">max_len")
  if (overflow_n > 0) {
    hist_tbl <- bind_rows(
      hist_tbl,
      tibble(bin_start = max_len, bin_end = NA_integer_, n = overflow_n)
    )
  }
  
  hist_tbl %>%
    mutate(
      bin_label = ifelse(is.na(bin_end),
                         paste0(">", max_len),
                         paste0("[", bin_start, ",", bin_end, ")")),
      prop  = n / sum(n),
      cprop = cumsum(prop)
    ) %>%
    select(bin_start, bin_end, bin_label, n, prop, cprop)
}