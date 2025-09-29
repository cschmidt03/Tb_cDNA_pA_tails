# This function analyses 3'-tails sequence that were extracted from soft-clipped 
# region adjactent to adapter. The region length default was 10 or less if softclip
# was shorter. Extraction is done by python function.
run_AT_tail_analysis <- function(df, seq_col = "extracted_rev", min_run = 5L) {
  stopifnot(seq_col %in% names(df))
  library(stringi)
  
  x <- toupper(df[[seq_col]])
  x[is.na(x)] <- ""
  
  # pure leading run of given base (for the FIRST run)
  count_pure <- function(s, base) {
    m <- stri_extract_first_regex(s, paste0("^", base, "+"))
    ifelse(is.na(m), 0L, nchar(m))
  }
  
  # FIRST run (pure, no mismatches)
  first <- stri_sub(x, 1L, 1L)
  T_count <- ifelse(first == "T", count_pure(x, "T"), 0L)
  G_count <- ifelse(first == "G", count_pure(x, "G"), 0L)
  C_count <- ifelse(first == "C", count_pure(x, "C"), 0L)
  A_count <- ifelse(first == "A", count_pure(x, "A"), 0L)
  
  # strip the first run
  offset <- pmax(T_count, G_count, C_count, A_count)
  rest <- ifelse(nchar(x) > offset, stri_sub(x, offset + 1L), "")
  
  # SECOND stretch: allow ONE mismatch **anywhere** in the stretch.
  # Return the *length of the matched prefix*; require #target>=min_run inside it.
  count_second_one_anywhere <- function(s, base, min_run) {
    if (is.na(s) || s == "") return(0L)
    ch <- strsplit(s, "", fixed = TRUE)[[1]]
    mism <- 0L; len <- 0L; nbase <- 0L
    for (c in ch) {
      if (c != base) mism <- mism + 1L
      if (mism > 1L) break
      len <- len + 1L
      if (c == base) nbase <- nbase + 1L
    }
    if (nbase >= min_run) len else 0L
  }
  
  A_count2 <- vapply(rest, count_second_one_anywhere, integer(1), base = "A", min_run = min_run)
  T_count2 <- vapply(rest, count_second_one_anywhere, integer(1), base = "T", min_run = min_run)
  
  transform(
    df,
    A_count  = as.integer(A_count),
    T_count  = as.integer(T_count),
    G_count  = as.integer(G_count),
    C_count  = as.integer(C_count),
    A_count2 = as.integer(A_count2),
    T_count2 = as.integer(T_count2)
  )
}
