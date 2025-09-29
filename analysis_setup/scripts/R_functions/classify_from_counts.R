# this function creates tails classification based on the extracted sequences and
# numerircal analysis of tail contents
classify_from_counts <- function(df, seq_col = "extracted_rev", min_run = 5L, mixed_thresh = 0.8) {
  stopifnot(all(c("A_count","T_count","G_count","C_count","A_count2","T_count2") %in% names(df)))
  stopifnot(seq_col %in% names(df))
  library(stringi); library(dplyr)
  
  x <- toupper(replace(df[[seq_col]], is.na(df[[seq_col]]), ""))
  tail_len <- nchar(x)
  
  A_total <- tail_len - nchar(stri_replace_all_fixed(x, "A", ""))
  T_total <- tail_len - nchar(stri_replace_all_fixed(x, "T", ""))
  G_total <- tail_len - nchar(stri_replace_all_fixed(x, "G", ""))
  C_total <- tail_len - nchar(stri_replace_all_fixed(x, "C", ""))
  
  first_base <- case_when(
    df$A_count > 0 ~ "A",
    df$T_count > 0 ~ "T",
    df$G_count > 0 ~ "G",
    df$C_count > 0 ~ "C",
    TRUE           ~ "None"
  )
  
  none_tail <- (df$A_count==0 & df$T_count==0 & df$G_count==0 & df$C_count==0 &
                  df$A_count2==0 & df$T_count2==0)
  
  A_pure   <- (first_base == "A" & df$A_count >= min_run)
  
  T_then_A <- (first_base == "T" & df$T_count >= 1L & df$A_count2 >= min_run)
  G_then_A <- (first_base == "G" & df$G_count >= 1L & df$A_count2 >= min_run)
  C_then_A <- (first_base == "C" & df$C_count >= 1L & df$A_count2 >= min_run)
  A_then_A <- (first_base == "A" & df$A_count >= 1L & df$A_count2 >= min_run)
  
  T_pure   <- (df$T_count >= min_run)
  G_pure   <- (df$G_count >= min_run)
  C_pure   <- (df$C_count >= min_run)
  
  A_then_T <- (first_base == "A" & df$A_count >= 1L & df$T_count2 >= min_run)
  G_then_T <- (first_base == "G" & df$G_count >= 1L & df$T_count2 >= min_run)
  C_then_T <- (first_base == "C" & df$C_count >= 1L & df$T_count2 >= min_run)
  
  just_A <- (first_base == "A")
  just_T <- (first_base == "T")
  just_G <- (first_base == "G")
  just_C <- (first_base == "C")
  
  mixed_A <- just_A & (A_total < mixed_thresh * tail_len) & (tail_len > 0)
  mixed_T <- just_T & (T_total < mixed_thresh * tail_len) & (tail_len > 0)
  mixed_G <- just_G & (G_total < mixed_thresh * tail_len) & (tail_len > 0)
  mixed_C <- just_C & (C_total < mixed_thresh * tail_len) & (tail_len > 0)
  
  class <- case_when(
    none_tail ~ "None",
    A_pure    ~ "A_pure",
    
    T_then_A  ~ "T_then_A",
    G_then_A  ~ "G_then_A",
    C_then_A  ~ "C_then_A",
    A_then_A  ~ "A_then_A",
    
    T_pure    ~ "T_pure",
    G_pure    ~ "G_pure",
    C_pure    ~ "C_pure",
    
    A_then_T  ~ "A_then_T",
    G_then_T  ~ "G_then_T",
    C_then_T  ~ "C_then_T",
    
    mixed_A | mixed_T | mixed_G | mixed_C ~ "mixed",
    
    just_A    ~ "A",
    just_T    ~ "T",
    just_G    ~ "G",
    just_C    ~ "C",
    TRUE      ~ "other"
  )
  
  mutate(df,
         class    = class,
         tail_len = tail_len,
         A_total  = A_total,
         T_total  = T_total,
         G_total  = G_total,
         C_total  = C_total
  )
}
