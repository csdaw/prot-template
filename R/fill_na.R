fill_na <- function(df, na_col, fill_col) {
  df[[na_col]][is.na(df[[na_col]])] <- ""
  df[[na_col]] <- as.character(ifelse(df[[na_col]] == "",
                                      df[[fill_col]],
                                      df[[na_col]]))
  return(df)
}
