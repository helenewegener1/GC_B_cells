function (input.data, samples = NULL, ID = NULL, chain = "both", 
          sequence = "nt", dist_type = "levenshtein", dist_mat = "BLOSUM80", 
          normalize = "length", gap_open = -10, gap_extend = -1, call.related.clones = TRUE, 
          group.by = NULL, threshold = 0.85, cluster.method = "components", 
          use.V = TRUE, use.J = TRUE, removeNA = FALSE, removeMulti = FALSE, 
          filterMulti = TRUE, filterNonproductive = TRUE) 
{
  processed_list <- input.data %>% .checkList() %>% .checkContigs() %>% 
    unname() %>% purrr::imap(function(x, i) {
      x <- subset(x, chain %in% c("IGH", "IGK", "IGL"))
      if (!is.null(ID)) 
        x$ID <- ID[i]
      if (filterNonproductive && "productive" %in% colnames(x)) {
        x <- subset(x, tolower(productive) == "true")
      }
      if (filterMulti) {
        x$save_chain <- x$chain
        x$chain <- ifelse(x$chain == "IGH", "IGH", "IGLC")
        x <- .filteringMulti(x)
        x$chain <- x$save_chain
        x$save_chain <- NULL
      }
      x
    }) %>% (function(x) {
      if (!is.null(samples)) {
        .modifyBarcodes(x, samples, ID)
      }
      else {
        x
      }
    }) %>% lapply(function(x) {
      data2 <- data.frame(x)
      data2 <- .makeGenes(cellType = "B", data2)
      unique_df <- unique(data2$barcode)
      Con.df <- data.frame(matrix(NA, length(unique_df), 9))
      colnames(Con.df) <- c("barcode", heavy_lines, light_lines)
      Con.df$barcode <- unique_df
      Con.df <- .parseBCR(Con.df, unique_df, data2)
      Con.df <- .assignCT(cellType = "B", Con.df)
      if (!is.null(group.by)) {
        Con.df[[group.by]] <- data2[[group.by]][1]
      }
      Con.df %>% mutate(length1 = nchar(cdr3_nt1)) %>% mutate(length2 = nchar(cdr3_nt2))
    })
  if (call.related.clones) {
    clusters <- clonalCluster(processed_list, sequence = sequence, 
                              chain = chain, threshold = threshold, group.by = group.by, 
                              use.V = use.V, use.J = use.J, cluster.method = cluster.method, 
                              dist_type = dist_type, dist_mat = dist_mat, normalize = normalize, 
                              gap_open = gap_open, gap_extend = gap_extend)
  }
  list_names <- if (!is.null(samples)) {
    if (is.null(ID)) 
      samples
    else paste0(samples, "_", ID)
  }
  else {
    paste0("S", seq_along(processed_list))
  }
  final_list <- purrr::map2(processed_list, seq_along(processed_list), 
                            function(df, i) {
                              if (call.related.clones) {
                                cluster_col <- clusters[[i]][, ncol(clusters[[i]])]
                                seq_col <- ifelse(sequence == "aa", "cdr3_aa", 
                                                  "cdr3_nt")
                                heavy_seq_col <- paste0(seq_col, "1")
                                light_seq_col <- paste0(seq_col, "2")
                                heavy_unique <- ifelse(!is.na(df[, "vgene1"]) & 
                                                         !is.na(df[, heavy_seq_col]), paste0(df[, "vgene1"], 
                                                                                             ".", df[, heavy_seq_col]), "NA")
                                light_unique <- ifelse(!is.na(df[, "vgene2"]) & 
                                                         !is.na(df[, light_seq_col]), paste0(df[, "vgene2"], 
                                                                                             ".", df[, light_seq_col]), "NA")
                                if (chain == "both") {
                                  heavy_part <- ifelse(is.na(cluster_col), heavy_unique, 
                                                       cluster_col)
                                  light_part <- ifelse(is.na(cluster_col), light_unique, 
                                                       cluster_col)
                                }
                                else if (chain == "IGH") {
                                  heavy_part <- ifelse(is.na(cluster_col), heavy_unique, 
                                                       cluster_col)
                                  light_part <- light_unique
                                }
                                else if (chain %in% c("IGL", "IGK", "Light")) {
                                  heavy_part <- heavy_unique
                                  light_part <- ifelse(is.na(cluster_col), light_unique, 
                                                       cluster_col)
                                }
                                else {
                                  heavy_part <- heavy_unique
                                  light_part <- light_unique
                                }
                                df[, "CTstrict"] <- paste0(heavy_part, "_", light_part)
                              }
                              else {
                                df[, "CTstrict"] <- paste0(df[, "vgene1"], ".", 
                                                           df[, "cdr3_aa1"], "_", df[, "vgene2"], ".", 
                                                           df[, "cdr3_aa2"])
                              }
                              if (!is.null(samples)) 
                                df$sample <- samples[i]
                              if (!is.null(ID)) 
                                df$ID <- ID[i]
                              df[df == "NA_NA" | df == "NA.NA_NA.NA" | df == "NA;NA_NA;NA" | 
                                   df == "NA"] <- NA
                              col_selection <- c("barcode", "sample", "ID", heavy_lines[c(1, 
                                                                                          2, 3)], light_lines[c(1, 2, 3)], CT_lines)
                              col_selection <- col_selection[col_selection %in% 
                                                               names(df)]
                              df <- df[, col_selection]
                              df <- df[!duplicated(df$barcode), ]
                              df <- df[rowSums(is.na(df)) < (ncol(df) - 1), ]
                              return(df)
                            })
  names(final_list) <- list_names
  if (removeNA) 
    final_list <- .removingNA(final_list)
  if (removeMulti) 
    final_list <- .removingMulti(final_list)
  return(final_list)
}
<bytecode: 0x84c3a5578>
  <environment: namespace:scRepertoire>