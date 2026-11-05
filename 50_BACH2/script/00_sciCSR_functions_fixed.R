# Corrected copy of sciCSR::getIGHreadType (stray ignore.case = "" removed from gsub calls)
getIGHreadType_fixed <- function(tb) {
  # fixed
  read_info <- tb[, c("CB", "UB")]
  tb <- as.matrix(tb[, -which(colnames(tb) %in% c("CB", "UB"))])
  columns <- colnames(tb)
  c_col <- which(grepl("C$", columns))
  i_col <- which(grepl("I$", columns))
  j_col <- which(columns == "VDJ")
  out <- apply(tb, MARGIN = 1, function(x) {
    o <- c("isotype" = NA, "transcript_type" = NA)
    c_reads <- x[c_col]; names(c_reads) <- columns[c_col]
    i_reads <- x[i_col]; names(i_reads) <- columns[i_col]
    all_c <- stringr::str_extract_all(names(c_reads[c_reads > 0]), "I[Gg][Hh][MmDdEeGgAaTtYyXxWw]")
    all_i <- stringr::str_extract_all(names(i_reads[i_reads > 0]), "I[Gg][Hh][MmDdEeGgAaTtYyXxWw]")
    all_c <- unique(unlist(all_c))
    all_i <- unique(unlist(all_i))
    if (length(all_c) > 1 | length(all_i) > 1) {
      return(paste(o, collapse = "_"))
    }
    if (all(i_reads == 0)) {
      c_isotypes <- gsub("_C$", "", columns[c_col[which(c_reads == max(c_reads))]])
      if (length(c_isotypes) == 1) {
        if (x[j_col] > 0) o <- c("isotype" = c_isotypes, "transcript_type" = "P")
        else o <- c("isotype" = c_isotypes, "transcript_type" = "C")
      }
      return(paste(o, collapse = "_"))
    } else {
      i_isotypes <- gsub("_I$", "", columns[i_col[which(i_reads == max(i_reads))]])
      if (all(c_reads == 0)) {
        if (length(i_isotypes) == 1) {
          if (x[j_col] > 0) o <- c("isotype" = i_isotypes, "transcript_type" = "P")
          else o <- c("isotype" = i_isotypes, "transcript_type" = "S")
        }
        return(paste(o, collapse = "_"))
      }
      c_isotypes <- gsub("_C$", "", columns[c_col[which(c_reads == max(c_reads))]])
      if (length(i_isotypes) == 1) {
        if (length(c_isotypes) >= 1 && i_isotypes %in% c_isotypes) {
          if (x[j_col] > 0) o <- c("isotype" = i_isotypes, "transcript_type" = "P")
          else o <- c("isotype" = i_isotypes, "transcript_type" = "S")
        }
        return(paste(o, collapse = "_"))
      } else {
        if (length(c_isotypes) == 1 && c_isotypes %in% i_isotypes) {
          if (x[j_col] > 0) o <- c("isotype" = c_isotypes, "transcript_type" = "P")
          else o <- c("isotype" = c_isotypes, "transcript_type" = "S")
        }
        return(paste(o, collapse = "_"))
      }
    }
  })
  out <- data.frame(read_info, out, stringsAsFactors = FALSE)
  colnames(out) <- c("CB", "UB", "anno")
  out <- out[which(out[, 3] != "NA_NA"), ]
  out
}


write_h5ad_for_sciCSR <- function(obj, file,
                                  harmony = "subset_RNA_harmony",
                                  umap = "subset_umap.harmony") {
  ad <- reticulate::import("anndata")
  sc <- reticulate::import("scanpy")
  
  # Expression: scaled data for the variable features (what sciCSR would write)
  X <- t(as.matrix(LayerData(obj, assay = "RNA", layer = "scale.data")))
  
  # Metadata: factors and glue to character; anndata cannot write NA in string columns
  obs <- obj[[]]
  for (col in colnames(obs)) {
    if (is.factor(obs[[col]]) || inherits(obs[[col]], "glue")) obs[[col]] <- as.character(obs[[col]])
    if (is.character(obs[[col]])) obs[[col]][is.na(obs[[col]])] <- "NA"
  }
  obs <- obs[rownames(X), , drop = FALSE]
  
  adata <- ad$AnnData(X = X, obs = obs,
                      var = data.frame(row.names = colnames(X)))
  reticulate::py_set_item(adata$obsm, "X_harmony", Embeddings(obj, harmony)[rownames(X), ])
  reticulate::py_set_item(adata$obsm, "X_umap",    Embeddings(obj, umap)[rownames(X), ])
  
  # Neighbour graph from the batch-corrected embedding
  sc$pp$neighbors(adata, use_rep = "X_harmony")
  adata$write_h5ad(file)
  invisible(file)
}

