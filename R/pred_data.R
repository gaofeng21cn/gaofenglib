#' @export
clean_dat <- function(gene.exp, keytype = "ENTREZID", column = "SYMBOL") {
  if(any(grepl("///", colnames(gene.exp)))) {
    ids <- strsplit(colnames(gene.exp), split = " /// ")
    lt.ids <- sapply(ids, length)

    gene.exp <- do.call(cbind, lapply(1:length(ids), function(i) sapply(1:lt.ids[i], function(x) gene.exp[, i])))
    colnames(gene.exp) <- unlist(ids)
  }

  library(org.Hs.eg.db)
  colnames(gene.exp) <- mapIds(org.Hs.eg.db, keys=colnames(gene.exp), keytype = keytype, column = column)
  gene.exp <- gene.exp[, !is.na(colnames(gene.exp))]
  m <- apply(gene.exp, 2, mad)
  i <- tapply(m, colnames(gene.exp), which.max)
  ind <- tapply(1:ncol(gene.exp), colnames(gene.exp), function(x) x)
  gene.exp[, mapply(function(a, b) a[b], ind, i)]
}
