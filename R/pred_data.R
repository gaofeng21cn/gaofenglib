#' @export
clean_dat <- function(gene.exp, keytype = "ENTREZID", column = "SYMBOL") {
  gene.exp <- gene.exp[, !(is.na(colnames(gene.exp)) | colnames(gene.exp) ==
                             "")]
  gene.exp <- gene.exp[, !sapply(1:ncol(gene.exp), function(i) sum(!is.na(gene.exp[,i])) == 0)]
  if (any(grepl("///", colnames(gene.exp)))) {
    ids <- strsplit(colnames(gene.exp), split = " /// ")
    lt.ids <- sapply(ids, length)
    gene.exp <- do.call(cbind, lapply(1:length(ids), function(i) sapply(1:lt.ids[i],
                                                                        function(x) gene.exp[, i])))
    colnames(gene.exp) <- unlist(ids)
  }
  if (keytype != column) {
    library(org.Hs.eg.db)
    colnames(gene.exp) <- mapIds(org.Hs.eg.db, keys = colnames(gene.exp),
                                 keytype = keytype, column = column)
  }
  gene.exp <- gene.exp[, !is.na(colnames(gene.exp))]
  m <- apply(gene.exp, 2, function(x) mad(x, na.rm = T))
  i <- tapply(m, colnames(gene.exp), which.max)
  ind <- tapply(1:ncol(gene.exp), colnames(gene.exp), function(x) x)
  gene.exp[, mapply(function(a, b) a[b], ind, i)]

}

#' @export
calc_oncotypedx_crc <- function(dat) {
    dat <- t(dat)
    stroma <- c("BGN", "FAP", "INHBA")
    CC <- c("MKI67", "MYC", "MYBL2")
    Individual <- "GADD45B"
    RF <- c("ATP5E", "GPX1", "PGK1", "VDAC2", "UBB")
    ONDX <- c(stroma, CC, Individual, RF)
    ge.onco <- dat[na.omit(match(ONDX, rownames(dat))), ]

    stroma_sc <- apply(ge.onco[na.omit(match(stroma, rownames(ge.onco))),
        ], 2, mean)
    CC_sc <- apply(ge.onco[na.omit(match(CC, rownames(ge.onco))), ],
        2, mean)
    Individual_sc <- as.numeric(as.vector(ge.onco[na.omit(match(Individual,
        rownames(ge.onco))), ]))
    RF_sc <- apply(ge.onco[na.omit(match(RF, rownames(ge.onco))), ],
        2, mean)
    corrected_stroma <- stroma_sc - RF_sc + 10
    corrected_CC <- CC_sc - RF_sc + 10
    corrected_individual <- Individual_sc - RF_sc + 10
    RS_score <- 0.15 * corrected_stroma - 0.3 * corrected_CC + 0.15 *
        corrected_individual
    oncotypeDX_score <- 44 * (RS_score + 0.82)
    oncotypeDX_class <- cut(oncotypeDX_score, breaks = c(-Inf, 30, 41,
        Inf), labels = c("Low", "Intermediate", "High"))
    DX <- data.frame(stroma = stroma_sc, CC = CC_sc, Individual = Individual_sc,
        RF = RF_sc, corrected_stroma = corrected_stroma, corrected_CC = corrected_CC,
        corrected_individual = corrected_individual, RS_score = RS_score,
        oncotypeDX_score = oncotypeDX_score, oncotypeDX_class = oncotypeDX_class)
    DX
}


#' @export
#' @import foreach
calc_resamp_cox <- function(data, rfs, resamp.ratio = 0.8, resamp.time = 1000, cores = 40, verbose=F) {
  doParallel::registerDoParallel(cores)
  set.seed(100)

  foreach(i = 1:resamp.time, .combine=cbind) %dopar% {
    ind <- sample(nrow(data), nrow(data) * resamp.ratio)
    if(verbose) cat(paste("Resampling:", i))
    apply(data, 2, function(gene) tryCatch({summary(survival::coxph(rfs[ind, ] ~ gene[ind]))$coefficients[5]},
                                           error= function(e) {return(NA)}))
  }
}
