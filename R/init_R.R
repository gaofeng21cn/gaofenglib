init_R <- function() {
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("knitr", "dplyr", "reshape2", "survival", "data.table", "DT", "shiny", "xtable",
             "devtools", "pROC", "survcomp", "extrafont", "extrafontdb", "RColorBrewer", "caret", "GEOquery", "Hmisc"))

  biocLite(c("org.Hs.eg.db"))

  devtools::install_github("gaofeng21cn/gaofenglib")
  devtools::install_github("renozao/NMF")
}