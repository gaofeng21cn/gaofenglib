#' @export
init_R <- function() {
  if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
  BiocManager::install(c("knitr", "dplyr", "reshape2", "survival", "data.table", "DT", "shiny", "xtable",
             "devtools", "pROC", "survcomp", "extrafont", "extrafontdb", "RColorBrewer", "caret", "GEOquery", "Hmisc"))

  BiocManager::install(c("org.Hs.eg.db"))

  devtools::install_github("gaofeng21cn/gaofenglib")
  devtools::install_github("gaofeng21cn/gfplot")
  devtools::install_github("rstudio/keras")
  devtools::install_github("renozao/NMF")
}

#' @export
#'
update_all <- function(auth_token = NULL) {
    devtools::install_github("gaofeng21cn/gaofenglib")
    devtools::install_github("gaofeng21cn/gfplot")
    devtools::install_github("CityUHK-CompBio/DeepCC@keras")
    devtools::install_github("gaofeng21cn/curatedClinicalData", auth_token = auth_token)
}
