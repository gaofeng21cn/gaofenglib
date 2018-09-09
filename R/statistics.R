#' @export
#' @import pROC
calc_logit <- function(value, label) {
  set.seed(100)

  roc2 <- pROC::roc(label, value)

  power.roc <- c(power.roc.test(roc2, alternative = "one.side")$power, NA, NA)

  auc <- pROC::ci(roc2)[c(2,1,3)]
  names(auc) <- c("auc", "95% ci(lo)", "95% ci(hi)")

  rets <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv",
            "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv")
  others <- pROC::coords(roc2, "best", ret=rets, best.policy="omit")

  #if(!is.null(dim(others))) others <- others[, which.max(others["sensitivity", ])]

  YhatFac <- cut(value, breaks=c(-Inf, others["threshold"], Inf), labels=c("lo", "hi"))
  odds <- vcd::oddsratio(table(label, YhatFac) , log=F)
  or <- c(exp(odds$coefficients), confint(odds)[1], confint(odds)[2])

  df <- data.frame(ci.coords(roc2, x="best", input = "threshold", ret=rets, progress="none", best.policy="omit", parallel=TRUE))
  df <- rbind(auc=auc, power=power.roc, or=or, data.frame(others, df[, -2]))

  colnames(df) <- c("best", "95% ci(lo)", "95% ci(hi)")

  df
}

#' @export
#' @import survivalROC

calc_cutoff_survivalroc <- function(rfs, rs, limit = 60) {
  p <- survivalROC::survivalROC(Stime = rfs[, 1], status = rfs[, 2], marker = rs, predict.time = limit, method="KM")
  idx <- with(p, which.min(1-TP+ FP))
  rs_cut <- p$cut.values[idx]
  rs_cut
}

#' @export
#' @import survcomp survival
#'
factor_analysis <- function(clin_factors, rfs) {
  res_single <- sapply(1:ncol(clin_factors), function(i) {
    idx <- !is.na(clin_factors[, i]) & !is.na(rfs)
    hr <- survcomp::hazard.ratio(as.numeric(clin_factors[idx, i]), rfs[idx, 1], rfs[idx, 2])
    c(HR=hr$hazard.ratio, CI95lo=hr$lower, CI95hi=hr$upper, P=hr$p.value)
  })
  colnames(res_single) <- colnames(clin_factors)
  res_single <- t(res_single)

  ind <- which(res_single[, 4] < 0.05)

  icpi_model <- survival::coxph(rfs ~ ., data=clin_factors[, ind])
  res_mul <- data.frame(summary(icpi_model)$conf.int[, -2], summary(icpi_model)$coefficients[, 5])
  colnames(res_mul) <- colnames(res_single)
  rownames(res_mul) <- rownames(res_single)[ind]

  res <- t(dplyr::full_join(as.data.frame(t(res_single)),
                            as.data.frame(t(res_mul))))

  colnames(res) <- rep(colnames(res_single), 2)
  res
}
