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
