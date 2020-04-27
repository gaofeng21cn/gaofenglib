#' @export
#' @import pROC
calc_logit <- function(value, label, string = F) {
    set.seed(100)
    roc2 <- pROC::roc(label, value)
    power.roc <- c(power.roc.test(roc2, alternative = "one.side")$power, 
        NA, NA)
    auc <- pROC::ci(roc2)[c(2, 1, 3)]
    names(auc) <- c("auc", "95% ci(lo)", "95% ci(hi)")
    rets <- c("threshold", "specificity", "sensitivity", "accuracy", 
        "tn", "tp", "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", 
        "1-accuracy", "1-npv", "1-ppv")
    others <- pROC::coords(roc2, "best", ret = rets, best.policy = "omit")
    YhatFac <- cut(value, breaks = c(-Inf, others["threshold"], Inf), 
        labels = c("lo", "hi"))
    odds <- vcd::oddsratio(table(label, YhatFac), log = F)
    or <- c(exp(odds$coefficients), confint(odds)[1], confint(odds)[2])
    df <- data.frame(ci.coords(roc2, x = "best", input = "threshold", 
        ret = rets, progress = "none", best.policy = "omit", parallel = TRUE))
    df <- rbind(auc = auc, power = power.roc, or = or, data.frame(others, 
        df[, -2]))
    colnames(df) <- c("best", "95% ci(lo)", "95% ci(hi)")
    if (string) {
        tmp <- sprintf("%3.2f (%3.2f - %3.2f)", df[, 1], df[, 2], df[, 
            3])
        names(tmp) <- rownames(df)
        tmp[2] <- sprintf("%3.2f", df["power", "best"])
        tmp[8:11] <- sprintf("%d (%d - %d)", as.integer(df[, 1]), as.integer(df[, 
            2]), as.integer(df[, 3]))[8:11]
        df <- tmp
        names(df) <- c("AUC", "Power", "Odds Ratio", "Threshold", "Specificity", 
            "Sensitivity", "Accuracy", "TN", "TP", "FN", "FP", "NPV", 
            "PPV", "1-Specificity", "1-Sensitivity", "1-Accuracy", "1-NPV", 
            "1-PPV")
    }
    
    df
}



#' @export
#' @import survivalROC

calc_cutoff_survivalroc <- function(rfs, rs, limit = 60) {
    p <- survivalROC::survivalROC(Stime = rfs[, 1], status = rfs[, 2], 
        marker = rs, predict.time = limit, method = "KM")
    idx <- with(p, which.min(1 - TP + FP))
    rs_cut <- p$cut.values[idx]
    rs_cut
}

#' @export
#' @import pROC

calc_or <- function(value, label) {
    roc2 <- pROC::roc(label, value)
    cutoff <- pROC::coords(roc2, "best", ret = "threshold", best.policy = "omit")
    YhatFac <- cut(value, breaks = c(-Inf, cutoff, Inf), labels = c("lo", 
        "hi"))
    odds <- vcd::oddsratio(table(label, YhatFac), log = F)
    or <- c(exp(odds$coefficients), confint(odds)[1], confint(odds)[2])
    as.numeric(or)
}

