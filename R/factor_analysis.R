#' @export
#' @import survcomp survival
#'
factor_analysis_cox <- function(clin_factors, rfs, limit = NULL, string = F, 
    ignore.mul.auto = T) {
    
    time <- rfs[, 1]
    event <- rfs[, 2] == 1
    if (!is.null(limit)) {
        event[time > limit] <- F
        time[time > limit] <- limit
    }
    rfs <- survival::Surv(time, event)
    
    
    res_single <- sapply(1:ncol(clin_factors), function(i) {
        idx <- !is.na(clin_factors[, i]) & !is.na(rfs)
        hr <- tryCatch(survcomp::hazard.ratio(as.numeric(clin_factors[idx, 
            i]), rfs[idx, 1], rfs[idx, 2]), error = function(e) {
            return(NULL)
        })
        c(HR = hr$hazard.ratio, CI95lo = hr$lower, CI95hi = hr$upper, 
            P = hr$p.value)
    })
    colnames(res_single) <- colnames(clin_factors)
    res_single <- t(res_single)
    
    ind <- which(res_single[, "P"] < 0.05)
    
    if (string) {
        res_single <- data.frame(HR = sprintf("%3.2f (%3.2f - %3.2f)", 
            res_single[, 1], res_single[, 2], res_single[, 3]), P = signif(res_single[, 
            4], 2))
    }
    
    if (length(ind) < 2) {
        if (ignore.mul.auto) {
            res <- res_single
        } else {
            res_mul <- res_single
            res_mul[-ind, ] <- NA
            res <- t(dplyr::full_join(as.data.frame(t(res_single)), as.data.frame(t(res_mul))))
            colnames(res) <- rep(colnames(res_single), 2)
        }
        
    } else {
        icpi_model <- survival::coxph(rfs ~ ., data = (clin_factors[, 
            ind]))
        res_mul <- cbind(summary(icpi_model)$conf.int[, -2], summary(icpi_model)$coefficients[, 
            5])
        
        rownames(res_mul) <- rownames(res_single)[ind]
        
        if (string) {
            res_mul <- data.frame(HR = sprintf("%3.2f (%3.2f - %3.2f)", 
                res_mul[, 1], res_mul[, 2], res_mul[, 3]), P = signif(res_mul[, 
                4], 2))
        } else {
            colnames(res_mul) <- colnames(res_single)
        }
        
        res <- t(dplyr::full_join(as.data.frame(t(res_single)), as.data.frame(t(res_mul))))
        colnames(res) <- rep(colnames(res_single), 2)
    }
    
    data.frame(res)
    
}

#' @export
factor_analysis_logit <- function(clin_factors, label, string = F) {
    
    
    res_single <- sapply(1:ncol(clin_factors), function(i) {
        value <- clin_factors[, i]
        or <- tryCatch({
            calc_or(value, label)
        }, error = function(e) {
            return(NULL)
        })
        mylogit <- glm(label ~ value, family = binomial)
        p <- coef(summary(mylogit))[2, 4]
        
        c(OR = or[1], CI95lo = or[2], CI95hi = or[3], P = p)
    })
    
    colnames(res_single) <- colnames(clin_factors)
    res_single <- t(res_single)
    ind <- which(res_single[, "P"] < 0.05)
    if (string) {
        res_single <- data.frame(OR = sprintf("%3.2f (%3.2f - %3.2f)", 
            res_single[, 1], res_single[, 2], res_single[, 3]), P = signif(res_single[, 
            4], 2))
    }
    
    mylogit <- glm(label ~ ., family = binomial, data = clin_factors)
    lreg.or <- exp(cbind(OR = coef(mylogit), confint(mylogit)))
    p <- coef(summary(mylogit))[, 4]
    res_mul <- cbind(lreg.or, p)[-1, ]
    
    
    
    if (string) {
        res_mul <- data.frame(OR = sprintf("%3.2f (%3.2f - %3.2f)", res_mul[, 
            1], res_mul[, 2], res_mul[, 3]), P = signif(res_mul[, 4], 
            2))
    } else {
        colnames(res_mul) <- colnames(res_single)
    }
    
    res <- data.frame(res_single, res_mul)
    res
}
