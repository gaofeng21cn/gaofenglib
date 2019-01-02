#' @export
#'
eval_logit <- function(dat1, dat2, index = 1:ncol(dat1$data)) {
    
    dat1$data <- data.frame(dat1$data)
    dat2$data <- data.frame(dat2$data)
    
    data <- data.frame(Event = dat1$label, dat1$data[, index])
    model <- glm(Event ~ ., data = data, family = binomial(logit))
    rs1 <- predict(model, dat1$data)
    rs2 <- predict(model, dat2$data)
    
    aucs1 <- c(dat1 = as.numeric(pROC::roc(dat1$label, rs1)$auc), dat2 = as.numeric(pROC::roc(dat2$label, 
        rs2)$auc), switch = 0, z = 0)
    
    data <- data.frame(Event = dat2$label, dat2$data[, index])
    model <- glm(Event ~ ., data = data, family = binomial(logit))
    rs1 <- predict(model, dat1$data)
    rs2 <- predict(model, dat2$data)
    
    aucs2 <- c(dat1 = as.numeric(pROC::roc(dat1$label, rs1)$auc), dat2 = as.numeric(pROC::roc(dat2$label, 
        rs2)$auc), switch = 1, z = 0)
    
    # z-score
    dat1$data <- data.frame(scale(dat1$data))
    dat2$data <- data.frame(scale(dat2$data))
    
    data <- data.frame(Event = dat1$label, dat1$data[, index])
    model <- glm(Event ~ ., data = data, family = binomial(logit))
    rs1 <- predict(model, dat1$data)
    rs2 <- predict(model, dat2$data)
    
    aucs3 <- c(dat1 = as.numeric(pROC::roc(dat1$label, rs1)$auc), dat2 = as.numeric(pROC::roc(dat2$label, 
        rs2)$auc), switch = 0, z = 1)
    
    data <- data.frame(Event = dat2$label, dat2$data[, index])
    model <- glm(Event ~ ., data = data, family = binomial(logit))
    rs1 <- predict(model, dat1$data)
    rs2 <- predict(model, dat2$data)
    
    aucs4 <- c(dat1 = as.numeric(pROC::roc(dat1$label, rs1)$auc), dat2 = as.numeric(pROC::roc(dat2$label, 
        rs2)$auc), switch = 1, z = 1)
    
    rbind(aucs1, aucs2, aucs3, aucs4)
}
