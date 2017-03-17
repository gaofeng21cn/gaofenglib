#' Plot Kaplan-Meier Curve
#'
#' This function plots Kaplan-Meier Curve for survival analysis.
#'
#' @param clinical a survival object created by Surv function in survival package
#' @param labels a vector containing subtyping labels of patients
#' @param annot a character indicating the annotation showed up in the plot
#' @param color a vector containing colors used for different subtypes
#' @param font a character indicating the font used in the plot (default: "Arial")
#' @param xlab a character indicating the label of x-axis (default: "Follow up (weeks)")
#' @param ylab a character indicating the label of y-axis (default: "DFS (prob.)")
#' @param lenged.pos a character indicating whether the legend is (default: "top")
#' @param risk.table a logical indicating whether show risk table or not (default: FALSE)
#' @return a ggplot2 object of the plot
#' @import ggplot2 cowplot
#' @export
#' @examples
#' clinical <- survival::Surv(t.rfs, e.rfs)
#' color.cms <- c("#E69E00","#0070B0","#CA78A6", "#009C73")
#' plot_KMCurve(clinical, labels, "GSE39582", color.cms)
plot_KMCurve <- function (clinical, labels, annot = NULL, color = NULL, font = "Arial",
                          xlab = "Follow up (weeks)", ylab = "DFS (prob.)", title = NULL, legend.pos = "top",
                          risk.table = F, palette = "nature")
{
  obj <- clinical ~ labels
  surv <- survival::survfit(obj)
  survstats <- survival::survdiff(obj)
  survstats$p.value <- 1 - pchisq(survstats$chisq, length(survstats$n) -
                                    1)
  if (!is.null(color)) {
    if (!is.null(names(color))) {
      labels <- factor(labels, levels = names(color))
    }
  } else {
    color <- "Set1"
    if(palette == "nature") color <- ggsci::pal_npg("nrc")(length(unique(labels)))
    if(palette == "lancet") color <- ggsci::pal_lancet("lanonc")(length(unique(labels)))
  }

  if (class(labels) == "factor") {
    legend.labs <- na.omit(levels(droplevels(labels)))
  } else {
    legend.labs <- na.omit(unique(labels))
    labels <- factor(labels, levels = legend.labs)
  }
  p <- survminer::ggsurvplot(surv, xlab = xlab, ylab = ylab, main = title,
                             palette = color, legend = legend.pos, legend.title = NULL,
                             legend.labs = legend.labs, risk.table = risk.table, risk.table.title = element_blank(),
                             risk.table.y.text = FALSE, ggtheme = theme(text = element_text(family = font)))
  p$plot <- p$plot + annotate("text", family = font, x = Inf,
                              y = Inf, label = ifelse(survstats$p.value == 0, "italic(P)<1%*%10^{-22}",
                                                      paste0("italic(P)==", fancy_scientific(survstats$p.value,
                                                                                             3))), hjust = 1.2, vjust = 2, parse = TRUE)
  if (!is.null(annot))
    p$plot <- p$plot + annotate("text", x = 0, y = 0, label = annot,
                                hjust = 0, vjust = 0)
  if (risk.table) {
    p$table <- p$table + theme(axis.title.y = element_blank())
    return(p)
  } else return(p$plot)
}

#' Transform scientific notation to expression form
fancy_scientific <- function(l, digits = 3) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE, digits = digits)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

#' @export
#' @import ggplot2 cowplot precrec
plot_ROC <-  function(scores, labels)
{
  sscurves <- precrec::evalmod(scores = scores, labels = labels)
  roc <- precrec::auc(sscurves)[1, 4]
  if (auc(sscurves)$aucs[1] < 0.5) {
    sscurves <- precrec::evalmod(scores = -scores, labels = labels)
    roc <- precrec::auc(sscurves)[1, 4]
  }
  autoplot(sscurves, curvetype = "ROC") + cowplot::theme_cowplot() +
    theme(legend.position = "none", text = element_text(family = "Arial")) +
    annotate("text", x = 0.2, y = 0.8, label = paste0("AUC, ", round(roc, 3)), size = 4) +
    ggsci::scale_color_npg() + xlab("False Positive") + ylab("True Positive") +   scale_y_continuous(labels=percent) + scale_x_continuous(labels=percent)
}



#' @export
#' @import ggplot2 cowplot
plot_RiskScore <- function(rs, event) {
  if(is.logical(event)) event <- factor(event, levels = c(T, F), labels = c("Dead/Recurrence", "Disease free"))

  df <- data.frame(pt=names(rs), rs=rs, event=event)
  df <- df %>% arrange(rs)
  df$pt <- factor(df$pt, levels = as.character(df$pt))

  ggplot(df, aes(pt, rs, fill=event)) + geom_bar(stat="identity", alpha=0.7) +
    scale_fill_brewer(palette = "Set1") + ylab("Risk score") +
    theme(axis.text.x=element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.8,0.2), legend.key.width = unit(1, "cm"))
}



