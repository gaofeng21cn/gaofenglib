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
plot_KMCurve <- function (time, event, labels, limit = NULL, annot = NULL, color = NULL, font = "Arial",
                           xlab = "Follow up (weeks)", ylab = "DFS (prob.)", title = NULL,
                           legend.pos = "top", risk.table = T, palette = "nature")
{
  if(!is.null(limit)) {
    clinical[clinical[, 1] > limit, 2] <- F
    clinical[clinical[, 1] > limit, 1] <- limit
  }

  df <- data.frame(futime=time, fustat=event, group=labels)
  surv <- survival::survfit(survival::Surv(futime, fustat) ~ group, data = df)

  survstats <- survival::survdiff(survival::Surv(futime, fustat) ~ group, data = df)
  survstats$p.value <- 1 - pchisq(survstats$chisq, length(survstats$n) -
                                    1)

  if (!is.null(color)) {
    if (!is.null(names(color))) {
      labels <- factor(labels, levels = names(color))
    }
  }
  else {
    color <- "Set1"
    if (palette == "nature")
      color <- (ggsci::pal_npg("nrc"))(length(unique(labels)))
    if (palette == "lancet")
      color <- (ggsci::pal_lancet("lanonc"))(length(unique(labels)))
    if (palette == "jco")
      color <- (ggsci::pal_jco("default"))(length(unique(labels)))
    if (palette == "jama")
      color <- c("#164870", "#10B4F3", "#FAA935", "#2D292A",
                 "#87AAB9", "#CAC27E", "#818282")[1:length(unique(labels))]
    if (palette == "jama_raju")
      color <- c("#3676BB", "#DDBB1B", "#858585", "#606060")[1:length(unique(labels))]
  }
  if (class(labels) == "factor") {
    legend.labs <- na.omit(levels(droplevels(labels[!is.na(clinical)])))
  }
  else if (class(labels) == "logical") {
    labels <- factor(labels, levels = c(F, T))
    legend.labs <- na.omit(levels(droplevels(labels)))
  }
  else {
    legend.labs <- na.omit(unique(labels))
    labels <- factor(labels, levels = legend.labs)
  }

  fancy_scientific <- function(l, dig=3) {
    # turn in to character string in scientific notation
    l <- format(l, digits = dig, scientific = TRUE)
    # quote the part before the exponent to keep all the digits
    l <- gsub("^(.*)e", "'\\1'e", l)
    # turn the 'e+' into plotmath format
    l <- gsub("e", "%*%10^", l)
    # return this as an expression
    parse(text=l)
  }

  p <- survminer::ggsurvplot(surv, xlab = xlab, ylab = ylab,
                             palette = color, legend = legend.pos,
                             legend.title = NULL, legend.labs = legend.labs, risk.table = risk.table,
                             risk.table.title = element_blank(), risk.table.y.text = FALSE,
                             ggtheme = theme(text = element_text(family = font)))
  p$plot <- p$plot + ggtitle(title) + annotate("text", family = font, x = 0,
                                               y = 0, label = ifelse(survstats$p.value == 0, "italic(P)<1%*%10^{-22}",
                                                                     paste0("italic(P)==", fancy_scientific(survstats$p.value, 3))), hjust = 0, vjust = -2, parse = TRUE)


  # HR
  if(length(legend.labs) == 2) {
    hr <- survcomp::hazard.ratio(labels[!is.na(clinical)], clinical[!is.na(clinical), 1], clinical[!is.na(clinical), 2])

    p$plot <- p$plot + annotate("text", x = 0, y = 0,
                                label = sprintf("HR = %3.2f (%3.2f - %3.2f)", hr$hazard.ratio, hr$lower, hr$upper),
                                hjust = 0, vjust = -1, parse = F)
  }

  if (!is.null(annot))
    p$plot <- p$plot + annotate("text", x = 0, y = 0, label = annot,
                                hjust = 0, vjust = 0)
  if (risk.table) {
    p$table <- p$table + theme(axis.title.y = element_blank())
    pp <- plot_grid(plotlist = list(p$plot + theme(axis.title.x = element_blank()), p$table + labs(x=xlab)), labels = "", ncol = 1, align = "hv", rel_heights = c(2.5,1))
    return(pp)
  }
  else return(p$plot)
}



switch_fill_color <- function(p, palette) {
  switch(palette, jco = {
    p + ggsci::scale_fill_jco()
  }, lancet = {
    p + ggsci::scale_fill_lancet()
  }, jama = {
    p + scale_fill_manual(values = c("#164870", "#10B4F3",
                                     "#FAA935", "#2D292A", "#87AAB9", "#CAC27E", "#818282"))
  }, lancet = {
    p + ggsci::scale_fill_npg()
  },
  p + scale_fill_brewer(palette = "Set1"))
}


#' @export
#' @import ggplot2 cowplot precrec scales
plot_ROC <-  function(scores, labels, fontsize=18, palette = "nature")
{
  sscurves <- precrec::evalmod(scores = scores, labels = labels)
  roc <- precrec::auc(sscurves)[1, 4]
  if (auc(sscurves)$aucs[1] < 0.5) {
    sscurves <- precrec::evalmod(scores = -scores, labels = labels)
    roc <- precrec::auc(sscurves)[1, 4]
  }
  p <- autoplot(sscurves, curvetype = "ROC") + cowplot::theme_cowplot(font_size = fontsize, font_family = "Arial", line_size = 1) +
    theme(legend.position = "none") +
    annotate("text", x = 0.7, y = 0.1, label = paste0("AUC, ", format(round(roc, 3), nsmall =3)), size=fontsize/3)
  switch(palette,
         "jco"= {
           p + ggsci::scale_color_jco() +
            xlab("False Positive") + ylab("True Positive") +   scale_y_continuous(labels=percent) + scale_x_continuous(labels=percent)
         },
         "lancet"= {
           p + ggsci::scale_color_lancet()
         },
         "jama"= {
           p + scale_color_manual(values = c("#164870", "#10B4F3", "#FAA935", "#2D292A", "#87AAB9", "#CAC27E", "#818282"))
         }, p + ggsci::scale_color_npg() )
}

#' @export
#' @import ggplot2 cowplot precrec
plot_multiROC <- function(scores, groups, fontsize=16, palette="nature") {
  msmdat1 <- precrec::mmdata(scores , groups, modnames = colnames(scores))
  mmcurves <- precrec::evalmod(msmdat1)

  inds <- subset(precrec::auc(mmcurves), curvetypes=="ROC")$auc < 0.5

  if(any(inds)) {
    scores[, inds] <- -scores[, inds]
    msmdat1 <- precrec::mmdata(scores , groups, modnames = colnames(scores))
    mmcurves <- precrec::evalmod(msmdat1)
  }

  p <- autoplot(mmcurves, "ROC")+ cowplot::theme_cowplot(font_size = fontsize, font_family = "Arial", line_size = 1)  +
    theme(legend.position = c(0.8, 0.2),
          legend.title = element_blank())
  switch(palette,
         "jco"= {
           p + ggsci::scale_color_jco() +
             xlab("False Positive") + ylab("True Positive") +   scale_y_continuous(labels=percent) + scale_x_continuous(labels=percent)
         },
         "lancet"= {
           p + ggsci::scale_color_lancet()
         },
         "jama"= {
           p + scale_color_manual(values = c("#164870", "#10B4F3", "#FAA935", "#2D292A", "#87AAB9", "#CAC27E", "#818282"))
         }, p + ggsci::scale_color_npg() )
}



#' @export
#' @import ggplot2 cowplot
plot_RiskScore <- function(rs, event, fontsize = 16, legend.position = c(0.2, 0.8), palette = "nature") {
  if(is.logical(event)) event <- factor(event, levels = c(T, F), labels = c("Dead/Recurrence", "Disease free"))

  if(is.null(names(rs))) names(rs) <- 1:length(rs)
  df <- data.frame(pt=names(rs), rs=rs, event=event)
  df <- df %>% arrange(rs)
  df$pt <- factor(df$pt, levels = as.character(df$pt))

  p <- ggplot(df, aes(pt, rs, fill=event)) + geom_bar(stat="identity", alpha=0.7) +
    cowplot::theme_cowplot(font_size = fontsize, font_family = "Arial") +
    ylab("Risk score") +
    theme(axis.text.x=element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          legend.position = legend.position, legend.key.width = unit(1, "cm"))

  switch(palette, jco = {
    p + ggsci::scale_fill_jco()
  }, lancet = {
    p + ggsci::scale_fill_lancet()
  }, jama = {
    p + scale_fill_manual(values = c("#164870", "#10B4F3",
                                     "#FAA935", "#2D292A", "#87AAB9", "#CAC27E", "#818282"))
  }, lancet = {
    p + ggsci::scale_fill_npg()
  },
  p + scale_fill_brewer(palette = "Set1"))
}

#' @export
#' @import ggplot2 cowplot
plot_Boxplot <- function(value, label, palette = "nature", fontsize = 18) {
  p <- qplot(x= label, y= value, geom= "boxplot", color= label) + cowplot::theme_cowplot(font_size = fontsize, font_family = "Arial", line_size = 1) + theme(legend.position = "none", axis.title = element_blank())

  switch(palette,
         "jco"= {
           p + ggsci::scale_color_jco()
         },
         "lancet"= {
           p + ggsci::scale_color_lancet()
         },
         "jama"= {
           p + scale_color_manual(values = c("#164870", "#10B4F3", "#FAA935", "#2D292A", "#87AAB9", "#CAC27E", "#818282"))
         }, p + ggsci::scale_color_npg() )
}

#' @export
#' @import ggfortify
plot_PCA <- function(data, labs, title="Evaluate the batch effect between groups") {
  library(ggfortify)

  df <- data.frame(group=labs, data, check.names = T)

  autoplot(prcomp(df[, -1]), data=df, colour="group") +
    theme(legend.title=element_blank()) +
    scale_color_brewer(palette = "Set1") +
    ggtitle(title)
}
