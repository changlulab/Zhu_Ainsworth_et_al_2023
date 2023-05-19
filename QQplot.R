rm(list = ls())
options(stringsAsFactors = F)
library(ggplot2)
library(stats)
library(car)
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/ChIP-seq_differential/enhancer_count/differential_peaks/")
SCZ_1 = read.csv('All_trt_SCZ_Control_2_neuron_enhancer.csv')
# SCZ_1 = all_matrix
pvalues = data.frame(SCZ_1$pvalue)
# qqnorm(new_pvalues$pvalue)
# shapiro.test(new_pvalues$pvalue)
pvalues$BH_adj_p = p.adjust(pvalues$SCZ_1.pvalue,method = 'BH')
new_pvalues = data.frame(pvalues[order(pvalues$BH_adj_p),2])
colnames(new_pvalues) = 'pvalue'
# n = nrow(new_pvalues)                         
# new_pvalues$expected = (1:n -0.5)/n

setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/Cell reports_review/figures for revision/QQ plot/new/")

pdf('QQ_plot_glia_enhancer_SCZ vs Control.pdf', width = 6, height = 5)
# qqPlot(new_pvalues$pvalue, main = "Normal QQ Plot")
# abline(a = mean(new_pvalues$pvalue), b = sd(new_pvalues$pvalue), col = "red")
# legend("topleft", legend = paste0("lambda = ", round(sd(new_pvalues$pvalue), 2)), col = "red", lty = 1, bty = "n")

inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}
inflation(new_pvalues$pvalue)
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}
gg_qqplot(new_pvalues$pvalue) +
  theme_bw(base_size = 24) +
  annotate(
    geom = "text",
    x = -Inf,
    y = Inf,
    hjust = -0.15,
    vjust = 1 + 0.15 * 3,
    label = sprintf("lambda value = %.2f", inflation(new_pvalues$pvalue)),
    size = 8
  ) +
  theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank()
    # panel.grid = element_line(size = 0.5, color = "grey80")
  )
dev.off()

ggplot(data.frame(new_pvalues$pvalue)) +
  geom_histogram(aes(x = new_pvalues$pvalue), bins = 25, color = "white", size = 0.3, boundary = 0.5) +
  theme_minimal(base_size = 20) +
  labs(x = NULL, y = NULL, title = "Histogram of p-values") +
  scale_y_continuous(expand = c(0.02, 0)) +
  theme(
    axis.line.x = element_line(size = 0.5),
    axis.ticks.x = element_line(size = 0.5),
    panel.grid.major.y = element_line(size = 0.5),
    panel.grid.minor.y = element_line(size = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
length(new_pvalues$pvalue[which(new_pvalues$pvalue < 0.05)])

