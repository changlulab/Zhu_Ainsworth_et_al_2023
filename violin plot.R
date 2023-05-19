library(ggplot2)

rm(list = ls())
options(stringsAsFactors = F)

QC_K27 = read.csv('RNA_QC.csv')

colnames(QC_K27)
pdf("RNA_deduplication.rate.pdf",height = 5,width = 8)
ggplot(QC_K27,aes(x=condition,y=deduplication.rate,fill=condition)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width=0.2,fill="white") +
  scale_color_brewer(palette="Dark2") +
  labs(title="Plot of total reads by Groups",x="Groups", y = "number of total reads (million)") +
  coord_cartesian(ylim = c(40,80)) +
  theme_minimal()

dev.off()

