library(ggplot2)
require(dplyr)
library(scales)

setwd("C:/Users/funky/Dropbox/document/Homework/PublishingStackedBoundaries/Fig4/data/")
setwd("C:/Users/Jojo/Dropbox/document/Homework/PublishingStackedBoundaries/Fig4/data/")

## -------------------------------------------------------
dfR_big = as.data.frame(matrix(nrow = 0, ncol = 7))
df_big =  as.data.frame(matrix(nrow = 0, ncol = 8))
colnames(dfR_big) = c("Assay","limb","RR","lCI","rCI","N","param")
colnames(df_big) = c("Assay","limb","p","pep","N","pstd","pepstd","param")
for (d in 1:6){
  Pfile = dir(pattern=paste0("o3modelPercent_pooled_", d,"Model_v4"))
  Rfile = dir(pattern=paste0("o3modelRR_pooled_", d,"Model_v4"))
  dfR = read.csv(Rfile, sep = '\t', header = FALSE, col.names = c("Assay", "limb", "RR", "lCI", "rCI", "N"))
  
  dfR$Assay[dfR$Assay==1] = "Merge"
  dfR$Assay[dfR$Assay==2] = "Stack"
  dfR$Assay[dfR$Assay==3] = "Neither"
  dfR$Assay = factor(dfR$Assay,levels = c("Neither", "Stack", "Merge"))
  dfR$Assay = factor(dfR$Assay,levels = c("Merge", "Stack", "Neither"))
  dfR$limb[dfR$limb==1] = "h"
  dfR$limb[dfR$limb==2] = "f"
  dfR$param = d
  
  dfR_big = rbind(dfR_big, dfR)
  
  ## -------------------------------------------------------
  # Plotting 3 models' percentage of cells for pooled data
  df = read.csv(Pfile, sep = '\t', header = FALSE, col.names = c("Assay", "limb", "p", "pep", "N"))
  
  df$Assay[df$Assay==1] = "Merge"
  df$Assay[df$Assay==2] = "Stack"
  df$Assay[df$Assay==3] = "Other"
  df$Assay[df$Assay==4] = "All EP-contact"
  df$Assay = factor(df$Assay,levels = c("All EP-contact","Merge", "Stack", "Other"))
  
  df$limb[df$limb==1] = "h"
  df$limb[df$limb==2] = "f"
  
  df$pstd = (df$p*(1-df$p)/df$N)^0.5
  df$pepstd = (df$pep*(1-df$pep)/df$N)^0.5
  df$param = d
  
  df_big = rbind(df_big, df)  
}

ggplot(dfR_big[dfR_big$limb=="h",], aes(x=Assay, y=RR, ymin = lCI, ymax = rCI, label=paste("n =",N))) + 
  geom_pointrange(position = position_dodge(0.5), color = "#aa2029", size=1) + 
  scale_y_log10()+
  geom_hline(yintercept=1, lty=2) + 
  #geom_text(angle=90, position=position_dodge(0.5), hjust=0.5, vjust = 1) + 
  ylab("Risk Ratio") + 
  xlab("") + 
  theme_bw() + 
  facet_wrap(~param, nrow = 1) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 30),
        axis.text.x = element_text(angle = 90),
        strip.background  = element_blank(),
        strip.text.x = element_blank())# +

#coord_flip()
ggsave("20230526_RR.eps", dpi = 300, units = "mm")

dodge = 0.6
ggplot(df_big[df_big$Assay!="All EP-contact",], aes(x=Assay, y=p, ymin = p-pstd, ymax = p+pstd, fill = limb, label=paste("n =",N))) + 
  geom_col(position = position_dodge(dodge), width = 0.5) + 
  geom_errorbar(position = position_dodge(dodge), width = 0.25) + 
  scale_fill_manual(values=c("#005680", "#aa2029"),labels = c("forelimb", "hindlimb")) + 
  ylab("Fraction of all cells") + 
  xlab("") + 
  scale_y_continuous(limits=c(0.1,0.5),oob = rescale_none) + 
  guides(fill=guide_legend("")) +
  theme_bw() + 
  facet_wrap(~param, nrow = 1) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 20),
        legend.position = "none",
        axis.text.x = element_text(angle = 90),
        strip.background  = element_blank(),
        strip.text.x = element_blank())# +)
ggsave("20230526_proportion_mai.eps", dpi = 300)
