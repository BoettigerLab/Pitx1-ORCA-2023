library(ggplot2)
require(dplyr)
library(scales)

setwd("C:/Users/funky/Dropbox/document/Homework/PublishingStackedBoundaries/20230418_E12E11/data/")
setwd("C:/Users/Jojo/Dropbox/document/Homework/PublishingStackedBoundaries/20230418_E12E11/data/")

## -------------------------------------------------------
# Plotting 3 models' RR for pooled data
#df = read.csv("n3modelRR_pooled.txt", sep = '\t', header = FALSE, col.names = c("Assay", "limb", "RR", "lCI", "rCI", "N"))
#df = read.csv("o3modelRR_pooled_sep1.txt", sep = '\t', header = FALSE, col.names = c("Assay", "limb", "RR", "lCI", "rCI", "N"))
#df = read.csv("sim3modelRR_pooled.txt", sep = '\t', header = FALSE, col.names = c("Assay", "limb", "RR", "lCI", "rCI", "N"))
#df = read.csv("sim3modelRR_pooled_soph.txt", sep = '\t', header = FALSE, col.names = c("Assay", "limb", "RR", "lCI", "rCI", "N"))
#df = read.csv("o3modelRR_pooled_20230526_E11_szMaxtriu.txt", sep = '\t', header = FALSE, col.names = c("Assay", "limb", "RR", "lCI", "rCI", "N"))
df = read.csv("o3modelRR_pooled_20230711.txt", sep = '\t', header = FALSE, col.names = c("Assay", "limb", "RR", "lCI", "rCI", "N"))

df$Assay[df$Assay==1] = "Merge"
df$Assay[df$Assay==2] = "Stack"
df$Assay[df$Assay==3] = "Neither"
df$Assay = factor(df$Assay,levels = c("Neither", "Stack", "Merge"))
df$Assay = factor(df$Assay,levels = c("Merge", "Stack", "Neither"))
df$limb[df$limb==1] = "h"
df$limb[df$limb==2] = "f"

ggplot(df[df$limb=="h",], aes(x=Assay, y=RR, ymin = lCI, ymax = rCI, label=paste("n =",N))) + 
  geom_pointrange(position = position_dodge(0.5), color = "#aa2029", size=1.5) + 
  scale_y_log10(breaks = c(0.6, 1, 1.5))+
  geom_hline(yintercept=1, lty=2) + 
  #geom_text(angle=90, position=position_dodge(0.5), hjust=0.5, vjust = 1) + 
  ylab("Risk Ratio") + 
  xlab("") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 30))# +
  #coord_flip()
ggsave("20230711_RR.eps", dpi = 300, units = "mm")
## -------------------------------------------------------
# Plotting 3 models' percentage of cells for pooled data
#df = read.csv("n3modelPercent_pooled.txt", sep = '\t', header = FALSE, col.names = c("Assay", "limb", "p", "pep", "N"))
#df = read.csv("o3modelPercent_pooled_sep1.txt", sep = '\t', header = FALSE, col.names = c("Assay", "limb", "p", "pep", "N"))
#df = read.csv("sim3modelPercent_pooled.txt", sep = '\t', header = FALSE, col.names = c("Assay", "limb", "p", "pep", "N"))
#df = read.csv("sim3modelPercent_pooled_soph.txt", sep = '\t', header = FALSE, col.names = c("Assay", "limb", "p", "pep", "N"))
#df = read.csv("o3modelPercent_pooled_20230526_E11_szMaxtriu.txt", sep = '\t', header = FALSE, col.names = c("Assay", "limb", "p", "pep", "N"))
df = read.csv("o3modelPercent_pooled_20230711.txt", sep = '\t', header = FALSE, col.names = c("Assay", "limb", "p", "pep", "N"))

df$Assay[df$Assay==1] = "Merge"
df$Assay[df$Assay==2] = "Stack"
df$Assay[df$Assay==3] = "Neither"
df$Assay[df$Assay==4] = "All EP-contact"
df$Assay = factor(df$Assay,levels = c("All EP-contact","Merge", "Stack", "Neither"))

df$limb[df$limb==1] = "h"
df$limb[df$limb==2] = "f"

df$pstd = (df$p*(1-df$p)/df$N)^0.5
df$pepstd = (df$pep*(1-df$pep)/df$N)^0.5
dodge = 0.6
ggplot(df[df$Assay!="All EP-contact",], aes(x=Assay, y=p, ymin = p-pstd, ymax = p+pstd, fill = limb, label=paste("n =",N))) + 
  geom_col(position = position_dodge(dodge), width = 0.5) + 
  geom_errorbar(position = position_dodge(dodge), width = 0.25) + 
  scale_fill_manual(values=c("#005680", "#aa2029"),labels = c("forelimb", "hindlimb")) + 
  ylab("Fraction of all cells") + 
  xlab("") + 
  scale_y_continuous(limits=c(0.25,0.45),oob = rescale_none) + 
  guides(fill=guide_legend("")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 20),
        legend.position = "none")
ggsave("20230526_proportion.eps", dpi = 300)

ggplot(df[df$Assay!="All EP-contact",], aes(x=Assay, y=pep, ymin = pep-pepstd, ymax = pep+pepstd, fill = limb, label=paste("n =",N))) + 
  geom_col(position = position_dodge(dodge), width = 0.5) + 
  geom_errorbar(position = position_dodge(dodge), width = 0.25) + 
  scale_fill_manual(values=c("#005680", "#aa2029"), labels = c("FL", "HL")) + 
  ylab("Fraction exhibiting \nconformation & E-P contact") + 
  xlab("") + 
  theme_bw() + 
  #ylim(0,0.4) + 
  guides(fill=guide_legend("")) +
  geom_hline(yintercept=df$pep[which(df$Assay=="All EP-contact" & df$limb=="f")], lty=2, color = "#005680") + 
  geom_hline(yintercept=df$pep[which(df$Assay=="All EP-contact" & df$limb=="h")], lty=2,  color = "#aa2029") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size = 30),
        axis.title = element_text(size = 25))

ggsave("20230711_EPproportion1.eps", dpi = 300)


df$pepPor = df$pep/df$pep[7]
df$pepPor[df$limb=="f"] = df$pep[df$limb=="f"]/df$pep[8]
ggplot(df[df$Assay!="All EP-contact",], aes(x=limb, y=pep, ymin = pep-pepstd, ymax = pep+pepstd, fill = Assay, label=paste("n =",N))) + 
  geom_bar(position = 'stack', stat='identity') + 
  #geom_errorbar(position = position_dodge(dodge), width = 0.25) + 
  scale_fill_manual(values=c("#92CB6D", "#ec1e24","#999999"), labels = c("Merge", "Stack","Other")) + 
  ylab("Fraction exhibiting \nconformation & E-P contact") + 
  xlab("") + 
  theme_bw() + 
  #ylim(0,0.4) + 
  #guides(fill=guide_legend("")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size = 30),
        axis.title = element_text(size = 25))

ggsave("20230711_EPproportion2.eps", dpi = 300)


ggplot(df[df$Assay!="All EP-contact",], aes(x=limb, y=pepPor, fill = Assay, label=paste("n =",N))) + 
  geom_bar(position = 'stack', stat='identity') + 
  #geom_errorbar(position = position_dodge(dodge), width = 0.25) + 
  scale_fill_manual(values=c("#92CB6D", "#ec1e24","#999999"), labels = c("Merge", "Stack","Other")) + 
  ylab("Fraction exhibiting \nconformation & E-P contact") + 
  xlab("") + 
  theme_bw() + 
  #ylim(0,0.4) + 
  #guides(fill=guide_legend("")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size = 30),
        axis.title = element_text(size = 25))

ggsave("20230711_EPproportion3.eps", dpi = 300)

for (l in c("h", "f")){
  df$pep_percent[df$limb == l] = df$pep[df$limb == l]/df$pep[df$limb == l & df$Assay == "All EP-contact"]
  df$pepstd_percent[df$limb == l] = df$pepstd[df$limb == l]/df$pep[df$limb == l & df$Assay == "All EP-contact"]
}


ggplot(df[df$Assay!="All EP-contact",], aes(x=Assay, y=pep_percent, ymin = pep_percent-pepstd_percent, ymax = pep_percent+pepstd_percent, fill = limb, label=paste("n =",N))) + 
  geom_col(position = position_dodge(dodge), width = 0.5) + 
  geom_errorbar(position = position_dodge(dodge), width = 0.25) + 
  scale_fill_manual(values=c("#005680", "#aa2029"), labels = c("FL", "HL")) + 
  ylab("Fraction of traces exhibiting \nconformation & E-P contact") + 
  xlab("") + 
  theme_bw() + 
  #ylim(0,0.4) + 
  guides(fill=guide_legend("")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size = 30),
        axis.title = element_text(size = 25))

ggsave("20230510_EPproportion_perc.eps", dpi = 300)
