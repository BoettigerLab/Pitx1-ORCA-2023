library(ggplot2)
require(dplyr)
library(scales)
library(tidyr)
library(hash)

setwd("C:/Users/funky/Dropbox/document/Homework/PublishingStackedBoundaries/20230407_thresholds/")
setwd("C:/Users/Jojo/Dropbox/document/Homework/PublishingStackedBoundaries/20230407_thresholds/")
data_folder = "data/"

FH_colors = c("#005680","#aa2029")

#thresholds = c(1,10,25,50,100,150,200,300,400,600,800,1000)
thresholds = c(25,50,100,150,200,300,400,600,800,1000)

# Plot hi-C correlation
hic_corr = as.data.frame(matrix(nrow = length(thresholds), ncol = 3))
colnames(hic_corr) = c("th", "HL_corr", "FL_corr")
hic_corr$th = thresholds
for (th in thresholds){
  th_str = sprintf("%04.f", th)
  df = read.csv(paste0(data_folder, th_str, "_hicCorr.txt"), sep = "\t", header = FALSE)
  hic_corr$HL_corr[hic_corr$th == th] = df[1,1]
  hic_corr$FL_corr[hic_corr$th == th] = df[1,3]
}
hic_corr = pivot_longer(hic_corr, cols = c(2,3), names_to = "limb")
ggplot(hic_corr, aes(x = th, y = value, color = limb)) + 
  geom_line() + 
  geom_point() + 
  scale_x_continuous(name = "threshold (nm)", trans='log2', 
                   breaks = thresholds, minor_breaks = NULL) +
  scale_y_continuous(minor_breaks = NULL) +
  scale_color_manual(values = FH_colors) +
  ylab("correlation with hi-C") + 
  theme_light()
ggsave("20230628_01_corr.eps", dpi = 300, units = "mm")


# Plot the enhancer distances to pitx
loci = c("PelB", "PDE", "RA4", "Pen", "Neurog1")
en_dist = as.data.frame(matrix(nrow = length(thresholds) * length(loci), ncol = 6))
colnames(en_dist) = c("th", "HL_EP", "HL_N", "FL_EP", "FL_N", "locus")
en_dist$locus = rep(loci, length(thresholds))
en_dist$th = rep(thresholds, each = length(loci))
for (th in thresholds){
  th_str = sprintf("%04.f", th)
  df = read.csv(paste0(data_folder, th_str, "_enhancerDiffs.txt"), sep = "\t", header = FALSE)
  en_dist$HL_EP[en_dist$th == th] = df[1,]
  en_dist$HL_N[en_dist$th == th] = df[2,]
  en_dist$FL_EP[en_dist$th == th] = df[3,]
  en_dist$FL_N[en_dist$th == th] = df[4,]
}
en_dist$p = 1
# Perform two-sample z test
for (i in 1:nrow(en_dist)){
  res = prop.test(x = c(en_dist$HL_EP[[i]], en_dist$FL_EP[[i]]), 
                  n = c(en_dist$HL_N[[i]], en_dist$FL_N[[i]]))
  en_dist$p[[i]] = res$p.value
}

#en_dist = pivot_longer(en_dist, cols = c(2,3), names_to = "limb")
en_dist$labs = ""
en_dist$labs[en_dist$p < 0.05] = "*"
#en_dist$value = unlist(en_dist$value)
en_dist$locus = factor(en_dist$locus, levels = loci)
#en_dist$diff = unlist(en_dist$H) - unlist(en_dist$F)
#en_dist$norm_H = unlist(en_dist$HL_EP)/unlist(en_dist$F)
en_dist$norm_H = unlist(en_dist$HL_EP)*unlist(en_dist$FL_N)/(unlist(en_dist$FL_EP) * unlist(en_dist$HL_N))
en_dist$color = "black"
en_dist$color[en_dist$norm_H > 1] = FH_colors[2]
en_dist$color[en_dist$norm_H < 1] = FH_colors[1]

ggplot(subset(en_dist, th != 1), aes(x = th, y = log2(norm_H)))+#, color = limb))+#)) + 
  geom_line() + 
  geom_point(aes(color = color)) + 
  geom_text(aes(label = labs), vjust = -0.2) +
  scale_x_continuous(name = "threshold (nm)",# trans='log2',
                     breaks = thresholds, minor_breaks = NULL) +
  scale_y_continuous(minor_breaks = NULL) +
  scale_color_manual(values = append(FH_colors, "black")) +
  ylab("Pitx1 contact frequency HL/FL") + 
  facet_grid(vars(locus)) +
  theme_light()
ggsave("20230713_02_FHenhancers.eps", dpi = 300, units = "mm")



# Plot proportion of HL and FL for 3 models.
loci = c("PelB", "PDE", "RA4", "Pen", "Neurog1")
model_dict <- hash(1:4, c("Merge", "Stack", "Other", "Any"))
limb_dict = hash(1:2, c("HL", "FL"))

Propor = as.data.frame(matrix(ncol = 6))
colnames(Propor) = c("model", "limb", "proportion", "EP_proportion", "N", "th")
for (th in thresholds){
  th_str = sprintf("%04.f", th)
  df = read.csv(paste0(data_folder, th_str, "_3modelPercent.txt"), sep = "\t", header = FALSE)
  df$th = th
  colnames(df) = colnames(Propor)
  Propor = rbind(Propor, df)
}
Propor = Propor[2:nrow(Propor),]
Propor$model = sapply(Propor$model, function(x){return(model_dict[[as.character(x)]])})
Propor$limb = sapply(Propor$limb, function(x){return(limb_dict[[as.character(x)]])})

ggplot(subset(Propor, model != "Any"), aes(x = th, y = proportion, color = limb)) + 
  geom_line() + 
  geom_point() + 
  scale_x_continuous(name = "threshold (nm)",trans='log2',
                     breaks = thresholds, minor_breaks = NULL) +
  scale_y_continuous(minor_breaks = NULL) +
  scale_color_manual(values = FH_colors) +
  ylab("Proportion of traces fitting a model") + 
  facet_grid(vars(model)) +
  theme_light()
ggsave("20230628_03_3modelProportion.eps", dpi = 300, units = "mm")

# Need more work on this plot: other ways of representations?
ggplot(subset(Propor, model != "Any"), aes(x = th, y = EP_proportion, color = limb)) + 
  geom_line() + 
  geom_point() + 
  scale_x_continuous(name = "threshold (nm)",trans='log2',
                     breaks = thresholds, minor_breaks = NULL) +
  scale_y_continuous(minor_breaks = NULL) +
  scale_color_manual(values = FH_colors) +
  ylab("Proportion of traces with EP_contact fitting a model") + 
  facet_grid(vars(model)) +
  theme_light()
ggsave("20230628_03_3modelProportionEP.eps", dpi = 300, units = "mm")

# Plot RR.
loci = c("PelB", "PDE", "RA4", "Pen", "Neurog1")
model_dict <- hash(1:3, c("Merge", "Stack", "Other"))
limb_dict = hash(1:2, c("HL", "FL"))

RR_df = as.data.frame(matrix(ncol = 7))
colnames(RR_df) = c("model", "limb", "RR", "CI_L", "CI_R", "N", "th")
for (th in thresholds){
  th_str = sprintf("%04.f", th)
  df = read.csv(paste0(data_folder, th_str, "_RR.txt"), sep = "\t", header = FALSE)
  df$th = th
  colnames(df) = colnames(RR_df)
  RR_df = rbind(RR_df, df)
}
RR_df = RR_df[2:nrow(RR_df),]
RR_df$model = sapply(RR_df$model, function(x){return(model_dict[[as.character(x)]])})
RR_df$limb = sapply(RR_df$limb, function(x){return(limb_dict[[as.character(x)]])})

ggplot(subset(RR_df, (limb == "HL") & (th > 1)), aes(x = th, y = RR, ymin = CI_L, ymax = CI_R)) + 
  geom_line() + 
  geom_pointrange() + 
  scale_x_continuous(name = "threshold (nm)", trans='log2',
                     breaks = thresholds, minor_breaks = NULL) +
  scale_y_continuous(minor_breaks = NULL, trans = 'log2') +
  geom_hline(yintercept = 1) + 
  #scale_color_manual(values = FH_colors) +
  ylab("Risk ratio") + 
  coord_cartesian(ylim=c(1/3, 3)) +
  facet_grid(vars(model)) +
  theme_light()
  
  #geom_vline(xintercept = 200)
  
  # vertical line for 200.
  # color scheme for HL FL
  # white background.

ggsave("20230710_05_RR.eps", dpi = 300, units = "mm")



# Plot the enhancer distances to pitx
loci = c("PelB", "PDE", "RA4", "Pen", "Neurog1")
en_dist = as.data.frame(matrix(nrow = length(thresholds) * length(loci), ncol = 5))
colnames(en_dist) = c("th", "H", "F", "RS_p", "locus")
en_dist$locus = rep(loci, length(thresholds))
en_dist$th = rep(thresholds, each = length(loci))
for (th in thresholds){
  th_str = sprintf("%04.f", th)
  df = read.csv(paste0(data_folder, th_str, "_enhancerDiffs.txt"), sep = "\t", header = FALSE)
  en_dist$H[en_dist$th == th] = df[1,]
  en_dist$F[en_dist$th == th] = df[2,]
  en_dist$RS_p[en_dist$th == th] = df[3,]
}
#en_dist = pivot_longer(en_dist, cols = c(2,3), names_to = "limb")
en_dist$labs = ""
en_dist$labs[en_dist$RS_p < 0.01] = "*"
#en_dist$value = unlist(en_dist$value)
en_dist$locus = factor(en_dist$locus, levels = loci)
#en_dist$diff = unlist(en_dist$H) - unlist(en_dist$F)
en_dist$norm_H = unlist(en_dist$H)/unlist(en_dist$F)
en_dist$color = "black"
en_dist$color[en_dist$norm_H > 1] = FH_colors[2]
en_dist$color[en_dist$norm_H < 1] = FH_colors[1]

ggplot(subset(en_dist, th != 1), aes(x = th, y = log2(norm_H)))+#, color = limb))+#)) + 
  geom_line() + 
  geom_point(aes(color = color)) + 
  geom_text(aes(label = labs), vjust = -0.2) +
  scale_x_continuous(name = "threshold (nm)",# trans='log2',
                     breaks = thresholds, minor_breaks = NULL) +
  scale_y_continuous(minor_breaks = NULL) +
  scale_color_manual(values = append(FH_colors, "black")) +
  ylab("Pitx1 contact frequency HL/FL") + 
  facet_grid(vars(locus)) +
  theme_light()
ggsave("20230628_02_FHenhancers.eps", dpi = 300, units = "mm")