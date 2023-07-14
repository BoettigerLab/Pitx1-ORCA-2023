#Plot distance to center plot in Fig 3
library(ggplot2)

setwd("C:/Users/funky/Dropbox/document/Homework/PublishingStackedBoundaries/FigSupp7/")
df = read.csv("SuppFig7_CenterDist.txt", header = FALSE)

groups = c("high","low")
newdf = df[,1:3]
colnames(newdf) = c("lower","median","upper")
newdf$group = "high"
newdf$pos = 1:200

for (i in 2){
  df2 = df[,(i*3-2):(i*3)]
  colnames(df2) = c("lower","median","upper")
  df2$group = groups[i]
  df2$pos = 1:200
  newdf = rbind(newdf, df2)
}
FH_colors = c("#ec1e24","#204497")
fill_colors = c("#f2766d","#8989c2")
ggplot(newdf, aes(x = pos, y = median, color = group)) +
  geom_linerange(aes(ymin = lower, ymax = upper, color = group)) + 
  geom_line() + 
  scale_x_continuous(limits = c(1,200)) +
  scale_color_manual(values = FH_colors) +
  scale_fill_manual(values = fill_colors) +
  theme_bw() 


# MultiContact ===================================================================
df = read.csv("SuppFig7_MultiContact.txt", header = FALSE)
groups = c("high","low")
newdf = df[,1:3]
colnames(newdf) = c("lower","median","upper")
newdf$group = "high"
newdf$pos = 1:200

for (i in 2){
  df2 = df[,(i*3-2):(i*3)]
  colnames(df2) = c("lower","median","upper")
  df2$group = groups[i]
  df2$pos = 1:200
  newdf = rbind(newdf, df2)
}
FH_colors = c("#ec1e24","#204497")
fill_colors = c("#f2766d","#8989c2")
ggplot(newdf, aes(x = pos, y = median, color = group)) +
  #geom_ribbon(aes(ymin = lower, ymax = upper, fill = group, color = NA)) + 
  geom_line() + 
  scale_x_continuous(limits = c(1,200)) +
  scale_color_manual(values = FH_colors) +
  scale_fill_manual(values = fill_colors) +
  theme_bw() 
