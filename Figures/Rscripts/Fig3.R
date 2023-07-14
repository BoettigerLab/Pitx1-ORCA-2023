#Plot distance to center plot in Fig 3
library(ggplot2)

setwd("C:/Users/funky/Dropbox/document/Homework/PublishingStackedBoundaries/Fig3")
df = read.csv("Fig3_CenterDist.txt", header = FALSE)
df = read.csv("SuppFig7_CenterDist.txt", header = FALSE)

groups = c("stack", "non-stack", "all")
newdf = df[,1:3]
colnames(newdf) = c("lower","median","upper")
newdf$group = "stack"
newdf$pos = 1:75

for (i in 2:3){
  df2 = df[,(i*3-2):(i*3)]
  colnames(df2) = c("lower","median","upper")
  df2$group = groups[i]
  df2$pos = 1:75
  newdf = rbind(newdf, df2)
}
FH_colors = c("#999999","#000000","#aa2029")
fill_colors = c("#dddddd","#444444","#ddaabb")
ggplot(newdf, aes(x = pos, y = median, color = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = group, color = NA)) + 
  geom_line() + 
  scale_x_continuous(limits = c(1,75)) +
  scale_color_manual(values = FH_colors) +
  scale_fill_manual(values = fill_colors) +
  theme_bw() 

# PB2 stuf ===================================================================
df = read.csv("Fig3_PB2.txt", header = FALSE)
df = read.csv("Fig3_th100_PB2.txt", header = FALSE)
df$pos = 1:75
df1 = df[,c(1,3)]
df2 = df[,c(2,3)]
colnames(df1) = c("median", "pos")
df1$group = "PB2"
df1$N = 2981
#df1$N = 570 #for th100

colnames(df2) = c("median", "pos")
df2$group = "all"
df2$N = 15057

newdf = rbind(df1, df2)
newdf$stdev = ((1 - newdf$median)*newdf$median/newdf$N)^0.5

FH_colors = c("#000000","#ec1e24")
fill_colors = c("#999999","#f2766d")
ggplot(newdf, aes(x = pos, y = median, color = group)) +
  geom_ribbon(aes(ymin = median-stdev, ymax = median+stdev, fill = group, color = NA)) + 
  geom_line() + 
  scale_x_continuous(limits = c(1,75)) +
  scale_color_manual(values = FH_colors) +
  scale_fill_manual(values = fill_colors) +
  theme_bw() 
