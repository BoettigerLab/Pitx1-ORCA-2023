library(ggplot2)


setwd("C:/Users/funky/Dropbox/document/Homework/PublishingStackedBoundaries/Fig1f3d/")
setwd("C:/Users/Jojo/Dropbox/document/Homework/PublishingStackedBoundaries/Fig1f3d/")

# Fig 1f ===================================================================
df = read.csv("Fig1f.txt", header = FALSE)
output = as.data.frame(matrix(nrow = nrow(df) * 2, ncol = 4))
colnames(output) = c("N1","N2","enhancer","limb")
df$enhancer = c("PelB", "RA3", "RA4", "Pen", "Neurog1")
output[1:5,1:2] = df[1:5, 1:2]
output[6:10,1:2] = df[1:5, 3:4]
output$enhancer = df$enhancer
output$enhancer = factor(output$enhancer, levels = df$enhancer)
output$limb[1:5] = "hindlimb"
output$limb[6:10] = "forelimb"
output$proportion = output$N1/output$N2
output$stdev = ((output$proportion)*(1 - output$proportion)/output$N2)^0.5

for (i in 1:5){
  result = prop.test(x = c(df$V1[i], df$V3[i]), n = c(df$V2[i], df$V4[i]), alternative = "two.sided")
  print(result$p.value)
}

ggplot(data = output, aes(x = enhancer, y = proportion, color = limb)) +
  geom_linerange(aes(ymin = proportion - stdev, ymax = proportion + stdev)) + 
  theme_bw() + 
  scale_color_manual(values=c("#005680", "#aa2029"),labels = c("forelimb", "hindlimb"))
ggsave("fig1f.eps", dpi = 300, units = "mm")

# Fig 3e ===================================================================
df1 = read.csv("Fig3d_20230616_pop.txt", sep = ",", header = FALSE)
df2 = read.csv("Fig3d_20230616_PB2.txt", sep = ",", header = FALSE)

output1 = as.data.frame(matrix(nrow = ncol(df1), ncol = 4))
output2 = as.data.frame(matrix(nrow = ncol(df2), ncol = 4))
colnames(output1) = c("mean", "stdev", "group", "distance")
colnames(output2) = c("mean", "stdev", "group", "distance")
output1[,1] = colMeans(df1)
output1[,2] = sapply(df1, sd)
output2[,1] = colMeans(df2)
output2[,2] = sapply(df2, sd)
output1$group = "pop"
output2$group = "PB2"
output1$distance = 0:20
output2$distance = 0:20

output = rbind(output1, output2)


ggplot(data = output, aes(x = distance, y = mean, color = group)) + 
  geom_ribbon(aes(ymin = mean-stdev, ymax = mean + stdev, fill = group)) + 
  geom_line() + 
  #ylim(0,3) + 
  theme_bw()

ggsave("20230614_fig3d.eps", dpi = 300, units = "mm")

# Fig 3f ===================================================================
df1 = read.csv("Fig3f_hasData.txt", sep = ",", header = FALSE)
df2 = read.csv("Fig3f_isStack.txt", sep = ",", header = FALSE)
colnames(df1) = c("mean", "stdev", "lb", "ub")
colnames(df2) = c("mean", "stdev", "lb", "ub")
df1$pop = "whole"
df2$pop = "stack"
df1$positions = 1:75
df2$positions = 1:75
df = rbind(df1, df2)

ggplot(data = df, aes(x = positions, y = mean, color = pop)) + 
  geom_ribbon(aes(ymin = lb, ymax = ub, fill = pop)) +
  geom_line() + 
  theme_bw()
