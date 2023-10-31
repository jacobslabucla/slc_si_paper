library(ggplot2)
library(dplyr)
library(cowplot)
library(here)
library(tidyr)
library(ggbeeswarm)


## Histo Scores --
histo<- read.csv("/home/julianne/Documents/slcsipaper/Ileum Spontaneous ICP-MS.csv")
histo$Genotype <- factor(histo$Genotype, levels=c("WT", "HET", "MUT"))

ggplot(data=histo,aes(x=Genotype,y=Grade, color=Genotype)) + 
  stat_summary(aes(x=Genotype, y=Grade), fun=median, geom="crossbar", colour="black")+
  geom_beeswarm(cex = 3,priority = "density",size=3)+
  scale_color_viridis_d(option = "D")+
  theme_cowplot(16) +
  theme(legend.position = "none")+
  ggtitle("Histology")+
  ylab("Score")+
  xlab("")+
  ylim(0,2.5)+
  theme(plot.title = element_text(hjust = 0.5))

wilcox.test(Grade~Genotype,data=histo)