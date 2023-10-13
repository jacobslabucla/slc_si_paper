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
### Accompanying statistics ---

## Rotarod -- 
rotarod_bw <- merge(rotarod, bw,by="MouseID")
rotarod_bw$SLC_Genotype.x <- factor(rotarod_bw$SLC_Genotype.x, levels = c("WT", "HET", "MUT"))
lm_day1 <- lm(Average_Latency ~  Weight + SLC_Genotype.x + Sex.x, data = subset(rotarod_bw, Day == "One"))
summary(lm_day1)
lm_day2 <- lm(Average_Latency ~  Weight + SLC_Genotype.x + Sex.x, data = subset(rotarod_bw, Day == "Two"))
summary(lm_day2)
lm_day3 <- lm(Average_Latency ~  Weight + SLC_Genotype.x + Sex.x, data = subset(rotarod_bw, Day == "Three"))
summary(lm_day3)

## OLM -- 
olm_bw <- merge(bw, olm,by="MouseID")
olm_bw$SLC_Genotype.x <- factor(olm_bw$SLC_Genotype.x, levels = c("WT", "HET", "MUT"))
lm_day1 <- lm(Percent_Time_Investigated~  SLC_Genotype.x + Sex.x, data = subset(olm_bw, Day == "Training"))
summary(lm_day1)
lm_day2 <- lm(Percent_Time_Investigated~  SLC_Genotype.x + Sex.x, data = subset(olm_bw, Day == "Testing"))
summary(lm_day2)

## Body Weight -- 
bw$SLC_Genotype <- factor(bw$SLC_Genotype, levels=c("WT","HET", "MUT"))
lm_day1 <- lm(Weight ~ Sex + SLC_Genotype, data = bw)
summary(lm_day1)

## Open Field --
of$SLC_Genotype <- factor(of$SLC_Genotype, levels=c("WT","HET", "MUT"))
lm_day1 <- lm(Center_Time ~ Sex + SLC_Genotype, data = of)
summary(lm_day1)

of_bw <- merge(bw, of,by="MouseID")
of_bw$SLC_Genotype.x <- factor(of_bw$SLC_Genotype.x, levels = c("WT", "HET", "MUT"))
lm_day1 <- lm(Distance ~  Weight+SLC_Genotype.x + Sex.x, data = of_bw)
summary(lm_day1)


