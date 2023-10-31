###### SLC ICP-MS Data Aggregated ---
### Date: 10.26.2022
### Contents: Includes data from 10/5/21 and 10/25/22 Core submissions
###### End ---

library(here)
library(ggplot2)
library(rlang)
library(cowplot)
library(ggpubr)
library(dplyr)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/SLC Spontaneous/ICPMS/")
here::i_am("ICP_MS_Analysis.R")

### Data Preprocessing ---
df<- read.csv("Analysis_ICP_MS.csv", header=TRUE, row.names=1)

# replace all n/a and declare all element columns as numerical
df[df=="n/a"]<-0
vector <- names(df)
elements <- vector[1:7]
df <- df %>% mutate_at(c(elements), as.numeric)
str(df)
df$Genotype_Batch <- paste0(df$Genotype, "_",df$Batch)
df$Genotype_Sex <- paste0(df$Genotype,"_",df$Sex)

# Subset by SampleType - with outliers
df_fp_col <- df %>% filter(SampleType=="FP-COL")
df_fp_si <- df %>% filter(SampleType=="FP-SI")
df_muc_col <- df %>% filter(SampleType=="MUC-COL")
df_muc_si <- df %>% filter(SampleType=="MUC-SI")
df_ts_col <- df %>% filter(SampleType=="TS-COL")
df_ts_si <- df %>% filter(SampleType=="TS-SI")


### Figures ---
generate_violin_plots <- function (input_data, column_index, X) {
  # read in file
  data<-as.data.frame(input_data)
  
  #Ensure correct ordering of levels 
  data$Genotype <- factor(data$Genotype, levels = c("WT","MUT"))
  data$Genotype_Batch <- factor(data$Genotype_Batch, levels=c("WT_One","WT_Two", "WT_Three", "MUT_One", "MUT_Two", "MUT_Three"))
  data$Genotype_Sex <- factor(data$Genotype_Sex, levels=c("WT_Male","WT_Female", "MUT_Male", "MUT_Female"))
  
  #Get correct y-axis label
  ylabel <- c(elements[column_index])
  
  ggplot(data=data,aes(x={{X}},y=data[, column_index], fill=Genotype)) + 
    geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+
    #scale_shape_manual(values=c(16,10))+
    scale_fill_viridis_d()+
    geom_point(size=1,position=position_jitter(width=0.25),alpha=0.8)+
    theme_cowplot(16) +
    theme(legend.position = "none") +
    labs(x="",y=ylabel)
    #ylim(min,max) +
    
}


element_plots <- list()
fp_col_plots <- list()
fp_si_plots <- list()
muc_col_plots <- list()
muc_si_plots <- list()
ts_col_plots <- list()
ts_si_plots <- list()
compare_vector<- c("WT","MUT")

# Loop through all elements - Genotype as X variable
for (int in 1:7){
    print(int)

    fp_col <- generate_violin_plots(df_fp_col, int, Genotype) +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle("FP Col")+
      stat_compare_means(comparisons = compare_vector,
                         method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
    fp_si <- generate_violin_plots(df_fp_si, int, Genotype)+
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle("FP SI")+
      stat_compare_means(comparisons = compare_vector,
                         method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
    muc_si <- generate_violin_plots(df_muc_si, int, Genotype) +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle("MUC SI")+
      stat_compare_means(comparisons = compare_vector,
                         method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
    muc_col <- generate_violin_plots(df_muc_col, int, Genotype)+
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle("MUC Col")+
      stat_compare_means(comparisons = compare_vector,
                         method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
    ts_si <- generate_violin_plots(df_ts_si, int, Genotype)+
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle("TS SI")+
      stat_compare_means(comparisons = compare_vector,
                         method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
    ts_col <- generate_violin_plots(df_ts_col, int, Genotype)+
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle("TS Col")+
      stat_compare_means(comparisons = compare_vector,
                         method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

    element_plots[[int]] <- cowplot::plot_grid(fp_col, fp_si,muc_col,muc_si, ts_col,ts_si, 
                                         rows = 3,cols=2)
    fp_col_plots[[int]] <- fp_col
    fp_si_plots[[int]] <- fp_si
    muc_col_plots[[int]] <- muc_col
    muc_si_plots[[int]] <- muc_si
    ts_col_plots[[int]] <- ts_col
    ts_si_plots[[int]] <- ts_si
}

# Loop through all elements - Genotype_Batch as X variable
element_plots <- list()
fp_col_plots <- list()
fp_si_plots <- list()
muc_col_plots <- list()
muc_si_plots <- list()
ts_col_plots <- list()
ts_si_plots <- list()
compare_vector<- c("WT","MUT")

for (int in 1:7){
  print(int)
  
  fp_col <- generate_violin_plots(df_fp_col, int, Genotype_Batch) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("FP Col")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  fp_si <- generate_violin_plots(df_fp_si, int, Genotype_Batch)+
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("FP SI")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  muc_si <- generate_violin_plots(df_muc_si, int, Genotype_Batch) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("MUC SI")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  muc_col <- generate_violin_plots(df_muc_col, int, Genotype_Batch)+
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("MUC Col")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ts_si <- generate_violin_plots(df_ts_si, int, Genotype_Batch)+
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("TS SI")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ts_col <- generate_violin_plots(df_ts_col, int, Genotype_Batch)+
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("TS Col")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  element_plots[[int]] <- cowplot::plot_grid(fp_col, fp_si,muc_col,muc_si, ts_col,ts_si, 
                                             rows = 3,cols=2)
  fp_col_plots[[int]] <- fp_col
  fp_si_plots[[int]] <- fp_si
  muc_col_plots[[int]] <- muc_col
  muc_si_plots[[int]] <- muc_si
  ts_col_plots[[int]] <- ts_col
  ts_si_plots[[int]] <- ts_si
}

# Loop through all elements - Genotype_Sex as X variable
element_plots <- list()
fp_col_plots <- list()
fp_si_plots <- list()
muc_col_plots <- list()
muc_si_plots <- list()
ts_col_plots <- list()
ts_si_plots <- list()

for (int in 1:7){
  print(int)
  
  fp_col <- generate_violin_plots(df_fp_col, int, Genotype_Sex) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("FP Col")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  fp_si <- generate_violin_plots(df_fp_si, int, Genotype_Sex)+
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("FP SI")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  muc_si <- generate_violin_plots(df_muc_si, int, Genotype_Sex) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("MUC SI")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  muc_col <- generate_violin_plots(df_muc_col, int, Genotype_Sex)+
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("MUC Col")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ts_si <- generate_violin_plots(df_ts_si, int, Genotype_Sex)+
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("TS SI")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ts_col <- generate_violin_plots(df_ts_col, int, Genotype_Sex)+
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("TS Col")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  element_plots[[int]] <- cowplot::plot_grid(fp_col, fp_si,muc_col,muc_si, ts_col,ts_si, 
                                             rows = 3,cols=2)
  fp_col_plots[[int]] <- fp_col
  fp_si_plots[[int]] <- fp_si
  muc_col_plots[[int]] <- muc_col
  muc_si_plots[[int]] <- muc_si
  ts_col_plots[[int]] <- ts_col
  ts_si_plots[[int]] <- ts_si
}



#FP COL
dev.new(width=15, height=10)
plot_grid(fp_col_plots[[1]],fp_col_plots[[2]],
          fp_col_plots[[3]],fp_col_plots[[4]],
          fp_col_plots[[5]],fp_col_plots[[6]],
          fp_col_plots[[7]],nrow = 2, ncol=4)

#FP SI
dev.new(width=15, height=10)
plot_grid(fp_si_plots[[1]],fp_si_plots[[2]],
          fp_si_plots[[3]],fp_si_plots[[4]],
          fp_si_plots[[5]],fp_si_plots[[6]],
          fp_si_plots[[7]],nrow = 2, ncol=4)

#MUC COL
dev.new(width=15, height=10)
plot_grid(muc_col_plots[[1]],muc_col_plots[[2]],
          muc_col_plots[[3]],muc_col_plots[[4]],
          muc_col_plots[[5]],muc_col_plots[[6]],
          muc_col_plots[[7]],nrow = 2, ncol=4)

#MUC SI
dev.new(width=15, height=10)
plot_grid(muc_si_plots[[1]],muc_si_plots[[2]],
          muc_si_plots[[3]],muc_si_plots[[4]],
          muc_si_plots[[5]],muc_si_plots[[6]],
          muc_si_plots[[7]],nrow = 2, ncol=4)

#TS COL
dev.new(width=15, height=10)
plot_grid(ts_col_plots[[1]],ts_col_plots[[2]],
          ts_col_plots[[3]],ts_col_plots[[4]],
          ts_col_plots[[5]],ts_col_plots[[6]],
          ts_col_plots[[7]],nrow = 2, ncol=4)

#TS SI
dev.new(width=15, height=10)
plot_grid(ts_si_plots[[1]],ts_si_plots[[2]],
          ts_si_plots[[3]],ts_si_plots[[4]],
          ts_si_plots[[5]],ts_si_plots[[6]],
          ts_si_plots[[7]],nrow = 2, ncol=4)

### Statistics --- 
element_stats_para <- list()
element_stats_nonpara <-list()


for (int in 1:7){
  print(int)
ts_si_para <- t.test(df_ts_si[,int]~Genotype,df_ts_si)
ts_si_nonpara <- wilcox.test(df_ts_si[,int]~Genotype,df_ts_si)
ts_col_para <- t.test(df_ts_col[,int]~Genotype,df_ts_col)
muc_si_para <- t.test(df_muc_si[,int]~Genotype,df_muc_si)
muc_si_nonpara <- wilcox.test(df_muc_si[,int]~Genotype,df_muc_si)
ts_col_nonpara <- wilcox.test(df_ts_col[,int]~Genotype,df_ts_col)
muc_col_para <- t.test(df_muc_col[,int]~Genotype,df_muc_col)
muc_col_nonpara <- wilcox.test(df_muc_col[,int]~Genotype,df_muc_col)
fp_si_para <- t.test(df_fp_si[,int]~Genotype,df_fp_si)
fp_si_nonpara <- wilcox.test(df_fp_si[,int]~Genotype,df_fp_si)
fp_col_para <- t.test(df_fp_col[,int]~Genotype,df_fp_col)
fp_col_nonpara <- wilcox.test(df_fp_col[,int]~Genotype,df_fp_col)

element_stats_para[[int]] <-list(print(fp_col_para), print(fp_si_para),print(muc_col_para),print(muc_si_para), print(ts_col_para),print(ts_si_para))
element_stats_nonpara[[int]] <-list(print(fp_col_nonpara), print(fp_si_nonpara),print(muc_col_nonpara),print(muc_si_nonpara), print(ts_col_nonpara),print(ts_si_nonpara))

}

names(df_ts_col)
# Iron 
element_stats_para[[1]]
element_stats_nonpara[[1]]

# Cobalt 
element_stats_para[[2]]
element_stats_nonpara[[2]]

# Copper
element_stats_para[[3]]
element_stats_nonpara[[3]]

# Zinc
element_stats_para[[4]]
element_stats_nonpara[[4]]

# Cadmium
element_stats_para[[5]]
element_stats_nonpara[[5]]

# Manganese
element_stats_para[[6]]
element_stats_nonpara[[6]]

# Selenium
element_stats_para[[7]]
element_stats_nonpara[[7]]