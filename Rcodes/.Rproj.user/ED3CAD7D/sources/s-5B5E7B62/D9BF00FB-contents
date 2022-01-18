# Created by Hyungseok Kim (hskimm@mit.edu), last modified on 1/18/22
# Continued from create_df.R
library("ggplot2")
library(tidyverse)

# subset df by each treatment
trt = '+1R1'
df_trt <- df[df$Treatment==trt,]
df_trt$"Voltage" <- as.factor(df_trt$"Voltage")
df_trt$"Replicate" <- as.factor(df_trt$"Replicate")

# create a boxplot, grouped by voltage
ggplot(df_trt, aes(x=Voltage, y=Velocity, colour=Replicate)) + 
  geom_boxplot() +
  # geom_point(position=position_jitterdodge())
  labs(title=trt, x="Voltage (V)", y = "Velocity (µm/s)")

setwd(paste(path,"figs",sep="/"))
ggsave(paste(trt,".eps", sep=""), width = 5, height = 4, units = "in")

# stat test
# i = 3; j = 4;
# g1 = df$Value[df[,1]==trt[i]]
# g2 = df$Value[df[,1]==trt[j]]
# res <- t.test(g1, g2, alternative='two.sided', paired=FALSE)
# sprintf("t test %s vs %s",trt[i],trt[j])
# print(res)