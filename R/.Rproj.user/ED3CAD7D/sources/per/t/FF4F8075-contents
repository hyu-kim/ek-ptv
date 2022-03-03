# Created by Hyungseok Kim (hskimm@mit.edu), last modified on 2/18/22
# Continued from create_df.R
library("ggplot2")
library(tidyverse)

# subset df by each treatment
trt = 'ARW1R1'
df_p <- df[(df$Treatment==trt),]
df_p$"Voltage" <- as.factor(df_p$"Voltage")
df_p$"Replicate" <- as.factor(df_p$"Replicate")

# create a boxplot, grouped by voltage
ggplot(df_p, aes(x=Voltage, y=Velocity, colour=Channel)) + 
  geom_boxplot() +
  # geom_point(position=position_jitterdodge())
  labs(title=paste("Bacterial strain",trt,sep=" "), x="Voltage (V)", y = "Velocity (µm/s)")

# setwd(paste(path,"figs",sep="/"))
ggsave(paste(trt,".eps", sep=""), width = 5, height = 4, units = "in")