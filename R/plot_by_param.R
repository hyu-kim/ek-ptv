# Created by Hyungseok Kim (hskimm@mit.edu), last run 2/11/22
# Continued from create_df.R
library("ggplot2")
library(tidyverse)

# subset df by each treatment
trt = 'Ax'
ph = phs[6]
df_p <- df[(df$Treatment==trt)&(df$pH==ph),]
df_p$"Voltage" <- as.factor(df_p$"Voltage")
df_p$"Replicate" <- as.factor(df_p$"Replicate")

# create a boxplot, grouped by voltage
ggplot(df_p, aes(x=Voltage, y=Velocity, colour=Replicate)) + 
  geom_boxplot() +
  # geom_point(position=position_jitterdodge())
  labs(title=paste("pH",ph,sep=" "), x="Voltage (V)", y = "Velocity (µm/s)")

# setwd(paste(path,"figs",sep="/"))
ggsave(paste("pH_",ph,".eps", sep=""), width = 5, height = 4, units = "in")