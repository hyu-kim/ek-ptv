# Created by Hyungseok Kim (hskimm@mit.edu), last run on 2/18/22
# 2021-12-28: copied from DEP-LPS project
library("ggplot2")
library(tidyverse)

rm(list=ls())

## import and clean
date = '2022-02-16'
path <- paste("/Users/hk/Desktop/LEMI/SFA/Electrokinetics/", date, " bact ek", sep="")
setwd(path)
filename = paste("info_",date,".txt",sep="")
rd_info = read.delim(filename,sep = ",",header=TRUE,dec = ".")

setwd(paste(path,"vy",sep="/"))
df = data.frame(Treatment=factor(), Replicate=factor(), Voltage=double(), Channel=factor(), Velocity=double()) # use only when beginning from scratch
for (ind in c(1:dim(rd_info)[1])) {
  # if (rd_info$'channel'[ind]>14) next # consider QW#1 device only this time
  str = sprintf("%s_R%d_Ch%02d_%s_%02dV_20X_001.ome.csv", rd_info$'cond'[ind], 
                rd_info$'rep'[ind], rd_info$'channel'[ind],
                rd_info$'light'[ind], rd_info$'voltage'[ind])
  rd = read.delim(str, header=TRUE, dec = ".")
  for (i in 1:dim(rd)[1]) {
    df_temp <- data.frame(
      Treatment = rd_info$'cond'[ind],
      Replicate = rd_info$'rep'[ind],
      Voltage = rd_info$'voltage'[ind],
      Channel = rd_info$'channel'[ind],
      Velocity = rd$'velocity'[i])
    df <- rbind(df, df_temp)
  }
}
df$"Replicate" <- as.factor(df$"Replicate")
df$"Channel" <- as.factor(df$"Channel")