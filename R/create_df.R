# Created by Hyungseok Kim (hskimm@mit.edu), last run on 1/18/22
# 2021-12-28: copied from DEP-LPS project
library("ggplot2")
library(tidyverse)

rm(list=ls())

## import and clean
date = '2021-12-31'
path <- paste("/Users/hyungseokkim/Desktop/LEMI/SFA/Electrokinetics/", date, " Pt mobility 4", sep="")
setwd(path)
filename = paste("info_",date,".txt",sep="")
rd_info = read.delim(filename,sep = ",",header=TRUE,dec = ".")

setwd(paste(path,"vy",sep="/"))
df = data.frame(Date=factor(), Treatment=factor(), Replicate=factor(), Voltage=double(), Velocity=double()) # use only when beginning from scratch
for (ind in c(1:dim(rd_info)[1])) {
  if (rd_info$'channel'[ind]>14) next # consider QW#1 device only this time
  str = sprintf("%s_R%d_Ch%02d_TR_%02dV_10X_001.ome.csv", rd_info$'cond'[ind], 
                rd_info$'rep'[ind], rd_info$'channel'[ind], 
                rd_info$'voltage'[ind])
  rd = read.delim(str, header=TRUE, dec = ".")
  for (i in 1:dim(rd)[1]) {
    df_temp <- data.frame(
      Date = rd_info$'date'[ind],
      Treatment = rd_info$'cond'[ind],
      Replicate = rd_info$'rep'[ind],
      Voltage = rd_info$'voltage'[ind],
      Velocity = rd$'velocity'[i])
    df <- rbind(df, df_temp)
  }
}
df$"Replicate" <- as.factor(df$"Replicate")