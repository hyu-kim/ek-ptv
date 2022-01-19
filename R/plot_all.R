# Created by Hyungseok Kim (hskimm@mit.edu), last modified on 1/18/22
# Continued from "create_df.R" THEN "lm.R"
# Plots mean+-sd of EK mobility by each treatment
library("ggplot2")
library(tidyverse)
library(Rmisc)

# mu_df_SE <- summarySE(
#   mu_df, 
#   measurevar = 'Mobility', groupvars = 'Treatment'
#   ) # summarizes df by each treatment and replicate
# 
# mu_df_SE %>%
#   mutate(Treatment = factor(Treatment, levels=trts)) %>%
#   ggplot( aes(x=Treatment, y=Mobility * 1e9)) +
#     geom_bar(stat="identity", position=position_dodge(), fill=NA, colour='black') + 
#     # geom_errorbar(aes(ymin=(Mobility-sd)*1e9, ymax=(Mobility+sd)*1e9), width=.2) + 
#     geom_point(stroke = 0.2, size = 1, pch = 21, position = position_jitterdodge(0.3)) +
#     labs(x='Strain', y=expression(paste('Mobility (','10'^{-9}, ' m'^2,'/ V-s)')))

mu_df %>%
  mutate(Treatment = factor(Treatment, levels=trts)) %>%
  mutate(Replicate = factor(Replicate, levels=reps)) %>%
  # ggplot(aes(x=pH, y=Mobility * 1e9, fill=Replicate)) + 
  ggplot(aes(x=pH, y=Mobility * 1e9)) + 
  # geom_bar(stat = "summary", fun = mean, size = 0.2, width = 0.75, fill=NA, colour='black') +
  geom_line() + geom_point() +
  # geom_point(stroke = 0.2, size = 1, pch = 21, position = position_jitterdodge(0.3)) +
  labs(x='pH', y=expression(paste('Mobility (','10'^{-9}, ' m'^2,'/ V-s)')))

setwd(paste(path,"figs",sep="/"))
ggsave('mobility.eps', width = 4, height = 3.5, units = "in")

### add a column with a date
mu_df$date <- date
write.csv(mu_df, paste(path,"/",date,"_mu.csv",sep=""),row.names=FALSE)

### combine two dates and plot
setwd("/Users/hk/Desktop/LEMI/SFA/Electrokinetics/codes")
path <- getwd()
dates = c('2021-12-22', '2021-12-31')
mu_df_all <- data.frame(matrix(rep(0,4), nrow=1))
names(mu_df_all) <- colnames(mu_df)
for (date in dates){
  mu_df_temp = read.csv(paste(path,"/",date,"_mu.csv",sep=""))
  mu_df_all <- rbind(mu_df_all, mu_df_temp)
}

mu_df_all <- mu_df_all[-1,]
mu_df_all$Treatment[mu_df_all$Treatment=='Ax'] <- 'Axenic'
mu_df_all$Treatment[mu_df_all$Treatment=='+1R1'] <- '+ARW1R1'
mu_df_all$date[mu_df_all$date=='2021-12-22'] <- 'log (5 d)'
mu_df_all$date[mu_df_all$date=='2021-12-31'] <- 'stat (14 d)'

mu_df_all %>%
  mutate(Treatment = factor(Treatment, levels=c('Axenic', '+13C1', '+ARW1R1'))) %>%
  mutate(Replicate = factor(Replicate, levels=reps)) %>%
  ggplot() + 
  geom_bar(stat = "summary", fun = mean, size = 0.2, width = 0.75, 
           position = position_dodge(),
           aes(x=Treatment, y=Mobility * 1e9, fill=date)) +
  geom_point(stroke = 0.2, size = 1, pch = 21, position = position_jitterdodge(0.1),
             aes(x=Treatment, y=Mobility * 1e9, fill=date)) +
  labs(x='Culture', y=expression(paste('Mobility (','10'^{-9}, ' m'^2,'/ V-s)')), fill='Growth stage')

ggsave('mobility_all.eps', width = 6, height = 3.5, units = "in")