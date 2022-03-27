# Created by Hyungseok Kim (hskimm@mit.edu), last modified on 2/18/22
# Continued from "create_df.R" THEN "lm.R"
# Plots EK mobility by each treatment
library("ggplot2")
library(tidyverse)
library(Rmisc)

mu_df2 <- mu_df[!(mu_df$Channel=='16'),]

mu_df2 %>%
  mutate(Strain = factor(Strain, levels=trts)) %>%
  mutate(Channel = factor(Channel, levels=chls)) %>%
  ggplot(aes(x=Strain, y=Mobility * 1e9)) + 
  geom_bar(stat = "summary", fun = mean, size = 0.2, width = 0.75, fill=NA, colour='black') +
  geom_point(size = 2, position = position_jitterdodge(0.1),
             aes(colour=Channel)) +
  scale_color_manual(values=c("#FF0000", "#00FF0C", "#0812F7", "#B32222", "#06920D", "#0D7DA7", "#EE7230", "#81EF86", "#95DEEC")) +
  # geom_point(stroke = 0.2, size = 1, pch = 21, position = position_jitterdodge(0.3)) +
  labs(x='Bacterial strain', y=expression(paste('Mobility (','- 10'^{-9}, ' m'^2,'/ V-s)'))) +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=0.15))

ggsave('mobility.eps', width = 3.5, height = 3.5, units = "in")

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