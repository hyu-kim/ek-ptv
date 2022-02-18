# Created by Hyungseok Kim (hskimm@mit.edu), last modified on 2/18/22
# Continued from "create_df.R" THEN "lm.R"
# Conducts statistical t-test to compare mobility between treatements

mu_df2 <- mu_df[!(mu_df$Channel=='16'),]

trt_test <- c('ARW1R1', '3-2')
mu_test1 <- mu_df2[mu_df2$Strain==trt_test[1],]$Mobility
mu_test2 <- mu_df2[mu_df2$Strain==trt_test[2],]$Mobility

str <- sprintf("t_%s_vs_%s.txt", trt_test[1], trt_test[2])
res2 <- t.test(mu_test1, mu_test2, alternative='two.sided', paired=FALSE)
sink(str)
print(paste(trt_test[1],'vs', trt_test[2]))
print(res2)
sink()