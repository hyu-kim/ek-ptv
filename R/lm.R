# Created by Hyungseok Kim (hskimm@mit.edu), last modified on 2/18/22
# Continued from create_df.R
# Generates linear model to obtain regression coefficient for mobility
# Parameters = Voltage (float), Replicate (categorical), (Treatment (categorical))

eps_r = 80; # relative permitivity []
eps_0 = 8.854e-12; # vacuum permitivity [F/m]
eta = 8.66e-4; # viscosity [Pa-s]
l = 10e-3; # channel length [m]

# subset df and loop
trts = c('13A', '4BL', 'ARW1Y1', 'ARW7G5W')
reps = c('1','2','3')
# chls = c('2','4','5','7','9','10','11','13')

setwd(paste(path,'figs',sep='/'))

mu_df <- data.frame(Treatment=factor(), Replicate=factor(), Channel=factor(), Mobility=double())
for (trt in trts){ # for passing into "plot_all.R"
  for (rep in reps){
    # trt <- rd_info$cond[rd_info$channel==chl][1]
    # rep <- rd_info$rep[rd_info$channel==chl][1]
    chl <- rd_info$rep[(rd_info$cond==trt)&(rd_info$rep==rep)][1]
    df_stat <- df[(df$Treatment==trt)&(df$Replicate==rep),]
    str <- sprintf("lm_%s_R%s_Ch%02d.txt", trt, rep, as.numeric(chl))
    res <- lm(Velocity ~ Voltage, data = df_stat)
    summary(res)
    sink(str)
    print(summary(res))
    sink()
    coef <- res$coefficients[[2]] # [µm / V-s] 
    mu <- coef * l * 1e-6 # [m2 / V-s]
    if (summary(res)$coefficients[2,4] < 1e-10){ # set a threshold using p-value
      mu_df_temp <- data.frame(
        Date = date,
        Strain = trt,
        Replicate = rep,
        Channel = chl,
        Mobility = mu
      )
      mu_df <- rbind(mu_df, mu_df_temp)
    }
  }
}

write.table(mu_df, "summ_ek.txt", sep=',', row.names = FALSE)