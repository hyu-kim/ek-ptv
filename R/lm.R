# Created by Hyungseok Kim (hskimm@mit.edu), last modified on 2/11/22
# Continued from create_df.R
# Generates linear model to obtain regression coefficient for mobility
# Parameters = Voltage (float), Replicate (categorical), (Treatment (categorical))

df$"Replicate" <- as.factor(df$"Replicate") # Replicate into categorical variables

eps_r = 80; # relative permitivity []
eps_0 = 8.854e-12; # vacuum permitivity [F/m]
eta = 8.66e-4; # viscosity [Pa-s]
l = 10e-3; # channel length [m]

# subset df and loop
trt <- "Ax"
reps <- c("1","2","3")
phs <- c(8.00, 8.30, 8.58, 8.90, 9.21, 9.48)

setwd(paste(path,"figs",sep="/"))
mu_df <- data.frame(Treatment=factor(), Replicate=factor(), pH=double(), Mobility=double())
for (ph in phs){ # for passing into "plot_all.R"
  for (rep in reps){
    df_params <- df[(df$pH==ph)&(df$Replicate==rep),]
    str <- sprintf("lm_%s_R%s_pH%g.txt", trt, rep, ph)
    res <- lm(Velocity ~ Voltage, data = df_params)
    summary(res)
    sink(str)
    print(summary(res))
    sink()
    coef <- res$coefficients[[2]] # [µm / V-s] 
    mu <- coef * l * 1e-6 # [m2 / V-s]
    mu_df_temp <- data.frame(
      Treatment = trt,
      Replicate = rep,
      pH = ph,
      Mobility = mu
    )
    mu_df <- rbind(mu_df, mu_df_temp)
  }
}

write.table(mu_df, "mu_ek.txt", sep=',', row.names = FALSE)