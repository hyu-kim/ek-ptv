# Created by Hyungseok Kim (hskimm@mit.edu), last modified on 12/28/21
# Continued from create_df.R
# Generates linear model to obtain regression coefficient for mobility
# Parameters = Voltage (float), Replicate (categorical), (Treatment (categorical))

df$"Replicate" <- as.factor(df$"Replicate") # Replicate into categorical variables

eps_r = 80; # relative permitivity []
eps_0 = 8.854e-12; # vacuum permitivity [F/m]
eta = 8.66e-4; # viscosity [Pa-s]
l = 10e-3; # channel length [m]

# subset df and loop
trts = c('Ax')
reps = c(1,2,3)
phs = c(7.76, 8, 8.21, 9.16, 9.43, 8.57, 9.69, 9.94)

mu_df <- data.frame(Treatment=factor(), Replicate=factor(), pH=double(), Mobility=double())
for (trt in trts){ # for passing into "plot_all.R"
  # for (rep in reps){
    for (ph in phs){
      rep <- rd_info$rep[rd_info$ph==ph][1]
      df_params <- df[(df$Treatment==trt)&(df$pH==ph),]
      res <- lm(Velocity ~ Voltage, data = df_params)
      print(summary(res))
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
  # }
}