# Created by Hyungseok Kim (hskimm@mit.edu), last modified on 12/28/21
# Continued from create_df.R
# Generates linear model to obtain regression coefficient for mobility
# Parameters = Voltage (float), Replicate (categorical), (Treatment (categorical))

df$"Replicate" <- as.factor(df$"Replicate") # Replicate into categorical variables
df$Date <- as.factor(df$Date)
df$Treatment <- as.factor(df$Treatment)

eps_r = 80; # relative permitivity []
eps_0 = 8.854e-12; # vacuum permitivity [F/m]
eta = 8.66e-4; # viscosity [Pa-s]
l = 10e-3; # channel length [m]

# subset df and loop
dates = c('2021-12-22', '2021-12-31')
trts = c('Ax', '+13C1', '+1R1')
reps = c(1,2,3)

mu_df <- data.frame(Date=factor(), Treatment=factor(), Replicate=factor(), Mobility=double())
for (trt in trts){ # for passing into "plot_all.R"
  for (d in dates){
    for (rep in reps){
      df_params <- df[(df$Treatment==trt)&(df$Date==d)&(df$Replicate==rep),]
      res <- lm(Velocity ~ Voltage, data = df_params)
      print(summary(res))
      coef <- res$coefficients[[2]] # [µm / V-s] 
      mu <- coef * l * 1e-6 # [m2 / V-s]
      mu_df_temp <- data.frame(
        Date = d,
        Treatment = trt,
        Replicate = rep,
        Mobility = mu
        )
      mu_df <- rbind(mu_df, mu_df_temp)
    }
  }
}

# Export mu_df into csv
