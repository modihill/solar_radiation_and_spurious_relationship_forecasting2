##### Importing the necessary libraries

library(readr)
library(dplyr)
library(tidyr)
library(knitr)
library(TSA)
library(tseries)
library(car)
library(dynlm)
library(Hmisc)
library(forecast)
library(xts)
library(ggplot2)
library(AER)
library(x12)
library(dLagM)
library(kableExtra)


##### Function Definations 

# Descriptive Analysis function
descriptive_analysis <- function(ts, object)
{
  plot(ts,
       ylab = c(paste0(toString(object))),
       main = c(paste0("Time Series Plot of ",toString(object))),
       type="o")
  points(y=ts,x=time(ts), pch=as.vector(season(ts)))
  
  acf(ts,
      lag.max = 48,
      main = c(paste0("ACF plot of ",toString(object))))
  
  print(adf.test(ts))
}


# Function for Decomposition
decom <- function(ts, ts_series)
{
  decom.x12 = x12(ts)
  plot(decom.x12 , sa=TRUE , trend=TRUE,
       main = c(paste0("Monthly ", toString(ts_series)), " X12 Decomposed Series"))
  
  plotSeasFac(decom.x12)
  
  
  decomposition <- stl(ts, t.window=15, s.window="periodic", robust=TRUE)
  plot(decomposition,
       main = c(paste0("Monthly ", toString(ts_series)), " STL Decomposed Series"))
}


# Function for Summary and Residual Analysis
summary_residual_analysis <-  function(m)
{
  summary(m, diagnostics = TRUE)
  checkresiduals(m$model)
  print(bgtest(m$model))
  print(vif(m$model))
  print(shapiro.test(m$model$residuals))
}


# Task-1 : Analysis and Prediction of Monthly Average Solar Radiation

## Data Preparation
{r, message=FALSE}
task_1 <- read_csv("data.x.csv")
solar_radiation <- read_csv("data1.csv")

head(solar_radiation)
class(solar_radiation)

After importing necessary libraries and dataset, We must convert each variable into time-series object for Time Series Analysis.


solar_radiation_TS <- ts(solar_radiation, start = c(1960,1), frequency = 12)
solar_TS <- ts(solar_radiation$solar, start = c(1960,1), frequency = 12)
ppt_TS <- ts(solar_radiation$ppt, start = c(1960,1), frequency = 12)

solar_radiation_TS %>% head()
class(solar_radiation_TS)

##### 1. Solar Radiation

head(solar_TS)
class(solar_TS)

##### 2. Precipitation

head(ppt_TS)
class(ppt_TS)


## Descriptive Analysis

#### 1. Solar Radiation
descriptive_analysis(solar_TS, "Monthly Average Solar Radiation")


#### 2. Precipitation
descriptive_analysis(ppt_TS, "Monthly Precipitation")


#### 3. Combined Scaled Time Series Plot
combined <- scale(solar_radiation_TS)

plot(combined, 
     plot.type="s", 
     col = c("#05386b","#f01b1d"), 
     main = "Scaled Time Series Plot of Solar Radiation and Precipitation")

legend("topleft", 
       lty=1, 
       col = c("#05386b","#f01b1d"), 
       c("Solar Radiation", "Precipitation"))


#Correlation
  
cor(solar_TS, ppt_TS) %>% round(3)


## Decomposition

#### 1. Solar Radiation
decom(solar_TS, "Average Solar Radiatoin")


#### 2. Precipitation
decom(ppt_TS, "Precipitation")


### Finite Distributed Lag Model
for (i in 1:12)
{
  model_dlm <- dlm(formula = solar ~ ppt, data = data.frame(solar_radiation), q = i)
  cat("q = ", i, 
      "AIC = ", AIC(model_dlm$model), 
      "BIC = ", BIC(model_dlm$model), 
      "MASE = ", MASE(model_dlm)$MASE, "\n")
}


model1 <- dlm(formula = solar ~ ppt, 
              data = data.frame(solar_radiation), 
              q = 12)

summary_residual_analysis(model1)



### PolyNomial Distributed Lag Model
model2 = polyDlm(x = as.vector(solar_radiation$ppt), 
                 y = as.vector(solar_radiation$solar), 
                 q = 12, 
                 k = 2, 
                 show.beta = TRUE)

summary_residual_analysis(model2)


### Koyck Distributed Lag Model
model3 = koyckDlm(x = as.vector(solar_radiation$ppt), 
                 y = as.vector(solar_radiation$solar))
summary_residual_analysis(model3)


### Autoregressive Distributed Lag Model
for (i in 1:5)
{
  for(j in 1:5)
  { 
    model_ardlm <- ardlDlm(formula = solar ~ ppt, 
                           data = data.frame(solar_radiation), 
                           p = i,
                           q = j)
    cat("p =", i, 
        "q =" , j,
        "AIC =", AIC(model_ardlm$model), 
        "BIC =", BIC(model_ardlm$model),
        "MASE =", MASE(model_ardlm)$MASE, "\n")
  }
}


#ardlDLM(3,5)
model4_1 = ardlDlm(formula = solar ~ ppt, 
        data = data.frame(solar_radiation), 
        p = 3,
        q = 5)
        
summary_residual_analysis(model4_1)

#ardlDLM(4,5)
model4_2 = ardlDlm(formula = solar ~ ppt, 
        data = data.frame(solar_radiation), 
        p = 4,
        q = 5)
        
summary_residual_analysis(model4_2)

#ardlDLM(5,5)
model4_3 = ardlDlm(formula = solar ~ ppt, 
        data = data.frame(solar_radiation), 
        p = 5,
        q = 5)
        
summary_residual_analysis(model4_3)


### Dynamic Models

#model5_1
model5_1 = dynlm(solar_TS ~ L(solar_TS , k = 1 ) + season(solar_TS))
summary(model5_1)
checkresiduals(model5_1)

#model5_2
model5_2 = dynlm(solar_TS ~ L(solar_TS , k = 1 ) + trend(solar_TS) + season(solar_TS))
summary(model5_2)
checkresiduals(model5_2)

#model5_3
model5_3 = dynlm(solar_TS ~ L(solar_TS , k = 1 ) + L(solar_TS , k = 2 ) + trend(solar_TS) + season(solar_TS))
summary(model5_3)
checkresiduals(model5_3)

#model5_4
model5_4 = dynlm(solar_TS ~ L(solar_TS , k = 1 ) + L(solar_TS , k = 2 ) + L(solar_TS , k = 3 ) +  season(solar_TS))
summary(model5_4)
checkresiduals(model5_4)

#model5_5
model5_5 = dynlm(solar_TS ~ L(solar_TS , k = 1 ) + L(solar_TS , k = 2 ) + L(solar_TS , k = 3 ) + trend(solar_TS) + season(solar_TS))
summary(model5_5)
checkresiduals(model5_5)


#Model Comparison
#The following table display the each fitted Time-series regression models with `MASE()`, `AIC()` and `BIC()`. 


attr(model3$model,"class") = "lm"
models <- c("Finite DLM", "Poly DLM", "Koyck", "ARDL_3_5", "ARDL_4_5", "ARDL_5_5", "dynlm_1", "dynlm_2", "dynlm_3", "dynlm_4", "dynlm_5")

aic <- AIC(model1$model, model2$model, model3$model, model4_1$model, model4_2$model, model4_3$model, model5_1,model5_2,model5_3, model5_4, model5_5)$AIC

bic <- BIC(model1$model, model2$model, model3$model, model4_1$model, model4_2$model, model4_3$model, model5_1,model5_2,model5_3, model5_4, model5_5)$BIC

mase <- MASE(model1$model, model2$model, model3$model, model4_1, model4_2, model4_3, lm(model5_1), lm(model5_2), lm(model5_3), lm(model5_4), lm(model5_5))$MASE

Model_Comparison <- data.frame(models, mase, aic, bic)
colnames(Model_Comparison) <- c("Model","MASE","AIC", "BIC")



## Exponential smoothing methods
HW_models = c("Holt_Winter additive method",
              "Holt_Winter multiplicative method with exponential trend",
              "Holt_Winter multiplicative method",
              "Holt_Winter additive method",
              "Holt_Winter multiplicative method with exponential trend",
              "Holt_Winter multiplicative method")

exponential = c(TRUE,FALSE)
seasonality = c("additive","multiplicative")
damped = c(TRUE,FALSE)
exponential_models <- expand.grid(exponential, seasonality, damped)
exponential_models <- exponential_models[-c(1,5),]

HW_AIC <- array(NA, 6)
HW_BIC <- array(NA, 6)
HW_MASE <- array(NA, 6)
levels <- array(NA, dim=c(6,3))

for (i in 1:6)
{
  HW_model <- hw(solar_TS,
                 exponential = exponential_models[i,1],
                 seasonal = toString(exponential_models[i,2],
                                     damped = exponential_models[i,3]))
  HW_AIC[i] <- HW_model$model$aic
  HW_BIC[i] <- HW_model$model$bic
  HW_MASE[i] <- accuracy(HW_model)[6]
  levels[i,1] <- exponential_models[i,1]
  levels[i,2] <- toString(exponential_models[i,2])
  levels[i,3] <- exponential_models[i,3]
  summary(HW_model)
  checkresiduals(HW_model)
  print(shapiro.test(HW_model$model$residuals))
}

results_HW = data.frame(HW_models, levels, HW_MASE, HW_AIC, HW_BIC)
colnames(results_HW) = c("Model", "Exponential","Seasonality","Damped","MASE","AIC", "BIC")

kbl(results_HW) %>% kable_paper()


## State-space Models
ets_models = c("AAA", "MAA", "MAM", "MMM")
damped = c(TRUE,FALSE)
ETS_models <- expand.grid(ets_models, damped)

ETS_AIC <- array(NA, 8)
ETS_BIC <- array(NA, 8)
ETS_MASE <- array(NA, 8)
levels <- array(NA, dim=c(8,2))

for (i in 1:8)
{
  ETS <- ets(solar_TS,
             model = toString(ETS_models[i, 1]), damped = ETS_models[i,2])
  ETS_AIC[i] <- ETS$aic
  ETS_BIC[i] <- ETS$bic
  ETS_MASE[i] <- accuracy(ETS)[6]
  levels[i,1] <- toString(ETS_models[i,1])
  levels[i,2] <- ETS_models[i,2]
  summary(ETS)
  checkresiduals(ETS)
  print(shapiro.test(ETS$residuals))
}

results_ETS = data.frame(levels, ETS_MASE, ETS_AIC, ETS_BIC)
colnames(results_ETS) = c("Model","Damped","MASE","AIC", "BIC")

kbl(results_ETS) %>% kable_paper()


### Auto-Fit State-Space Model
ETS_auto_model = ets(solar_TS,model="ZZZ")
ETS_auto_model$method


## Model Comparison
Formating results_ETS table 


results_ETS$Damped <- factor(results_ETS$Damped,
                             levels = c(TRUE, FALSE),
                             labels = c("Dumped"," "))

results_ETS <- unite(results_ETS,
                     "Model", c("Model","Damped"), sep = "_")

kbl(results_ETS) %>% kable_paper()

Formating results_HW table 


results_HW$Damped <- factor(results_HW$Damped,
                             levels = c(TRUE, FALSE),
                             labels = c("Dumped"," "))

results_HW <- unite(results_HW,
                     "Model", c("Model","Damped"), sep = "_")

results_HW <- results_HW[,-c(2,3)]
kbl(results_HW) %>% kable_paper()

#Merge Model Comparison Table
Model_Comparison <- rbind(Model_Comparison, results_ETS, results_HW)

sorted_MASE <- Model_Comparison %>% arrange(MASE)
kbl(sorted_MASE) %>%
  kable_paper()



## Forecasting
prediction_1 <- hw(solar_TS,
                   seasonal = "multiplicative",
                   dumped = TRUE,
                   h = 2*frequency(solar_TS))

prediction_2 <- hw(solar_TS,
                   seasonal = "multiplicative",
                   dumped = TRUE,
                   exponential = TRUE,
                   h = 2*frequency(solar_TS))


prediction_3 <- ets(solar_TS,
                    model="AAA",
                    damped = T)
prediction_3 <- forecast(prediction_3)

plot(prediction_3,
     main = "Next Two years Forecast of Solar Radiation)",
     ylab = "Solar Radiation",
     fcol = "#e8d31d")

lines(fitted(prediction_3), col = "#e8d31d")

lines(fitted(prediction_2), col = "#039fbe")
lines(prediction_2$mean, col = "#039fbe", lwd = 2)

lines(fitted(prediction_1), col = "#b20238")
lines(prediction_1$mean, col = "#b20238", lwd = 2)


legend("bottomleft",
       lty = 1,
       col = c("black", "#b20238", "#039fbe", "#e8d31d"),
       c("Data", "Holt-Winters' Multiplicative_Damped", "Holt-Winters' Multiplicative Exponential_Damped", "ETS(A,Ad,A)"))


# final prediction
plot(prediction_1, fcol = "#b20238", 
     main = "Forecasting of Solar Radiation in 2015 and 2106",
     ylab = "Radiation")
lines(fitted(prediction_1), col = "#b20238")

legend("topleft", 
       lty = 1, 
       col = c("black", "#b20238"), 
       c("Data", "Prediction"))




# confidence interval
kbl(prediction_1) %>% kable_paper()



# Task - 2 : Demonstrate whether correlation between Residential PPI and Population Change is Spurious or not.

# data reading
task_2 <- read_csv("data2.csv")

head(task_2)
class(task_2)


## Data Preparation

# time series object
task_2_TS <- ts(task_2[,2:3], start = c(2003,3), frequency = 4)
price_TS <- ts(task_2$price, start = c(2003,3), frequency = 4)
change_TS <- ts(task_2$change, start = c(2003,3), frequency = 4)


task_2_TS %>% head()
class(task_2_TS)

##### 1. Residential Property Price Index

head(price_TS)
class(price_TS)

##### 2. Population Change

head(change_TS)
class(change_TS)


## Descriptive Analysis
plot(task_2_TS,
     main = "Time Series plot Residential PPI and Population Change",
     yax.flip = TRUE, 
     type = 'o')


## Cross-Covariance Function
ccf(as.vector(price_TS), 
    as.vector(change_TS), 
    ylab = "Cross-Covariance Function", 
    main = "CCF plot of Residential PPI and Population Change")


## Cheking the Stationarity

#### 1. Residential Property Price Index(PPI)
descriptive_analysis(price_TS, "Residential Property Price Index")


#### 2. Population Change
descriptive_analysis(change_TS, "Population Change")


## Transformation to Stationary Series

#### 1. Residential Property Price Index(Differenced)

price_TS_diff <- diff(diff(price_TS),4)
descriptive_analysis(price_TS_diff, "Residential PPI (Differenced)")


#### 2. Population Change(Differenced)

change_TS_diff <- diff(diff(change_TS),4)
descriptive_analysis(change_TS_diff, "Population Change(Differenced)")


## Prewhitening
prewhiten(as.vector(price_TS_diff), 
          as.vector(change_TS_diff), 
          ylab='Cross-Covariance Function', 
          main = "CCF plot of Residential PPI and Population Change after prewhitening")