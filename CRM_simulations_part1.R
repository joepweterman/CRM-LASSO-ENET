# Cellwise Robust M-regression (CRM) - Simulation Study -------------------------------------------



rm(list = ls())
# Setup -------------------------------------------------------------------------------------------

# simulation setting - - - - - - - - - - - - - - -
n <- 400               # number of cases
p <- 50                # number of predictor variables
pct_case_out <- 0.05   # percentage of casewise outliers
pct_cell_out <- 0.10   # percentage of cellwise outliers for each casewise outlier
n_sims <- 100          # number of simulations


# CRM input parameters - - - - - - - - - - - - - -
maxiter             <- 100
tolerance           <- 0.01
outlyingness.factor <- 1.5
spadieta            <- seq(0.9, 0.1, -0.1)



# Load packages -----------------------------------------------------------------------------------
library(crmReg)
library(plyr)
library(FNN)
library(glmnet)
library(MASS)
library(cellWise)
library(lubridate)
library(robustbase)
source("cellwiseheatmap.R")
source("crm.R")
source("daprpr.R")
source("HampelWeightFunction.R")
source("impute_outlying_cells.R")
source("predict.crm.R")
source("scaleResidualsByMAD.R")
source("spadimo.R")
source("spadimo_lasso.R")
source("spadimo_enet.R")
source("crm_lasso.R")
source("crm_enet.R")

# Start simulation procedure ----------------------------------------------------------------------
results_MAE       <- data.frame(matrix(nrow = n_sims, ncol = 3))
results_MSEP      <- data.frame(matrix(nrow = n_sims, ncol = 3))
results_RMSEI     <- data.frame(matrix(nrow = n_sims, ncol = 3))
results_precision <- data.frame(Precision = rep(NA, n_sims))
results_recall    <- data.frame(Recall = rep(NA, n_sims))
results_time      <- data.frame(Time = rep(NA, n_sims))
results_time_l      <- data.frame(Time = rep(NA, n_sims))
results_time_e      <- data.frame(Time = rep(NA, n_sims))

names(results_MAE)   <- c("CRM", "CRM-LASSO", "CRM-ElasticNet")
names(results_MSEP)  <- c("CRM", "CRM-LASSO", "CRM-ElasticNet")
names(results_RMSEI) <- c("CRM", "CRM-LASSO", "CRM-ElasticNet")


t_start <- proc.time()
set.seed(2019)
cat(paste("\n*", n_sims ,"simulations - started at", format(Sys.time(), "%X"),
          "============================================================\n\n"))


for (j_sim in 1:n_sims) {


  # Generate clean sample -------------------------------------------------------------------------
  mu <- rep(0, p)
  Sigma <- diag(p)
  Sigma[(row(Sigma) - col(Sigma)) == 1 | row(Sigma) - col(Sigma) == -1] <- 0.5
  X <- mvrnorm(n = n, mu = mu, Sigma = Sigma)

  slopes <- rnorm(n = p, mean = 0, sd = 1)
  slopes <- 10 * slopes / sqrt(sum(slopes^2))
  intercept <- 10
  noise <- rnorm(n = n, mean = 0, sd = 0.5)

  y <- intercept + X %*% slopes + noise
  betas <- c(intercept, slopes)



  # Add contamination in design matrix ------------------------------------------------------------
  Xc <- X
  contamination <- colMeans(X) + 6 * apply(X, 2, sd)

  case_outliers <- sample(n, size = n * pct_case_out)
  outliers_mat_flag <- matrix(FALSE, nrow = n, ncol = p)

  for (i in case_outliers) {
    cell_outliers <- sample(p, size = p * pct_cell_out)
    Xc[i, cell_outliers] <- contamination[cell_outliers] + rnorm(length(cell_outliers))
    outliers_mat_flag[i, cell_outliers] <- TRUE
  }



  # Apply DDC to contaminated predictor data ------------------------------------------------------
  cat(" DDC:")
  DDC_Ximp <- DDC(Xc)$Ximp



  # Collect data samples --------------------------------------------------------------------------
  data_clean        <- cbind.data.frame(y, X)
  data_contaminated <- cbind.data.frame(y, Xc)
  data_DDC          <- cbind.data.frame(y, DDC_Ximp)
  names(data_clean) <- names(data_contaminated) <- names(data_DDC) <- c("y", paste0("X", 1:p))



  # Fit regression models -------------------------------------------------------------------------
  crm_fit <- suppressWarnings(crm(formula   = y ~ .,
                                  data      = data_contaminated,
                                  maxiter   = maxiter,
                                  tolerance = tolerance,
                                  outlyingness.factor = outlyingness.factor,
                                  spadieta  = spadieta,
                                  center    = "median",
                                  scale     = "Qn",
                                  regtype   = "MM",
                                  verbose   = FALSE))

  crm_lts_fit <- suppressWarnings(crm(formula   = y ~ .,
                                      data      = data_contaminated,
                                      maxiter   = maxiter,
                                      tolerance = tolerance,
                                      outlyingness.factor = outlyingness.factor,
                                      spadieta  = spadieta,
                                      center    = "median",
                                      scale     = "Qn",
                                      regtype   = "LTS",
                                      alphaLTS  = 0.5,
                                      verbose   = FALSE))
  crm_fit_l <- suppressWarnings(crm_lasso(formula   = y ~ .,
                                  data      = data_contaminated,
                                  maxiter   = maxiter,
                                  tolerance = tolerance,
                                  outlyingness.factor = outlyingness.factor,
                                  spadieta  = spadieta,
                                  center    = "median",
                                  scale     = "Qn",
                                  regtype   = "MM",
                                  verbose   = FALSE))

  crm_lts_fit_l <- suppressWarnings(crm_lasso(formula   = y ~ .,
                                      data      = data_contaminated,
                                      maxiter   = maxiter,
                                      tolerance = tolerance,
                                      outlyingness.factor = outlyingness.factor,
                                      spadieta  = spadieta,
                                      center    = "median",
                                      scale     = "Qn",
                                      regtype   = "LTS",
                                      alphaLTS  = 0.5,
                                      verbose   = FALSE))

  crm_fit_e <- suppressWarnings(crm_enet(formula   = y ~ .,
                                  data      = data_contaminated,
                                  maxiter   = maxiter,
                                  tolerance = tolerance,
                                  outlyingness.factor = outlyingness.factor,
                                  spadieta  = spadieta,
                                  center    = "median",
                                  scale     = "Qn",
                                  regtype   = "MM",
                                  verbose   = FALSE))

  crm_lts_fit_e <- suppressWarnings(crm_enet(formula   = y ~ .,
                                      data      = data_contaminated,
                                      maxiter   = maxiter,
                                      tolerance = tolerance,
                                      outlyingness.factor = outlyingness.factor,
                                      spadieta  = spadieta,
                                      center    = "median",
                                      scale     = "Qn",
                                      regtype   = "LTS",
                                      alphaLTS  = 0.5,
                                      verbose   = FALSE))

  # Evaluate performance --------------------------------------------------------------------------

  # Mean Absolute Error - - - - - - - - - - - - -
  results_MAE$CRM[j_sim]               <- mean(abs(    crm_fit$coefficients - betas))
  results_MAE$`CRM-LASSO`[j_sim]       <- mean(abs(    crm_fit_l$coefficients - betas))
  results_MAE$`CRM-ElasticNet`[j_sim]  <- mean(abs(    crm_fit_e$coefficients - betas))


  # Mean Squared Error of Prediction - - - - - - -
  results_MSEP$CRM[j_sim]               <- mean((    crm_fit$fitted.values - data_clean$y)[-case_outliers]^2)
  results_MSEP$`CRM-LASSO`[j_sim]       <- mean((     crm_fit_l$fitted.values - data_clean$y)[-case_outliers]^2)
  results_MSEP$`CRM-ElasticNet`[j_sim]  <- mean(( crm_fit_e$fitted.values - data_clean$y)[-case_outliers]^2)


  # Root Mean Squared Error of Imputation - - - -
  results_RMSEI$CRM[j_sim]               <- sqrt(mean((as.matrix(data_clean[, -1] - crm_fit$data.imputed[, -1]))^2))
  results_RMSEI$`CRM-LASSO`[j_sim]       <- sqrt(mean((as.matrix(data_clean[, -1] - crm_fit_l$data.imputed[, -1]))^2))
  results_RMSEI$`CRM-ElasticNet`[j_sim]  <- sqrt(mean((as.matrix(data_clean[, -1] - crm_fit_e$data.imputed[, -1]))^2))


  # Execution time - - - - - - - - - - - - - - - -
  results_time$Time[j_sim] <- crm_fit$time
  results_time_l$Time[j_sim] <- crm_fit_l$time
  results_time_e$Time[j_sim] <- crm_fit_e$time


  # Elapsed time - - - - - - - - - - - - - - - - -
  t_end <- seconds_to_period(round((proc.time() - t_start)[3]))
  cat(paste0(" - ", round(100 * j_sim / n_sims, 2),
             "% - ", sprintf("%02d:%02d:%02d", t_end@hour, minute(t_end), second(t_end)), "\n\n"))
}


t_end <- seconds_to_period(round((proc.time() - t_start)[3]))
cat(paste("\nTime elapsed:",
          sprintf("%02d:%02d:%02d", t_end@hour, minute(t_end), second(t_end)), "\n\n"))



# Study results -----------------------------------------------------------------------------------

# Mean Absolute Error - - - - - - - - - - - - - -
boxplot(results_MAE, ylab = "MAE", col = "lightblue",
        ylim = c(0, max(results_MAE)), las = 1, cex.axis = 1.2, cex.lab = 1.3)
points(colMeans(results_MAE), col = "red", pch = 18, cex = 1.7)
text(rep(0, 8), labels = round(colMeans(results_MAE), 4), cex = 1.4,
     font = ifelse(1:3 == which.min(colMeans(results_MAE)), 2, 1))


# Mean Squared Error of Prediction - - - - - - - -
boxplot(results_MSEP, ylab = "MSEP", col = "lightblue",
        ylim = c(-0.5, max(results_MSEP)), las = 1, cex.axis = 1.2, cex.lab = 1.3)
points(colMeans(results_MSEP), col = "red", pch = 18, cex = 1.7)
text(rep(-0.5, 8), labels = round(colMeans(results_MSEP), 4), cex = 1.4,
     font = ifelse(1:3 == which.min(colMeans(results_MSEP)), 2, 1))


# Root Mean Squared Error of Imputation - - - - -
#par(mfrow = c(1, 3))
boxplot(results_RMSEI, ylab = "RMSEI", col = "lightblue",
        ylim = c(0, max(results_RMSEI)), las = 1, cex.axis = 1.2, cex.lab = 1.3)
points(colMeans(results_RMSEI), col = "red", pch = 18, cex = 1.7)
text(rep(0, 3), labels = round(colMeans(results_RMSEI), 4), cex = 1.4,
     font = ifelse(1:3 == which.min(colMeans(results_RMSEI)), 2, 1))


# Execution time - - - - - - - - - - - - - - - - -
cat("CRM average execution time:", round(mean(results_time$Time), 1), "seconds\n\n")
cat("CRM-LASSO average execution time:", round(mean(results_time_l$Time), 1), "seconds\n\n")
cat("CRM-ElasticNet average execution time:", round(mean(results_time_e$Time), 1), "seconds\n\n")

