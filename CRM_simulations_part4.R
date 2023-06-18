# Cellwise Robust M-regression (CRM) - Simulation Study -------------------------------------------



rm(list = ls())
# Setup -------------------------------------------------------------------------------------------

# simulation setting - - - - - - - - - - - - - - -
n <- 400                                   # number of cases
p <- 50                                    # number of predictor variables
pct_case_out <- 0.05                       # percentage of casewise outliers
pct_cell_out_seq <- seq(0.05, 0.5, 0.05)   # sequences of percentages of cellwise outliers for each casewise outlier
n_sims <- 10                               # number of simulations for each percentage of cellwise outlier


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
library(pcaPP)
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
results_MAE       <- list()
results_MSEP      <- list()
results_RMSEI     <- list()
results_precision <- list()
results_recall    <- list()
results_time      <- list()

t_start <- proc.time()
set.seed(2020)
cat(paste("\n* Simulations started at", format(Sys.time(), "%X"),
          "============================================================\n"))


for (pct_cell_out_value in pct_cell_out_seq) {

  pct_cell_out_index <- which(pct_cell_out_seq == pct_cell_out_value)
  cat(paste("\n*", pct_cell_out_index, "out of", length(pct_cell_out_seq), "========================================\n"))

  results_MAE_k       <- data.frame(matrix(nrow = n_sims, ncol = 3))
  results_MSEP_k      <- data.frame(matrix(nrow = n_sims, ncol = 3))
  results_RMSEI_k     <- data.frame(matrix(nrow = n_sims, ncol = 3))
  results_precision_k <- data.frame(Precision = rep(NA, n_sims))
  results_recall_k    <- data.frame(Recall = rep(NA, n_sims))
  results_time_k      <- data.frame(Time = rep(NA, n_sims))

  names(results_MAE_k)   <- c("CRM", "CRM-LASSO", "CRM-ElasticNet")
  names(results_MSEP_k)  <- c("CRM", "CRM-LASSO", "CRM-ElasticNet")
  names(results_RMSEI_k) <- c("CRM", "CRM-LASSO", "CRM-ElasticNet")


  for (j_sim in 1:n_sims) {


    # Generate clean sample -----------------------------------------------------------------------
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



    # Add contamination in design matrix ----------------------------------------------------------
    Xc <- X
    contamination <- colMeans(X) + 6 * apply(X, 2, sd)

    case_outliers <- sample(n, size = n * pct_case_out)
    outliers_mat_flag <- matrix(FALSE, nrow = n, ncol = p)

    for (i in case_outliers) {
      cell_outliers <- sample(p, size = p * pct_cell_out_value)
      Xc[i, cell_outliers] <- contamination[cell_outliers] + rnorm(length(cell_outliers))
      outliers_mat_flag[i, cell_outliers] <- TRUE
    }
    outliers_mat_index <- which(outliers_mat_flag, arr.ind = TRUE)



    # Apply DDC to contaminated predictor data ----------------------------------------------------
    cat(" DDC:")
    DDC_Ximp <- DDC(Xc)$Ximp



    # Collect data samples ------------------------------------------------------------------------
    data_clean        <- cbind.data.frame(y, X)
    data_contaminated <- cbind.data.frame(y, Xc)
    data_DDC          <- cbind.data.frame(y, DDC_Ximp)
    names(data_clean) <- names(data_contaminated) <- names(data_DDC) <- c("y", paste0("X", 1:p))



    # Fit regression models -----------------------------------------------------------------------
    crm_fit <- suppressWarnings(crm(formula   = y ~ .,
                                    data      = data_contaminated,
                                    maxiter   = maxiter,
                                    tolerance = tolerance,
                                    outlyingness.factor = outlyingness.factor,
                                    spadieta  = spadieta,
                                    center    = "median",
                                    scale     = "qn",
                                    regtype   = "MM",
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

    # Evaluate performance ------------------------------------------------------------------------

    # Mean Absolute Error - - - - - - - - - - - -
    results_MAE_k$`CRM`[j_sim]             <- mean(abs(    crm_fit$coefficients - betas))
    results_MAE_k$`CRM-LASSO`[j_sim]       <- mean(abs(     crm_fit_l$coefficients - betas))
    results_MAE_k$`CRM-ElasticNet`[j_sim]  <- mean(abs( crm_fit_e$coefficients - betas))


    # Mean Squared Error of Prediction - - - - - -
    results_MSEP_k$`CRM`[j_sim]             <- mean((    crm_fit$fitted.values - data_clean$y)[-case_outliers]^2)
    results_MSEP_k$`CRM-LASSO`[j_sim]       <- mean((     crm_fit_l$fitted.values - data_clean$y)[-case_outliers]^2)
    results_MSEP_k$`CRM-ElasticNet`[j_sim]  <- mean(( crm_fit_e$fitted.values - data_clean$y)[-case_outliers]^2)



    # Root Mean Squared Error of Imputation - - -
    results_RMSEI_k$`CRM`[j_sim]            <- sqrt(mean((as.matrix(data_clean[, -1] - crm_fit$data.imputed[, -1]))^2))
    results_RMSEI_k$`CRM-LASSO`[j_sim]      <- sqrt(mean((as.matrix(data_clean[, -1] - crm_fit_l$data.imputed[, -1]))^2))
    results_RMSEI_k$`CRM-ElasticNet`[j_sim] <- sqrt(mean((as.matrix(data_clean[, -1] - crm_fit_e$data.imputed[, -1]))^2))

    # Execution time - - - - - - - - - - - - - - -
    results_time_k$Time[j_sim] <- crm_fit$time


    # Elapsed time - - - - - - - - - - - - - - - -
    t_end <- seconds_to_period(round((proc.time() - t_start)[3]))
    cat(paste0(" - ", round(100 * (j_sim + (pct_cell_out_index - 1) * n_sims) / (n_sims * length(pct_cell_out_seq)), 2),
               "% - ", sprintf("%02d:%02d:%02d", t_end@hour, minute(t_end), second(t_end)), "\n\n"))


  } # end of forloop "j_sim in 1:n_sims"


  results_MAE[[pct_cell_out_index]]       <- results_MAE_k
  results_MSEP[[pct_cell_out_index]]      <- results_MSEP_k
  results_RMSEI[[pct_cell_out_index]]     <- results_RMSEI_k
  results_time[[pct_cell_out_index]]      <- results_time_k


} # end of forloop "pct_cell_out_value in pct_cell_out_seq"


t_end <- seconds_to_period(round((proc.time() - t_start)[3]))
cat(paste("\nTime elapsed:",
          sprintf("%02d:%02d:%02d", t_end@hour, minute(t_end), second(t_end)), "\n\n"))



# Study results -----------------------------------------------------------------------------------
results_MAE       <- cbind.data.frame(do.call(rbind, results_MAE),       pct_cell_out = sort(rep(pct_cell_out_seq, n_sims)))
results_MSEP      <- cbind.data.frame(do.call(rbind, results_MSEP),      pct_cell_out = sort(rep(pct_cell_out_seq, n_sims)))
results_RMSEI     <- cbind.data.frame(do.call(rbind, results_RMSEI),     pct_cell_out = sort(rep(pct_cell_out_seq, n_sims)))
results_time      <- cbind.data.frame(do.call(rbind, results_time),      pct_cell_out = sort(rep(pct_cell_out_seq, n_sims)))

results_MAE_mean       <- aggregate(. ~ pct_cell_out, data = results_MAE,       FUN = mean)
results_MSEP_mean      <- aggregate(. ~ pct_cell_out, data = results_MSEP,      FUN = mean)
results_RMSEI_mean     <- aggregate(. ~ pct_cell_out, data = results_RMSEI,     FUN = mean)
results_time_mean      <- aggregate(. ~ pct_cell_out, data = results_time,      FUN = mean)


# Mean Absolute Error - - - - - - - - - - - - - -
plot(results_MAE_mean$pct_cell_out, results_MAE_mean$`CRM`, type = "l", xlab = "% cellwise outliers", ylab = "Average MAE",
     ylim = range(results_MAE_mean[, -1]), las = 1, cex.axis = 1.2, cex.lab = 1.3, lwd = 2)
lines(results_MAE_mean$pct_cell_out, results_MAE_mean$`CRM-LASSO`,      lwd = 2, lty = 2)
lines(results_MAE_mean$pct_cell_out, results_MAE_mean$`CRM-ElasticNet`,  lwd = 2, lty = 3)
legend("topright", legend = names(results_MAE_mean)[-1], lty = 1:3, lwd = 2)


# Mean Squared Error of Prediction - - - - - - - -
plot(results_MSEP_mean$pct_cell_out, results_MSEP_mean$`CRM`, type = "l", xlab = "% cellwise outliers", ylab = "Average MSEP",
     ylim = range(results_MSEP_mean[, -1]), las = 1, cex.axis = 1.2, cex.lab = 1.3, lwd = 2)
lines(results_MSEP_mean$pct_cell_out, results_MSEP_mean$`CRM-LASSO`,      lwd = 2, lty = 2)
lines(results_MSEP_mean$pct_cell_out, results_MSEP_mean$`CRM-ElasticNet`,  lwd = 2, lty = 3)
legend("topright", legend = names(results_MSEP_mean)[-1], lty = 1:3, lwd = 2)


# Root Mean Squared Error of Imputation - - - - -
plot(results_RMSEI_mean$pct_cell_out, results_RMSEI_mean$`CRM`, type = "l", xlab = "% cellwise outliers", ylab = "Average RMSEI",
     ylim = range(results_RMSEI_mean[, -1]), las = 1, cex.axis = 1.2, cex.lab = 1.3, lwd = 2)
lines(results_RMSEI_mean$pct_cell_out, results_RMSEI_mean$`CRM-LASSO`, lwd = 2, lty = 2)
lines(results_RMSEI_mean$pct_cell_out, results_RMSEI_mean$`CRM-ElasticNet`, lwd = 2, lty = 3)
legend("topleft", legend = names(results_RMSEI_mean)[-1], lty = 1:3, lwd = 2)

# Execution time - - - - - - - - - - - - - - - - -
cat("CRM average execution time:", round(mean(results_time$Time), 1), "seconds")

