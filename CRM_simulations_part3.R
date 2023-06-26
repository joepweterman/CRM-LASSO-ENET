# Cellwise Robust M-regression (CRM) - Simulation Study -------------------------------------------



rm(list = ls())
# Setup -------------------------------------------------------------------------------------------

# simulation setting - - - - - - - - - - - - - - -
n <- 400                                # number of cases
p <- 50                                 # number of predictor variables
pct_case_out_seq <- seq(0, 0.5, 0.05)   # sequence of percentages of casewise outliers
pct_cell_out <- 0.10                    # percentage of cellwise outliers for each casewise outlier
n_sims <- 10                            # number of simulations for each value of pct_case_out


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
results_time      <- list()

t_start <- proc.time()
set.seed(2019)
cat(paste("\n* Simulations started at", format(Sys.time(), "%X"),
          "============================================================\n"))


for (j_sim in 1:n_sims) {
  cat(paste0("\n* Simulation ", j_sim, "/", n_sims, " =======================================\n"))

  n_pct_case_out <- length(pct_case_out_seq)

  results_MAE_j_sim       <- data.frame(matrix(nrow = n_pct_case_out, ncol = 3))
  results_MSEP_j_sim      <- data.frame(matrix(nrow = n_pct_case_out, ncol = 3))
  results_RMSEI_j_sim     <- data.frame(matrix(nrow = n_pct_case_out, ncol = 3))
  results_time_j_sim      <- data.frame(Time = rep(NA, n_pct_case_out))

  names(results_MAE_j_sim)   <- c("CRM", "CRM-LASSO", "CRM-ElasticNet")
  names(results_MSEP_j_sim)  <- c("CRM", "CRM-LASSO", "CRM-ElasticNet")
  names(results_RMSEI_j_sim) <- c("CRM", "CRM-LASSO", "CRM-ElasticNet")



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



  # Setup for contamination -----------------------------------------------------------------------
  Xc <- X
  contamination <- colMeans(X) + 6 * apply(X, 2, sd)
  case_outliers <- c()
  uncontaminated_rows <- 1:n
  outliers_mat_flag <- matrix(FALSE, nrow = n, ncol = p)



  for (case_pct_index in 1:n_pct_case_out) {

    # Add contamination in design matrix ----------------------------------------------------------
    extra_pct_case_out <- ifelse(case_pct_index == 1,
                                 pct_case_out_seq[1],
                                 pct_case_out_seq[case_pct_index] - pct_case_out_seq[case_pct_index - 1])

    new_case_outliers   <- sample(uncontaminated_rows, size = round(n * extra_pct_case_out))
    case_outliers       <- c(case_outliers, new_case_outliers)
    uncontaminated_rows <- setdiff(uncontaminated_rows, case_outliers)

    cell_outliers <- NULL
    for (i in new_case_outliers) {
      cell_outliers <- sample(p, size = p * pct_cell_out)
      Xc[i, cell_outliers] <- contamination[cell_outliers] + rnorm(length(cell_outliers))
      outliers_mat_flag[i, cell_outliers] <- TRUE
    }

    cat(paste0("\n # casewise outliers = ", length(case_outliers),
               " (",100 * pct_case_out_seq[case_pct_index], "%)\n"))



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
    results_MAE_j_sim$`CRM`[case_pct_index]             <- mean(abs( crm_fit$coefficients - betas))
    results_MAE_j_sim$`CRM-LASSO`[case_pct_index]       <- mean(abs( crm_fit_l$coefficients - betas))
    results_MAE_j_sim$`CRM-ElasticNet`[case_pct_index]  <- mean(abs( crm_fit_e$coefficients - betas))

    # Mean Squared Error of Prediction - - - - - -
    if (length(case_outliers) == 0) {
      results_MSEP_j_sim$`CRM`[case_pct_index]            <- mean(( crm_fit$fitted.values - data_clean$y)^2)
      results_MSEP_j_sim$`CRM-LASSO`[case_pct_index]      <- mean(( crm_fit_l$fitted.values - data_clean$y)^2)
      results_MSEP_j_sim$`CRM-ElasticNet`[case_pct_index] <- mean(( crm_fit_e$fitted.values - data_clean$y)^2)
    } else {
      results_MSEP_j_sim$`CRM`[case_pct_index]             <- mean(( crm_fit$fitted.values - data_clean$y)[-case_outliers]^2)
      results_MSEP_j_sim$`CRM-LASSO`[case_pct_index]       <- mean(( crm_fit_l$fitted.values - data_clean$y)[-case_outliers]^2)
      results_MSEP_j_sim$`CRM-ElasticNet`[case_pct_index]  <- mean(( crm_fit_e$fitted.values - data_clean$y)[-case_outliers]^2)
    }


    # Root Mean Squared Error of Imputation - - -
    results_RMSEI_j_sim$CRM[case_pct_index]            <- sqrt(mean((as.matrix(data_clean[, -1] - crm_fit$data.imputed[, -1]))^2))
    results_RMSEI_j_sim$`CRM-LASSO`[case_pct_index]      <- sqrt(mean((as.matrix(data_clean[, -1] - crm_fit_l$data.imputed[, -1]))^2))
    results_RMSEI_j_sim$`CRM-ElasticNet`[case_pct_index] <- sqrt(mean((as.matrix(data_clean[, -1] - crm_fit_e$data.imputed[, -1]))^2))

    # Execution time - - - - - - - - - - - - - - -
    results_time_j_sim$Time[case_pct_index] <- crm_fit$time


    # Elapsed time - - - - - - - - - - - - - - - -
    t_end <- seconds_to_period(round((proc.time() - t_start)[3]))
    cat(paste0(" - ", round(100 * (case_pct_index + (j_sim - 1) * n_pct_case_out) /
                              (n_sims * n_pct_case_out), 2),
               "% - ", sprintf("%02d:%02d:%02d", t_end@hour, minute(t_end), second(t_end)), "\n\n"))


  } # end of forloop "case_pct_index in 1:n_pct_case_out"


  results_MAE[[j_sim]]       <- results_MAE_j_sim
  results_MSEP[[j_sim]]      <- results_MSEP_j_sim
  results_RMSEI[[j_sim]]     <- results_RMSEI_j_sim
  results_time[[j_sim]]      <- results_time_j_sim


} # end of forloop "j_sim in 1:n_sims"


t_end <- seconds_to_period(round((proc.time() - t_start)[3]))
cat(paste("\nTime elapsed:",
          sprintf("%02d:%02d:%02d", t_end@hour, minute(t_end), second(t_end)), "\n\n"))



# Study results -----------------------------------------------------------------------------------
results_MAE       <- cbind.data.frame(do.call(rbind, results_MAE),       pct_case_out = rep(pct_case_out_seq, n_sims))
results_MSEP      <- cbind.data.frame(do.call(rbind, results_MSEP),      pct_case_out = rep(pct_case_out_seq, n_sims))
results_RMSEI     <- cbind.data.frame(do.call(rbind, results_RMSEI),     pct_case_out = rep(pct_case_out_seq, n_sims))
results_time      <- cbind.data.frame(do.call(rbind, results_time),      pct_case_out = rep(pct_case_out_seq, n_sims))

results_MAE_mean       <- aggregate(. ~ pct_case_out, data = results_MAE,       FUN = mean)
results_MSEP_mean      <- aggregate(. ~ pct_case_out, data = results_MSEP,      FUN = mean)
results_RMSEI_mean     <- aggregate(. ~ pct_case_out, data = results_RMSEI,     FUN = mean)
results_time_mean      <- aggregate(. ~ pct_case_out, data = results_time,      FUN = mean)


# Mean Absolute Error - - - - - - - - - - - - - -
plot(100 * results_MAE_mean$pct_case_out, results_MAE_mean$`CRM`, type = "l", xlab = "% casewise outliers",
     ylab = "Average MAE", ylim = range(results_MAE_mean[, -1]), las = 1, cex.axis = 1.2, cex.lab = 1.3, lwd = 2, col = 2)
lines(100 * results_MAE_mean$pct_case_out, results_MAE_mean$`CRM-LASSO`,      lwd = 2, lty = 1, col = 3)
lines(100 * results_MAE_mean$pct_case_out, results_MAE_mean$`CRM-ElasticNet`,  lwd = 2, lty = 1, col = 4)
points(100 * results_MAE_mean$pct_case_out, results_MAE_mean$`CRM`, pch = 19, col = 2)
points(100 * results_MAE_mean$pct_case_out, results_MAE_mean$`CRM-LASSO`, pch = 22, col = 3)
points(100 * results_MAE_mean$pct_case_out, results_MAE_mean$`CRM-ElasticNet`, pch = 23, col = 4)
legend("topleft", legend = names(results_MAE_mean)[-1], col = 2:4, lwd = 2)


# Mean Squared Error of Prediction - - - - - - - -
plot(100 * results_MSEP_mean$pct_case_out, results_MSEP_mean$`CRM`, type = "l", xlab = "% casewise outliers",
     ylab = "Average MSEP", ylim = range(results_MSEP_mean[, -1]), las = 1, cex.axis = 1.2, cex.lab = 1.3, lwd = 2, col = 2)
lines(100 * results_MSEP_mean$pct_case_out, results_MSEP_mean$`CRM-LASSO`,      lwd = 2, lty = 1, col = 3)
lines(100 * results_MSEP_mean$pct_case_out, results_MSEP_mean$`CRM-ElasticNet`,  lwd = 2, lty = 1, col = 4)
points(100 * results_MSEP_mean$pct_case_out, results_MSEP_mean$`CRM`, pch = 19, col = 2)
points(100 * results_MSEP_mean$pct_case_out, results_MSEP_mean$`CRM-LASSO`, pch = 22, col = 3)
points(100 * results_MSEP_mean$pct_case_out, results_MSEP_mean$`CRM-ElasticNet`, pch = 23, col = 4)
legend("topleft", legend = names(results_MSEP_mean)[-1], col = 2:4, lwd = 2)


# Root Mean Squared Error of Imputation - - - - -
plot(100 * results_RMSEI_mean$pct_case_out, results_RMSEI_mean$CRM, type = "l", xlab = "% casewise outliers",
     ylab = "Average RMSEI", ylim = range(results_RMSEI_mean[, -1]), las = 1, cex.axis = 1.2, cex.lab = 1.3, lwd = 2, col = 2)
lines(100 * results_RMSEI_mean$pct_case_out, results_RMSEI_mean$`CRM-LASSO`, lwd = 2, lty = 1, col = 3)
lines(100 * results_RMSEI_mean$pct_case_out, results_RMSEI_mean$`CRM-ElasticNet`, lwd = 2, lty = 2, col = 4)
points(100 * results_RMSEI_mean$pct_case_out, results_RMSEI_mean$`CRM`, pch = 19, col = 2)
points(100 * results_RMSEI_mean$pct_case_out, results_RMSEI_mean$`CRM-LASSO`, pch = 22, col = 3)
points(100 * results_RMSEI_mean$pct_case_out, results_RMSEI_mean$`CRM-ElasticNet`, pch = 23, col = 4)
legend("topleft", legend = names(results_RMSEI_mean)[-1], col = 2:4, lwd = 2)

# Execution time - - - - - - - - - - - - - - - - -
cat("CRM average execution time:", round(mean(results_time$Time), 1), "seconds")


