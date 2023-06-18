spadimo_e <- function(data, weights, obs,
                    control = list(scaleFun  = Qn,
                                   nlatent   = 1,
                                   etas      = NULL,
                                   csqcritv  = 0.975,
                                   stopearly = FALSE,
                                   trace     = FALSE,
                                   plot      = TRUE))
{

  # SPADIMO : SPArse DIrections of Maximal Outlyingness
  # -----------------------------------------------------------------------------------------------------------------------
  # Inputs: data as data frame
  #         weights as numeric vector: case weights from robust estimator
  #         obs as integer: the case number under consideration
  #         control: a list with control parameters
  #                     nlatent: integer number of latent variables for sparse PLS regression (via SNIPLS) (default is 1)
  #                     etas: vector of decreasing sparsity parameters (NULL is default)
  #                     csqcritv: probability level for internal chi-squared quantile (used when n > p)
  #                     stopearly: if TRUE, method stops as soon as the reduced case is no longer outlying,
  #                                else it loops through all values of eta
  #                     trace: print intermediate results
  #                     plot: show heatmaps and graph of the results
  # -----------------------------------------------------------------------------------------------------------------------
  # Outputs: outlvars: vector containing individual variable names contributing most to obs's outlyingness
  #          outlvarslist: list of variables contributing to obs's outlyingness for different values of eta
  #          a as vector: sparse direction of maximal outlyingness
  #          alist: list of sparse directions of maximal outlyingness for different values of eta
  #          o.before: outlyingness of original case (n < p) or
  #                    PCA outlier flag (n >= p) before removing outlying variables
  #          o.after: outlyingness of reduced case (n > p) or
  #                   PCA outlier flag (n >= p) after removing outlying variables
  #          eta: cutoff where obs is no longer outlying
  #          time: time to execute the SPADIMO algorithm
  #          control: a list with control parameters that are used
  # -----------------------------------------------------------------------------------------------------------------------
  # Written by Sven Serneels, BASF Corp., 4/2016 - 7/2016 & Sebastiaan Höppner, KU Leuven, 4/2017 - 7/2017
  # -----------------------------------------------------------------------------------------------------------------------

  if (missing(data)) {
    stop("Argument 'data' is missing, with no default.")
  }
  if (missing(weights)) {
    stop("Argument 'weights' is missing, with no default.")
  }
  if (missing(obs)) {
    stop("Argument 'obs' is missing, with no default.")
  }
  if (missing(control)) {
    control <- list(scaleFun  = Qn,
                    nlatent   = 1,
                    etas      = NULL,
                    csqcritv  = 0.975,
                    stopearly = FALSE,
                    trace     = FALSE,
                    plot      = TRUE)
  }

  # (1) Initiate some values
  starttimer <- proc.time()
  obs <- as.integer(obs)
  x <- as.matrix(data)
  n <- nrow(x)
  p <- ncol(x)
  w <- weights

  if (is.null(control$scaleFun)) {
    control$scaleFun <- Qn
  }
  if (is.null(control$nlatent)) {
    control$nlatent <- 1
  }
  if (is.null(control$etas)) {
    if (n > p) {
      control$etas <- seq(0.9, 0.1, -0.05)
    } else if (n <= p) {
      control$etas <- seq(0.6, 0.1, -0.05)
    }
  }
  if (is.null(control$csqcritv)) {
    control$csqcritv <- 0.975
  }
  if (is.null(control$stopearly)) {
    control$stopearly <- FALSE
  }
  if (is.null(control$trace)) {
    control$trace <- FALSE
  }
  if (is.null(control$plot)) {
    control$plot <- TRUE
  }

  etas <- sort(control$etas, decreasing = TRUE)
  stopcrit <- FALSE
  a.list <- list()
  outlvars.list <- list()


  # (2) Robust standardization of data
  robLoc <- apply(x, 2, weighted.mean, w) #robLoc <- apply(x, 2, median)
  robScale <- apply(x, 2, control$scaleFun)
  if (any(robScale < .Machine$double.eps)) {
    stop("The following variables have 0 scale:\n ",
         paste(names(which(robScale < .Machine$double.eps)), collapse = ", "))
  }
  z <- scale(x, center = robLoc, scale = robScale)


  # (3) If  n > p: compute robust Mahalanobis distance
  #     If n <= p: apply robust PCA method (ROBPCA)
  outlyingness.before <- NA
  PCA.outlflag.before <- NA
  if (n > p) {
    Sigmaw <- (t(z) %*% diag(w) %*% z) / (sum(w)-1)
    outlyingness.before <- sqrt(t(z[obs, ]) %*% chol2inv(chol(Sigmaw)) %*% z[obs, ])
    if (control$trace) printer(type = 1, n = n, p = p, control = control, outlyingness.before = outlyingness.before)
  } else if (n <= p) {
    set.seed(2017)
    PCA.outlflag.before <- !(PcaHubert(z, alpha = 0.75, k = 0, kmax = 10, maxdir = 250, scale = FALSE)@flag[obs])
    if (control$trace) printer(type = 1, n = n, p = p, PCA.outlflag.before = PCA.outlflag.before)
  }


  # (4) Loop through set of etas to be screened and
  #     stop if chi-squared criterion is satisfied or case is not outlying according to ROBPCA
  for (i in 1:length(etas)) {

    # (4.1) Estimate sparse direction of maximal outlyingness and set of contributing variables
    reg <- spadimo.exs(Z = z,
                       w = w,
                       obs = obs,
                       nlatent = control$nlatent,
                       eta = etas[i])
    outlvars <- reg$outlvars
    a.list[[i]] <- reg$a
    outlvars.list[[i]] <- outlvars

    if (control$trace) printer(type = 2, etas = etas, i = i, reg = reg)

    if (length(outlvars) == p) {
      if (stopcrit == FALSE) {
        stopcrit <- TRUE
        outlvars.crit <- outlvars
        eta.crit <- reg$eta
        a.crit <- reg$a
        outlyingness.after.crit <- NA
        PCA.outlflag.after.crit <- NA
      }
      warning('All variables were flagged!', call. = FALSE)
      break
    }


    # (4.2) Remove contributing variables from the data set
    z.reduced <- as.matrix(z[,-outlvars])
    df <- p - length(outlvars) # No. remaining variables = ncol(z.reduced)
    csqcrit.after <- qchisq(control$csqcritv, df)


    # (4.3) If  n > p: compute robust Mahalanobis distance of reduced case (without outlying cells)
    #       If n <= p: apply robust PCA method (ROBPCA) on reduced case
    outlyingness.after <- 0
    PCA.outlflag.after <- FALSE
    if (n > p) {
      Sigmaw <- (t(z.reduced) %*% diag(w) %*% z.reduced) / (sum(w)-1)
      outlyingness.after <- sqrt(t(z.reduced[obs, ]) %*% chol2inv(chol(Sigmaw)) %*% z.reduced[obs, ])
      if (control$trace) printer(type = 3, n = n, p = p, control = control, outlyingness.after = outlyingness.after, csqcrit.after = csqcrit.after, df = df)
    } else if (n <= p) {
      set.seed(2017)
      PCA.outlflag.after <- !(PcaHubert(z.reduced, alpha = 0.75, k = 0, kmax = 10, maxdir = 250, scale = FALSE)@flag[obs])
      if (control$trace) printer(type = 3, n = n, p = p, PCA.outlflag.after = PCA.outlflag.after)
    }


    # (4.4) Check if reduced case is still outlying
    if (n > p & outlyingness.after^2 < csqcrit.after & stopcrit == FALSE) {
      stopcrit <- TRUE
      outlvars.crit <- outlvars
      eta.crit <- reg$eta
      a.crit <- reg$a
      outlyingness.after.crit <- outlyingness.after
      PCA.outlflag.after.crit <- NA
      if (control$stopearly) break
    } else if (n <= p & PCA.outlflag.after == FALSE & stopcrit == FALSE) {
      stopcrit <- TRUE
      outlvars.crit <- outlvars
      eta.crit <- reg$eta
      a.crit <- reg$a
      outlyingness.after.crit <- NA
      PCA.outlflag.after.crit <- PCA.outlflag.after
      if (control$stopearly) break
    } else if (i == length(etas) & stopcrit == FALSE) {
      warning('Algorithm did not converge; reduced case remains outlying.', call. = FALSE)
      outlvars.crit <- outlvars
      eta.crit <- reg$eta
      a.crit <- reg$a
      outlyingness.after.crit <- NA
      PCA.outlflag.after.crit <- PCA.outlflag.after
    }

  } #end of for loop over etas

  endtimer <- proc.time() - starttimer

  if (control$trace) printer(type = 4, endtimer = endtimer)

  if (control$plot) printer(type = 5, n = n, p = p, x = x, i = i, obs = obs, control = control, etas = etas, eta.crit = eta.crit,
                            outlvars.crit = outlvars.crit, rownamesData = rownames(data), colnamesData = colnames(data),
                            a.list = a.list, outlvars.list = outlvars.list)

  return(list(outlvars     = outlvars.crit,
              outlvarslist = outlvars.list,
              a            = a.crit,
              alist        = a.list,
              eta          = eta.crit,
              o.before     = ifelse(n > p, as.numeric(outlyingness.before), PCA.outlflag.before),
              o.after      = ifelse(n > p, as.numeric(outlyingness.after.crit), PCA.outlflag.after.crit),
              time         = round(endtimer[3], 4),
              control      = control))
}




spadimo_elasticnet <- function(data, weights, obs, control = list(etas = NULL, alpha = 0.5)) {
  x <- as.matrix(data)
  n <- nrow(x)
  p <- ncol(x)
  w <- weights

  if (is.null(control$etas)) {
    if (n > p) {
      control$etas <- seq(0.9, 0.1, -0.05)
    } else if (n <= p) {
      control$etas <- seq(0.6, 0.1, -0.05)
    }
  }

  etas <- sort(control$etas, decreasing = TRUE)
  alpha <- control$alpha

  outlvars_list <- list()
  a_list <- list()

  for (i in 1:length(etas)) {
    elasticnet_model <- glmnet(x, w, family = "gaussian", alpha = alpha, lambda = etas[i])
    outlvars <- which(coef(elasticnet_model) != 0)[-1]
    a <- coef(elasticnet_model)[-1]
    outlvars_list[[i]] <- outlvars
    a_list[[i]] <- a
  }

  outlvars_crit <- outlvars_list[[length(etas)]]
  eta_crit <- etas[length(etas)]
  a_crit <- a_list[[length(etas)]]

  return(list(outlvars = outlvars_crit, outlvarslist = outlvars_list, a = a_crit, alist = a_list,
              eta = eta_crit, control = control))
}

# Function for printing (intermediate) results
#----------------------------------------------
printer <- function(type,
                    x = NULL,
                    n = NULL,
                    p = NULL,
                    i = NULL,
                    df = NULL,
                    obs = NULL,
                    reg = NULL,
                    etas = NULL,
                    a.list = NULL,
                    control = NULL,
                    endtimer = NULL,
                    eta.crit = NULL,
                    rownamesData = NULL,
                    colnamesData = NULL,
                    csqcrit.after = NULL,
                    outlvars.crit = NULL,
                    outlvars.list = NULL,
                    outlyingness.after = NULL,
                    PCA.outlflag.after = NULL,
                    outlyingness.before = NULL,
                    PCA.outlflag.before = NULL,
                    nOutliers = NULL,
                    nOutliers.crit = NULL)
{
  # require(gplots)
  # require(ggplot2)

  # type = 1: if  n > p: print robust Mahalanobis distance and chi-squared quantile
  #           if n <= p: print outlier flag from robust PCA method (ROBPCA)
  if (type == 1) {
    if (n > p) {
      cat(paste0('\n --------------------------------------------------------------------------',
                 '\n squared outlyingness of original case = ', round(outlyingness.before^2, 4),
                 '\n qchisq(', control$csqcritv, ', df = ', p, ') = ', round(qchisq(control$csqcritv, p), 4)))
    } else if (n <= p) {
      cat(paste('\n --------------------------------------------------------------------------',
                '\n orignal case outlying according to ROBPCA :', PCA.outlflag.before))
    }
  }


  # type = 2: print the flagged variables and the corresponding value for eta
  if (type == 2) {
    cat(paste0('\n\n --------------------------------------------------------------------------',
               '\n eta = ', etas[i], ' (SNIPLS, nlatent = ', reg$nlatent, ')',
               '\n ', length(reg$outlvars), ' variables retained: ', '\n'))
    print(reg$outlvars)
  }


  # type = 3: if  n > p: print robust Mahalanobis distance of reduced case (without outlying cells) and chi-squared quantile
  #           if n <= p: print outlier flag from robust PCA method (ROBPCA) on reduced case
  if (type == 3) {
    if (n > p) {
      cat(paste0('\n squared outlyingness of reduced case = ', round(outlyingness.after^2, 4),
                 '\n qchisq(' , control$csqcritv, ', df = ', df, ') = ', round(csqcrit.after, 4)))
    } else if (n <= p) {
      cat(paste('\n reduced case outlying according to ROBPCA :', PCA.outlflag.after))
    }
  }


  # type = 4: print computation time
  if (type == 4) {
    cat(paste('\n\n computation time:', round(endtimer[3], 4), 'sec.\n'))

  }


  # type = 5: show
  #           - graph of No. flagged variables versus eta values and identify the point when reduced case is no longer outlying
  #           - heatmap of sparse directions of maximal outlyingness for different values of eta
  #           - heatmap of the observations in which the identified variables are indicated for different values of eta
  if (type == 5) {
    # Heatmaps
    if (i > 1) { # only draw heatmaps and screeplot if i > 1
      heatmap <- matrix(0, nrow = i, ncol = p)
      cellnotes.a <- matrix(0, nrow = i, ncol = p)
      cellnotes.obs <- matrix(0, nrow = i, ncol = p)
      for (j in 1:i) {
        heatmap[j, which(a.list[[j]] > 0)] <- 1
        heatmap[j, which(a.list[[j]] < 0)] <- -1
        if (p <= 50) {
          cellnotes.a[j, ] <- round(a.list[[j]], 2)
          cellnotes.obs[j, ] <- round(x[obs, ], 2)
        } else {
          cellnotes.a[j, ] <- rep(NA, p)
          cellnotes.obs[j, ] <- rep(NA, p)
        }
      }
      rownames(heatmap) <- round(etas[1:i], 2)
      if (is.null(colnamesData)) {
        colnames(heatmap) <- 1:p
      } else {
        colnames(heatmap) <- colnamesData
      }

      Color <- c()
      if(-1 %in% c(heatmap)) {
        Color <- c(Color, 'blue')
      }
      if (0 %in% c(heatmap)) {
        Color <- c(Color, 'lightgray')
      }
      if (1 %in% c(heatmap)) {
        Color <- c(Color, 'red')
      }

      if (p <= 50) {
        heatmap.2(heatmap,
                  scale        = 'none',
                  col          = Color,
                  cellnote     = cellnotes.a,
                  notecol      = ifelse(as.vector(t(heatmap[i:1, ])) != 0, 'white', 'black'),
                  xlab         = 'sparse direction of maximal outlyingness',
                  ylab         = expression(eta),
                  notecex      = 1.0,
                  Rowv         = FALSE,
                  Colv         = FALSE,
                  dendrogram   = 'none',
                  density.info = 'none',
                  trace        = 'none',
                  key          = FALSE,
                  margins      = c(5, 4.5),
                  lhei         = c(0.5, 12),
                  lwid         = c(0.5, 10),
                  colsep       = 1:p,
                  rowsep       = 1:i,
                  sepcolor     = 'white',
                  sepwidth     = c(0.01, 0.01))
      }

      if (p <= 500) {
        heatmap.2(heatmap,
                  scale        = 'none',
                  col          = Color,
                  cellnote     = cellnotes.obs,
                  notecol      = ifelse(as.vector(t(heatmap[i:1, ])) != 0, 'white', 'black'),
                  xlab         = 'case',
                  ylab         = expression(eta),
                  notecex      = 0.8,
                  Rowv         = FALSE,
                  Colv         = FALSE,
                  dendrogram   = 'none',
                  density.info = 'none',
                  trace        = 'none',
                  key          = FALSE,
                  margins      = c(5, 4.5),
                  lhei         = c(0.5, 12),
                  lwid         = c(0.5, 10),
                  colsep       = 1:p,
                  rowsep       = 1:i,
                  sepcolor     = 'white',
                  sepwidth     = c(0.01, 0.01))
      }


      # Screeplot
      screeplot <- ggplot(data = data.frame(etas = etas[1:i], nOutliers = sapply(outlvars.list, length)), aes(x = etas, y = nOutliers)) +
        geom_point(size = 3.5) +
        geom_line() +
        geom_point(data = data.frame(eta.crit = eta.crit, nOutliers.crit = length(outlvars.crit)), aes(x = eta.crit, y = nOutliers.crit),
                   shape = 18, size = 7, col = 'red') +
        labs(x = expression(eta),
             y = 'No. flagged variables') +
        theme(plot.title = element_text(size = 15),
              text = element_text(size = 20),
              axis.title = element_text(size = 15)) +
        ggtitle(substitute(paste(eta, ' = ', eta.crit, ' , ', nflags, ' variables flagged'),
                           list(eta.crit = eta.crit, nflags = length(outlvars.crit))))
      plot(screeplot)
      # screeplot <- ggplot(data = data.frame(etas = etas[1:i], nOutliers = sapply(outlvars.list, length)), aes(x = etas, y = nOutliers)) +
      #   geom_point(size = 4.2) +
      #   geom_line(size = 0.8) +
      #   geom_point(data = data.frame(eta.crit = eta.crit, nOutliers.crit = length(outlvars.crit)), aes(x = eta.crit, y = nOutliers.crit),
      #              shape = 18, size = 8, col = 'red') +
      #   coord_cartesian(ylim = c(0, 12)) +
      #   scale_x_continuous(breaks = seq(min(etas), max(etas), 0.1), labels = seq(min(etas), max(etas), 0.1)) +
      #   scale_y_continuous(breaks = seq(0, 12, 4), labels = seq(0, 12, 4)) +
      #   labs(x = expression(eta),
      #        y = 'no. flagged variables') +
      #   theme(plot.title = element_text(size = 22),
      #         text = element_text(size = 24),
      #         axis.title = element_text(size = 24)) +
      #   ggtitle(substitute(paste(eta, ' = ', eta.crit, ' , ', nflags, ' variable flagged'),
      #                      list(eta.crit = eta.crit, nflags = length(outlvars.crit))))
      # plot(screeplot)

    } # end of drawing heatmaps and screeplot

  } # end of type = 5

}
