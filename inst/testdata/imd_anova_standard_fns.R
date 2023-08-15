# ANOVA and G-test functions ---------------------------------------------------

# Function to construct the null space projection matrix.
proj_mat <- function(X, ngroups){
  #If the X matrix has atleast two rows, find projection matrix
  #into the null space corresponding to X
  if(!is.null(nrow(X))){
    Imat <- diag(1,nrow(X))
    
    # Zero out the first ngroups cols of X because we only want to remove the
    # effect of the covariates (which appear in the last n_covariates cols of
    # X).
    X[, 1:ngroups] <- 0
    
    Px <- MASS::ginv(t(X)%*%X)%*%t(X)
    
    Px <- X%*%Px
    return(Imat-Px)
  }

  #If only one datapoint is left then just return 1
  return(1)
}


# Project each row of the data.matrix into X's null space.
project_to_null <- function(data_mat, Xmatrix, ngroups){

  data_no_x <- data_mat

  for(i in 1:nrow(data_mat)){

    to_rem <- which(is.na(data_mat[i,]))
    if(length(to_rem)>0){
      roi <-  data_mat[i,-to_rem]
      IPxi <- proj_mat(Xmatrix[-to_rem,], ngroups)
      data_no_x[i,-to_rem] <- IPxi%*%matrix(roi,ncol=1)
    }else{
      roi <- data_mat[i,]
      IPxi <- proj_mat(Xmatrix, ngroups)
      data_no_x[i,] <- IPxi%*%matrix(roi,ncol=1)
    }

  }

  return(data.frame(data_no_x))

}

# This function and the following function, two_factor_anova_r, are used for two
# main effects (with or without covariates) when there is a significant
# interaction between the main effects for some of the biomolecules (rows of
# e_data). These functions are not used for examples when there is no
# significant interaction between the main effects, for any of the biomolecules,
# because the mean and variance computations are the same for every row.
run_two_factor <- function (data, gpData, red_df) {

  #Create design matrix for reduced model, i.e., no interaction effect
  colnames(gpData)[-c(1,2)] <- c("Factor1","Factor2")
  gpData <- cbind(gpData,y=1:nrow(gpData))

  #Create design matrix for the full model, i.e., all first order and
  #interaction effects
  Xred <- unname(model.matrix(lm(y ~ Factor1 + Factor2 - 1, data = gpData)))
  Xfull <- unname(model.matrix(lm(y ~ Factor1*Factor2 - 1, data = gpData)))

  res <- two_factor_anova_r(
    y = data,
    X_full = Xfull,
    X_red = Xred,
    red_df = red_df,
    group_ids = as.numeric(factor(gpData$Group,
                                  levels = unique(gpData$Group)))
  )

  #Get the unique group levels to translate ANOVA parameters into group means
  red_gpData <- dplyr::distinct(gpData, Group, .keep_all=TRUE)
  red_gpData <- dplyr::arrange(red_gpData, y)

  res$par_estimates <- res$par_estimates[, 1:length(red_gpData$Group)]
  colnames(res$par_estimates) <- red_gpData$Group

  return (res)

}

two_factor_anova_r <- function (y,
                                X_full,
                                X_red,
                                red_df,
                                group_ids) {

  n <- nrow(y)
  p_red <- ncol(X_red)
  p_full <- ncol(X_full)
  pval <- vector(mode = "numeric", length = n)
  Fstat <- vector(mode = "numeric", length = n)
  sig_est <- vector(mode = "numeric", length = n)
  par_ests <- vector(mode = "numeric", length = p_full)

  parmat <- matrix(nrow = n,
                   ncol = p_full)
  group_sizes <- matrix(nrow = n,
                        ncol = p_full)

  #Loop over rows in y
  for (i in 1:n) {

    yrowi <- y[i, ]
    to_remove <- which(!is.finite(yrowi))

    yrowi_nona <- yrowi
    X_red_nona <- X_red
    X_full_nona <- X_full
    num_to_remove <- length(to_remove)
    group_ids_nona <- group_ids # -1 #Make group_ids zero indexed
    # Indices in R start at 1. There is no need to subtract 1 from group_ids to
    # correctly subset rows/columns by group.

    if (num_to_remove > 0) {

      # Subset matrices by the to_remove object.
      yrowi_nona <- yrowi_nona[-to_remove]
      X_red_nona <- X_red_nona[-to_remove, ]
      X_full_nona <- X_full_nona[-to_remove, ]
      group_ids_nona <- group_ids_nona[-to_remove]

    }

    # Remove completely empty columns.
    csums <- colSums(X_full_nona)
    zero_cols <- which(csums == 0)
    non_zero_cols <- X_full_nona[, -zero_cols]

    #Subtract off df spent elsewhere, e.g., on covariates
    df_red <- (nrow(X_red_nona) -
                 Matrix::rankMatrix(X_red_nona)[[1]] -
                 red_df[i, ])
    df_full <- (nrow(X_full_nona) -
                  Matrix::rankMatrix(X_full_nona)[[1]] -
                  red_df[i, ])

    PxRed <- X_red_nona %*% MASS::ginv(
      t(X_red_nona) %*% X_red_nona
    ) %*% t(X_red_nona)

    # Create an identity matrix with the same dimension as PxRed
    diag_mat <- diag(nrow(PxRed))

    sigma2_red <- yrowi_nona %*% (diag_mat - PxRed) %*% yrowi_nona / df_red


    if ((df_red - df_full) <= 0) {

      #If interaction can't be estimated, automatically select smaller model
      pval[[i]] <- 1
      Fstat[[i]] <- 0

    } else {

      PxFull <- (X_full_nona %*% MASS::ginv(t(X_full_nona) %*% X_full_nona) %*%
                   t(X_full_nona))

      # Create an identity matrix with the same dimension as PxFull
      diag_mat <- diag(nrow(PxFull))

      sigma2_full <- (yrowi_nona %*% (diag_mat - PxFull) %*% yrowi_nona /
                        df_full)

      Fstat[[i]] <- (sigma2_red  * df_red - sigma2_full * df_full) / sigma2_full
      pval[[i]] <- pf(q = Fstat[[i]],
                      df1 = df_red - df_full,
                      df2 = df_full,
                      lower.tail = FALSE,
                      log.p = FALSE)

      if (!is.finite(pval[[i]])) {

        pval[[i]] <- 1

      }

    }

    if (pval[[i]] < 0.05) {

      #Reject null hypothesis that reduced model is good enough, use full model
      XFinal <- X_full_nona
      sig_est[[i]] <- sigma2_full

    } else {

      XFinal <- X_full_nona

      # Zero out the interaction terms if they're insignificant.
      XFinal[, (p_red + 1):p_full] <- 0
      sig_est[[i]] <- sigma2_red

    }

    #"Parameter estimates" are the group means: Xbeta_hat=X(XpX)^{-1}XpY
    par_ests_temp <- (XFinal %*% MASS::ginv(t(XFinal) %*% XFinal) %*%
                        t(XFinal) %*% yrowi_nona)

    #Find groups that had at least one non-missing value
    group_ids_nona_unq <- unique(group_ids_nona)

    # Fill in the par_ests vector, put NaN if group was missing or the average
    # effect in groups with no missing data
    if (length(group_ids_nona_unq) < p_full) {

      for (j in 1:p_full) {

        missing_gp <- which(group_ids_nona_unq == j)

        if (length(missing_gp) > 0) {

          par_ests[[j]] <- mean(par_ests_temp[which(group_ids_nona == j)])

        } else {

          par_ests[[j]] <- NA

        }

      }

    } else {

      for (k in 1:p_full) {

        par_ests[[k]] <- mean(par_ests_temp[which(group_ids_nona == k)])

      }

    }

    parmat[i, ] <- par_ests

    gsizes <- rep(0, p_full)

    # Compute group sizes after accounting for NaNs
    for (j in 1:p_full) {

      size_j <- which(group_ids_nona == j)

      if (length(size_j) > 0) {

        gsizes[[j]] <- length(size_j)

      }

    }

    group_sizes[i, ] <- gsizes # For now don't return the interaction groups

  }

  return (list(par_estimates = parmat,
               group_sizes = group_sizes,
               Sigma2 = sig_est,
               Fstats = Fstat,
               pvalue = pval))

}

a <- function (data, groups) {

  # nas <- !is.na(data)

  nova <- vector(length = nrow(data))

  for (e in 1:nrow(data)) {

    # nova[[e]] <- tryCatch (
    #   summary(
    #     aov(as.numeric(data[e, nas[e, ]]) ~ groups[nas[e, ]])
    #   )[[1]]$`Pr(>F)`[[1]],
    #   error = function (e) {NaN}
    # )

    nova[[e]] <- tryCatch (
      anova(
        lm(as.numeric(data[e, ]) ~ groups)
      )$`Pr(>F)`[[1]],
      error = function (e) {NaN}
    )

  }

  return (nova)

}

# A function for calculating the g-test statistic for two groups.
g <- function (obs1, obs2, abs1, abs2, n1, n2, n_total) {

  samurai_obs <- obs1 * log(obs1 / ((obs1 + obs2) * n1 / n_total))
  samurai_obs[is.na(samurai_obs)] <- 0
  samurai_abs <- abs1 * log(abs1 / ((abs1 + abs2) * n1 / n_total))
  samurai_abs[is.na(samurai_abs)] <- 0

  samurai <- samurai_obs + samurai_abs

  katana_obs <- obs2 * log(obs2 / ((obs1 + obs2) * n2 / n_total))
  katana_obs[is.na(katana_obs)] <- 0
  katana_abs <- abs2 * log(abs2 / ((abs1 + abs2) * n2 / n_total))
  katana_abs[is.na(katana_abs)] <- 0

  katana <- katana_obs + katana_abs

  fierce <- unname(samurai + katana)

  return (2 * fierce)

}

gflag <- function (obs1, obs2, abs1, abs2, pvals, cutoff) {

  maiden <- obs1 / (obs1 + abs1)

  dragon <- obs2 / (obs2 + abs2)

  knight <- sign(maiden - dragon)

  # If any p-values are NA change them to 1.
  pvals[is.na(pvals)] <- 1

  knight[which(pvals >= cutoff)] <- 0

  return (unname(knight))

}

aflag <- function (grp1, grp2, pvals, cutoff) {

  ninja <- rep(0, length(pvals))

  star <- which(pvals < cutoff)

  ninja[star] <- sign(grp1 - grp2)[star]

  return (ninja)

}

aflag_diff <- function (diff, pvals, cutoff) {

  sorceress <- rep(0, length(pvals))

  spell <- which(pvals < cutoff)

  sorceress[spell] <- sign(diff)[spell]

  return (sorceress)

}
