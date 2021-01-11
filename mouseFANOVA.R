normdata <- function (x) {
  m <- mean(x)
  s <- sd(x)
  y <- (x-mean(x)) / s
  return(list(y=y, mean=m, sd=s))
}


summary_func <- function(x, col) { 
    c(mean=mean(x[[col]], na.rm=TRUE), 
      se=sd(x[[col]], na.rm=TRUE)/sqrt(sum(!is.na(x[[col]]))))
}


rotate_trajectory <- function(trajectory, initial_pos, final_pos) { 
  theta <- -atan((final_pos[2] - initial_pos[2]) / 
                 (final_pos[1] - initial_pos[1]))
  rotmat <- rbind(c(cos(theta), -sin(theta)), c(sin(theta), cos(theta)))

  a = rotmat %*% trajectory
  return(a)
}


deviance_matern <- function(b, X, residuals) { 
  require(GPfit)
  dev <- 0
  for (i in 1:dim(residuals)[1]) { 
    if (any(is.na(residuals[i,]))) { next } 
    dev <- dev + GP_deviance(b, X, residuals[i,], 
			     corr=list(type='matern', nu=5/2))
  }
  return(dev)
}

fit_beta <- function(out_t, residuals) { 
  X = out_t/ max(out_t)
  #opt <- optimize(err_beta, X=X, residuals=residuals, interval=c(0,10))
  opt <- optimize(deviance_matern, X=X, residuals=residuals, interval=c(0,10))
  return(opt$minimum)
}


collapse_to_cell <- function(dat) { 
  in_t <- dat$times
  in_t_norm <- normdata(in_t)

  n_data <- dim(dat$dat_tibble)[1]

  conditions <- levels(dat$dat_tibble$Condition)
  n_conditions <- length(conditions)

  subjects <- levels(dat$dat_tibble$Subject)
  n_subjects <- length(subjects)

  responses <- levels(dat$dat_tibble$Response)
  n_responses <- length(responses)

  all_subj <- rep(subjects[1:n_subjects], each=n_conditions * n_responses)
  all_cond <- rep(rep(conditions[1:n_conditions], each=n_responses),
                  n_subjects)
  all_resp <- rep(responses[1:n_responses], n_subjects * n_conditions)
  all_traj <- vector("list")
  all_aad <- c()
  all_mad <- c()
  all_mad_time <- c()
  all_xflips <- c()
  all_oflips <- c()
  for (sn in 1:n_subjects) { 
    for (cn in 1:n_conditions) { 
      for (rn in 1:n_responses) { 
        which_trials <- with(dat$dat_tibble, Response==responses[rn] &
                                             Condition==conditions[cn] & 
                                             Subject==subjects[sn])

        samples <- with(dat$dat_tibble, CleanTrajectory[which_trials])
                               
                                               
	                    		       
  	if (length(samples) > 0) { 
     	  x <- rep(in_t_norm$y, length(samples))
          y <- unlist(samples)
          y <- y[grepl('Perpendicular', names(y))]
          sp <- smooth.spline(x, y)$y
	} else { 
      	  x_safe <- rep(in_t_norm$y, each=1)
      	  sp <- rep(NA, length(in_t_norm$y))
        }
	all_traj <- append(all_traj, list(sp))
        all_aad <- c(all_aad, with(dat$dat_tibble, mean(AAD[which_trials])))
        all_mad <- c(all_mad, with(dat$dat_tibble, mean(MAD[which_trials])))
        all_mad_time <- c(all_mad_time, 
                          with(dat$dat_tibble, mean(MADtime[which_trials])))
        all_xflips <- c(all_xflips, 
                        with(dat$dat_tibble, mean(Xflips[which_trials])))
        all_oflips <- c(all_oflips, 
                        with(dat$dat_tibble, mean(Oflips[which_trials])))
                                   
      }
    }
  }

  all_subj <- factor(all_subj)
  all_cond <- factor(all_cond)
  all_resp <- factor(all_resp)
  cell_tibble <- tibble(Subject=all_subj, Condition=all_cond, 
			   Response=all_resp, AAD=all_aad, MAD=all_mad,
                           MADtime=all_mad_time, Xflips=all_xflips,
                           Oflips=all_oflips, Trajectory=all_traj)
  return(cell_tibble)
}





anova_model <- function(out_t, collapsed, factors, heldout=NA) { 
  if (!is.na(heldout)) { 
    for (x in heldout) { 
      collapsed <- subset(collapsed,Subject != x)
    }
  }

  n_subjects <- length(levels(collapsed$Subject))
  n_time_bins <- length(collapsed$Trajectory[[1]])
  n_samples <- dim(collapsed)[1]
  n_factors <- length(factors)


  gm <- apply(matrix(unlist(collapsed$Trajectory), n_samples, n_time_bins,
		     byrow=TRUE), 2, mean, na.rm=TRUE)

  # Main effects model
  alpha <- vector("list", n_factors)
  for (fn in 1:n_factors) { 
    f <- factors[fn]
    f_levels <- levels(collapsed[[f]])
    n_levels <- length(f_levels)
    alpha[[fn]] <- vector("list", n_levels)

    for (l in 1:n_levels) { 
      samples <- collapsed$Trajectory[collapsed[f]==f_levels[l]]

      gm_mat <- matrix(rep(gm, each=length(samples)), nrow=length(samples))
      samples_mat <- matrix(unlist(samples), nrow=length(samples), 
                            byrow=TRUE)
      alpha[[fn]][[l]] <- apply(samples_mat - gm_mat, 2, mean, na.rm=TRUE)
                                             
    }
  }

  f_levels_1 <- levels(collapsed[[factors[1]]])
  n_levels_1 <- length(f_levels_1)
  if (n_factors == 2) { 
    # Interaction model (note this only works for 2 factors!)
    f_levels_2 <- levels(collapsed[[factors[2]]])
    n_levels_2 <- length(f_levels_2)
    alpha_beta <- vector("list", n_levels_1)
    for (l_1 in 1:n_levels_1) { 
      alpha_beta[[l_1]] <- vector("list", n_levels_2)
      for (l_2 in 1:n_levels_2) { 
        samples <- collapsed$Trajectory[
                        collapsed[factors[1]]==f_levels_1[l_1] & 
  		        collapsed[factors[2]]==f_levels_2[l_2]]
        gm_mat <- matrix(rep(gm, each=length(samples)), 
                         nrow=length(samples))
        a1_mat <- matrix(rep(alpha[[1]][[l_1]], each=length(samples)), 
  		         nrow=length(samples))
        a2_mat <- matrix(rep(alpha[[2]][[l_2]], each=length(samples)), 
  		         nrow=length(samples))
        samples_mat <- matrix(unlist(samples), nrow=length(samples), 
                              byrow=TRUE)
        alpha_beta[[l_1]][[l_2]] <- 
          apply(samples_mat - a1_mat - a2_mat - gm_mat, 2, mean, na.rm=TRUE)
      }
    }
  }

  residuals <- matrix(NA, n_samples, n_time_bins)
  if (n_factors ==1) { 
    for (r in 1:n_samples) { 
      l_1 <- which(f_levels_1 == collapsed[[factors[1]]][r] )
      residuals[r,] <- collapsed$Trajectory[[r]] - alpha[[1]][[l_1]] - gm
    }
    return(list(gm=gm, alpha=alpha, residuals=residuals))
  } else { 
    for (r in 1:n_samples) { 
      l_1 <- which(f_levels_1 == collapsed[[factors[1]]][r] )
      l_2 <- which(f_levels_2 == collapsed[[factors[2]]][r] )
      residuals[r,] <- collapsed$Trajectory[[r]] - gm - (alpha[[1]][[l_1]] 
          					 + alpha[[2]][[l_2]]
          					 + alpha_beta[[l_1]][[l_2]])
    }
    return(list(gm=gm, alpha=alpha, alpha_beta=alpha_beta, 
                residuals=residuals))
  }
}


loocv_anova_md <- function(dat, collapsed) { 
  sampling_rate = 60;
  max_time = 40;
  in_t = seq(0, max_time, by = (1/sampling_rate))
  subjects <- levels(collapsed$Subject)
  n_subjects <- length(subjects)
  conditions <- levels(collapsed$Condition)
  n_conditions <- length(conditions)
  responses <- levels(collapsed$Response)
  n_responses <- length(responses)

  deviance_null <- rep(NA, n_subjects)
  deviance_condition <- rep(NA, n_subjects)
  deviance_choice <- rep(NA, n_subjects)
  deviance_main <- rep(NA, n_subjects)
  deviance_full <- rep(NA, n_subjects)

  for (sn in 1:n_subjects) { 
    anova_full <- anova_model(dat$times, collapsed, 
			c("Condition", "Response"), heldout=subjects[sn])
    b <- fit_beta(dat$times, anova_full$residuals)

    sdat <- subset(collapsed, Subject==subjects[sn])
    n_trials <- dim(sdat)[1]

    residual_null <- matrix(NA, n_trials, length(dat$times))
    residual_condition <- matrix(NA, n_trials, length(dat$times))
    residual_choice <- matrix(NA, n_trials, length(dat$times))
    residual_main <- matrix(NA, n_trials, length(dat$times))
    residual_full <- matrix(NA, n_trials, length(dat$times))

    for (tr in 1:dim(sdat)[1]) {
      cn <- which(conditions == sdat$Condition[[tr]])
      rn <- which(responses == sdat$Response[[tr]])
      y <- sdat$Trajectory[[tr]]
      residual_null[tr,] <- y - anova_full$gm
      residual_condition[tr,] <- y - anova_full$gm - anova_full$alpha[[1]][[cn]]
      residual_choice[tr,] <- y - anova_full$gm - anova_full$alpha[[2]][[rn]]
      residual_main[tr,] <- (y - anova_full$gm 
			       - anova_full$alpha[[1]][[cn]]
                               - anova_full$alpha[[2]][[rn]])
      residual_full[tr,] <- (y - anova_full$gm 
			       - anova_full$alpha[[1]][[cn]]
                               - anova_full$alpha[[2]][[rn]]
                               - anova_full$alpha_beta[[cn]][[rn]])


    }
    deviance_null[sn] <- deviance_matern(b, dat$times, residual_null)
    deviance_condition[sn] <- deviance_matern(b, dat$times, residual_condition)
    deviance_choice[sn] <- deviance_matern(b, dat$times, residual_choice)
    deviance_main[sn] <- deviance_matern(b, dat$times, residual_main)
    deviance_full[sn] <- deviance_matern(b, dat$times, residual_full)
  }
  return(list(null = sum(deviance_null),
	      condition = sum(deviance_condition),
	      choice = sum(deviance_choice),
	      main = sum(deviance_main),
	      full = sum(deviance_full)))
}

loocv_anova <- function(dat, collapsed) { 
  #sampling_rate = 60;
  #max_time = 40;
  #in_t = seq(0, max_time, by = (1/sampling_rate))
  #subjects <- levels(collapsed$Subject)
  #n_subjects <- length(subjects)

  #f <- rep(NA, n_subjects)
  #n <- rep(NA, n_subjects)

  sampling_rate = 60;
  max_time = 40;
  in_t = seq(0, max_time, by = (1/sampling_rate))
  subjects <- levels(collapsed$Subject)
  n_subjects <- length(subjects)
  conditions <- levels(collapsed$Condition)
  n_conditions <- length(conditions)
  responses <- levels(collapsed$Response)
  n_responses <- length(responses)

  deviance_null <- rep(NA, n_subjects)
  deviance_choice <- rep(NA, n_subjects)

  for (sn in 1:n_subjects) { 
    anova_full <- anova_model(dat$times, collapsed, 
			      c("Response"), heldout=subjects[sn])
    b <- fit_beta(dat$times, anova_full$residuals)

    sdat <- subset(collapsed, Subject==subjects[sn])
    n_trials <- dim(sdat)[1]

    residual_null <- matrix(NA, n_trials, length(dat$times))
    residual_choice <- matrix(NA, n_trials, length(dat$times))

    for (tr in 1:dim(sdat)[1]) {
      rn <- which(responses == sdat$Response[[tr]])
      y <- sdat$Trajectory[[tr]]
      residual_null[tr,] <- y - anova_full$gm
      residual_choice[tr,]<- y - anova_full$gm - anova_full$alpha[[1]][[rn]]

    }
    deviance_null[sn] <- deviance_matern(b, dat$times, residual_null)
    deviance_choice[sn] <- deviance_matern(b, dat$times, residual_choice)
  }
  return(list(null = sum(deviance_null),
	      choice = sum(deviance_choice)))
}
