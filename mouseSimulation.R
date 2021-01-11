initialize_sim <- function(id, simulation_parameters) {
  attach(simulation_parameters)
  param_list <- vector("list", length=n_subjects * length(risky.dev) * 
                                      length(risky.ctl))
  i <- 0
  for (dev in risky.dev) {
    for (ctl in risky.ctl) {
      for (subj in 1:n_subjects) { 
        i <- i + 1
        x <- list(sim=id, subj=subj, dev=dev, ctl=ctl)
        param_list[[i]] <- vector("list")
        param_list[[i]] <- x
      }
    }
  }
  detach(simulation_parameters)
  return(param_list)
}

run_sim_subj <- function(id, subj_params, fixed_parameters) { 
  attach(fixed_parameters)
  all_tr <- 1:(2*n_trials)
  risky.dev <- subj_params$risky.dev
  risky.ctl <- subj_params$risky.ctl
  n_trials <- subj_params$n_trials
  subj <- subj_params$subj
    
  safe <- rtrajectory(n_trials, initial_acceleration, max_acceleration, 
                      initial_pos, final_pos, sampling_freq, max_t,
                      control_noise = safe.ctl,
                      first_deviation = safe.dev)
  risky <- rtrajectory(n_trials, initial_acceleration, max_acceleration, 
                       initial_pos, final_pos, sampling_freq, max_t,
                       control_noise = risky.ctl,
                       first_deviation = risky.dev)

  chose_safe <- rep(c(TRUE, FALSE), each=n_trials)
  all_aad <- rep(NA, 2*n_trials)
  all_mad <- rep(NA, 2*n_trials)
  all_mad_time <- rep(NA, 2*n_trials)
  all_xflips <- rep(NA, 2*n_trials)
  all_oflips <- rep(NA, 2*n_trials)
  all_raw <- vector("list", 2*n_trials)
  all_clean <- vector("list", 2*n_trials)
  all_rt <- rep(NA, 2*n_trials)
  for (tr in 1:(2*n_trials)) { 
    if (tr <= n_trials) { 
      all_raw[[tr]] <- list(x=safe[tr,1,], y=safe[tr,2,])
      a <- rotate_trajectory(safe[tr,,], initial_pos, final_pos)
      all_rt[tr] <- t_vec[length(t_vec) - 
         max(which(rev(apply(apply(safe[tr,,], 1, diff) ==0, 1, any)[-1])))]
    } else { 
      tx <- tr - n_trials
      all_raw[[tr]] <- list(x=risky[tx,1,], y=risky[tx,2,])
      a <- rotate_trajectory(risky[tx,,], initial_pos, final_pos)
      all_rt[tr] <- t_vec[length(t_vec) - 
        max(which(rev(apply(apply(risky[tx,,], 1, diff) ==0, 1, any)[-1])))]
    }
    all_clean[[tr]] <- list(Direct=a[1,], Perpendicular=a[2,])

    all_aad[tr] <- mean(abs(a[2,]))
    all_mad[tr] <- max(abs(a[2,]))
    all_mad_time[tr] <- min(t_vec[abs(a[2,]) == all_mad[tr]])

    tmp <- all_raw[[tr]]$x
    all_xflips[tr] <- sum(diff(sign(tmp[tmp!=0])) !=0)
    all_oflips[tr] <- sum(diff(sign(tmp)) != 0) - 1

  }

  detach(fixed_parameters)
  all_responses <- factor(chose_safe, levels=c(FALSE, TRUE), 
                         labels=c("Risky", "Safe"))
  subj_tibble <- tibble(Subject=id, Condition="simulated", 
                        Trial=1:(2*n_trials),
                        Response=all_responses, AAD=all_aad, MAD=all_mad, 
                        MADtime=all_mad_time,
                        Xflips=all_xflips, Oflips=all_oflips, 
                        RT=all_rt,
                        RawTrajectory=all_raw, CleanTrajectory=all_clean)
  return(subj_tibble)
}
