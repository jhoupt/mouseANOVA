rm(list=ls())
require("parallel")
require("nlme")
require("geoR")
require("tibble")
require("nlme")
gitdir <- "~/Documents/GitHub/mouseANOVA/"
source(paste(gitdir, "mouseFANOVA.R", sep=""))
source(paste(gitdir, "SimulatedTrajectory.R", sep=""))
source(paste(gitdir, "mouseSimulation.R", sep=""))
options(mc.cores = parallel::detectCores()-1)
cl <- makeCluster(15, type="FORK")

fixed_parameters <- list(
  sampling_freq = 60,
  max_t = 3,
  initial_acceleration = 100,
  max_acceleration = 1000,
  initial_pos = c(0,0), 
  final_pos = c(180, 184),
  ttest_alpha = 0.05,
  safe.ctl = .01,
  safe.dev = .4,
  t_vec = c()
)
fixed_parameters$t_vec <- with(fixed_parameters, 
			       seq(0, max_t, by=1/sampling_freq))

experiment_parameters <- list(
  n_trials = 15,
  n_subjects = 20,
  risky.ctl = 0.01,
  risky.dev = .4
  #risky.dev = c(.425, .45, .475),
  #risky.ctl = c(.03, .04, .05)
)

#exportCluster(cl=cl, )

run_sim_experiment <- function(id=1, experiment_parameters, fixed_parameters, load_from_file=FALSE) { 
  set.seed(id)
  attach(experiment_parameters)
  subj_params = list(subj=NA, dev=NA, ctl=NA)

  if (load_from_file) { 
    pattern <- paste('ctl-', risky.ctl, '_dev-', risky.dev, sep="")
    fname <- list.files(pattern=pattern) 
    tryCatch({
      load(fname)
      dat_tibble <- X[[id]]$dat
    }, error = function(e) { print("Unable to load file."); load_from_file <- FALSE })
  }
    
  if (!load_from_file) { 
    dat_tibble <- NULL
    #all_tr <- vector("list", length=n_subjects)
    all_tr <- 1:(2*n_trials)
    all_tr_dat <- c()
    all_tr_mat <- c()
    for (s in 1:n_subjects) { 
      source(paste(gitdir, "mouseSimulation.R", sep=""))
      subj_params$risky.dev <- risky.dev
      subj_params$risky.ctl <- risky.ctl
      subj_params$n_trials <- n_trials
      subj_params$subj <- s

      sj <- run_sim_subj(s, subj_params, fixed_parameters)
      while(any(is.na(sj$RT))) { 
        sj <- run_sim_subj(s,subj_params, fixed_parameters)
      }

      dat_tibble <- rbind(dat_tibble, sj)
    }
    dat_tibble$Subject <- factor(dat_tibble$Subject)
    dat_tibble$Condition <- factor(dat_tibble$Condition)
  }

  dat <- list(times = fixed_parameters$t_vec, dat_tibble=dat_tibble)

  n_times <- length(fixed_parameters$t_vec)

  perpendicular <- matrix(NA, n_subjects * n_trials * 2, 101)
                          

  for (i in 1:length(dat_tibble$CleanTrajectory)) { 
    x <- dat_tibble$CleanTrajectory[[i]]$Perpendicular
    rtx <- which(fixed_parameters$t_vec==dat_tibble$RT[i])
    x <- smooth.spline(fixed_parameters$t_vec[1:rtx], x[1:rtx])
    perpendicular[i,] <- predict(x, seq(0,dat_tibble$RT[i], length.out=101))$y
  }

  ttest.result.l <- rep(NA, 100)
  for (tx in 2:101) { 
    dat_tibble$x <- perpendicular[,tx]
    #lines(density(x), col=tx)
    ttest.result.l[tx-1] <- tryCatch({
          anova(lme(x ~ Response, random=~1|Subject, dat_tibble))$"p-value"[2]

    }, error = function(e) { 1 })
  }
  ttest.result <- mean(ttest.result.l < fixed_parameters$ttest_alpha)

  mad.result <- anova(lme(MAD ~ Response, 
			  random = ~1|Subject, dat_tibble))$"p-value"[2]

  mad_t.result <- anova(lme(MADtime~ Response, 
			  random = ~1|Subject, dat_tibble))$"p-value"[2]

  aad.result <- anova(lme(AAD ~ Response, 
			  random = ~1|Subject, dat_tibble))$"p-value"[2]

  if (load_from_file) { 
    fanova.result <- X[[id]]$result['fanova']
  } else { 
    collapsed <- collapse_to_cell(dat)
    x <- loocv_anova(dat, collapsed)
    fanova.result <- x$choice - x$null
  }

  result <- c(ttest.result, mad.result, mad_t.result, aad.result, 
              fanova.result)
  names(result) <- c("sequential", "mad", "mad_t", "aad", "fanova")
  detach(experiment_parameters)
  return(list(dat=dat_tibble, result=result))
}


plot.sim <- function(ctl, dev) { 
  require(ggplot2)
  require(gridExtra)
  all_dat <- c()
  res <- c()
  for (ct in ctl) { 
  for (dv in dev) { 
    pattern <- paste('ctl-', ct, '_dev-', dv, "_", sep="")
    fname <- list.files(pattern=pattern) 
    tryCatch({
      load(fname)
      for ( i in 1:length(X) ) { 
        all_dat <- rbind(all_dat, c(ct, dv, X[[i]]$result))
      }
    }, error = function(e) { print(paste(pattern, "not found!")) })
  }
  }
  all_dat <- data.frame(all_dat)
  names(all_dat) <- c("Control_Noise", "Initial_Deviation", "t_test", "MAD",
                       "MAD_Time", "AAD", "fANOVA")
  all_dat$Initial_Deviation <- factor(all_dat$Initial_Deviation)
  all_dat$Control_Noise <- factor(all_dat$Control_Noise)

  p1 <- ggplot(subset(all_dat, Control_Noise==0.01), 
               aes(x=Initial_Deviation, y=fANOVA)) +
	           geom_violin(scale="width")
  p2 <- ggplot(subset(all_dat, Control_Noise==0.01), 
               aes(x=Initial_Deviation, y=AAD)) +
	  geom_violin(scale="width")
  p3 <- ggplot(subset(all_dat, Control_Noise==0.01), 
               aes(x=Initial_Deviation, y=MAD)) +
	  geom_violin(scale="width")
  p4 <- ggplot(subset(all_dat, Control_Noise==0.01), 
               aes(x=Initial_Deviation, y=MAD_Time)) +
	  geom_violin(scale="width")
  p5 <- ggplot(subset(all_dat, Control_Noise==0.01), 
               aes(x=Initial_Deviation, y=t_test)) +
	  geom_violin(scale="width")
  grid.arrange(p1, p2, p3, p4, p5, nrow=5)

  p1 <- ggplot(subset(all_dat, Initial_Deviation==0.4), 
               aes(x=Control_Noise, y=fANOVA)) +
	  geom_violin(scale="width")
  p2 <- ggplot(subset(all_dat, Initial_Deviation==0.4), 
               aes(x=Control_Noise, y=AAD)) +
	  geom_violin(scale="width")
  p3 <- ggplot(subset(all_dat, Initial_Deviation==0.4), 
               aes(x=Control_Noise, y=MAD)) +
	  geom_violin(scale="width")
  p4 <- ggplot(subset(all_dat, Initial_Deviation==0.4), 
               aes(x=Control_Noise, y=MAD_Time)) +
	  geom_violin(scale="width")
  p5 <- ggplot(subset(all_dat, Initial_Deviation==0.4), 
               aes(x=Control_Noise, y=t_test)) +
	  geom_violin(scale="width")
  grid.arrange(p1, p2, p3, p4, p5, nrow=5)
  return(all_dat)
}


plot.roc <- function(all_dat, param=c("Initial_Deviation", "Control_Noise")) { 
  test_names <- c("t_test", "MAD", "MAD_Time", "AAD", "fANOVA")
  setEPS()
  postscript(paste(param, "roc.eps", sep="_"), width=7.45, height = 2.5)
  par(mfrow=c(1,3), mar=c(4, 3, 3, 2)+.1)

  if (param == "Initial_Deviation") { 
    base_targ <- .4
    dat <- subset(all_dat, Control_Noise == 0.01)
    lvls = c(.425, .45, .475)
  } else { 
    base_targ <- .01
    dat <- subset(all_dat, Initial_Deviation == 0.4)
    lvls = c(.02, .03, .04)
  }
  dat$t_test <- 1-dat$t_test
  for( comp in lvls) { 
    plot(c(0,0), c(1,1), xlim=c(0,1), ylim=c(0,1), 
         type='n', xlab="", ylab="",
         main="")
    if(param == "Initial_Deviation") { 
      title(main=paste(comp, "ms"), line=.8)
    } else { 
      title(main=comp, line=.8)
    }
    mtext(side=1, text="False Positives", line=2.1, cex=.7)
    mtext(side=2, text="True Positives", line=2, cex=.7)
    abline(0,1,lty=2)

    j <- 1
    for (test.name in test_names) { 
      tmp_cr <- seq(.9*min(dat[,test.name]), 1.1*max(dat[,test.name]), length.out=300)
      tmp_x <- rep(NA, length(tmp_cr))
      tmp_y <- rep(NA, length(tmp_cr))
      for (i in 1:length(tmp_cr)) { 
        tmp_x[i] <- mean(dat[dat[param]==base_targ, test.name] < tmp_cr[i])
                         
        tmp_y[i] <- mean(dat[dat[param]==comp, test.name] <= tmp_cr[i])
      }
      tmp_roc <- splinefun(tmp_x, tmp_y, method="monoH.FC")
      #tmp_roc <- approxfun(tmp_x, tmp_y)
      curve(tmp_roc(x), 0,1, pch=j, col=j, add=TRUE)
      x.seq <- seq(0,1,length.out=8)
      points(x.seq, tmp_roc(x.seq), pch=j, col=j)
      j <- j + 1
    }
  }
  legend.names <- test_names
  legend.names[1] <- "t-Test"
  legend.names[3] <- "MAD Time"
  legend("bottomright", legend=legend.names, col=1:5, pch=1:5, cex=.9)
  dev.off()
}



#for (dev in c(.425, .45, .475)) { 
#  experiment_parameters$risky.dev <- dev
#  X <- parLapply(cl, 1:150, run_sim_experiment, experiment_parameters=experiment_parameters, fixed_parameters=fixed_parameters, load_from_file=FALSE)
#  with(experiment_parameters, 
#       save(X, file=paste("sim_ctl-", risky.ctl, "_dev-", risky.dev, 
#                          "_", date(), ".Rdata", sep="")))
#}

