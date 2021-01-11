require("plyr")
require("tibble")
require("cowplot")
require("BayesFactor")
require("ggplot2")
require("xlsx")

#setwd("~/Documents/Research/ResponseDynamicAnalysis/fANOVA/Analysis")
gitdir <- "~/Documents/GitHub/mouseANOVA/"
source(paste(gitdir, "mouseFANOVA.R", sep=""))


loadMultipleDevicesSubject_fill <- function(n_time_bins) { 
    data_path = '~/Documents/Research/ResponseDynamicAnalysis/fANOVA/Data/MultipleDevices/'
    dat_file <-paste(data_path, "alldat.Rdata", sep="")
    if (file.exists(dat_file)) { 
      fname <- load(dat_file)
      return(dat)
    } else { 
      files = dir(path=data_path, 
		  pattern='S[[:digit:]]+_O[[:digit:]]+_[[:alpha:]]+.csv')
      #files <- files[files != "8.xlsx"]
      
      # Desired Times
      max_time = 40;
      sampling_rate = 60;
      time_bins = rep(0,n_time_bins);
      tmp = 3*n_time_bins;
      while(length(time_bins) > (n_time_bins - 1)) {
          time_bins = unique(10^(seq(0,log10(800), length.out=tmp)))
          tmp = tmp-1;
      }
      time_bins = c(time_bins, floor(max_time*sampling_rate));
      out_t = time_bins / sampling_rate;
      
      subj_re <- regexpr("S[[:digit:]]+_", files)
      attr(subj_re, "match.length") <- attr(subj_re, "match.length") -1
      subjects <- unique(substr(files, 2, attr(subj_re, "match.length")))
      n_subjects <- length(subjects)

      cond_re <- regexpr("_[[:alpha:]]+.csv", files)
      attr(cond_re, "match.length") <- attr(cond_re, "match.length") - 5
      conditions <- unique(substr(files, cond_re + 1, cond_re + 
				             attr(cond_re, "match.length")))
      n_conditions <- length(conditions)

      n_files <- length(files)
      dat_tibble <- NULL

      for (fn in 1:n_files) { 
          print(paste("File ", fn, " (", files[fn], ")", sep=""))
          file_name <- paste(data_path, files[fn], sep='')
          dat_tibble <- rbind(dat_tibble, read.subject.md(file_name, out_t))
      }

      dat_tibble$Subject <- factor(dat_tibble$Subject)
      dat_tibble$Condition <- factor(dat_tibble$Condition)
       
      dat <- list(times=out_t, dat_tibble=dat_tibble)
      save(dat, file=dat_file)
      return(dat)
    }
}


read.subject.md <- function(file_name, out_t) { 
  raw <- read.table(file_name, sep=",", header=TRUE, fill=TRUE)
  max_samps <- (dim(raw)[2] - 11) / 2
  n_trials <- dim(raw)[1]

  pos <- raw[,12:dim(raw)[2]]

  # Most complicated way to to this? 
  subj_re <- regexpr("S[[:digit:]]+_", file_name)
  attr(subj_re, "match.length") <- attr(subj_re, "match.length") - 2
  subject <- substr(file_name, subj_re, 
		    subj_re + attr(subj_re, "match.length"))
                  
  cond_re <- regexpr("_[[:alpha:]]+.csv", file_name)
  attr(cond_re, "match.length") <- attr(cond_re, "match.length") - 5
  condition <- substr(file_name, cond_re + 1, cond_re + 
			         attr(cond_re, "match.length"))

                  

  sampling_rate = 60;
  max_time = 60;
  in_t = seq(0, max_time, by = (1/sampling_rate))

  targLoc1 = c(140, 100);
  targLoc2 = c(500, 100);
  start_loc = c(320,384);
  
  targLoc1 = targLoc1 - start_loc;
  targLoc2 = targLoc2 - start_loc;
  
  theta = -acos(targLoc2 %*% c(1,0) / sqrt(targLoc2[1]^2 + targLoc2[2]^2) )
  #theta = -acos(targLoc1 %*% targLoc2/(sqrt(targLoc1[1]^2 + targLoc1[2]^2) *
  #                                     sqrt(targLoc2[1]^2 + targLoc2[2]^2)))
  rotmat = rbind(c(cos(theta), -sin(theta)), c(sin(theta), cos(theta)))


  in_t_norm <- normdata(in_t)
  out_t_norm <- (out_t - in_t_norm$mean) / in_t_norm$sd

  x_idx <- 2 * 1:(max_samps) -1
  X <- pos[,x_idx]
  Y <- pos[,x_idx+1]

  X <- X - start_loc[1]
  Y <- start_loc[2] - Y # Assumes flipped Y in data (e.g., screen)

  trial_n <- raw[,1]

  trials <-  unique(trial_n)
  t_vec <- 1:dim(X)[2]


  prob_left <- as.numeric(substr(raw[,4], 1, 2)) / 100
  prob_right <- as.numeric(substr(raw[,6], 1, 2)) / 100
  choice <- raw[,10]


  # 1: Left column is safer; 2: right is safer; 0: equal probability
  which_safe <- (prob_right < prob_left) + 2 * (prob_left < prob_right)
  chose_safe <- (which_safe == choice)

  all_tr <- 1:length(trials)

  all_raw <- vector("list", length=length(trials))
  all_clean <- vector("list", length=length(trials))
  all_risky <- c()
  all_safe <- c()
  all_aad <- c()
  all_mad <- c()
  all_mad_time <- c()
  all_xflips <- c()
  all_oflips <- c()

  for (tr in all_tr) { 
    #trial_dat <- c()
    #trial_dat$trial_n <- trial_n[tr]
    #trial_dat$chose_safe <- chose_safe[tr]

    trial_t <- t_vec[!is.na(X[tr,])]

    # If chose left, reflect to right
    Xtr <- unname(unlist(X[tr,trial_t]))
    if (choice[tr]==1) {
        Xtr = -Xtr
    }
    Ytr <- unname(unlist(Y[tr,trial_t]))

    all_raw[[tr]] <- list(x=Xtr, y=Ytr)

    # a[1,:] = Distance along direct path over time
    # a[2,:] = Deviations from direct path over time
    a = rotmat %*% rbind( Xtr, Ytr )
    #trial_dat$rotated_pos <- list(x=a[1,], y=a[2,])

    all_aad[tr] <- mean(abs(a[2,]))
    all_mad[tr] <- max(abs(a[2,]))
    all_mad_time[tr] <- min(trial_t[abs(a[2,])==all_mad[tr]])
    tmp <- diff(Xtr)
    all_xflips[tr] <- sum(diff(sign(tmp[tmp!=0])) !=0)
    all_oflips[tr] <- sum(diff(sign(Xtr)) != 0) - 1
    if(Xtr[1] == 1) {all_oflips[tr] <- all_oflips[tr] + 1}

    # Take the values of the deviations from path and pad with zeros
    tmp <- max(which(!is.na(unlist(X[tr,]))))
    pad_length <- sum(out_t > in_t[tmp+1])
    a1 <- c(a[1,], rep(0, pad_length))
    a2 <- c(a[2,], rep(0, pad_length))

    in_t_tr_norm <- c(in_t_norm$y[1:tmp], out_t_norm[out_t > in_t[tmp+1]])

    a1_smooth = smooth.spline(in_t_tr_norm,a1)
    a1_smooth = predict(a1_smooth, out_t_norm)$y

    a2_smooth = smooth.spline(in_t_tr_norm,a2)
    a2_smooth = predict(a2_smooth, out_t_norm)$y

    all_clean[[tr]] <- list(Direct=a1_smooth, Perpendicular=a2_smooth)
  }

  all_responses <- factor(chose_safe, levels=c(FALSE, TRUE), 
                         labels=c("Risky", "Safe"))
  subj_tibble <- tibble(Subject=subject, Condition=condition, Trial=all_tr,
                        Response=all_responses, AAD=all_aad, MAD=all_mad, 
                        MADtime=all_mad_time,
                        Xflips=all_xflips, Oflips=all_oflips, 
                        RawTrajectory=all_raw, CleanTrajectory=all_clean)
  return(subj_tibble)
}


analyze.md <- function() { 
  n_time_bins <- 100
  dat <- loadMultipleDevicesSubject_fill(n_time_bins)
  collapsed <- collapse_to_cell(dat)

  dat_frame <- data.frame(dat$dat_tibble[,1:9])
  dat_frame_c <- data.frame(collapsed[,1:8])


  ## BF ANOVA
  if (FALSE) { 
    madt_aov <- anovaBF(MADtime ~ Condition * Response + Subject, 
                        data=dat_frame, whichRandom="Subject")

    mad_aov <- anovaBF(MAD ~ Condition * Response + Subject, 
                        data=dat_frame, whichRandom="Subject")

    aad_aov <- anovaBF(AAD ~ Condition * Response + Subject, 
                        data=dat_frame, whichRandom="Subject")
  } else { 
    madt_aov <- anovaBF(MADtime ~ Condition * Response + Subject, 
                        data=subset(dat_frame_c, !is.na(AAD)), 
                        whichRandom="Subject")
    mad_aov <- anovaBF(MAD ~ Condition * Response + Subject, 
                        data=subset(dat_frame_c, !is.na(AAD)), 
                        whichRandom="Subject")
    aad_aov <- anovaBF(AAD ~ Condition * Response + Subject, 
                        data=subset(dat_frame_c, !is.na(AAD)), 
                        whichRandom="Subject")
    xflips_aov <- anovaBF(Xflips ~ Condition * Response + Subject, 
                          data=subset(dat_frame_c, !is.na(AAD)), 
                          whichRandom="Subject")
    oflips_aov <- anovaBF(Oflips ~ Condition * Response + Subject, 
                          data=subset(dat_frame_c, !is.na(AAD)), 
                          whichRandom="Subject")


    dvnames <- c("AAD", "MAD", "MADtime", "Xflips", "Oflips")
    dvlabs <- c("Average Absolute Deviation", 
                "Mean Absolute Deviation", 
                "MAD Time", "Horizontal Direction Changes", 
                "Origin Crossings")
    ylabs <- c(rep("Deviation (pixels)", 2), "Time (ms)", 
               rep("Count", 2))
    plots <- vector("list", length(dvnames))
    for (dn in 1:length(dvnames)) { 
      dv <- dvnames[dn]
      sumdat <- ddply(dat_frame_c, c("Condition", "Response"), 
                       .fun=summary_func, dv)
      #plots[[dn]] <- ggplot(sumdat, 
      #          aes(x=Condition, y=mean, fill=Response)) + 
      #          geom_bar(position=position_dodge(), stat="identity") + 
      #          geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
      #                        width=.2,
      #                        position=position_dodge(.9)) + 
      plots[[dn]] <- ggplot(dat_frame_c, 
                aes(x=Condition, y=.data[[dv]], fill=Response)) + 
                #geom_violin() + 
                geom_boxplot() +
                ylim(0,max(dat_frame_c[[dv]])) +
                labs(title=dvlabs[dn], y=ylabs[dn], 
                     breaks=c("risky", "safe"),
                     labels=c("Risky", "Safe"))
      if (dn == length(dvnames)) { 
        legend <- get_legend(plots[[dn]])
      }
      plots[[dn]] <- plots[[dn]] + theme(legend.position="none")
    }
    toprow <- plot_grid(plots[[1]], plots[[2]], plots[[3]], 
                        align='v', axis='l', nrow=1)
    bottomrow <- plot_grid(plots[[4]], plots[[5]], legend, nrow=1)
    setEPS()
    postscript("multiple_devices_standard_metrics.eps", width=9.31, 
                height=6.2)
    plot_grid(toprow, bottomrow, ncol=1)
    dev.off()
  }

  x <- loocv_anova_md(dat, collapsed)
  anova_full <- anova_model(dat$times, collapsed, c("Condition","Response"))
  sd_t <- sqrt(apply((anova_full$residuals)^2, 2, mean, na.rm=TRUE)) / 
            sqrt(dim(anova_full$residuals)[1])
  

  df_mu <- data.frame(Time=dat$times, Trajectory= anova_full$gm)
  df_dev <- data.frame(Time=rep(dat$times,3), 
                       Condition=rep(c("Mouse", "Trackpad", "Wii"), each=n_time_bins),
                       Trajectory= with(anova_full, 
                                        c(alpha[[1]][[1]], alpha[[1]][[2]], alpha[[1]][[3]])))
                                         
  df_choice <- data.frame(Time=rep(dat$times,2), 
                       Response=rep(c("Risky", "Safe"), each=n_time_bins),
                       Trajectory= with(anova_full, 
                                        c(alpha[[2]][[1]], alpha[[2]][[2]])))
                                         
  df_interaction <- data.frame(Time=rep(dat$times, 3*2), 
        Condition=rep(c("Mouse", "Trackpad", "Wii"), each=n_time_bins * 2),
        Response=rep(rep(c("Risky","Safe"), each=n_time_bins),3), 
        Trajectory=with(anova_full, c(alpha_beta[[1]][[1]], alpha_beta[[1]][[2]],
                                      alpha_beta[[2]][[1]], alpha_beta[[2]][[2]],
                                      alpha_beta[[3]][[1]], alpha_beta[[3]][[2]])))
              
  df_interaction$CxR <- paste(df_interaction$Condition, df_interaction$Response, sep="-")

  n_residuals <- dim(anova_full$residuals)[1]
  df_residual <- data.frame(Time=rep(dat$times,n_residuals), 
                            Trajectory = c(t(anova_full$residuals)),
                            idx=rep(1:n_residuals, each=n_time_bins))

  df_full <- data.frame(Time=rep(dat$times, 3*2), 
        Condition=rep(c("Mouse", "Trackpad", "Wii"), each=n_time_bins * 2),
        Response=rep(rep(c("Risky","Safe"), each=n_time_bins),3), 
        Deviation=with(anova_full, 
            c(gm + alpha[[1]][[1]] + alpha[[2]][[1]] + alpha_beta[[1]][[1]],
              gm + alpha[[1]][[1]] + alpha[[2]][[2]] + alpha_beta[[1]][[2]],
              gm + alpha[[1]][[2]] + alpha[[2]][[1]] + alpha_beta[[2]][[1]],
              gm + alpha[[1]][[2]] + alpha[[2]][[2]] + alpha_beta[[2]][[2]],
              gm + alpha[[1]][[3]] + alpha[[2]][[1]] + alpha_beta[[3]][[1]],
              gm + alpha[[1]][[3]] + alpha[[2]][[2]] + alpha_beta[[3]][[2]]
            )))
  df_full$Deviation.lwr <- df_full$Deviation - rep(sd_t, 2)
  df_full$Deviation.upr <- df_full$Deviation + rep(sd_t, 2)
  hues = seq(15, 375, length = 2 + 1)
  fill_cols <- hcl(h = hues, l = 80, c = 40)[1:2]
  line_cols <- hcl(h = hues, l = 65, c = 100)[1:2]
  line_types <- c("solid", "dashed", "dotted")

  
  plots_gp <- vector("list", 5)
  for (cn in 1:3) { 
    plots_gp[[cn]] <- ggplot(data=subset(df_full,Condition==conditions[cn]),
                             aes(x=Time, y=Deviation, 
                             group=Response)) + 
                      geom_ribbon(aes(ymin=Deviation.lwr, 
                                      ymax=Deviation.upr, fill=Response)) +
                      scale_fill_manual(values=hclvals) + 
                      geom_line(aes(color=Response),
                                    linetype=line_types[cn]) + 
                      xlim(0,12) + 
                      labs(title=conditions[[cn]], y="Deviation (pixels)",
                           breaks=c("risky", "safe"),
                           labels=c("Risky", "Safe"))

  }
  postscript("multiple_devices_gp_condition.eps", width=9.31, height=9.1)
  plot_grid(plots_gp[[1]], plots_gp[[2]], plots_gp[[3]], ncol=1)
  dev.off()
  responses <- levels(df_full$Response)
  for (rn in 1:2) { 
    plots_gp[[3+rn]] <- ggplot(data=subset(df_full,Response==responses[rn]),
                             aes(x=Time, y=Deviation, 
                             group=Condition)) + 
                      geom_ribbon(aes(ymin=Deviation.lwr, 
                                      ymax=Deviation.upr), 
                                      fill=fill_cols[rn]) +
                      #scale_fill_manual(values=hclvals) + 
                      geom_line(aes(linetype=Condition), 
                                color=line_cols[rn]) + 
                      xlim(0,12) + 
                      labs(title=responses[[rn]], y="Deviation (pixels)")
  }
  postscript("multiple_devices_gp_response.eps", width=9.31, height=6.1)
  plot_grid(plots_gp[[4]], plots_gp[[5]], ncol=1)
  dev.off()




  plots_aov <- vector("list", 5)
  plots_aov[[1]] <- ggplot(data=df_mu, aes(x=Time, y=Trajectory)) + geom_line() + 
                    labs(title="Grand Mean") +
                    xlim(0,12) + theme(legend.position="none")

  plots_aov[[2]] <- ggplot(data=df_dev, aes(x=Time, y=Trajectory)) + 
                    geom_line(aes(linetype=Condition)) + 
                    labs(title="Condition") +
                    xlim(0,12) + ylim(-22, 22) + theme(legend.position="none")
  plots_aov[[3]] <- ggplot(data=df_choice, aes(x=Time, y=Trajectory)) + 
                    geom_line(aes(color=Response)) + 
                    labs(title="Choice") +
                    xlim(0,12) + ylim(-22, 22) + theme(legend.position="none")
  plots_aov[[4]] <- ggplot(data=df_interaction, aes(x=Time, y=Trajectory, group=CxR)) + 
                    geom_line(aes(color=Response, linetype=Condition)) + 
                    labs(title="Interaction") +
                    xlim(0,12) + ylim(-22, 22) + theme(legend.position="none")
  plots_aov[[5]] <- ggplot(data=df_residual, aes(x=Time, y=Trajectory, group=idx)) + 
                    geom_line(aes(color=idx)) +scale_fill_grey(start = 0, end = .9) +  
                    labs(title="Residual") +
                    xlim(0,12) + ylim(-160, 160) + theme(legend.position="none")
  setEPS()
  postscript("multiple_devices_aov.eps", width=9.31, height=3.1)
  plot_grid(plots_aov[[1]], plots_aov[[2]], plots_aov[[3]], plots_aov[[4]], plots_aov[[5]], nrow=1)
  dev.off()


}

