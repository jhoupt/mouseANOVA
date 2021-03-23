require("plyr")
require("tibble")
require("cowplot")
require("BayesFactor")
require("ggplot2")
require("xlsx")
gitdir <- "~/Documents/GitHub/mouseANOVA/"
#setwd("~/Documents/Research/ResponseDynamicAnalysis/fANOVA/Analysis")
source(paste(gitdir, "mouseFANOVA.R", sep=""))

loadKJ3Subject_fill <- function(n_time_bins) { 
    data_path = '~/Documents/Research/ResponseDynamicAnalysis/fANOVA/Data/KJExperiment3/';
    dat_file <-paste(data_path, "alldat.Rdata", sep="")
    if (file.exists(dat_file)) { 
      fname <- load(dat_file)
      return(dat)
    } else { 
      files = dir(path=data_path, pattern='[[:digit:]]+.xlsx')
      files <- files[files != "8.xlsx"]
      
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
      
      n_files <- length(files)
      dat_tibble <- NULL

      for (fn in 1:n_files) { 
          print(paste("File ", fn, " (", files[fn], ")", sep=""))
          file_name <- paste(data_path, files[fn], sep='')
          dat_tibble<- rbind(dat_tibble, read.subject.kj3(file_name, out_t))
          #all_data[[fn]] <- read.subject.kj3(file_name, out_t, 
          #                                   collapse_to_cell=FALSE)
          #                                   #collapse_to_cell=TRUE)
      }
      dat_tibble$Subject <- factor(dat_tibble$Subject)
      dat_tibble$Condition <- factor(dat_tibble$Condition)
      dat <- list(times=out_t, dat_tibble=dat_tibble)
      
      save(dat, file=dat_file)
      return(dat)
    }
}


read.subject.kj3 <- function(file_name, out_t) { 
  library(xlsx)
  raw = read.xlsx(file_name, 1, header=TRUE, startRow=2);
  raw = raw[2:dim(raw)[1],]

  # Most complicated way to to this? 
  subject <- gsub("[.]", "", regmatches(file_name, 
                  regexpr("[[:digit:]]+[.]", file_name))[1])

  sampling_rate = 60;
  max_time = 40;
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

  X <- raw[,grepl('XTracked', names(raw))]
  x_idx <- as.numeric(sub('XTracked', '', 
                      names(raw)[grepl('XTracked', names(raw))]))
  x_idx <- sort(x_idx, index.return=TRUE)$ix
  X <- X[,x_idx]
  X <- X - start_loc[1]

  Y <- raw[,grepl('YTracked', names(raw))]
  y_idx <- as.numeric(sub('YTracked', '', 
                      names(raw)[grepl('YTracked', names(raw))]))
  y_idx <- sort(y_idx, index.return=TRUE)$ix
  Y <- Y[,y_idx]
  Y <- start_loc[2] - Y # Assumes flipped Y in data (e.g., screen)


  subject <- unique(raw[,"Subject"])
  trial_n <- raw[,'Trial']

  trials <-  unique(trial_n)
  t_vec <- 1:dim(X)[2]


  prob_left <- raw[,"LeftProb"]
  prob_right <- raw[,"RightProb"]
  choice <- raw[,"Stimulus.ACC"]


  # 1: Left column is safer; 2: right is safer; 0: equal probability
  which_safe <- (prob_right < prob_left) + 2 * (prob_left < prob_right)
  chose_safe <- (which_safe == choice)


  #all_trials <- vector("list", length=length(trials))
  #all_trial_mat <- matrix(NA, length(trials), length(out_t))
  #all_risky <- c()
  #all_safe <- c()

  all_tr <- 1:length(trials)
  all_raw <- vector("list", length=length(trials))
  all_clean <- vector("list", length=length(trials))
  all_aad <- rep(NA, length(trials))
  all_mad <- rep(NA, length(trials))
  all_mad_time <- rep(NA, length(trials))
  all_xflips <- rep(NA, length(trials))
  all_oflips <- rep(NA, length(trials))



  for (tr in all_tr) { 
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
    trial_dat$rotated_pos <- list(x=a[1,], y=a[2,])

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
  subj_tibble <- tibble(Subject=subject, Condition="Mouse", Trial=all_tr,
                        Response=all_responses, AAD=all_aad, MAD=all_mad, 
                        MADtime=all_mad_time,
                        Xflips=all_xflips, Oflips=all_oflips, 
                        RawTrajectory=all_raw, CleanTrajectory=all_clean)
  return(subj_tibble)
}


analyze.md <- function() { 
  dat <- loadKJ3Subject_fill(100)
  collapsed <- collapse_to_cell(dat)

  dat_frame <- data.frame(dat$dat_tibble[,1:9])
  dat_frame_c <- data.frame(collapsed[,1:8])


  ## BF ANOVA
  if (FALSE) { 
    madt_aov <- anovaBF(MADtime ~ Response + Subject, 
                        data=dat_frame, whichRandom="Subject")

    mad_aov <- anovaBF(MAD ~ Response + Subject, 
                        data=dat_frame, whichRandom="Subject")

    aad_aov <- anovaBF(AAD ~ Response + Subject, 
                        data=dat_frame, whichRandom="Subject")
  } else { 
    madt_aov <- anovaBF(MADtime ~ Response + Subject, 
                        data=subset(dat_frame_c, !is.na(AAD)), 
                        whichRandom="Subject")
    mad_aov <- anovaBF(MAD ~ Response + Subject, 
                        data=subset(dat_frame_c, !is.na(AAD)), 
                        whichRandom="Subject")
    aad_aov <- anovaBF(AAD ~ Response + Subject, 
                        data=subset(dat_frame_c, !is.na(AAD)), 
                        whichRandom="Subject")
    xflips_aov <- anovaBF(Xflips ~ Response + Subject, 
                          data=subset(dat_frame_c, !is.na(AAD)), 
                          whichRandom="Subject")
    oflips_aov <- anovaBF(Oflips ~ Response + Subject, 
                          data=subset(dat_frame_c, !is.na(AAD)), 
                          whichRandom="Subject")


    dvnames <- c("AAD", "MAD", "MADtime", "Xflips", "Oflips")
    dvlabs <- c("Mean Absolute Deviation", 
                "Maximum Absolute Deviation", 
                "MaxAD Time", "Horizontal Direction Changes", 
                "Origin Crossings")
    ylabs <- c(rep("Deviation (pixels)", 2), "Time (ms)", 
               rep("Count", 2))
    plots <- vector("list", length(dvnames))
    for (dn in 1:length(dvnames)) { 
      dv <- dvnames[dn]
      sumdat <- ddply(dat_frame_c, c("Response"), 
                       .fun=summary_func, dv)
      #plots[[dn]] <- ggplot(sumdat, 
      #          aes(x=Response, y=mean, fill=Response)) + 
      #          geom_bar(position=position_dodge(), stat="identity") + 
      #          geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
      #                        width=.2,
      #                        position=position_dodge(.9)) + 
      plots[[dn]] <- ggplot(dat_frame_c, 
                aes(x=Response, y=.data[[dv]], fill=Response)) + 
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
    postscript("kj3_standard_metrics.eps", width=9.31, height=6.2)
                
    plot_grid(toprow, bottomrow, ncol=1)
    dev.off()
    
  }
  x <- loocv_anova(dat, collapsed)
  anova_full <- anova_model(dat$times, collapsed, c("Response"))
  sd_t <- sqrt(apply((anova_full$residuals)^2, 2, mean)) / sqrt(52)
  
  df_full <- data.frame(Time=rep(dat$times, 2), 
        Response=rep(c("Risky","Safe"), each=length(dat$times)), 
        Deviation=c(anova_full$gm + anova_full$alpha[[1]][[1]],
                    anova_full$gm + anova_full$alpha[[1]][[2]]))
  df_full$Deviation.lwr <- df_full$Deviation - rep(sd_t, 2)
  df_full$Deviation.upr <- df_full$Deviation + rep(sd_t, 2)

  hues = seq(15, 375, length = 2 + 1)

  postscript("kj3_gp.eps", width=9.31, height=3.1)
  ggplot(data=df_full, aes(x=Time, y=Deviation, group=Response)) + 
     geom_ribbon(aes(ymin=Deviation.lwr,ymax=Deviation.upr,fill=Response)) +
     scale_fill_manual(values=hcl(h = hues, l = 80, c = 40)[1:2]) + 
     geom_line(aes(color=Response)) + 
     xlim(0,12) + 
     labs(title="Koop and Johnson (2013)", y="Deviation (pixels)",
          breaks=c("risky", "safe"),
          labels=c("Risky", "Safe"))
  dev.off()

}
