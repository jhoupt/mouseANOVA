require("fields")
require("ggplot2")
require("cowplot")

calcSigmaMatern <- function(X1,X2,nu=5/2, l=1) {
  library(fields)
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- Matern(sqrt( (X1[i]-X2[j])^2 ), range=l, nu=nu)
    }
  }
  return(Sigma)
}




tvec <- seq(0,12,length.out=1200)

mu <- sin(tvec) / tvec
df_mu_x <- data.frame(Time=tvec, Trajectory = mu)

a1 <- -(tvec-6)^2 / 72 + .5
#a1 <- Time^3 / 12^3
a2 <- -a1
df_dev_x <- data.frame(Time=rep(tvec, 2), 
                       Condition=rep(c("a1", "a2"), each=length(tvec)),
                       Trajectory = c(a1, a2))
                                                        
nsamps <- 10
Sigma <- calcSigmaMatern(tvec, tvec)
noise <- .3 * mvrnorm(nsamps, mu=rep(0, length(tvec)), Sigma = Sigma)
df_res_x <- data.frame(Time=rep(tvec, nsamps), 
                       Trajectory=c(t(noise)),
                       idx=rep(1:nsamps, each=length(tvec)))

dat <- matrix(rep(mu, nsamps), nsamps, length(tvec), byrow=TRUE) + 
       matrix(c(rep(a1, nsamps/2), rep(a2, nsamps/2)), nsamps, 
              length(Time), byrow=TRUE) + 
       noise
df_dat_x <- data.frame(Time=rep(tvec, nsamps), 
                       Trajectory=c(t(dat)),
                       Condition=rep(c("a1","a2"),
                                     each=length(tvec)*nsamps/2),
                       idx=rep(1:nsamps, each=length(tvec)))


plots_aov_x <- vector("list", 4)
plots_aov_x[[1]] <- ggplot(data=df_dat_x, aes(x=Time, y=Trajectory, 
                                              group=idx)) + 
                    geom_line(aes(color=idx, linetype=Condition)) +
                    scale_fill_grey(start = 0, end = .9) +  
                    labs(title="Data") +
                    xlim(0,12) + ylim(-1.5, 1.5) + 
                    theme(legend.position="none")

plots_aov_x[[2]] <- ggplot(data=df_mu_x, aes(x=Time, y=Trajectory)) + 
                    geom_line() + 
                    labs(title="Grand Mean") +
                    xlim(0,12) + ylim(-1.5,1.5) + 
                    theme(legend.position="none")

plots_aov_x[[3]] <- ggplot(data=df_dev_x, aes(x=Time, y=Trajectory)) + 
                    geom_line(aes(linetype=Condition)) + 
                    labs(title="Condition") +
                    xlim(0,12) + ylim(-1.5, 1.5) + 
                    theme(legend.position="none")

plots_aov_x[[4]] <- ggplot(data=df_res_x, aes(x=Time, y=Trajectory, 
                                              group=idx)) + 
                    geom_line(aes(color=idx)) +
                    scale_fill_grey(start = 0, end = .9) +  
                    labs(title="Residual") +
                    xlim(0,12) + ylim(-1.5, 1.5) + 
                    theme(legend.position="none")


setEPS()
postscript("demo_aov.eps", width=9.31, height=3.1)
plot_grid(plots_aov_x[[1]], plots_aov_x[[2]], plots_aov_x[[3]], 
          plots_aov_x[[4]], nrow=1)
dev.off()
