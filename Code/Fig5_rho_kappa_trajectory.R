
source("./solve_eqs.R") ## functions for integrating the ODEs
source("./plot.R") ## plotting functions

theme_set(theme_bw())
theme_update(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

br <- c(1.07, 2.01, 2.76, 5.18)


##set.seed(57321) ## set random seed for reproducibility
S <- 2 ## number of species
w <- 0.1 ## competition width
h2 <- 0.5 ## heritability
bshape <- "b_quadratic" ## shape of intrinsic growth function
sigma <- c(0.3, 0.3) ## vector of species trait standard deviations

tmax <- 1e3 ## time to integrate equations for
stepout <- tmax/5000 ## time step size for output
time <- seq(0, tmax, by=stepout) ## sampling points in time

trajectory<-function(th, K1, mu1, mu2){
  theta <- th ## width of intrinsic growth function
  K <- c(K1, 1)
  params <- list(w=w, h2=h2, bshape=bshape, theta=theta, sigma=sigma, K=K)
  muinit <- c(mu1, mu2) ## initial trait means
  dm <- outer(muinit, muinit, FUN="-") ## difference matrix of trait means
  sv <- outer(sigma^2, sigma^2, FUN="+") ## sum matrix of trait vars
  alpha <- exp(-dm^2/(2*sv+w^2))*w/sqrt(2*sv+w^2) ## alpha matrix
  b <- match.fun(bshape)(K, sigma, theta, muinit)$b
  ninit <- solve(alpha, b)
  if (min(ninit)<=0) ninit <- rep(1, S) ## initial species densities
  ic <- c(ninit, muinit) ## initial conditions coerced into a vector

  sol <- ode(func=eqs, y=ic, parms=params, times=time) %>% ## solve ODEs
    organize_results(params) ## put results in tidy table
  rk <- rhokappa(sol, params) %>% ## get rho and kappa ratio through time
    rename(rho=`Niche overlap`, kratio=`Competitive difference`) %>%
    mutate(D=sqrt((rho-lag(rho))^2+(log(kratio)-log(lag(kratio)))^2)) %>%
    filter((D>1e-3)|(time==0))

  kratiolimit <- max(2, abs(log(max(rk$kratio))))
  separatrix <- tibble(rho=seq(0,1,by=0.01)) %>%
    mutate(lower=log(rho), upper=log(rho)*(-1))

  output<-list(rk, separatrix,kratiolimit)
  return(output)
}

#PLOT FOR
th=1.554
output<-trajectory(th, .97, 0 , .3)
rk<-output[[1]]
separatrix<-output[[2]]
kratiolimit<-output[[3]]

output<-trajectory(th, 2.5, 0 , .3)
rk2<-output[[1]]

output<-trajectory(th, .9, -.6, .6)
rk3<-output[[1]]

output<-trajectory(th, .4, 0, .3)
rk4<-output[[1]]

output<-trajectory(th, 3, -.65 , .65)
rk5<-output[[1]]

output<-trajectory(th, .4, -.6, .6)
rk6<-output[[1]]

plot1<-ggplot(rk) +
	ggtitle(expression(paste(theta/sigma," = " ,5.18)))+
  geom_line(data=separatrix, aes(x=rho, y=upper), colour="#377eb8") +
  geom_line(data=separatrix, aes(x=rho, y=lower), colour="#377eb8") +
  geom_ribbon(data=separatrix, aes(x=rho, ymin=upper, ymax=kratiolimit),
              colour=NA, fill="#377eb8", alpha=0.1) +
  geom_ribbon(data=separatrix, aes(x=rho, ymax=lower, ymin=-kratiolimit),
              colour=NA, fill="#377eb8", alpha=0.1) +
  geom_hline(yintercept=0, colour="#377eb8") +
  geom_path(data=rk, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk, time==0), aes(x=rho, y=log(kratio))) +
  geom_path(data=rk2, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk2, time==0), aes(x=rho, y=log(kratio))) +
  geom_path(data=rk3, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk3, time==0), aes(x=rho, y=log(kratio))) +
  geom_path(data=rk4, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk4, time==0), aes(x=rho, y=log(kratio))) +
  geom_path(data=rk5, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk5, time==0), aes(x=rho, y=log(kratio))) +
    geom_path(data=rk6, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk6, time==0), aes(x=rho, y=log(kratio))) +
  scale_x_continuous(name=expression(paste(rho)), limits=c(0, 1),
                     expand=c(0, 0), labels=abbreviate) +
  scale_y_continuous(name=expression(paste(log(kappa[1]/kappa[2]))),
                     limits=c(-kratiolimit, kratiolimit), expand=c(0, 0))


#PLOT FOR
th=0.828
output<-trajectory(th, .97, 0 , .3)
rk<-output[[1]]
separatrix<-output[[2]]
kratiolimit<-output[[3]]

# output<-trajectory(th, 2.5, 0 , .3)
# rk2<-output[[1]]

output<-trajectory(th, 1, -.6, .6)
rk3<-output[[1]]

# output<-trajectory(th, .4, 0, .3)
# rk4<-output[[1]]

output<-trajectory(th, 1.7, -.65 , .65)
rk5<-output[[1]]

output<-trajectory(th, .77, -.5, .6)
rk6<-output[[1]]

plot2<-ggplot(rk) +
	ggtitle(expression(paste(theta/sigma," = " ,2.76)))+
  geom_line(data=separatrix, aes(x=rho, y=upper), colour="#377eb8") +
  geom_line(data=separatrix, aes(x=rho, y=lower), colour="#377eb8") +
  geom_ribbon(data=separatrix, aes(x=rho, ymin=upper, ymax=kratiolimit),
              colour=NA, fill="#377eb8", alpha=0.1) +
  geom_ribbon(data=separatrix, aes(x=rho, ymax=lower, ymin=-kratiolimit),
              colour=NA, fill="#377eb8", alpha=0.1) +
  geom_hline(yintercept=0, colour="#377eb8") +
  geom_path(data=rk, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk, time==0), aes(x=rho, y=log(kratio))) +
  # geom_path(data=rk2, aes(x=rho, y=log(kratio)),
            # arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  # geom_point(data=filter(rk2, time==0), aes(x=rho, y=log(kratio))) +
  geom_path(data=rk3, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk3, time==0), aes(x=rho, y=log(kratio))) +
  # geom_path(data=rk4, aes(x=rho, y=log(kratio)),
            # arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  # geom_point(data=filter(rk4, time==0), aes(x=rho, y=log(kratio))) +
  geom_path(data=rk5, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk5, time==0), aes(x=rho, y=log(kratio))) +
    geom_path(data=rk6, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk6, time==0), aes(x=rho, y=log(kratio))) +
  scale_x_continuous(name=expression(paste(rho)), limits=c(0, 1),
                     expand=c(0, 0), labels=abbreviate) +
  scale_y_continuous(name=expression(paste(log(kappa[1]/kappa[2]))),
                     limits=c(-kratiolimit, kratiolimit), expand=c(0, 0))



  #PLOT FOR
  th=0.603
output<-trajectory(th, 1.5, -.2 , .3)
rk<-output[[1]]
separatrix<-output[[2]]
kratiolimit<-output[[3]]


output<-trajectory(th, .8, -.20, .3)
rk3<-output[[1]]

output<-trajectory(th, 3, -.3 , .3)
rk5<-output[[1]]

output<-trajectory(th, .6, -.3, .3)
rk6<-output[[1]]

plot3<-ggplot(rk) +
	ggtitle(expression(paste(theta/sigma," = " ,2.01)))+
  geom_line(data=separatrix, aes(x=rho, y=upper), colour="#377eb8") +
  geom_line(data=separatrix, aes(x=rho, y=lower), colour="#377eb8") +
  geom_ribbon(data=separatrix, aes(x=rho, ymin=upper, ymax=kratiolimit),
              colour=NA, fill="#377eb8", alpha=0.1) +
  geom_ribbon(data=separatrix, aes(x=rho, ymax=lower, ymin=-kratiolimit),
              colour=NA, fill="#377eb8", alpha=0.1) +
  geom_hline(yintercept=0, colour="#377eb8") +
  geom_path(data=rk, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk, time==0), aes(x=rho, y=log(kratio))) +
  # geom_path(data=rk2, aes(x=rho, y=log(kratio)),
            # arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  # geom_point(data=filter(rk2, time==0), aes(x=rho, y=log(kratio))) +
  geom_path(data=rk3, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk3, time==0), aes(x=rho, y=log(kratio))) +
  geom_path(data=rk5, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk5, time==0), aes(x=rho, y=log(kratio))) +
    geom_path(data=rk6, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk6, time==0), aes(x=rho, y=log(kratio))) +
  scale_x_continuous(name=expression(paste(rho)), limits=c(0, 1),
                     expand=c(0, 0), labels=abbreviate) +
  scale_y_continuous(name=expression(paste(log(kappa[1]/kappa[2]))),
                     limits=c(-kratiolimit, kratiolimit), expand=c(0, 0))

#PLOT FOR
  th=0.321
output<-trajectory(th, .97, 0 , .1)
rk<-output[[1]]
separatrix<-output[[2]]
kratiolimit<-output[[3]]

output<-trajectory(th, 2.5, -.1 , .1)
rk2<-output[[1]]

output<-trajectory(th, .9, 0, .1)
rk3<-output[[1]]

output<-trajectory(th, .8, 0, .1)
rk4<-output[[1]]

output<-trajectory(th, 3, 0 , .1)
rk5<-output[[1]]

output<-trajectory(th, 1.2, 0, .1)
rk6<-output[[1]]

plot4<-ggplot(rk) +
	ggtitle(expression(paste(theta/sigma," = " ,1.07)))+
  geom_line(data=separatrix, aes(x=rho, y=upper), colour="#377eb8") +
  geom_line(data=separatrix, aes(x=rho, y=lower), colour="#377eb8") +
  geom_ribbon(data=separatrix, aes(x=rho, ymin=upper, ymax=kratiolimit),
              colour=NA, fill="#377eb8", alpha=0.1) +
  geom_ribbon(data=separatrix, aes(x=rho, ymax=lower, ymin=-kratiolimit),
              colour=NA, fill="#377eb8", alpha=0.1) +
  geom_hline(yintercept=0, colour="#377eb8") +
  geom_path(data=rk, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk, time==0), aes(x=rho, y=log(kratio))) +
  geom_path(data=rk2, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk2, time==0), aes(x=rho, y=log(kratio))) +
  geom_path(data=rk3, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk3, time==0), aes(x=rho, y=log(kratio))) +
  geom_path(data=rk4, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk4, time==0), aes(x=rho, y=log(kratio))) +
  geom_path(data=rk5, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk5, time==0), aes(x=rho, y=log(kratio))) +
    geom_path(data=rk6, aes(x=rho, y=log(kratio)),
            arrow=arrow(length=unit(0.04, "npc"), angle=15, ends="last", type="closed")) +
  geom_point(data=filter(rk6, time==0), aes(x=rho, y=log(kratio))) +
  scale_x_continuous(name=expression(paste(rho)), limits=c(0, 1),
                     expand=c(0, 0), labels=abbreviate) +
  scale_y_continuous(name=expression(paste(log(kappa[1]/kappa[2]))),
                     limits=c(-kratiolimit, kratiolimit), expand=c(0, 0))

plot_grid(plot4,plot3,plot2,plot1, nrow=1)
