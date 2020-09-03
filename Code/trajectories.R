## This script shows various dynamical trajectories of communities
## in the space spanned by niche overlap and competitive differences.
## It is Figure 5 in the manuscript.

source("./functions.R") ## for integrating ODEs
theme_set(theme_bw()) ## set default theme
theme_update(panel.grid.major=element_blank(), panel.grid.minor=element_blank())


h2 <- c(0.1, 0.1) ## heritability
sigma <- c(1, 1) ## vector of species trait standard deviations
tmax <- 1e4 ## time to integrate equations for
stepout <- tmax/500 ## time step size for output
time <- seq(0, tmax, by=stepout) ## sampling points in time
separatrix <- tibble(rho=seq(0, 1, by=0.01)) %>%
  mutate(lower=log(rho), upper=log(rho)*(-1))

param <- expand_grid( ## table for factorial numerical experiment
  w=c(1, 3, 5), ## competition width
  theta=c(5.32, 8.54), ## intrinsic growth width
  K1=c(1.0, 1.5), ## species 1 max intrinsic growth rate
  mu1=c(0, 4), ## species 1 initial trait mean
  mu2=c(-1, -3) ## species 2 initial trait mean
)

traj <- tibble() ## create empty tibble
for (i in 1:nrow(param)) {
  pars <- list(w=param$w[i], theta=param$theta[i], K=c(param$K1[i], 1),
               sigma=sigma, h2=h2)## set parameters
  muinit <- c(param$mu1[i], param$mu2[i]) ## initial trait means
  dm <- outer(muinit, muinit, FUN="-") ## difference matrix of trait means
  sv <- outer(sigma^2, sigma^2, FUN="+") ## sum matrix of trait vars
  alpha <- exp(-dm^2/(2*sv+pars$w^2))*pars$w/sqrt(2*sv+pars$w^2) ## alpha matrix
  b <- pars$K-(muinit^2+pars$sigma^2)/pars$theta^2 ## sp-level intrinsic growth
  ninit <- solve(alpha, b) ## initial densities fit to local equilibrium condition
  if (min(ninit)<=0) ninit <- c(1, 1) ## unless local coexistence is not possible
  ic <- c(ninit, muinit) ## initial conditions coerced into a vector
  sol <- ode(func=eqs, y=ic, parms=pars, times=time) %>% ## solve ODEs
    organize_results(pars) %>% ## put results in tidy table
    rhokappa(pars) %>% ## get rho and kappa ratio through time
    mutate(D=sqrt((rho-lag(rho))^2+(log(kapparatio)-log(lag(kapparatio)))^2)) %>%
    filter((D>1e-4)|(time==0)) %>%
    mutate(w=param$w[i],
           par=paste0("theta == ", param$theta[i], "~~K[1] == ", param$K1[i]),
           muinit=paste(param$mu1[i], param$mu2[i], sep="_"))
  traj <- bind_rows(traj, sol)
}

## create plot
ggplot(traj) +
  geom_line(data=separatrix, aes(x=rho, y=upper), colour="#0072B2") +
  geom_line(data=separatrix, aes(x=rho, y=lower), colour="#0072B2") +
  geom_ribbon(data=separatrix, aes(x=rho, ymin=upper, ymax=2),
              colour=NA, fill="#0072B2", alpha=0.1) +
  geom_ribbon(data=separatrix, aes(x=rho, ymax=lower, ymin=-2),
              colour=NA, fill="#0072B2", alpha=0.1) +
  geom_hline(yintercept=0, colour="#0072B2") +
  geom_path(data=traj, aes(x=rho, y=log(kapparatio), linetype=muinit, colour=par),
            na.rm=TRUE, alpha=0.8, arrow=arrow(length=unit(0.04, "npc"), angle=15,
                                               ends="last", type="closed")) +
  scale_x_continuous(name=expression(paste(rho)), limits=c(0, 1),
                     expand=c(0, 0), labels=abbreviate) +
  scale_y_continuous(name=expression(paste(log(kappa[1]/kappa[2]))),
                     limits=c(-2, 2), expand=c(0, 0)) +
  scale_colour_viridis_d(name=NULL, end=0.92, labels=scales::label_parse()) +
  scale_linetype_manual(values=rep(1, 4), guide=FALSE) +
  facet_wrap(~w, nrow=1, labeller=label_bquote(omega==.(w))) +
  theme(legend.text.align=0)
