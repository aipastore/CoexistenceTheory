
clargs <- commandArgs(trailingOnly=TRUE)
if (length(clargs)>0) { ## command-line arguments
    theta <- as.numeric(clargs[1]) ## width of intrinsic growth function
    K1 <- as.numeric(clargs[2]) ## species 1 max intrinsic growth rate
} else { ## sample parameters if no command line arguments are given
    theta <- 2
    K1 <- 1.01
}

source("./solve_eqs.R") ## code for solving the dynamical equations

fact <- expand.grid( ## table for factorial numerical experiment
    theta=theta, ## width of intrinsic growth function
    K1=K1, ## species 1 max intrinsic growth rate (command line parameter)
    mu1=seq(-4, 0, l=21), ## sp 1 init trait; positive ones omitted due to symmetry
    mu2=seq(-4, 4, l=201), ## species 2 initial trait mean
    h2=seq(0.1, 0.5, l=5), ## heritability
    w=0.1, ## competition width
    sigma=0.3, ## intraspecific trait std dev (for both species)
    tmax=10000, ## time to integrate dynamics for
    fext=0.001 ## functional extinction threshold
) %>% as_tibble %>% ## convert data frame to tibble
filter(mu2<=mu1) %>% ## species 1 always to the left of species 2 initially
## confine initial trait means to [-theta, theta] interval:
filter(mu1 >= -theta, mu1 <= theta, mu2 >= -theta, mu2 < theta)

ninit <- c(0.2, 0.2) ## initial species densities
write("theta,K1,mu1,mu2,h2,w,sigma,tmax,fext,nichei,fiti,nichef,fitf,t_ext",
      stdout()) ## write header

for (r in 1:nrow(fact)) { ## go through each row of experiments
    params <- list(w=fact$w[r], h2=fact$h2[r], bshape="b_quadratic", ## list of
                   theta=fact$theta[r], K=c(fact$K1[r], 1), ## model parameters
                   sigma=c(fact$sigma[r], fact$sigma[r]))
    muinit <- c(fact$mu1[r], fact$mu2[r]) ## set initial trait means
    ic <- c(ninit, muinit) ## initial conditions coerced into a vector
    time <- seq(0, fact$tmax[r], by=1) ## sampling points in time
    sol <- ode(func=eqs, y=ic, parms=params, times=time) %>% ## solve ODEs
        organize_results(params) ## put results in tidy table
    rkappa <- sol %>% ## extract niche overlap and fit diffs at beginning and end
        rhokappa(params) %>%
        filter(n1>fact$fext[r], n2>fact$fext[r])
    write(paste(paste(fact[r,], collapse=","),
                rkappa$`Niche overlap`[1], ## initial niche overlap
                rkappa$`Competitive difference`[1], ## initial fit diff
                rkappa$`Niche overlap`[nrow(rkappa)], ## final niche overlap
                rkappa$`Competitive difference`[nrow(rkappa)], ## final fit diff
                rkappa$time[nrow(rkappa)], ## time of 1st ext (= tmax if none)
                sep=","), stdout()) ## write row to stdout
}

