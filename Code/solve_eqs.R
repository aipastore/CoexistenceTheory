
require(deSolve)
require(tidyverse)


erf <- function(x) 2*pnorm(x*sqrt(2))-1 ## error function
lseq<-function(from, to, length.out) exp(seq(log(from), log(to), length.out = length.out))  ##sequence on a log scale


## Calculate parameters assuming quadratic intrinsic growth function
## Input & output: same as b_Gaussian
b_quadratic <- function(k, s, theta, m) {
    b <- k-(m^2+s^2)/theta^2
    g <- -2*m*s^2/theta^2
    return(list(b=b, g=g))
}


## Apply smoothed step function to array
## Input:
## - n: vector or arbitrary array of values
## - a: cutoff threshold
## Output:
## - array of values with smoothed step function applied to them
cutoff <- function(n) {
    return(ifelse(n<1, (1*(n>0))*(n*n*n*(10+n*(-15+6*n))), 1))
}


## Right-hand side of dynamical equations
## Input:
## - time: moment of time (this is here for compatibility reasons; the
##         system of equations does not actually depend on time explicitly)
## - state: the dynamical state variables (densities and trait means)
##          packaged into a single vector
## - pars: list of parameters, with the following elements:
##         $w: width of competition kernel
##         $theta: width of intrinsic growth function
##         $sigma: vector of species' intraspecific standard deviations
##         $h2: heritability of trait (assumed equal for all species)
##         $bshape: the name of the function (as a string) to be executed for
##                  calculating b and g, which depend on the shape of the
##                  intrinsic growth function
## Output:
## - time derivatives of the state variables coerced into a single vector
eqs <- function(time, state, pars) {
    S <- length(pars$sigma) ## number of species
    n <- state[1:S] ## species densities
    m <- state[(S+1):(2*S)] ## species trait means
    dm <- outer(m, m, FUN="-") ## difference matrix of trait means
    sv <- outer(pars$sigma^2, pars$sigma^2, FUN="+") ## sum matrix of trait vars
    alpha <- exp(-dm^2/(2*sv+pars$w^2))*pars$w/sqrt(2*sv+pars$w^2) ## alpha matrix
    beta <- alpha*2*pars$sigma^2*(-dm)/(2*sv+pars$w^2) ## beta matrix
    p <- match.fun(pars$bshape)(pars$K, pars$sigma, pars$theta, m) ## get b and g
    ## The dndt equations are multiplied by a cutoff function, to make
    ## behavior smooth as abundances approach 0. The cutoff function is
    ## a smoothed-out step function whose derivative exists everywhere.
    dndt <- n*(p$b-alpha%*%n)*cutoff(n/(1e-7)) ## equations for abundances
    dmdt <- pars$h2*(p$g-beta%*%n) ## equations for trait means
    ## return equations by first flattening them back into a single vector
    return(list(c(dndt, dmdt)))
}

## Organize simulation results into tidy tibble
## Input:
## - sol: output produced by the function ode
## - pars: list of parameters, with the following elements:
##         $w: width of competition kernel
##         $theta: width of intrinsic growth function
##         $sigma: vector of species' intraspecific standard deviations
##         $h2: heritability of trait (assumed equal for all species)
##         $bshape: the name of the function (as a string) to be executed for
##                  calculating b and g, which depend on the shape of the
##                  intrinsic growth function
## Output:
## - a tibble with columns: time, species, n (density), m (trait mean),
##   sigma, w, theta, h2, and bshape
organize_results <- function(sol, pars) {
    S <- length(pars$sigma) ## number of species
    dat <- sol %>% as.data.frame %>% as_tibble ## convert to tibble
    names(dat)[1] <- "time" ## name the first column "time"
    names(dat)[2:(S+1)] <- paste0("n_", 1:S) ## name abundance columns (n_k)
    names(dat)[(S+2):(2*S+1)] <- paste0("m_", 1:S) ## name trait mean columns
    dat <- dat %>%
        gather("variable", "v", 2:ncol(dat)) %>% ## long format
        separate(variable, c("type", "species"), sep="_") %>%
        spread(type, v) %>% ## separate columns for densities n and trait means m
        select(time, species, n, m) %>% ## rearrange columns
        mutate(species=as.integer(species), sigma=pars$sigma[species], w=pars$w,
               theta=pars$theta, h2=pars$h2, bshape=pars$bshape) ## add parameters
    return(as_tibble(dat))
}

## Obtain niche overlap rho and the kappa ratio through time
## Input:
## - sol: output produced by the function ode
## - pars: list of parameters, with the following elements:  
##         $w: width of competition kernel
##         $kappa:  difference in height of intrinsic growth curve   
##         $theta: width of intrinsic growth function
##         $sigma: vector of species' intraspecific standard deviations
##         $h2: heritability of trait (assumed equal for all species)
##         $bshape: the name of the function (as a string) to be executed for
##                  calculating b and g, which depend on the shape of the
##                  intrinsic growth function
## Output:
## - a tibble with columns: time, niche overlap, fitness ratio
rhokappa <- function(sol, pars) {
    m1 <- sol %>% filter(species==1) %>% pull(m)
    m2 <- sol %>% filter(species==2) %>% pull(m)
    n1 <- sol %>% filter(species==1) %>% pull(n)
    n2 <- sol %>% filter(species==2) %>% pull(n)
    sv <- outer(pars$sigma^2, pars$sigma^2, FUN="+")
    rho <- rep(0, length(m1))
    kapparatio <- rep(0, length(m1))
    for (i in 1:length(m1)) {
        mm <- c(m1[i], m2[i])
        dm <- outer(mm, mm, FUN="-")
        alpha <- exp(-dm^2/(2*sv+pars$w^2))*pars$w/sqrt(2*sv+pars$w^2)
        b <- match.fun(pars$bshape)(pars$K, pars$sigma, pars$theta, mm)$b   
        rho[i] <- sqrt(alpha[1,2]*alpha[2,1]/(alpha[1,1]*alpha[2,2]))
        kapparatio[i] <- (b[1]/b[2])*sqrt(alpha[2,1]*alpha[2,2]/
                                          (alpha[1,2]*alpha[1,1]))
    }
    return(tibble(time=time, `Niche overlap`=rho,
                  `Competitive difference`=kapparatio, n1=n1, n2=n2))
}

