## Functions aiding numerical integration of the ordinary differential equations
## (ODEs), obtaining niche overlap and competitive differences from the results,
## and organizing them into a tidy output format.


require(deSolve) ## for solving differential equations
require(tidyverse) ## for efficient data manipulation and plotting


## Apply smoothed step function to array
## Input:
## - n: vector or arbitrary array of values
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
##         $K: vector of species-specific intrinsic rates
##         $sigma: vector of species' intraspecific standard deviations
##         $h2: vector of heritabilities
## Output:
## - time derivatives of the state variables coerced into a single vector
eqs <- function(time, state, pars) {
  S <- length(pars$sigma) ## number of species
  n <- state[1:S] ## species densities
  m <- state[(S+1):(2*S)] ## species trait means
  b <- pars$K-(m^2+pars$sigma^2)/pars$theta^2 ## species-level intrinsic growth
  g <- -2*m*pars$sigma^2/pars$theta^2 ## species-level selection pressure
  dm <- outer(m, m, FUN="-") ## difference matrix of trait means
  sv <- outer(pars$sigma^2, pars$sigma^2, FUN="+") ## sum matrix of trait vars
  alpha <- exp(-dm^2/(2*sv+pars$w^2))*pars$w/sqrt(2*sv+pars$w^2) ## alpha matrix
  beta <- alpha*2*pars$sigma^2*(-dm)/(2*sv+pars$w^2) ## beta matrix
  ## The dndt equations are multiplied by a cutoff function, to make
  ## behavior smooth as abundances approach 0. The cutoff function is
  ## a smoothed-out step function whose derivative exists everywhere.
  dndt <- n*(b-alpha%*%n)*cutoff(n/(1e-7)) ## equations for abundances
  dmdt <- pars$h2*(g-beta%*%n) ## equations for trait means
  ## return equations by first flattening them back into a single vector
  return(list(c(dndt, dmdt)))
}


## Organize simulation results into tidy tibble
## Input:
## - sol: output produced by the function ode
## - pars: list of parameters, with the following elements:
##         $w: width of competition kernel
##         $theta: width of intrinsic growth function
##         $K: vector of species-specific intrinsic rates
##         $sigma: vector of species' intraspecific standard deviations
##         $h2: vector of heritabilities
## Output:
## - a tibble with columns: time, species, n (density), m (trait mean),
##   w, theta, K, sigma, and h2
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
    mutate(species=as.integer(species), w=pars$w, theta=pars$theta, ## add params
           K=pars$K[species], sigma=pars$sigma[species], h2=pars$h2[species])
  return(dat)
}


## Obtain niche overlap rho and the kappa ratio through time
## Input:
## - sol: output produced by the function ode
## - pars: list of parameters, with the following elements:  
##         $w: width of competition kernel
##         $theta: width of intrinsic growth function
##         $K: vector of species-specific intrinsic rates
##         $sigma: vector of species' intraspecific standard deviations
##         $h2: vector of heritabilities
## Output:
## - a tibble with columns: time, niche overlap, fitness ratio
rhokappa <- function(sol, pars) {
  m1 <- sol %>% filter(species==1) %>% pull(m) ## sp. 1 trait mean through time
  m2 <- sol %>% filter(species==2) %>% pull(m) ## sp. 2 trait mean through time
  n1 <- sol %>% filter(species==1) %>% pull(n) ## sp. 1 density through time
  n2 <- sol %>% filter(species==2) %>% pull(n) ## sp. 2 density through time
  sv <- outer(pars$sigma^2, pars$sigma^2, FUN="+") ## sum matrix of trait vars
  rho <- rep(0, length(m1)) ## vector to store niche overlap through time
  kapparatio <- rep(0, length(m1)) ## vector to store comp. diff. through time
  for (i in 1:length(m1)) { ## at each point in time:
    mm <- c(m1[i], m2[i]) ## create vector of trait means
    b <- pars$K-(mm^2+pars$sigma^2)/pars$theta^2 ## intrinsic growth rates
    dm <- outer(mm, mm, FUN="-") ## difference matrix of trait means
    alpha <- exp(-dm^2/(2*sv+pars$w^2))*pars$w/sqrt(2*sv+pars$w^2) ## alpha matrix
    rho[i] <- sqrt(alpha[1,2]*alpha[2,1]/(alpha[1,1]*alpha[2,2])) ## niche overlap
    kapparatio[i] <- (b[1]/b[2])*sqrt(alpha[2,1]*alpha[2,2]/(
      alpha[1,2]*alpha[1,1])) ## competitive difference
  } ## end of for loop
  ## put the results into columns of a tibble, and return this tibble:
  return(tibble(time, rho, kapparatio, m1, m2, n1, n2))
}
