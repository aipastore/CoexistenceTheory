#Script generates analytical equilibrium solutions for model see results folder for R output data

require(nleqslv)
require(tidyverse)

## equilibrium equations as a function of state vars (state) and parameters (pars)
eqb <- function(state, pars) {
  S <- length(pars$sigma) ## number of species
  n <- state[1:S] ## species densities
  m <- state[(S+1):(2*S)] ## species trait means
  dm <- outer(m, m, FUN="-") ## difference matrix of trait means
  sv <- outer(pars$sigma^2, pars$sigma^2, FUN="+") ## sum matrix of trait vars
  alpha <- exp(-dm^2/(2*sv+pars$w^2))*pars$w/sqrt(2*sv+pars$w^2) ## alpha matrix
  beta <- alpha*2*pars$sigma^2*(-dm)/(2*sv+pars$w^2) ## beta matrix
  b <- pars$K+(m+pars$sigma^2)/pars$theta^2 ## b vector
  g <- -2*m*pars$sigma^2/pars$theta^2 ## g vector
  return(c(b-alpha%*%n, g-beta%*%n)) ## equilibrium equations
}

## define table for factorial numerical experiment
fact <- expand.grid(
  theta=exp(seq(log(0.5), log(10), l=20)), ## width of intrinsic growth function
  K1=exp(seq(log(0.2), log(5), l=30)), ## species 1 max intrinsic growth rate
  mu1=seq(-10, 10, l=41), ## sp 1 init trait; positive ones omitted due to symmetry
  mu2=seq(-10, 10, l=41), ## species 2 initial trait mean
  w=c(0.5, 1, 3, 5) ## competition width
) %>%
  as_tibble %>% ## convert data frame to tibble
  filter(mu2<=mu1) %>% ## species 1 always to the left of species 2 initially
  ## confine initial trait means to [-theta, theta] interval:
  filter(mu1 >= -theta, mu1 <= theta, mu2 >= -theta, mu2 <= theta)

## write header to stdout (redirect output to write to file):
write(paste0("theta,K1,w,ninit1,ninit2,muinit1,muinit2,",
             "nfinal1,nfinal2,mufinal1,mufinal2,singular"), stdout())

for (r in 1:nrow(fact)) { ## go through each row of experiments
  params <- list(w=fact$w[r], theta=fact$theta[r], K=c(fact$K1[r], 1),
                 sigma=c(1, 1)) ## model parameters
  muinit <- c(fact$mu1[r], fact$mu2[r])
  dm <- outer(muinit, muinit, FUN="-") ## difference matrix of trait means
  sv <- outer(params$sigma^2, params$sigma^2, FUN="+") ## sum matrix of trait vars
  alpha <- exp(-dm^2/(2*sv+params$w^2))*params$w/sqrt(2*sv+params$w^2)
  b <- params$K+(muinit+params$sigma^2)/params$theta^2
  if (abs(det(alpha))<1e-10) ninit <- c(1, 0.9) else ninit <- solve(alpha, b)
  if (min(ninit)<=0) next
  ic <- c(ninit, muinit) ## initial conditions coerced into a vector
  sol <- nleqslv(c(ninit, muinit), eqb, pars=params, jacobian=TRUE) ## eqb eqs
  sing <- (max(Re(eigen(sol$jac, only.values=TRUE)$values))>=0) ## stable?
  if (sol$termcd %in% c(5, 6)) sing <- TRUE ## is Jacobian too ill-conditioned?
  sing <- ifelse(sing, "yes", "no")
  write(paste(fact$theta[r], fact$K1[r], fact$w[r], ninit[1], ninit[2], muinit[1],
              muinit[2], sol$x[1], sol$x[2], sol$x[3], sol$x[4], sing,
              sep=","), stdout()) ## write row of data to stdout
}

