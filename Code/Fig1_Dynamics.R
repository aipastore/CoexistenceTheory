#Generates an example of model dynamics

source("./solve_eqs.R") ## code for solving the dynamical equations
source("./Figplot.R") ## code for organizing and plotting results


mu1=-0.01 # intial trait value of species one - will loop from -theta to 0 with even steps
mu2=mu1+0.1# initial trait difference between species to change rho - will loop  from 0 to theta-mu1 with even steps
k1= 1.15 #this will be looped soon from like 1/10 to 10 with a logarithmic step size to change kappa ratio

S <- 2 ## number of species
w <- 0.1 ## competition width
h2 <- 0.5 ## heritability
bshape <- "b_quadratic" ## shape of intrinsic growth function
theta <-1 ## width of intrinsic growth function
sigma <- c(0.3, 0.301) ## vector of species trait standard deviations

K<-c(k1,1) # vector of maximum intrinsic growth rates.  Incurs fitness differences
params <- list(w=w, h2=h2, bshape=bshape, theta=theta, sigma=sigma, K=K) ## param list

muinit <- c(mu1, mu2) ## initial species trait means
dm <- outer(muinit, muinit, FUN="-") ## difference matrix of trait means
sv <- outer(sigma^2, sigma^2, FUN="+") ## sum matrix of trait vars
alpha <- exp(-dm^2/(2*sv+w^2))*w/sqrt(2*sv+w^2) ## alpha matrix
b <- match.fun(bshape)(K, sigma, theta, muinit)$b
ninit <- solve(alpha, b)
if (min(ninit)<=0) ninit <- rep(1, S) ## initial species densities
ic <- c(ninit, muinit) ## initial conditions coerced into a vector

tmax=100 ## time to integrate equations for
stepout <- tmax/1500 ## time step size for output
time <- seq(0, tmax, by=stepout) ## sampling points in time

sol <- ode(func=eqs, y=ic, parms=params, times=time) %>% ## solve ODEs
    organize_results(params) ## put results in tidy table

#choose data at functional extinction rather than max time
f.ext<-1e-3  ##declare a value of n that is considered functionally extinct

#extract niche overlap and fitness differences at 3 time points t=c(0, tmax, ext)
rk<-sol %>% rhokappa(params)
initialrk<-rk[time==0,] ##initial niche overlap and fitness differences

tmaxrk<-rk[time==tmax,]  ## rho and kapparatio at tmax

ext.rk.row<-min(c(which(rk$n1<f.ext),which( rk$n2 <f.ext),tmax/stepout+1))  #which row is the first occurence of a functionally extinct population
extrk<-rk[ext.rk.row,] #rho and kappa ratio at functional extinction


plot_all(sol,params, max(time), c(-1.1, 1.1, NA), 1001) ## plot results

output.names<-c("nichei", "fiti", "nichef","fitf", "x1","mu1","d", "k1", "w", "h2", "theta", "sigma" )
output<-c(initialrk$'niche overlap', initialrk$'fitness ratio', extrk$'niche overlap', extrk$'fitness ratio', mu1, d, w, h2, theta, sigma[1])

