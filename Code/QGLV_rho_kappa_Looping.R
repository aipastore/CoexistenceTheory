#Pastore modifying Barabas's codes

source("./solve_eqs.R") ## code for solving the dynamical equations
source("./plot.R") ## code for organizing and plotting results

#create an empty dataframe to add all the output data to
output<- data.frame("nichei"=double(),
												"fiti"=double(), 
												"nichef"=double(),
												"fitf"=double(),
												"mu1"=double(),
												"d"=double(),
												 "k1"=double(),
												  "w"=double(),
												   "h2"=double(), 
												   "theta"=double(),
												    "sigma"=double(),
												    "ext.time"=double(),
												    "f.ext"=double())
names(output)<-c("nichei", "fiti", "nichef","fitf","mu1","d", "k1", "w", "h2", "theta", "sigma" ,"ext.time","f.ext")

S <- 2 ## number of species
w <- 3 ## competition width
h2 <- 0.4 ## heritability
bshape <- "b_quadratic" ## shape of intrinsic growth function
sigma <- c(1, 1.001) ## vector of species trait standard deviations
ninit <- c(1, 0.9) ## initial species densities
w=1

tmax=3000 ## time to integrate equations for 
stepout <- tmax/1500 ## time step size for output
time <- seq(0, tmax, by=stepout) ## sampling points in time
f.ext<-1e-3  ##declare a value of n that is considered functionally extinct
		

			
				theta.vector <-  c(1.507637, 2.066557, 2.832682, 3.882830, 5.322295, 7.295406, 10.000000 )
					for(theta in theta.vector){ # start theta loop
						print(theta)
						
						k1.vector=exp(seq(log(0.2), log(5), l=30)) #initial max fitness for species 1 to allow for different kappa ratios.  Step are on a log scale.
							for(k1 in k1.vector){ #start k loop
								
								mu1.vector = seq(-theta, theta, l=11)
								#mu1.vector=seq(-theta, 0, length.out=loops.mu) # intial trait value of species one - will loop from -theta to 0 with even steps
									for(mu1 in mu1.vector){ #start mu loop
										
										mu2.vector=seq(-theta, theta, l=11) # initial trait difference between species to change rho - will loop  from 0 to theta-mu1 with even steps
											for(mu2 in mu2.vector){ #start d loop
							
													
													kappa<-c(k1,1) # vector of maximum intrinsic growth rates.  Incurs fitness differences
													params <- list(w=w, h2=h2, bshape=bshape, theta=theta, sigma=sigma, K=kappa) ## param list
													
													muinit <- c(mu1,mu2) ## initial species trait means
													ic <- c(ninit, muinit) ## initial conditions coerced into a vector
										
													sol <- ode(func=eqs, times=time, y=ic, parms=params) %>% ## solve ODEs
													    organize_results(params) ## put results in tidy table
													
													
													#extract niche overlap and fitness differences at 3 time points t=c(0, tmax, ext)
													rk<-sol %>% rhokappa(params)
													initialrk<-rk[time==0,] ##initial niche overlap and fitness differences
													
													tmaxrk<-rk[time==tmax,]  ## rho and kapparatio at tmax
													
													ext.rk.row<-min(c(which(rk$n1<f.ext)-1,which( rk$n2 <f.ext)-1,tmax/stepout+1))  #which row is the first occurence of a functionally extinct population
													extrk<-rk[ext.rk.row,] #rho and kappa ratio at functional extinction
													
													#plot_all(sol,params, max(time), c(-1.1, 1.1, NA), 1001) ## plot results
																			
													var.output<-data.frame(initialrk$'Niche overlap', initialrk$'Competitive difference', extrk$'Niche overlap', extrk$'Competitive difference', mu1, mu2, k1, w, h2, theta, sigma[1], extrk$time,f.ext)
													names(var.output)<-c("nichei", "fiti", "nichef","fitf","mu1","mu2", "k1", "w", "h2", "theta", "sigma" ,"ext.time", "f.ext")
													output<-rbind(output, var.output)
													
												}#end d loop
										}#end mu loop
								}#end k loop
						}#end theta loop
	
		
# and write out the results for posterity's sake
resultsdir <- "../results/"
write.table(
	output,
	file=paste0(resultsdir,"output.csv"),
	sep=",",
	quote=FALSE
)



