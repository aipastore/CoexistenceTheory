
require(shiny)
require(deSolve)
require(tidyverse)
require(cowplot)
theme_set(theme_cowplot())


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
##         $K: vector of species-specific intrinsic growth rates
##         $sigma: vector of species' intraspecific standard deviations
##         $h2: heritability of trait (assumed equal for all species)
## Output:
## - time derivatives of the state variables coerced into a single vector
eqs <- function(time, state, pars) {
  n <- state[1:2] ## species densities
  m <- state[3:4] ## species trait means
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
##         $K: vector of species-specific intrinsic growth rates
##         $sigma: vector of species' intraspecific standard deviations
##         $h2: heritability of trait (assumed equal for all species)
## Output:
## - a tibble with columns: time, species, n (density), m (trait mean),
##   sigma, w, theta, and h2 (heritabilities)
organize_results <- function(sol, pars) {
  dat <- sol %>% as.data.frame %>% as_tibble ## convert to tibble
  names(dat)[1] <- "time" ## name the first column "time"
  names(dat)[2:3] <- paste0("n_", 1:2) ## name abundance columns (n_k)
  names(dat)[4:5] <- paste0("m_", 1:2) ## name trait mean columns (m_k)
  dat <- dat %>%
    gather("variable", "v", 2:ncol(dat)) %>% ## long format
    separate(variable, c("type", "species"), sep="_") %>%
    spread(type, v) %>% ## separate columns for densities n and trait means m
    select(time, species, n, m) %>% ## rearrange columns
    mutate(species=as.integer(species), sigma=pars$sigma[species], w=pars$w,
           theta=pars$theta, h2=pars$h2[species]) %>% ## add parameters
    mutate(species=as.character(species))
  return(dat)
}


## Obtain niche overlap rho and the kappa ratio through time
## Input:
## - sol: output produced by the function ode
## - pars: list of parameters, with the following elements:
##         $w: width of competition kernel
##         $theta: width of intrinsic growth function
##         $K: vector of species-specific intrinsic growth rates
##         $sigma: vector of species' intraspecific standard deviations
##         $h2: heritability of trait (assumed equal for all species)
## Output:
## - a tibble with columns: time, niche overlap, fitness ratio
rhokappa <- function(sol, pars) {
  m1 <- sol %>% filter(species==1) %>% pull(m)
  m2 <- sol %>% filter(species==2) %>% pull(m)
  n1 <- sol %>% filter(species==1) %>% pull(n)
  n2 <- sol %>% filter(species==2) %>% pull(n)
  time <- sol %>% pull(time) %>% unique
  sv <- outer(pars$sigma^2, pars$sigma^2, FUN="+")
  rho <- rep(0, length(m1))
  kapparatio <- rep(0, length(m1))
  for (i in 1:length(m1)) {
    mm <- c(m1[i], m2[i])
    dm <- outer(mm, mm, FUN="-")
    alpha <- exp(-dm^2/(2*sv+pars$w^2))*pars$w/sqrt(2*sv+pars$w^2)
    b <- pars$K-(mm^2+pars$sigma^2)/pars$theta^2
    rho[i] <- sqrt(alpha[1,2]*alpha[2,1]/(alpha[1,1]*alpha[2,2]))
    kapparatio[i] <- (b[1]/b[2])*sqrt(alpha[2,1]*alpha[2,2]/
                                        (alpha[1,2]*alpha[1,1]))
  }
  return(tibble(time=time, `niche overlap`=rho,
                `log competitive difference`=log(kapparatio)))
}


## Plot niche overlap and competitive difference through time
## Input:
## - dat: data generated by organize_results, defined in solve_eqs.R
## - parms: list of parameters
## Output:
## - a ggplot2 plot
plot_rhokappa <- function(dat, parms){
  rhokappa(dat, parms) %>%
    gather(quantity, value, c(`niche overlap`, `log competitive difference`)) %>%
    mutate(quantity=factor(quantity, ordered=TRUE, levels=c(
      "niche overlap", "log competitive difference"))) %>%
    ggplot() +
    aes(x=time, y=value, colour=quantity) +
    geom_line(na.rm=TRUE) +
    ggtitle("Niche overlap and log competitive difference through time") +
    scale_colour_manual(values=c("#0072B2", "#E69F00")) %>%
    return
}


## Plot species densities through time
## Input:
## - dat: data generated by organize_results, defined in solve_eqs.R
## Output:
## - a ggplot2 plot
plot_density <- function(dat) {
  ggplot(dat) +
    geom_line(aes(x=time, y=n, colour=species)) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    ggtitle("Population densities through time") +
    scale_colour_manual(values=c("#0072B2", "#E69F00")) %>%
    return
}


## Plot species trait values through time
## Input:
## - dat: data generated by organize_results, defined in solve_eqs.R
## Output:
## - a ggplot2 plot
plot_trait <- function(dat) {
  ggplot(dat) +
    aes(x=time, y=m, ymin=m-sigma, ymax=m+sigma, colour=species, fill=species) +
    geom_ribbon(colour=NA, alpha=0.15) +
    geom_line() +
    ylab("trait value") +
    ggtitle("Trait evolution through time") +
    scale_colour_manual(values=c("#0072B2", "#E69F00")) +
    scale_fill_manual(values=c("#0072B2", "#E69F00")) %>%
    return
}


## Make a plot of the trait distributions at some given moment
## Input:
## - dat: data generated by organize_results, defined in solve_eqs.R
## - moment: at which time point should the community state be plotted
## Output:
## - a ggplot2 plot
plot_snapshot <- function(dat, moment=0) {
  sp <- dat %>% pull(species) %>% unique ## set of species
  theta <- dat$theta[1] ## environmental breadth
  traitaxis <- seq(-1.1*theta, 1.1*theta, l=501) ## trait axis
  snap <- dat %>% filter(time==moment) %>% select(-time) ## time = moment
  traits <- expand_grid(species=sp, trait=traitaxis) ## data frame with traits
  traits["density"] <- 0 ## add column for population densities
  for (i in sp) {
    v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits$density[(traits$species==i)] <- v$n* ## trait normally distributed,
      dnorm(traits$trait[(traits$species==i)], v$m, v$sigma) ## times density
  }
  landscape <- tibble(trait=traitaxis, r=0) ## for plotting intrinsic rates
  landscape$r <- 1-landscape$trait^2/theta^2
  landscape$r[landscape$r<0] <- NA
  traits$density[traits$density<1e-4] <- NA
  landscape$r <- landscape$r*max(traits$density, na.rm=TRUE) ## scale to plot
  ggplot(traits) + ## generate figure
    geom_line(aes(x=trait, y=density, colour=species), na.rm=TRUE) +
    geom_ribbon(aes(x=trait, ymin=0, ymax=density, fill=species), alpha=0.15) +
    geom_line(data=landscape, aes(x=trait, y=r), linetype="dashed",
              colour="#999999", alpha=0.75, na.rm=TRUE) +
    ggtitle("Snapshot of species' trait distributions at final state") +
    scale_x_continuous(name="trait value", limits=range(traitaxis)) +
    scale_y_continuous(name="density", limits=c(0, NA)) +
    scale_colour_manual(values=c("#0072B2", "#E69F00")) +
    scale_fill_manual(values=c("#0072B2", "#E69F00")) %>%
    return
}


makefig <- function(N1, N2, mu1, mu2, theta, sig1, sig2, h21, h22, K1, w, tmax) {
  params <- list( ## create parameter list
    w=w,
    theta=theta,
    K=c(K1, 1),
    sigma=c(sig1, sig2),
    h2=c(h21, h22)
  )
  ic <- c(N1, N2, mu1, mu2) ## initial conditions
  stepout <- tmax/1500 ## time step size for output
  time <- seq(0, tmax, by=stepout) ## sampling points in time
  sol <- ode(func=eqs, y=ic, parms=params, times=time) %>% ## solve ODEs
    organize_results(params) ## put results in tidy table
  return(plot_grid(plot_snapshot(sol, tmax),
                   plot_density(sol),
                   plot_trait(sol),
                   plot_rhokappa(sol, params),
                   ncol=1, align="v", axis="r"))
}


## ---------------------------- shiny app interface ----------------------------

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      withMathJax(),
      sliderInput("K1", "species 1 growth (dis)advantage",
                  min=0.2, max=5, value=1, step=0.01),
      sliderInput("w", "competition width",
                  min=0.1, max=5, value=2, step=0.01),
      sliderInput("theta", "environmental breadth",
                  min=1, max=10, value=5, step=0.01),
      sliderInput("mu1", "species 1 initial trait mean",
                  min=-10, max=10, value=-0.6, step=0.01),
      sliderInput("mu2", "species 2 initial trait mean",
                  min=-10, max=10, value=-0.7, step=0.01),
      sliderInput("tmax", "simulation time",
                  min=1, max=5000, value=1000, step=10),
      sliderInput("N1", "species 1 initial density",
                  min=0.01, max=1, value=0.5, step=0.01),
      sliderInput("N2", "species 2 initial density",
                  min=0.01, max=1, value=0.25, step=0.01),
      sliderInput("sigma1", "species 1 trait std dev",
                  min=0.1, max=5, value=1, step=0.01),
      sliderInput("sigma2", "species 2 trait std dev",
                  min=0.1, max=5, value=1, step=0.01),
      sliderInput("h21", "species 1 heritability",
                  min=0, max=1, value=0.5, step=0.01),
      sliderInput("h22", "species 2 heritability",
                  min=0, max=1, value=0.5, step=0.01)
    ),
    mainPanel(
      plotOutput("model_plot")
    )
  )
)


server <- function(input, output) {
  output$model_plot <- renderPlot({
    makefig(input$N1, input$N2, input$mu1, input$mu2, input$theta, input$sigma1,
            input$sigma2, input$h21, input$h22, input$K1, input$w, input$tmax)
  }, height=800, width=600)
}


shinyApp(ui=ui, server=server)

