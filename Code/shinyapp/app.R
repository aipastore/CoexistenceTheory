
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
## - snap: data generated by organize_results, defined in solve_eqs.R,
##         filtered for one particular moment in time
## - title: the plot label at the top
## Output:
## - a ggplot2 plot
plot_snapshot <- function(snap, title) {
  sp <- c("1", "2") ## the two species
  theta <- snap$theta[1] ## environmental breadth
  traitaxis <- seq(-1.1*theta, 1.1*theta, l=501) ## trait axis
  traits <- expand_grid(species=sp, trait=traitaxis) ## data frame with traits
  traits["density"] <- 0 ## add column for population densities
  for (i in sp) { ## for each species:
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
    ggtitle(title) +
    scale_x_continuous(name="trait value", limits=range(traitaxis)) +
    scale_y_continuous(name="density", limits=c(0, NA)) +
    scale_colour_manual(values=c("#0072B2", "#E69F00")) +
    scale_fill_manual(values=c("#0072B2", "#E69F00")) %>%
    return
}


## create figure that will be displayed in the shiny app
## Input:
## - N1, N2, mu1, mu2: initial densities and trait means for species 1 & 2
## - theta: environmental breadth
## - sig1, sig2: intraspecific trait standard deviations for species 1 & 2
## - h21, h22: trait heritabilities for species 1 & 2
## - w: competition width
## - tmax: time to integrate dynamics for
## Output:
## - a composite ggplot2 plot; if tmax==0, then only the trait distributions
##   are shown, otherwise the time series of the population densities, trait
##   values, and niche overlap and log competitive difference as well
makefig <- function(N1, N2, mu1, mu2, theta, sig1, sig2, h21, h22, K1, w, tmax) {
  params <- list( ## create parameter list
    w=w, ## competition width
    theta=theta, ## environmental breadth
    K=c(K1, 1), ## species' max growth capacities
    sigma=c(sig1, sig2), ## intraspecific trait standard deviations
    h2=c(h21, h22) ## heritabilities
  )
  ic <- c(N1, N2, mu1, mu2) ## initial conditions
  stepout <- tmax/1500 ## time step size for output
  time <- seq(0, tmax, by=stepout) ## sampling points in time
  if (tmax==0) { ## if no duration, only show initial trait distribution
    plt <- ode(func=eqs, y=ic, parms=params, times=c(0, 0.1)) %>% ## "solve" eqs
      organize_results(params) %>% ## put results in tidy table
      filter(time==0) %>% ## filter for initial moment
      plot_snapshot(title="Initial trait distributions") %>%
      plot_grid(nrow=4)
  } else { ## otherwise, show distributions and also time series
    sol <- ode(func=eqs, y=ic, parms=params, times=time) %>% ## solve ODEs
      organize_results(params) ## put results in tidy table
    plt <- plot_grid(plot_snapshot(filter(sol, time==tmax),
                                   title=paste("Snapshot of species' trait",
                                               "distributions at final state")),
                     plot_density(sol),
                     plot_trait(sol),
                     plot_rhokappa(sol, params), ## organize the four sub-figures
                     ncol=1, align="v", axis="r") ## in one column
  }
  return(plt)
}


## ---------------------------- shiny app interface ----------------------------

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      withMathJax(),
      tabsetPanel(
        id="tabs",
        type="tabs",
        tabPanel("Basic",
                 sliderInput("K1", "species 1 growth capacity, \\(K_1\\)",
                             min=0.2, max=5, value=1, step=0.01),
                 sliderInput("w", "competition width, \\(\\omega\\)",
                             min=0.1, max=5, value=2, step=0.01),
                 sliderInput("theta", "environmental breadth, \\(\\theta\\)",
                             min=1, max=10, value=5, step=0.01),
                 sliderInput("mu1", "species 1 initial trait mean, \\(\\mu_1\\)",
                             min=-10, max=10, value=-3, step=0.01),
                 sliderInput("mu2", "species 2 initial trait mean, \\(\\mu_2\\)",
                             min=-10, max=10, value=-4, step=0.01),
                 sliderInput("tmax", "total integration time, \\(t\\)",
                             min=0, max=5000, value=300, step=10)
        ),
        tabPanel("Advanced",
                 sliderInput("sigma1", "species 1 trait std dev, \\(\\sigma_1\\)",
                             min=0.1, max=5, value=1, step=0.01),
                 sliderInput("sigma2", "species 2 trait std dev, \\(\\sigma_2\\)",
                             min=0.1, max=5, value=1, step=0.01),
                 sliderInput("h21", "species 1 heritability, \\(h^2_1\\)",
                             min=0, max=1, value=0.5, step=0.01),
                 sliderInput("h22", "species 2 heritability, \\(h^2_2\\)",
                             min=0, max=1, value=0.5, step=0.01),
                 sliderInput("N1", "species 1 initial density, \\(N_1\\)",
                             min=0.01, max=1, value=0.5, step=0.01),
                 sliderInput("N2", "species 2 initial density, \\(N_2\\)",
                             min=0.01, max=1, value=0.25, step=0.01)
        )
      )
    ),
    mainPanel(
      plotOutput("plt")
    )
  )
)


server <- function(input, output) {
  output$plt <- renderPlot({
    makefig(input$N1, input$N2, input$mu1, input$mu2, input$theta, input$sigma1,
            input$sigma2, input$h21, input$h22, input$K1, input$w, input$tmax)
  }, height=800, width=600)
}


shinyApp(ui=ui, server=server)

