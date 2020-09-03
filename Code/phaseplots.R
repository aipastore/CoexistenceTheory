## This script shows final niche overlap and competitive difference
## values, as a funcion of model parameters and initial conditions
## (Figure 4 in the manuscript). The results are in fact independent
## of initial conditions and depend only on the parameters.

require(tidyverse) ## for efficient data manipulation and plotting
theme_set(theme_bw()) ## set default theme
theme_update(panel.grid.major=element_blank(), panel.grid.minor=element_blank())


dat <- readRDS("../Results/alldata.rds") %>% ## load data
  mutate(theta=round(theta, 2)) %>% ## round theta values
  filter(w>0.5) %>% ## drop lowest environmental breadth value
  filter(coexf) ## only those with coexistence at the end


separatrix <- tibble(rho=seq(0,1,by=0.01)) %>% ## for plotting coex. regions
  mutate(lower=log(rho), upper=log(rho)*(-1)) ## calculate upper & lower limits

## create plot
dat %>%
  ggplot +
  geom_line(data=separatrix, aes(x=rho, y=upper), colour="#0072B2") +
  geom_line(data=separatrix, aes(x=rho, y=lower), colour="#0072B2") +
  geom_ribbon(data=separatrix, aes(x=rho, ymin=upper, ymax=max(separatrix$upper)),
              colour=NA, fill="#0072B2", alpha=0.1) +
  geom_ribbon(data=separatrix, aes(x=rho, ymax=lower, ymin=min(separatrix$lower)),
              colour=NA, fill="#0072B2", alpha=0.1) +
  geom_hline(yintercept=0, colour="#0072B2") +
  geom_point(aes(x=nichef, y=log(fitf), colour=theta), size=0.9, na.rm=TRUE) +
  geom_path(aes(x=nichef, y=log(fitf), group=theta, colour=theta), alpha=0.3) +
  scale_x_continuous(name=expression(paste(rho)), limits=c(0, 1),
                     expand=c(0, 0), labels=abbreviate) +
  scale_y_continuous(name=expression(paste(log(kappa[1]/kappa[2]))),
                     limits=c(-3, 3), expand=c(0, 0)) +
  scale_colour_continuous(type="viridis", name=expression(paste(theta))) +
  facet_wrap(~w, nrow=1, labeller=label_bquote(omega==.(w)))
