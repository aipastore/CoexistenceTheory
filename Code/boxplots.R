## This script generates box plots of the change in niche overlap and
## competitive differences, as a function of the environmental breadth.
## It is Figure 2 in the manuscript.

require(tidyverse) ## for efficient data manipulation and plotting
require(cowplot) ## for joining and labeling plots in a grid
theme_set(theme_bw()) ## set default theme
theme_update(panel.grid.major=element_blank(), panel.grid.minor=element_blank())


dat <- readRDS("../Data/alldata.rds") %>% ## load data
  filter(w==3, coexf, fiti>0) ## only those with w=3 and coexistence at the end

## plot change in niche overlap against environmental breadth
p_rho <- dat %>%
  mutate(drho=nichef-nichei) %>%
  ggplot +
  aes(x=theta, group=theta, y=drho) +
  geom_boxplot(colour="#0072B2", fill="#0072B2", alpha=0.1, coef=2) +
  xlab(expression(paste("environmental breadth, ", theta))) +
  ylab(expression(paste("change in ", rho)))+
  geom_hline(yintercept=0, alpha=0.5)

## plot change in competitive differences against environmental breadth
p_kratio <- dat %>%
  mutate(dkappa=abs(log(fitf)-log(fiti))) %>%
  ggplot +
  aes(x=theta, group=theta, y=dkappa) +
  geom_boxplot(colour="#0072B2", fill="#0072B2", alpha=0.1, coef=2) +
  xlab(expression(paste("environmental breadth, ", theta))) +
  ylab(expression(paste("change in ", log(kappa[1]/kappa[2]))))

## join and label the two plots
plot_grid(p_rho, p_kratio, ncol=2, align="hv", labels=c("A", "B"))
