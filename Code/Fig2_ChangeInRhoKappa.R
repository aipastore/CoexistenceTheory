#generates figures 2
require(tidyverse)
require(cowplot)
theme_set(theme_bw())
cpal <- c("#999999","#E69F00","#56B4E9","#009E73","#0072B2","#CC79A7","#D55E00")


#figure 3

require(tidyverse)
require(cowplot)
theme_set(theme_bw())
theme_update(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

dat <- readRDS("../results/eqb_data.rds") %>% ## load data
  filter(singular=="no", w==3) %>%
  mutate(theta=round(theta, 2))

## figure 2
plot_change_niche <- dat %>%
  mutate(drho=nichef-nichei) %>%
  ggplot +
  aes(x=theta, group=theta, y=drho) +
  geom_boxplot(colour="#377eb8", fill="#377eb8", alpha=0.1) +
  xlab(expression(paste("Environmental breadth, ", theta/sigma))) +
  ylab(expression(paste("Change in ", rho)))+
  geom_hline(yintercept=0, alpha=0.5)

plot_change_comp_diff <- dat %>%
  mutate(dkappa=abs(log(fitf)-log(fiti))) %>%
  ggplot +
  aes(x=theta, group=theta, y=dkappa) +
  geom_boxplot(colour="#377eb8", fill="#377eb8", alpha=0.1) +
  xlab(expression(paste("Environmental breadth, ", theta/sigma))) +
  ylab(expression(paste("Change in ", log(kappa[1]/kappa[2]))))

plot_grid(plot_change_niche, plot_change_comp_diff, ncol=2, align="hv",
          labels=c("A", "B"))

