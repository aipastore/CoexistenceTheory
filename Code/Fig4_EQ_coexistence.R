
require(tidyverse)

separatrix <- tibble(rho=seq(0,1,by=0.01)) %>% ## for plotting coex. regions
  mutate(lower=log(rho), upper=log(rho)*(-1)) ## calculate upper & lower limits

readRDS("../results/eqb_data.rds") %>% ## load data
  mutate(theta=round(theta, 2))%>%
  filter(nfinal1>1e-8, nfinal2>1e-8) %>%
  filter(singular=="no") %>% ## only nonsingular outcomes
  filter(fitf>0) %>% ## drop meaningless nonpositive final fitness ratios
  filter(nichei<1)%>%
  filter(w>0.5)%>%

  ggplot + ## start plotting
  geom_line(data=separatrix, aes(x=rho, y=upper), colour="#377eb8") +
  geom_line(data=separatrix, aes(x=rho, y=lower), colour="#377eb8") +
  geom_ribbon(data=separatrix, aes(x=rho, ymin=upper, ymax=max(separatrix$upper)),
              colour=NA, fill="#377eb8", alpha=0.1) +
  geom_ribbon(data=separatrix, aes(x=rho, ymax=lower, ymin=min(separatrix$lower)),
              colour=NA, fill="#377eb8", alpha=0.1) +
  geom_hline(yintercept=0, colour="#377eb8") +
  geom_point(aes(x=nichef, y=log(fitf), colour=theta),
             alpha=0.5, size=0.75) + ## plot data points
  scale_x_continuous(name=expression(paste(rho)), limits=c(0, 1),
                     expand=c(0, 0), labels=abbreviate) +
  scale_y_continuous(name=expression(paste(log(kappa[1]/kappa[2]))),
                     limits=c(-3, 3), expand=c(0, 0)) +
  scale_colour_continuous(type="viridis", name=expression(paste(theta/sigma))) +
  facet_wrap(~w, nrow=1, labeller=label_bquote(omega==.(w))) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

## ggsave("../figures/Fig5_Eq_Coexistence.png", width=7.2, height=3.2, dpi=300)
