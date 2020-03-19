
require(tidyverse)
require(cowplot)
theme_set(theme_bw())
theme_update(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

all.dat <- read.csv("../results/output.csv") %>% ## load data
mutate(thetasig = theta/sigma) %>%
mutate(thetasig=round(thetasig,2)) 

all<-sort(unique(all.dat$thetasig)) 


eq.dat <- readRDS("../results/eqb_data.rds") %>% ## load data
 filter(singular=="no", w==3) %>%
 mutate(thetasig=round(theta, 2)) %>%
 
 eq<-sort(unique(eq.dat$thetasig))
 select.eq<-which(eq %in% all)
 

## figure 4
plot_change_niche_detail <- eq.dat %>%
  filter(thetasig %in% eq[ select.eq]) %>%
  ggplot +
  geom_point(aes(x=log(fiti), y=log(fitf), colour=nichei), alpha=0.3, size=0.5) +
   coord_cartesian(xlim =c(-4, 4), ylim = c(-7, 7))+
  geom_hline(yintercept=0, linetype="dashed", color = "grey")+
  geom_vline(xintercept=0, linetype="dashed", color = "grey")+
  geom_abline(intercept=0, slope=1)+
  facet_wrap(~theta, nrow=1, labeller=label_bquote(theta/sigma==.(theta))) +
  scale_x_continuous(name=expression(paste("initial ", log(kappa[1]/kappa[2]))),
                     labels=abbreviate) +
  scale_y_continuous(name=expression(paste("final ", log(kappa[1]/kappa[2]))),
                     labels=abbreviate) +
  scale_colour_continuous(type="viridis", name=expression(paste("initial ", rho)))

plot_change_comp_detail <- all.dat %>%
  filter(thetasig %in% all) %>%
  ggplot +
  geom_point(aes(x=log(fiti), y=log(fitf), colour=nichei), alpha=0.3, size=0.5) +
  coord_cartesian(xlim =c(-4, 4), ylim = c(-7, 7))+
  geom_hline(yintercept=0, linetype="dashed", color = "grey")+
  geom_vline(xintercept=0, linetype="dashed", color = "grey")+
  geom_abline(intercept=0, slope=1)+
  facet_wrap(~theta, nrow=1, labeller=label_bquote(theta==.(theta))) +
  scale_x_continuous(name=expression(paste("initial ", log(kappa[1]/kappa[2]))),
                     labels=abbreviate) +
  scale_y_continuous(name=expression(paste("final ", log(kappa[1]/kappa[2]))),
                     labels=abbreviate) +
  scale_colour_continuous(type="viridis", name=expression(paste("initial ", rho)))

plot_grid(plot_change_niche_detail, plot_change_comp_detail, ncol=1, align="hv",
          labels=c("A", "B"))

#ggsave("../figures/Stable_vs_All_initial_vs_final.png", width=7, height=4, dpi=300)

