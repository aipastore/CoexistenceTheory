## This script generates a heat map of the likelihood of coexistence,
## given the change in niche overlap (abscissa) and the change in
## competitive difference (ordinate) from the initial to the final state.
## It is Figure 3 in the manuscript.

require(tidyverse) ## for efficient data manipulation and plotting
theme_set(theme_bw()) ## set default theme
theme_update(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

dat <- readRDS("../Results/alldata.rds") ## load data

res <- 15 ## grid resolution of heat map
dat %>%
  filter(fiti>0, fitf>0) %>% ## only meaningful (positive) comp. differences
  ## convert final coexistence from logical to numeric; calculate change in
  ## niche overlap and log competitive diffs from initial to the final state
  mutate(coexf=1*coexf, dn=nichef-nichei, df=log(fitf)-log(fiti)) %>%
  ## classify niche overlap values into bins
  mutate(dn=cut(dn, breaks=seq(-1, 1, l=res+1),
                labels=seq(-1, 1, l=res),
                right=FALSE) %>% as.character %>% as.numeric) %>%
  ## classify competitive difference values into bins
  mutate(df=cut(df, breaks=seq(-6, 6, l=res+1),
                labels=seq(-6, 6, l=res),
                right=FALSE) %>% as.character %>% as.numeric) %>%
  ## determine the likelihood of coexistence in each bin
  group_by(dn, df) %>%
  summarise(coex=mean(coexf)) %>%
  ## create plot
  ggplot +
  aes(x=dn, y=df, fill=coex) +
  geom_raster(interpolate=TRUE, na.rm=TRUE) +
  scale_fill_continuous(type="viridis", name="probability of\ncoexistence") +
  scale_x_continuous(name=expression(paste("change in ", rho)), expand=c(0, 0),
                     label=abbreviate) +
  scale_y_continuous(name=expression(paste("change in ", log(kappa[1]/kappa[2]))),
                     expand=c(0, 0))
