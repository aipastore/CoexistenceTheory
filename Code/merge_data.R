
require(tidyverse)

clargs <- commandArgs(trailingOnly=TRUE)
if (length(clargs)>0) outfile <- clargs[1] else outfile <- "../results/dat.rds"
setwd(Sys.getenv("PWD"))
fs <- Sys.glob("../results/dat_*.csv")
dat <- read_csv(fs[1])
for (f in fs[-1]) {
    tab <- read_csv(f)
    dat <- bind_rows(dat, tab)
}
##write_csv(dat, path=outfile)
saveRDS(dat, file=outfile)

