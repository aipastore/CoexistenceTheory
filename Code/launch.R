
require(tidyverse)

tab <- expand.grid(K1=exp(seq(log(0.2), log(5), l=30))) %>% ## sp 1 max gr.rates
  as_tibble %>% ## convert to tibble
  mutate(out=sprintf("> ../results/dat_%02d.csv", row_number())) ## output files

for (r in 1:nrow(tab)) {
  jobname <- sprintf("rka%02d.sh", r)
  f <- file(jobname)
  writeLines(c(
    paste("#!/bin/bash"),
    paste("#SBATCH -J", jobname),
    paste("#SBATCH -A liu-2017-00089-2"),
    paste("#SBATCH -t 00:05:00"),
    paste("#SBATCH --mem=8000"),
    paste("#SBATCH -n 1"),
    paste("#"),
    paste("# Run single task in foreground"),
    paste("module add R/3.6.0-nsc1-gcc-7.3.0"),
    paste("cd", Sys.getenv("PWD")),
    paste("Rscript eqb_eqs.R", paste(tab[r,], collapse=" ")),
    paste("#"),
    paste("# Script ends here")
  ), f)
  close(f)
  system(paste("sbatch", jobname))
  Sys.sleep(0.1)
}

