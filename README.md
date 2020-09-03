# CoexistenceTheory

Code and data to replicate results for the manuscript "The evolution of niche overlap and competitive differences".

In the `Code` directory:
* `get_data.R`: The main script that runs the model over many different combinations of parameters and initial conditions, and generates the model's output. By default, this output is dumped on the screen; redirect the output to a file to save the results (see the comments within the script as well). The result of this run, compressed into an `.rds` file, can be found in the `Data` directory.
* `functions.R`: Functions aiding numerical integration of the ordinary differential equations (ODEs), obtaining niche overlap and competitive differences from the results, and organizing them into a tidy output format. It is used by `get_data.R`, `sample_dynamics.R`, and `trajectories.R`.
* `boxplots.R`, `phaseplots.R`, `predict_coex.R`, `sampe_dynamics.R`, and `trajectories.R`: These scripts generate the figures in the manuscript.

In the `Data` directory:
* `alldata.rds`: The output generated by `get_data.R`, compressed into `.rds` format to reduce size.
