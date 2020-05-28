# Hidden treat repo information

This repository contains code replicating the numerical simulations and figures in:

'Causal inference from treatment-control studies having an additional factor with unknown assignment mechanism' (URL TBD.)

There are two R files:

- simulation.R
- simulation_functions.R

Inputs:

The main file is simulation.R.  To run the simulations, change the base directory at the top of the file to point
to this repo locally, such as:

```
base.dir = '~/hidden_treat/'
```
Sourcing the file will then run the simulations. This file sources simulation_functions.R, which contains
helper functions for the simulations.

Outputs:

Sourcing the simulation file reproduces the plots in the paper and saves them into the 'plots' subdirectory.

