The scripts and data in this directory can be used to reproduce all figures
from the paper that show results of Clawpack runs.
Note that the default resolution (values of Nx, Ny) set here is much lower
than that used in the paper.  This is to reduce the computational cost
to something feasible on a good workstation in a few hours.  But for most
cases, this is insufficient to accurately capture the dynamics.  To
truly reproduce the figures in the paper, you should set Nx=32, Ny=128.
But for this you will need (as of 2013) at least a large cluster.

Figure names correspond to the arXiv version of the paper.

For example, to reproduce figure 3a:

   > import reproduce
   > reproduce.reproduce_figure('3a',Nx=32, Ny=128)

This puts the simulation results in a subdirectory of computed_data/,
and the figures in the figures/ directory.

You can also just plot the figure from the provided simulation results:

   > import reproduce
   > reproduce.figure_3a(use_saved_data=True)

This uses data stored in a subdirectory of saved_data/.  Note that the saved data
for Figure 14 is stored in .zip files that you must unzip before using.
