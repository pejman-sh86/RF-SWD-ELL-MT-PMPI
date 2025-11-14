# RF-SWD-ELL-MT-PMPI

This projrct presents probabilistic single- and multi-physics inversion with variable-complexity of one-dimensional receiver function, surface (fundamental mode Rayleigh) wave dispersion, rayleigh wave ellipticity, and magnetotelluric data as decribed in the following paper:

Shahsavari, P., Dettmer, J., Unsworth, M. J., & Schaeffer, A. (2025). Multiphysics inversion with variable complexity of receiver-function, surface-wave dispersion and magnetotelluric data reduces uncertainty for lithosphere structure. arXiv preprint arXiv:2510.15779. (Submitted to Geophysicsl Journal International)

## Installation
NVIDIA HPC SDK must be installed. https://developer.nvidia.com/hpc-sdk \
The code has been written and tested for Linux distributions. \
To install, in the command line, run \
$ git clone https://github.com/pejman-sh86/RF-SWD-ELL-MT-PMPI \
$ cd RF-SWD-ELL-MT-PMPI/src \
$ make clean \
$ make

If compiling and linking is done successfuly, a binary file called "prjmh_temper_rf" is created inside the bin derectorty.

### Run examples
For example to run the MT-SWD-RF PMPI for example 1: \
Go to the directory example1_partial_coupling \
$ cd <path/to/RF-SWD-ELL-MT-PMPI/example1_partial_coupling/> \
$ cd MT_SWD_RF \
$ mpirun -np 12 <path/to/RF-SWD-ELL-MT-PMPI/src/bin/prjmh_temper_r> \
Here 12 is the number of parallel MPI processes. All the results in the paper were produced by using 12 MPI processes. \
To terminate reversible-jump Markov chain Monte Carlo sampling, press CTRL + C. \
To plot the results, put <path/to/RF-SWD-ELL-MT-PMPI/plotting_scripts/> inside your MATLAB path \
Then inside the MATLAB command window run the following script \
rf_plot_rjhist_varpar2.m
