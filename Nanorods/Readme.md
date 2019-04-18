# NanoRods analysis

These files are a set of scripts to analyse data for NanoRods and Quantum Dots experiment.
* main.py goes through all the step
* runall_*.py files do each step separately for all files in a set of folders specified  by
given variable

There are several steps.
   1. Extracting ROIs
   2. Producing a movie with selected ROIs
   3. Analysis of blinking NRs
     3.1. Fourier transform and Spectrogram
     3.2. Selecting blinking periods - statistics ON state 
     3.3. ??
   4. Diffusion analysis
     4.1. Fits to 2d gaussians
     4.2. MSD fits
     4.3. Laterality//anisotropy
     4.4. ??


## Getting Started

First do the setup

python setup.py build_ext --inplace

Then, go to main.py, change variable wdirs to contain appropiate directories, with full address,
and run it

### Prerequisites

Some libraries are required:
 - scipy
 - cython
 - openmp
 - numpy
 - matplotlib
 - others

### Installing

## Built With
