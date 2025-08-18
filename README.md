# SPatial ECdna Intratumor Evolution Simulation (SPECIES)

This repository contains the code for simulating the spatial patterns of ecDNA content within expanding tumour cell populations. The program generates an output file tissue.csv containing a list of the (x,y)-coordinates of all cells in the system, and the number of ecDNA copies the carry. Tested using Apple clang version 15.0.0 and Python v3.11.0.

Execute a single simulation by running:

```
g++ -o ./spatial_ecDNA_patterns_2D ./spatial_ecDNA_patterns_2D.cpp
./spatial_ecDNA_patterns_2D [--verbose] [-n N] [-k K] [-s S] [-q Q] [-x X]
```

where\
&nbsp; --verbose &emsp;&emsp; verbose flag (optional)\
&nbsp; -n &emsp;&emsp; maximum tumour size (in number of cells)\
&nbsp; -q &emsp;&emsp; cell pushing strength\
&nbsp; -k &emsp;&emsp; ecDNA copy number in initial cell\
&nbsp; -s &emsp;&emsp; selection strength. Scalar multiplier for the ecDNA dependent birth rate function, with s=0 giving rise to neutral growth and s>0 giving rise to positive ecDNA copy number dependent selection\
&nbsp; -x &emsp;&emsp; random seed

Flags should be specified before numerical arguments.

Run multiple simulations in a batch by editing the variables in the batch_run_parameters.txt file, then running:

```
bash ./batch_run.sh
```

Plot the final spatial data in tissue.csv by running:

```
python3 ./plot_spatial_data.py [--path PATH] [--plotCopyNumber]
```

where [PATH] is the file path the relevant tissue.csv file, and [--plotCopyNumber] specifies colouring cells according to their ecDNA copy number (if flag not specified, then cells are coloured according to ecDNA+ or ecDNA- status).

To perform parameter fitting (approximate Bayesian computation with rejection sampling) of the SPECIES spatial model to multi-region ecDNA copy number dsitributions (GB-UK cohort, data found in /GB_UK_ecDNA_copy_number_distributions/), run the following command:

```
bash ./ABC_full_algorithm.sh
```

Uniform prior distributions for parameters s and k are assumed, and the list of q values to simulate must be explicitly specified. Parameters of the prior distributions can be edited in the ABC_input_parameters.txt file.

