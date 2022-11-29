# Spatial ecDNA pattern simulations, M J Haughey 2022

This repository contains the code for simulating the spatial patterns of ecDNA content within expanding tumour cell populations. The program generates an output file tissue.csv containing a list of the (x,y)-coordinates of all cells in the system, and the number of ecDNA copies the carry. Tested using Apple clang version 14.0.0 and Python v3.11.0.

Execute a simulation by running:

```
g++ -o ./spatial_ecDNA_patterns_2D ./spatial_ecDNA_patterns_2D.cpp
./spatial_ecDNA_patterns_2D [--verbose] [--hubs] [--copyNumberDependentSelection] [-q Q] [-n N] [-s S] [-x X]
```

where\
&nbsp; --verbose &emsp;&emsp; verbose flag (optional)\
&nbsp; --hubs &emsp;&emsp; ecDNA clustering flag (optional)\
&nbsp; --copyNumberDependentSelection &emsp;&emsp; flag to set ecDNA copy number dependent selection model. Default is flat selection based simply on ecDNA+ or ecDNA- status of cell (optional)\
&nbsp; -q &emsp;&emsp; cell pushing strength\
&nbsp; -n &emsp;&emsp; ecDNA copy number in initial cell\
&nbsp; -s &emsp;&emsp; selection strength. Scalar multiplier for the ecDNA dependent birth rate function, with s=0 giving rise to neutral growth and s>0 giving rise to positive ecDNA copy number dependent selection\
&nbsp; -x &emsp;&emsp; random seed

Flags should be specified before numerical arguments. Plot the final spatial data in tissue.csv by running:

```
python3 ./plot_spatial_data.py [--path PATH] [--plotCopyNumber]
```

where [PATH] is the file path the relevant tissue.csv file, and [--plotCopyNumber] specifies colouring cells according to their ecDNA copy number (if flag not specified, then cells are coloured according to ecDNA+ or ecDNA- status).
