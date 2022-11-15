# Spatial ecDNA pattern simulations, M J Haughey 2022

This repository contains the code for simulating the spatial patterns of ecDNA content within expanding tumour cell populations. The program generates an output file tissue.csv containing a list of the (x,y)-coordinates of all cells in the system, and the number of ecDNA copies the carry. Tested using Apple clang version 14.0.0 and Python v3.11.0.

Execute a simulation by running:

```
g++ -o ./spatial_ecDNA_patterns_2D ./spatial_ecDNA_patterns_2D.cpp
./spatial_ecDNA_patterns_2D [-v1] [-q Q] [-x X]
```

where\
&nbsp; -v1 &emsp;&emsp; verbose flag (optional)\
&nbsp; -q &emsp;&emsp; cell pushing strength\
&nbsp; -x &emsp;&emsp; random seed

Plot the final spatial data in tissue.csv by running:

```
python3 ./plot_spatial_data.py [--path PATH]
```

where [PATH] is the file path the relevant tissue.csv file.
