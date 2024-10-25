# SPatial ECdna Intratumor Evolution Simulation (SPECIES)

This repository contains the code for simulating the spatial patterns of ecDNA content within expanding tumour cell populations. The program generates an output file tissue.csv containing a list of the (x,y)-coordinates of all cells in the system, and the number of ecDNA copies the carry. Tested using Apple clang version 15.0.0 and Python v3.11.0.

Execute a simulation by running:

```
g++ -o ./spatial_ecDNA_patterns_2D ./spatial_ecDNA_patterns_2D.cpp -std=c++20
./spatial_ecDNA_patterns_2D [--verbose] [-q Q] [-n N] [-s S] [-x X]
```

where\
&nbsp; --verbose &emsp;&emsp; verbose flag (optional)\
&nbsp; -q &emsp;&emsp; cell pushing strength\
&nbsp; -n &emsp;&emsp; ecDNA copy number in initial cell\
&nbsp; -s &emsp;&emsp; selection strength. Scalar multiplier for the ecDNA dependent birth rate function, with s=0 giving rise to neutral growth and s>0 giving rise to positive ecDNA copy number dependent selection\
&nbsp; -x &emsp;&emsp; random seed

Flags should be specified before numerical arguments. Plot the final spatial data in tissue.csv by running:

```
python3 ./plot_spatial_data.py [--path PATH] [--plotCopyNumber]
```

where [PATH] is the file path the relevant tissue.csv file, and [--plotCopyNumber] specifies colouring cells according to their ecDNA copy number (if flag not specified, then cells are coloured according to ecDNA+ or ecDNA- status).


Developed for and implemented in Noorani & Haughey et al. (2024):

Imran Noorani<sup>*1,2,3</sup>, Magnus Haughey<sup>*4</sup>, Jens Luebeck<sup>5</sup>, Andrew Rowan<sup>1</sup>, Eva Grönroos<sup>1</sup>, Francesco Terenzi<sup>4</sup>, Ivy Tsz-Lo Wong<sup>6,7</sup>, Jeanette Kittel<sup>8</sup>, Chris Bailey<sup>1</sup>, Clare Weeden<sup>1</sup>, Donald Bell<sup>9</sup>, Eric Joo<sup>3</sup>, Vittorio Barbe<sup>1</sup>, Matthew G. Jones<sup>10</sup>, Emma Nye<sup>11</sup>, Mary Green<sup>11</sup>, Lucy Meader<sup>11</sup>, Emma Jane Norton<sup>12,13</sup>, Mark Fabian<sup>12</sup>, Nnennaya Kanu<sup>14</sup>, Mariam Jamal-Hanjani<sup>14, 15, 16</sup>, Thomas Santarius<sup>17</sup>, James Nicoll<sup>12,13</sup>, Delphine Boche<sup>13</sup>, Howard Y Chang<sup>18,19</sup>, Vineet Bafna<sup>20,21</sup>, Weini Huang<sup>22,23</sup>, Paul S Mischel+<sup>6,7</sup>, Charles Swanton+<sup>1,14,24</sup>, Benjamin Werner+<sup>4</sup>

1. Cancer Evolution and Genome Instability Laboratory, The Francis Crick Institute, London, United Kingdom
2. Department of Neurosurgery, National Hospital for Neurology and Neurosurgery, London, United Kingdom
3. Institute of Neurology, University College London, London, United Kingdom
4. Evolutionary Dynamics Group, Centre for Cancer Genomics and Computational Biology, Barts Cancer Centre, Queen Mary University of London, London, United Kingdom
5. Department of Computer Science and Engineering, University of California, San Diego, La Jolla, United States
6. Sarafan ChEM-H, Stanford University, Stanford, United States
7. Department of Pathology, Stanford University, Stanford, United States
8. UCL Cancer Institute, London, United Kingdom
9. Crick Advanced Light Microscopy, The Francis Crick Institute, London, United Kingdom
10. Center for Personal Dynamic Regulomes, Stanford University, Stanford, United States
11. Experimental Histopathology, the Francis Crick Institute, United Kingdom
12. Department of Cellular Pathology, University Hospital Southampton NHS Foundation Trust, Southampton, United Kingdom
13. Clinical Neurosciences, Clinical and Experimental Sciences, Faculty of Medicine, University of Southampton, Southampton, United Kingdom
14. Cancer Research UK Lung Cancer Centre of Excellence, University College London Cancer Institute, London, United Kingdom
15. Cancer Metastasis Laboratory, University College London Cancer Institute, London, United Kingdom
16. Department of Medical Oncology, University College London Hospitals, London, United Kingdom
17. Department of Neurosurgery, Cambridge University Hospital, Cambridge
18. Departments of Dermatology and Genetics, Stanford University School of Medicine, Stanford, United States
19. Howard Hughes Medical Institute, Stanford University, Stanford, United States
20. Department of Computer Science and Engineering, University of California, San Diego, La Jolla, United States
21. Halıcıoğlu Data Science Institute, University of California, San Diego, La Jolla, United States
22. Department of Mathematics, Queen Mary University of London, London, United Kingdom
23. Group of Theoretical Biology, The State Key Laboratory of Biocontrol, School of Life Science, Sun Yat-sen University, Guangzhou, China
24. Department of Oncology, University College London Hospitals, London, United Kingdom

* These authors contributed equally to this work.


## Abstract

Extrachromosomal DNA (ecDNA) oncogene amplification is associated with treatment resistance and shorter survival in cancer. Currently, the spatial dynamics of ecDNA, and their evolutionary impact, are poorly understood. Here, we investigate ecDNA spatial-temporal evolution by integrating computational modeling with samples from 94 treatment-naive human IDH-wildtype glioblastoma patients. Random ecDNA segregation combined with ecDNA-conferred fitness advantages induce predictable spatial ecDNA copy-number patterns which depend on ecDNA oncogenic makeup. EGFR-ecDNAs often reach high copy-number, confer strong fitness advantages and do not co-amplify other oncogenes on the same ecDNA. In contrast, PDGFRA-ecDNAs reach lower copy-number, confer weaker fitness advantages and co-amplify other oncogenes. EGFR structural variants occur exclusively on ecDNA, arise from and are intermixed with wild-type EGFR-ecDNAs. Modeling suggests wild-type and variant EGFR-ecDNAs often accumulate before clonal expansion, even in patients co-amplifying multiple ecDNA species. Our results suggest a potential time window in which early ecDNA detection may facilitate more effective intervention.

