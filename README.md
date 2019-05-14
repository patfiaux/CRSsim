# RELICS-performance
This github account contains code and instructions for simulating CRISPR regulatory screens and assessing performance of different methods used in CRISPR regulatory screens

# 1. Simulating CRISPR regulatory screen data
The simulations mirror the experimental steps by taking the following variables into account:

| Experimental step | Simulated step |
| ----------------- | -------------- |
| NA | Size and number of enhancers|
| NA | Strength of enhancer |
| NA | Guide efficiency |
| Transduction of cells with guide library | Generate guide distribution |
| Load sorter with cells* | Specify number of cells to sort |
| Sort cells** | Specify sorting probabilities for sorting using Dirichlet-Multinomial |
| Specify sequencing depth | Specify sequencing depth |
| Sequence pools | Specify whether PCR-duplicates can be detected or not |

\* If a selection screen is performed the user would specify the number of cells in the ‘before’ pool

** dropout rate

## 1.1 Installations and Setup
The simulations are done in R [R](https://cran.r-project.org/bin/windows/base/). Please make sure you have R version 3.5.1 or higher

Download source code to your desired location: `git clone https://github.com/patfiaux/RELICS.git`

To run RELCIS you need the packages below. If you don't have them, install them using the command after the '#'):
### R packages
```
MCMCpack # 
```




