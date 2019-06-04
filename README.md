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

## 1.1 Installations and Setup for Simulations
The simulations are done in R [R](https://cran.r-project.org/bin/windows/base/). Please make sure you have R version 3.5.1 or higher

Download source code to your desired location: `git clone https://github.com/patfiaux/RELICS.git`

To simulate data you need the packages below. If you don't have them, install them using the command after the '#'):
### R packages
```
MCMCpack # install.packages('MCMCpack')

transport # install.packages('transport')

```

### Bioconductor packages
```
if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")
    
IRanges # BiocManager::install("IRanges", version = "3.8")

GenomicRanges # BiocManager::install("GenomicRanges", version = "3.8")
```

## Simulation quickstart with example data
### 1. source the script
```
source('/path/to/script/RELICS_sim.r')
```

### 2. Setting up simulation flags. 
There are several different flags which have no defaults and need to be supplied by the user. Below is the outline on how to set the most important flags to get the simulation going.

Flags are set up in a list format

```
sim.flags <- list()
```

Set the output name of the simulation
```
sim.flags$simName <- 'Example_simulation'
```

Provide info about the guide targets. Either supply them directly, as is done here, or generate them (see details 'Advanced Simulations')
```
sim.flags$guides <- example.info
```

Provide the original guide distribution for each replicate. This can be supplied from an existing data set, as is done below. Another option is to generate the distributions using a zero inflated negative binomial distribution. See 'Advanced Simulations' for details.
```
sim.flags$inputGuideDistr <- cbind(before_1 = example.counts$before_repl1, 
  before_2 = example.counts$before_repl2)
```  

```
sim.flags$poolNames <- c('before', 'after')
sim.flags$exon <- example.gene
sim.flags$nrEnhancers <- 5
sim.flags$enhancerSize <- 50  # base pairs

sim.flags$crisprSystem <- 'CRISPRi'

sim.flags$seqDepth <- list(repl1 = c(16656607, 19431422),
    repl2 = c(20217155, 21585515))
sim.flags$pcrDupl <- 'duplicate'  # keep PCR duplicates in this case to reproduce the shape of Fulco's distribution

sim.flags$guideEfficiency <- 'high'
sim.flags$enhancerStrenth <- 'high'

sim.flags$selectionScreen <- 'yes'
sim.flags$selectionStrength <- 'high'
simulate_data(sim.flags)
```

## Advanced simulations
