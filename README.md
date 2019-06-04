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

## Simulation quickstart with example data (Selection screen)
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

If guide targets are provided then the gene of interest should also be provided. The assumption is that guides are also targeting the gene of interest and can serve as positive controls.
```
sim.flags$exon <- example.gene
```

Provide the original guide distribution for each replicate. This can be supplied from an existing data set, as is done below. Another option is to generate the distributions using a zero inflated negative binomial distribution. See 'Advanced Simulations' for details.
```
sim.flags$inputGuideDistr <- cbind(before_1 = example.counts$before_repl1, 
  before_2 = example.counts$before_repl2)
```  

Specify the screen type. Either a selection screen like here or a FACSscreen, where cells are sorted into different pools. See 'FACS screen simulation' section for an example simulation of a FACS screen.
```
sim.flags$selectionScreen <- 'yes'
```

Add names for the different pools in each replicate. In this case each replicate will have a before selection and after selection pool. The names will be before_repl1, after_repl1, before_repl2, ...
```
sim.flags$poolNames <- c('before', 'after')
```

Sepcify the crispry system used. For CRISPRi the default range of effect is assumed to be 1kb. However this can be changed. See the 'Advanced Simulations' section for more details and how to simulate other CRISPR systems.
```
sim.flags$crisprSystem <- 'CRISPRi'
```

Specify the number and the size of enhancers to be simulated. They will be placed at random positions within the screen.
```
sim.flags$nrEnhancers <- 5
sim.flags$enhancerSize <- 50  # base pairs
```

Specify the sequencing depth for each of the pools. In addition, the simulations can simulate data sets where PCR duplicates are either accounted for or not. If they are accounted for set the 'pcrDupl' flag to 'no'.
```
sim.flags$seqDepth <- list(repl1 = c(16656607, 19431422),
    repl2 = c(20217155, 21585515))
sim.flags$pcrDupl <- 'yes'
```

Specify:

    - the selection strength, how strong the effect of disrupting the gene of interest is (high, low). 
    
    - the guide efficiency, what proportion of guides have an effect (high, medium, low)
    
    - the enhancer strength, how strong the signal from the enhancers is (high, medium, low)
    
    
The parameters for all of these can also be manually set. See the 'Advanced Simulations' section for details.
```
sim.flags$selectionStrength <- 'high'
sim.flags$guideEfficiency <- 'high'
sim.flags$enhancerStrenth <- 'high'
```

Run the simulations
```
simulate_data(sim.flags)
```

## FACS screen simulation quickstart with example data
Outlined below are the main differences to the flags set above

Instead of a selection screen, the FACS screen flag is set.
```
sim.flags$FACSscreen <- 'yes'
```

As above, each pool per replicate has to be names. However, in this case there will be more than one pool. Here an example where wells are sorted form an input pool into a high, medium and low expression pool.
```
sim.flags$poolNames <- c('input', 'high', 'medium', 'low')
```

The sequencing depth has to be specified for each of the four pools for each of the replicates.
```
sim.flags$seqDepth <- list(repl1 = rep(18e6, 4), repl2 = rep(18e6, 4) )
```


## Advanced Simulations

```
sim.flags$guides <- example.info
```

and CRISPRa it is assumed that the area of effect is 1kb. For Cas9 it is assumed to be 20bp and for a dualCRISPR system the sange is specified already by the start and , however this can be changed
