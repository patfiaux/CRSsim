# CRSsim
This github account contains code and instructions for simulating CRISPR regulatory screens and assessing performance of different methods used in CRISPR regulatory screens. It contains the following sections:
[Simulations](1.-Simulating-CRISPR-regulatory-screen-data)

#1. Simulating CRISPR regulatory screen data
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
### 1.1.1 R packages
```
MCMCpack # install.packages('MCMCpack')

transport # install.packages('transport')

```

### 1.1.2 Bioconductor packages
```
if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")
    
IRanges # BiocManager::install("IRanges", version = "3.8")

GenomicRanges # BiocManager::install("GenomicRanges", version = "3.8")
```

## 1.2 Simulation quickstart with example data (Selection screen)
### 1.2.1. source the script
```
source('/path/to/script/RELICS_sim.r')
```

### 1.2.2. Setting up simulation flags. 
There are several different flags which have no defaults and need to be supplied by the user. Below is the outline on how to set the most important flags to get the simulation going.

Flags are set up in a list format

```
sim.flags <- list()
```

Set the output name of the simulation
```
sim.flags$simName <- 'Example_simulation'
```

Provide info about the guide targets. Either supply them directly, as is done here, or generate them (see details 'Advanced Simulations'). The input is a data frame with columns for chromosome, start position and end position: 'chrom', 'start', 'end'. Each row is a guide and details the target information for ech guide. For Cas9, CRISPRi and CRISPRa screens, the difference in start and end should be set to something small, such as: start = target site - 20, end = target site.
```
sim.flags$guides <- example.info
```

If guide targets are provided then the gene of interest should also be provided. The assumption is that guides are also targeting the gene of interest and can serve as positive controls. The input is a data frame with columns called: 'chrom', 'start', 'end'
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

Specify the sequencing depth for each of the pools. In addition, the simulations can simulate data sets where PCR duplicates are either accounted for or not. If they are accounted for set the 'pcrDupl' flag to 'no'. The sequencing depth is given in list format, where each list entry is a replicate with the corresponding sequencing depths.
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

## 1.2.3 FACS screen simulation quickstart with example data
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

## Keep in mind!
An average guide count of 15 vs 100 vs 500 has a major effect on detecting true signal when everything else is held constant. Make sure you adjust (seqDepth) when changing the number of guides used for simulating the data.


# 2. Analyzing simulated data and evaluate performance


## 3. Advanced Simulations

Guides and their targets can be simulated if not readily available. Both single-guide as well as dual guide screens can be simulated. For both types the number of guides (nrGuides) have to be specified, as well as the screen type (screenType) and the step size between guides (stepSize). In addition to that, if a dual CRISPR screen is chosen, the deletion size has to be specified (stepSize).

screenType: CRISPRi, CRISPRa, Cas9, dualCRISPR

If this option is chosen, all guides are abritrarily chosen to be located on chromosome 1 and ~5% of the guides will be selected to serve as positive control.

```
sim.flags$guides <- generate_guide_info(list(nrGuides = 10000, screenType = 'dualCRISPR', stepSize = 20, deletionSize = 1000))
```

The input count distribution for the different replicates can either be taken from an existing data set. However, it is also possible both generalize existing distributions using the zero-inflated negative binomial distribution (ZINB). The ZINB has both a mean (rate) and a dispersion parameter as well as a parameter indicating the fraction of the disctribution originating from the zero mass (eta). Below are the steps to both obtain as well as use the parameters from a ZINB
```
# obtain ZINB parameters which descibe the distribution
before.repl1.par <- obtain_ZINB_pars(example.counts$before_repl1)
example.rate <- before.repl1.par$rate              # rate = 76.5
example.dispersion <- before.repl1.par$dispersion  # dispersion = 2.6
example.eta <- before.repl1.par$eta                # eta = 1e-4

# to generate a ZINB distribution with 15000 guides
before.repl1.simulated <- create_ZINB_shape(15000, example.eta, example.rate, example.dispersion)
before.repl2.simulated <- create_ZINB_shape(15000, example.eta, example.rate, example.dispersion)

# combine the two simulated replicates and set them as input distributions
sim.flags$inputGuideDistr <- cbind(before_1 = before.repl1.simulated, before_2 = before.repl2.simulated)
```  

Currently, 4 different CRISPR systems can be simulated: CRISPRi, CRISPRa, Cas9, dualCRISPR

By default, CRISPRi and CRISPRa are assumed to have an effect range of 1kb and Cas9 of 20bp. However, it is also possible to manually set this size using the 'crisprEffectRange'
For dualCRISPR the effect range is the deletion itself. The deletion size introduced by two guides has to be represented by 'start' being the target site of guide 1 and 'end' the target site of guide 2.
```
# example for how to change the effect range of a CRISPR system used
sim.flags$crisprSystem <- 'CRISPRi'
sim.flags$crisprEffectRange <- 500
```

Both the guide efficiency and the enhancer strength are simulated from a beta distribution. The two parameters can be specified by setting 'guideEfficiency' and 'enhancerStrenth' to either 'high', 'medium' or 'low'. However, it is also possible to directly specify the two shape parameters of the beta distribution. As a general rule, if the shape parameters provided are large, the variance seen in the distribution is reduced. The larger shape 1 parameter is comapred to shape 2 parameter, the more the distribution will be skewed towards 1. To get an intuition you can also plot the historgram from randomly drawing from a beta and then subsequently change the shape parameters. The default parameters are:

    - high: enhancerShape1 = 7, enhancerShape2 = 2
    
    - medium: enhancerShape1 = 5, enhancerShape2 = 5
    
    - low: enhancerShape1 = 2, enhancerShape2 = 7
    
(same for guideShape1, guideShape2)

```
hist(rbeta(10000, shape1 = 8, shape2 = 1))  # randomly generate 10000 instances of the beta distribution

sim.flags$guideEfficiency <- 'high'
sim.flags$enhancerStrenth <- 'high'

# the above is equivalent to what's below
sim.flags$enhancerShape1 <- 7
sim.flags$enhancerShape2 <- 2
sim.flags$guideShape1 <- 7
sim.flags$guideShape2 <- 2
```


Selection strength: To understand the details of the selection strength it is helpful to have some understanding of the Dirichlet distribution. Both for the selection screen as well as for the FACS screen the prbability of each guide either being selected, or sorted into a given pool is given by a random Dirichlet deviate. As an example; 

Cells are sorted into three pools (high, medium and low gene expression) and the probability of a negative control to be sorted into any of the given pools at random is captured by the probabilities: 0.48, 0.48, 0.04 . However, for each guide this probability shifts slightly. The Dirichlet deviate provides this variation. All cells containing the negative control above are subsequently assigned to the pools with following probabilities:  `rdirichlet(1, c(48, 48, 4))`

Continuing the example above: assume the probability of a positive control to be sorted into any given pool to be represented by the following: 0.45, 0.45, 0.1. All cells containing this positive control above are subsequently assigned to the pools with following probabilities:  `rdirichlet(1, c(45, 45, 10))`

If the above does not seem like a big difference, view it as benig a 2.5 fold change in being sorted into the low pool. Thatchange is not insignificant!

To manually set the diriclet probabilities use the 'posSortingFrequency' and the 'negSortingFrequency' flags. To continue the example from above:
```
sim.flags$posSortingFrequency <- c(45, 45, 10)
sim.flags$negSortingFrequency <- c(48, 48, 4)
```

Note: 

1. The sum of the frequencies does not have to be 1. 

2. The larger the numbers chosen the less variable the sorting becomes.


The defaults for the 'high' flags were chosen due to their capability of accurately representing the sorting parameters for either a selection screen or a FACS screen. The default 'low' flags were chosen as an arbitrary fraction of the 'high' selection.

The default flags used for 'high':
```
# for a selection screen:
sim.flags$posSortingFrequency <- c(1)
sim.flags$negSortingFrequency <- c(5)

# in a FACS screen sorted into 3 pools: 
# '97' is repeated for all pools except the last one
sim.flags$posSortingFrequency <- c(97, 97, 13) * 0.5
sim.flags$negSortingFrequency <- c(97, 97, 3) * 0.5
```

The default flags used for 'low':
```
# for a selection screen:
sim.flags$posSortingFrequency <- c(4)
sim.flags$negSortingFrequency <- c(5)

# in a FACS screen sorted into 3 pools: 
# '97' is repeated for all pools except the last one
sim.flags$posSortingFrequency <- c(97, 97, 5) * 0.5
sim.flags$negSortingFrequency <- c(97, 97, 3) * 0.5
```

