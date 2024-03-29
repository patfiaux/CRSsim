suppressPackageStartupMessages(require(MCMCpack))
suppressPackageStartupMessages(require(transport))
suppressPackageStartupMessages(require(IRanges))
suppressPackageStartupMessages(require(GenomicRanges))
suppressPackageStartupMessages(require(splines))

total_wasserstein <- function(input.df1, input.df2){
  out.wasserstein <- 0
  per.col.wst <- c()

  for(i in 1:ncol(input.df1)){
    temp.wst <- wasserstein1d(input.df1[,i], input.df2[,i])
    out.wasserstein <- out.wasserstein + temp.wst
    per.col.wst <- c(per.col.wst, temp.wst)
  }

  return(list(total_wst = out.wasserstein, per_col_wst = per.col.wst))
}

counts_from_probs <- function(input.probs, input.distr){
  data.probs <- input.probs  # probabilities
  out.df <- data.frame(matrix(nrow = length(input.distr), ncol = ncol(input.probs)))
  for(i in 1:length(input.distr)){
    # for each guide, compute the counts that were sorted into each pool
    out.df[i, ] <- t(rmultinom(1, size = input.distr[i], prob = data.probs[i,] ))
  }
  return(out.df)
}

# input.counts: veotr, number of guides in pool per guide
# pcr.duplicate: boolean, if true sequence with, else without duplicates
# input.seq.depth: how deep to sequence
# guide.nr: total number of guides (could have 0 counts for some, need to maintain guide order)
sequence_pool <- function(input.counts, pcr.duplicate, input.seq.depth, guide.nr){

  out.seq <- rep(0, guide.nr)  # initialize empty vector
  temp.obs <- rep(c(1:length(input.counts)), input.counts)  # guide id the number of times a guide is present
    # 1,1,1,2,3,3,3,3,4,4,5,5,5,...

  if(pcr.duplicate){ # == 'duplicate' | pcr.type == 'yes'){
    temp.sequenced <- table(sample(temp.obs, input.seq.depth, replace=TRUE))
  } else {
    if(length(temp.obs) < input.seq.depth){
      temp.sequenced <- table(temp.obs)
    } else {
      temp.sequenced <- table(sample(temp.obs, input.seq.depth, replace=FALSE))
    }
  }

  # table output returns id as name, and number of observations as counts
    #   1   21   23   51   58   67
    # 2014   18  371   19    5    4

  seq.pos <- as.numeric(names(temp.sequenced))
  out.seq[seq.pos] <- as.numeric(temp.sequenced)
  return(out.seq)
}

# input.pools, data frame with counts per guide-id
# pcr.duplicate: boolean, if true sequence with, else without duplicates
# input.seq.depth: how deep to sequence
experiment_sequencing <- function(input.pools, pcr.duplicate, input.seq.depth){
  out.df <- c()

  for(i in 1:ncol(input.pools)){
    temp.seq <- sequence_pool(input.pools[,i], pcr.duplicate, input.seq.depth[i], nrow(input.pools))
    out.df <- cbind(out.df, temp.seq)
  }
  return(out.df)
}

set_default_flags <- function(input.list){

  out.list <- input.list
  input.list.names <- names(input.list)

  if('guideFile' %in% names(input.list)){
    
    out.list$guides <- read.csv(input.list$guideFile, stringsAsFactors = F)
    out.list$exon <- read.csv(input.list$exon, stringsAsFactors = F)
  } else {
    if(! 'nrGuides' %in% names(input.list) | ! 'crisprSystem' %in% names(input.list) | ! 'stepSize' %in% names(input.list)){
      print('Error, missing flags for generating guides. Please specify: nrGuides, crisprSystem, stepSize')
    } else if(input.list$crisprSystem == 'dualCRISPR'){
      if(! 'deletionSize' %in% names(input.list)){
        print('Error. deletionSize must be specified for dualCRISPR screens')
      } else {
        out.list$guides <- generate_guide_info(list(nrGuides = input.list$nrGuides,
          crisprSystem = input.list$crisprSystem, stepSize = input.list$stepSize,
          deletionSize = input.list$deletionSize))
      }
    } else {
      out.list$guides <- generate_guide_info(list(nrGuides = input.list$nrGuides,
        crisprSystem = input.list$crisprSystem, stepSize = input.list$stepSize))
    }

    max.dist <- max(out.list$guides$end) - min(out.list$guides$start)
    gene.start <-  min(out.list$guides$start) + round(max.dist/ 2)
    gene.end <- round(gene.start + 0.05 * max.dist)
    gene.exons <- data.frame(chrom = out.list$guides$chrom[1], start = gene.start,
      end = gene.end, stringsAsFactors = F)
    out.list$exon <- gene.exons

  }

  if(! 'outDir' %in% input.list.names){
    out.list$outDir <- ''
  }
  if(! 'randomizeRepl' %in% input.list.names){
    # in case the correlation between replicates is to be broken
    out.list$randomizeRepl <- 'no'
  }
  
  # is this to ensure that the length of the input guide distribution is equal to the length of the guide targets data?
  if(nrow(out.list$inputGuideDistr) != nrow(out.list$guides)){
    if(nrow(out.list$inputGuideDistr) >= nrow(out.list$guides)){
      out.list$inputGuideDistr <- out.list$inputGuideDistr[sample(1:nrow(out.list$inputGuideDistr), size = nrow(out.list$guides)),]
    } else {
      out.list$inputGuideDistr <- out.list$inputGuideDistr[sample(1:nrow(out.list$inputGuideDistr), size = nrow(out.list$guides), replace = T),]
    }
  }
  
  if(! 'nrNeg' %in% input.list.names){
    out.list$nrNeg <- 300
  }
  
  # Number of simulations to generate with specified parameters
  if(! 'nrSims' %in% input.list.names){
    out.list$nrSims <- 1
  }
  
  # Total number of cells going though FACS
  if(! 'inputCellNr' %in% input.list.names){
    out.list$inputCellNr <- rep(50e6, ncol(out.list$inputGuideDistr))
  }
  
  # Parameterization of Beta distribution for enhancers
  if('enhancerStrength' %in% input.list.names){
    if(out.list$enhancerStrength == 'high'){
      out.list$enhancerShape1 <- 7
      out.list$enhancerShape2 <- 2
    } else if(out.list$enhancerStrength == 'medium'){
      out.list$enhancerShape1 <- 5
      out.list$enhancerShape2 <- 5
    } else if(out.list$enhancerStrength == 'low'){
      out.list$enhancerShape1 <- 2
      out.list$enhancerShape2 <- 7
    } else if(out.list$enhancerStrength == 'random'){
      out.list$enhancerShape1 <- 1
      out.list$enhancerShape2 <- 1
    } else{
      print('Error: please specify valid enhancerStrength (high, medium, low)')
      break()
    }
  }
  
  # Parameterization of Beta distribution for guide efficiency
  if('guideEfficiency' %in% input.list.names){
    if(out.list$guideEfficiency == 'high'){
      out.list$guideShape1 <- 7
      out.list$guideShape2 <- 2
    } else if(out.list$guideEfficiency == 'medium'){
      out.list$guideShape1 <- 5
      out.list$guideShape2 <- 5
    } else if(out.list$guideEfficiency == 'low'){
      out.list$guideShape1 <- 2
      out.list$guideShape2 <- 7
    } else if(out.list$guideEfficiency == 'random'){
      out.list$guideShape1 <- 1
      out.list$guideShape2 <- 1
    } else if(out.list$guideEfficiency == 'high-low'){
      out.list$guideShape1 <- 0.1
      out.list$guideShape2 <- 0.1
    } else {
      print('Error: please specify valid guideEfficiency (high, medium, low)')
      break()
    }
  }

  if(! ('screenType' %in% names(out.list) )){ #| 'FACSscreen' %in% names(out.list)) ){
    print('No info about screen type. Set screenType to either selectionScreen or FACSscreen')
    break()
  }
  
  if(out.list$screenType == 'selectionScreen'){
    # add the second sequencing depth parameter
    
    # # set seqDepth to default 100
    temp.seq.depth <- c()
    
    # modify seqDepth to 100 if not given in flags
    if (! 'seqDepth' %in% input.list.names) {
      temp.seq.depth <- list(repl1 = rep(100 * length(out.list$inputGuideDistr), 2),
                             repl2 = rep(100 * length(out.list$inputGuideDistr),2) )
    } else {
      temp.seq.depth<- out.list$seqDepth
    }

    # does this adjust sequencing depths based on the input guide distribution?
    adj.seq.depth <- lapply(temp.seq.depth, function(x){
      return(c(x, x[length(x)]))
    })

    # Set default dirichlet parameters
    out.list$seqDepth <- adj.seq.depth
    if('selectionStrength' %in% input.list.names){
      if(out.list$selectionStrength == 'very-high'){
        out.list$posSortingFrequency <- c(1, 95)
        out.list$negSortingFrequency <- c(10, 95)
      } else if(out.list$selectionStrength == 'high'){
        out.list$posSortingFrequency <- c(1, 95)
        out.list$negSortingFrequency <- c(5, 95)
      } else if(out.list$selectionStrength == 'low'){
        out.list$posSortingFrequency <- c(4, 95)
        out.list$negSortingFrequency <- c(5, 95)
      } else {
        print('Specified selectionStrength not yet implemented. Pick either high or low or specify parameters manually')
        break()
      }
    } else if(! 'posSortingFrequency' %in% names(out.list)){
      print('strength of selectionScreen has not been specified. Either set selectionStrength to high, low or specify the parameters manually')
      break()
    }
  } else if(out.list$screenType == 'FACSscreen'){
      if('selectionStrength' %in% input.list.names){
        if(out.list$selectionStrength == 'very-high'){
          out.list$posSortingFrequency <- c(rep(97, (length(out.list$seqDepth[[1]]) - 2) ), 13*2)*0.5
          out.list$negSortingFrequency <- c(rep(97, (length(out.list$seqDepth[[1]]) - 2) ), 3)*0.5
        } else if(out.list$selectionStrength == 'high'){
          out.list$posSortingFrequency <- c(rep(97, (length(out.list$seqDepth[[1]]) - 2) ), 13)*0.5
          out.list$negSortingFrequency <- c(rep(97, (length(out.list$seqDepth[[1]]) - 2) ), 3)*0.5
        } else if(out.list$selectionStrength == 'low'){
          out.list$posSortingFrequency <- c(rep(97, (length(out.list$seqDepth[[1]]) - 2) ), 5)*0.5
          out.list$negSortingFrequency <- c(rep(97, (length(out.list$seqDepth[[1]]) - 2) ), 3)*0.5
        } else {
          print('Specified selectionStrength not yet implemented. Pick either high or low or specify parameters manually')
        }
      } else if(! 'posSortingFrequency' %in% names(out.list)){
        print('strength of FACSscreen has not been specified. Either set selectionStrength to high, low or specify the parameters manually')
      }
  } else {
    print('Error: no screenType specified. Select either selectionScreen or FACSscreen')
    break()
  }

  if(! 'enhancerShape1' %in% names(out.list)){
    print('Error: please specify valid enhancerStrength (high, medium, low) or set the parameters for enhancerShape1 manually')
    break()
  }
  if(! 'enhancerShape2' %in% names(out.list)){
    print('Error: please specify valid enhancerStrength (high, medium, low) or set the parameters for enhancerShape2 manually')
    break()
  }
  if(! 'guideShape1' %in% names(out.list)){
    print('Error: please specify valid guideEfficiency (high, medium, low) or set the parameters for guideShape1 manually')
    break()
  }
  if(! 'guideShape2' %in% names(out.list)){
    print('Error: please specify valid guideEfficiency (high, medium, low) or set the parameters for guideShape2 manually')
    break()
  }
  
  # type of area of effect. options are 'normal' (default), 'uniform', 'logistic' (need to implement)
  if(! 'areaOfEffect_type' %in% input.list.names){
    out.list$areaOfEffect_type <- 'normal'
  }
  
  if(out.list$areaOfEffect_type == 'uniform'){
    if(! 'crisprEffectRange' %in% input.list.names){
      if(out.list$crisprSystem == 'CRISPRi' | out.list$crisprSystem == 'CRISPRa'){
        out.list$crisprEffectRange <- 500
      } else if(out.list$crisprSystem == 'Cas9'){
        out.list$crisprEffectRange <- 10
      } else if(out.list$crisprSystem == 'dualCRISPR'){
        out.list$crisprEffectRange <- 10
        if(! 'deletionProb' %in% input.list.names){
          out.list$deletionProb <- 0.2
        }
      } else {
        print("Error: please specify a valid CRISPR system for the 'uniform' area of effect (Cas9, CRISPRi, CRISPRa, dualCRISPR)")
        break()
      }
    }
  }
  
  if(out.list$areaOfEffect_type == 'normal'){
    if(! 'normal_areaOfEffect_sd' %in% input.list.names){
      if(out.list$crisprSystem %in% c('CRISPRi', 'CRISPRa', 'dualCRISPRi', 'dualCRISPRa') ){
        
        out.list$normal_areaOfEffect_sd <- 170
        if(! 'crisprEffectRange' %in% input.list.names){
          out.list$crisprEffectRange <- 415
        }
      } else if(out.list$crisprSystem %in% c('Cas9', 'dualCRISPR') ){
        out.list$normal_areaOfEffect_sd <- 8.5
        if(! 'crisprEffectRange' %in% input.list.names){
          out.list$crisprEffectRange <- 21
        }
        if(out.list$crisprSystem == 'dualCRISPR'){
          if(! 'deletionProb' %in% input.list.names){
            out.list$deletionProb <- 0.2
          }
        }
      } else {
        print("Error: please specify a valid CRISPR system for the 'uniform' area of effect (Cas9, CRISPRi, CRISPRa, dualCRISPR)")
        break()
      }
    }
  }
  
  if(! 'dispersionType' %in% input.list.names){
    out.list$dispersionType <- 'exponential' # will probably have to be changed
  }
  if(out.list$dispersionType == 'spline'){
    if(! 'splineModelLoc' %in% input.list.names){
      print("Error: please specify a location for the spline model.")
      break()
      # out.list$splineModelLoc <- ''
    }
  }

  return(out.list)
}


#' @title generate sgRNA target coordinates
#' @param input.list: list containing all the parameters: $nrGuides, $stepSize, if dualCRISPR: $deletionSize
#' @return data.frame: chrom, start, end
#' @export generate_guide_info()

generate_guide_info <- function(input.list){

  guide.chrom <- rep('chr1', input.list$nrGuides)
  guide.start <- c()
  guide.end <- c()
  if(input.list$crisprSystem == 'dualCRISPR'){
    guide.start <- 2000 + c(1:input.list$nrGuides) * input.list$stepSize
    guide.end <- guide.start + input.list$deletionSize
  } else {
    guide.start <- 2000 + c(1:input.list$nrGuides) * input.list$stepSize - 10
    guide.end <- guide.start + 20
  }

  out.df <- data.frame(chrom = guide.chrom, start = guide.start, end = guide.end, stringsAsFactors = F)

  return(out.df)

}


#' @title main function for simulating tiling CRISPR screen data
#' @param input.list: list containing all the necessary parameters: $nrGuides, $stepSize, if dualCRISPR: $deletionSize
#' @param input.parameter.list: RELICS analysis parameters in list format. Default = NULL
#' @param record.all.fs: logical, if information of all intermediate FS should be recorded, in addition to the final setof FS
#' @return $nrSims data.frames, each simulation returns data with info and counts as well as enhancers
#' @export simulate_data()

simulate_data <- function(input.list){
  
  input.list <- set_default_flags(input.list)
  
  updated.info <- format_guide_info(input.list$guides, input.list)
  
  dir.create(file.path(paste0(input.list$outDir, input.list$simName)))
  
  for(sim in 1:input.list$nrSims){
    print(paste0('Starting simulation ', sim, ' out of ', input.list$nrSims))
    input.list$inputPools <- list()
    input.list$enhancer <- generate_enhancers(input.list, updated.info, input.list$exon)
    
    print('Generating initial pool')
    for(repl in 1:length(input.list$inputCellNr)){
      if(input.list$randomizeRepl == 'yes'){
        randomized.counts <- sample(input.list$inputGuideDistr[, repl])
        input.list$inputPools[[paste0('repl',repl)]] <- rmultinom(1, input.list$inputCellNr[repl], 
                                                                  randomized.counts)
      } else {
        
        # generating multinomial distribution based on the input guide distribution and number of input cells for given replicate
        input.list$inputPools[[paste0('repl',repl)]] <- rmultinom(1, input.list$inputCellNr[repl], 
                                                                  input.list$inputGuideDistr[, repl])
      }
    }
    
    temp.sim <- c()
      
    if(input.list$crisprSystem %in% c('CRISPRi', 'CRISPRa', 'Cas9')){
      temp.sim <- single_guide_replicate_simulation(input.list, updated.info, sim)
    } else if(input.list$crisprSystem %in% c('dualCRISPR')){
      temp.sim <- paired_guide_replicate_simulation(input.list, updated.info, sim)
    }
    
    pool.prop <- as.data.frame(t(temp.sim$pool_proportions) )
    colnames(pool.prop) <- colnames(temp.sim$counts)
    
    write.csv(pool.prop, file = paste0(input.list$outDir, input.list$simName,'/',
                                             input.list$simName, '_sim', sim, '_poolGuideProportion.csv'), row.names = F)
    
    write.csv(temp.sim$counts, file = paste0(input.list$outDir, input.list$simName,'/',
                                             input.list$simName, '_sim', sim, '_counts.csv'), row.names = F)
    
    input.list$enhancer$enhancerStrength <- temp.sim$enhancerStrength
    ordered.enhancers <- input.list$enhancer[order(input.list$enhancer$start),]
    
    write.csv(ordered.enhancers, file = paste0(input.list$outDir, input.list$simName,'/',
                                               input.list$simName, '_sim', sim, '_enhancers.csv'), row.names = F)
    
    write.table(ordered.enhancers, file = paste0(input.list$outDir, input.list$simName,'/',
                                               input.list$simName, '_sim', sim, '_enhancers.bed'), 
                row.names = F, col.names = F, sep = '\t', quote = F)
    
    
    init.labels <- rep('chr', nrow(updated.info))
    final.info <- updated.info
    final.info$label <- init.labels
    exon.label.pos <- sim_enhancer_overlaps(updated.info, input.list$exon) # also works for exons
    final.info$label[exon.label.pos] <- 'exon'
    enh.label.pos <- sim_enhancer_overlaps(updated.info, ordered.enhancers)
    final.info$label[enh.label.pos] <- 'pos'
    random.neg <- sample(which(final.info$label == 'chr'), input.list$nrNeg)
    for(i in 1:length(random.neg)){
      final.info$chrom[random.neg[i]] <- NA
      final.info$start[random.neg[i]] <- NA
      final.info$end[random.neg[i]] <- NA
      final.info$label[random.neg[i]] <- 'neg'
    }
    
    if(input.list$crisprSystem %in% c('CRISPRi', 'CRISPRa', 'Cas9')){
      final.info$guideEfficiency <- temp.sim$guide_efficiencies
    } else if(input.list$crisprSystem %in% c('dualCRISPR')){
      final.info$guide1_Efficiency <- temp.sim$g1_efficiency
      final.info$guide2_Efficiency <- temp.sim$g2_efficiency
    }
    
    write.csv(final.info, file = paste0(input.list$outDir, input.list$simName,'/',
                                        input.list$simName, '_sim', sim, '_info.csv'), row.names = F)
  }
}


#' @title adjust the info file (accounting for single or paired guide design)
#' @param input.list: list containing all the necessary parameters: $crisprEffectRange, $crisprSystem, 
#' @param input.info: info file
#' @param sim.nr: which simulation number
#' @return list: counts, enhancerStrength, guide_efficiencies
#' @export format_guide_info()

format_guide_info <- function(input.info, input.list){
  
  updated.info <- input.info
  
 # account for single or paired guide design
  # paired guides are already spaced apart (ex. 1kb)
  if(input.list$crisprSystem %in% c('dualCRISPR')){
    
    if(input.list$areaOfEffect_type == 'normal'){
      if(! 'sgTarget_1' %in% names(updated.info)){
        updated.info$sgTarget_1 <- updated.info$start
      }
      if(! 'sgTarget_2' %in% names(updated.info)){
        updated.info$sgTarget_2 <- updated.info$end
      }
      updated.info$sgStart_1 <- updated.info$sgTarget_1 - input.list$crisprEffectRange
      updated.info$sgEnd_1 <- updated.info$sgTarget_1 + input.list$crisprEffectRange
      updated.info$sgStart_2 <- updated.info$sgTarget_2 - input.list$crisprEffectRange
      updated.info$sgEnd_2 <- updated.info$sgTarget_2 + input.list$crisprEffectRange
      updated.info$start <- updated.info$sgStart_1
      updated.info$end <- updated.info$sgEnd_2
      
    } else {
      updated.info$start <- updated.info$start - input.list$crisprEffectRange
      updated.info$end <- updated.info$end + input.list$crisprEffectRange
    }
    
  } else {
    
    if(! 'sgTarget' %in% names(updated.info)){
      updated.info$sgTarget <- round((updated.info$start + updated.info$end)/2)
    }
    
    updated.info$start <- updated.info$sgTarget - input.list$crisprEffectRange
    updated.info$end <- updated.info$sgTarget + input.list$crisprEffectRange

  }
  
  return(updated.info)
}


#' @title give guides different sorting probabilities based on overlaps with functional sequences
#' @param input.frame: list containing all the necessary parameters: $posSortingFrequency, $negSortingFrequency, if dualCRISPR: $deletionSize
#' @param input.info: info file
#' @param sim.nr: which simulation number
#' @return list: counts, enhancerStrength, guide_efficiencies
#' @export paired_guide_replicate_simulation()

# give different genomic properties different sorting probabilities
# this method samples from either the null sorting distrubuton or samples
# from an enhancer-strength distribution
# x% of guides have an enhancer effect, remaining guides have no effect
# this implementation also records the guide efficiency
paired_guide_replicate_simulation <- function(input.frame, input.info, sim.nr){
  
  out.sim.data <- c()
  out.repl.sort.fractions.filt <- c()

  effect.diff <- c()
  
  # Divide Negative Sorting Frequency to Obtain Proportions
  input.frame$negProp <- input.frame$negSortingFrequency / sum(input.frame$negSortingFrequency)
  input.frame$posProp <- input.frame$posSortingFrequency / sum(input.frame$posSortingFrequency)
  
  # does this lead to vectors like (4, 0)?
  if(input.frame$screenType == 'selectionScreen'){
    effect.diff <- input.frame$negProp - input.frame$posProp
  } else {
    effect.diff <- input.frame$posProp - input.frame$negProp
  }
  
  sort.factor <- rbeta(nrow(input.frame$enhancer), shape1 = input.frame$enhancerShape1, shape2 = input.frame$enhancerShape2)
  guide.1.efficiencies <- rbeta(nrow(input.info), shape1 = input.frame$guideShape1, shape2 = input.frame$guideShape2)
  guide.2.efficiencies <- rbeta(nrow(input.info), shape1 = input.frame$guideShape1, shape2 = input.frame$guideShape2)
  
  # for the number of replicates
  for(i in 1:length(input.frame$inputPools)){

    input.distr <- input.frame$inputPools[[i]]
    
    # Obtain Dispersion for Each Guide and Save in Vector
    guide.disp <- c()
    if(input.frame$dispersionType == 'independent'){
      guide.disp <- rep(sum(input.frame$negSortingFrequency), length(input.distr))
    } else if(input.frame$dispersionType == 'exponential'){

      # exponential fit without linear term (From MYC Repl 2)
      guide.disp <- -89.17911 + 17.78130 * log(input.distr)
      
      # set floor value of 3 for dispersion to avoid negative probabilities
      guide.disp[guide.disp < 3] <- 3
      
    } else if(input.frame$dispersionType == 'spline'){
      disp.model <- readRDS(input.frame$splineModelLoc)
      total.counts.df <- data.frame(counts = input.distr, stringsAsFactors = F)
      spline.predict <- suppressWarnings(predict(disp.model, total.counts.df))
      guide.disp <- spline.predict
      guide.disp[guide.disp < 1] <- 1
    }
    
    # compute the per-guide dirichlet parameters by including the dispersion
    temp.adj.bkg.freq <- t(t(guide.disp)) %*% input.frame$negProp 
    temp.adj.fs.freq <- t(t(guide.disp)) %*% input.frame$posProp
    
    
    print(paste0('Generating replicate ', i))
    # add counts, assuming negative sorting probabilities
    temp.sort.prob <- matrix(0, nrow = length(input.distr), ncol = length(input.frame$negSortingFrequency))
    for (j in 1:length(input.distr)) {
      temp.sort.prob[j, ] <- rdirichlet(1, temp.adj.bkg.freq[j, ])
    }
    
    if(length(which(is.na(temp.sort.prob))) > 0){
      temp.nan.rows <- which(is.na(temp.sort.prob[,1]))
      temp.fill.values <- rep(1 / ncol(temp.sort.prob), ncol(temp.sort.prob))
      temp.sort.prob[which(is.na(temp.sort.prob[,1])),] <- temp.fill.values
    }
    
    repl.counts <- counts_from_probs(temp.sort.prob, input.distr)
    
    # combine both uniform and normal into one function
    repl.counts <- compute_pairedGuide_AoE(input.frame, input.info, effect.diff,
                                           repl.counts, sort.factor, guide.1.efficiencies, 
                                           guide.2.efficiencies, input.distr,
                                           temp.adj.bkg.freq, temp.adj.fs.freq, guide.disp)

    repl.sequenced <- experiment_sequencing(cbind(input.distr, repl.counts),
                                            input.frame$pcrDupl, input.frame$seqDepth[[i]])
    repl.sort.fractions <- colSums(cbind(input.frame$inputPools[[i]], repl.counts)) / sum(colSums(cbind(input.frame$inputPools[[i]], repl.counts)))
    
    repl.sequenced.filtered <- repl.sequenced
    repl.sort.fractions.filt <- repl.sort.fractions
    if(input.frame$screenType == 'selectionScreen'){ 
      repl.sequenced.filtered <- repl.sequenced[,c(1,2)]
      repl.sort.fractions.filt <- repl.sort.fractions[c(1,2)]
    }
    
    colnames(repl.sequenced.filtered) <- paste(paste0('sim', sim.nr, '_repl', i), input.frame$poolNames, sep = '_')
    out.sim.data <- cbind(out.sim.data, repl.sequenced.filtered)
    
    out.repl.sort.fractions.filt <- c(out.repl.sort.fractions.filt, repl.sort.fractions.filt)
    
  }
  return(list(counts = out.sim.data, 
              enhancerStrength = sort.factor, 
              g1_efficiency = guide.1.efficiencies, 
              g2_efficiency = guide.2.efficiencies,
              pool_proportions = out.repl.sort.fractions.filt))

}


#' @title Compute the area of effect using a normal distribution for paired guide simulations
#' @param input.frame: list containing all the necessary parameters: $posSortingFrequency, $negSortingFrequency, if dualCRISPR: $deletionSize
#' @param input.info: info data frame, $chrom, $start, $end
#' @param effect.diff: difference in dirichlet probabilities between null and alternative
#' @param repl.counts: data frame wih the null counts of the replicate
#' @param sort.factor: beta distributed fractions of signal strength relative to exons
#' @param guide.1.efficiencies: guide 1 efficiency
#' @param guide.2.efficiencies: guide 2 efficiency
#' @param input.pool.counts: counts per cell in the input pool
#' @param adj.bkg.freq: dirichlet background parameters (including dispersion)
#' @param adj.fs.freq: dirichlet FS parameters (including dispersion)
#' @param guide.disp: vector, contains the dispersion for each guide
#' @return data frame: repl.counts
#' @export compute_pairedGuide_AoE()

compute_pairedGuide_AoE <- function(input.frame, input.info, effect.diff, 
                                    repl.counts, sort.factor, guide.1.efficiencies, 
                                    guide.2.efficiencies, input.pool.counts,
                                    adj.bkg.freq, adj.fs.freq, guide.disp){
  
  out.repl.counts <- repl.counts
  
  # compute overlaps of full guide-pair range with regulatory regions
  guide.regions <- GRanges(seqnames = input.info$chrom, 
                           ranges = IRanges(input.info$start, input.info$end))
  
  # combine exons and enhancers into one data frame
  all.enhancers <- input.frame$enhancer
  all.enhancers$sortFactor <- sort.factor
  
  all.exons <- input.frame$exon
  all.exons$sortFactor <- 1
  
  functional.sequences <- rbind(all.exons, all.enhancers)
  functional.sequences.ranges <- GRanges(seqnames = functional.sequences$chrom,
                                         ranges = IRanges(functional.sequences$start, functional.sequences$end))
  
  guide.overlaps <- as.data.frame(findOverlaps(guide.regions, functional.sequences.ranges, type = 'any'))
  guide.overlap.list <- split(guide.overlaps, guide.overlaps$queryHits)
  
  # for each overlapping guide
  #   set up variable for keeping track of all FS counts
  #   compute the bp-density / proportion of cells affected at each position (dependent on uniform or normal distribution)
  #   order the regulatory regions overlapped according to their strength
  #   set up variable for storing unique cell indexes
  for(i in 1:length(guide.overlap.list)){
    
    temp.fs.counts <- c()
    
    temp.overl <- guide.overlap.list[[i]]
    temp.overl.fs <- functional.sequences[unique(temp.overl$subjectHits),, drop = F]
    temp.overl.fs.ordered <- temp.overl.fs[order(temp.overl.fs$sortFactor),]
    
    temp.guide.idx <- temp.overl$queryHits[1]
    temp.guide <- input.info[temp.guide.idx,,drop = F]
    temp.g1.efficiency <- guide.1.efficiencies[temp.guide.idx]
    temp.g2.efficiency <- guide.2.efficiencies[temp.guide.idx]
    temp.bp.prop <- compute_paired_bp_sort_prop(temp.guide, input.frame, temp.g1.efficiency, temp.g2.efficiency)
    temp.bp.prop.ranges <- GRanges(seqnames = temp.bp.prop$chrom,
                                   ranges = IRanges(temp.bp.prop$start, temp.bp.prop$end))
    
    temp.unique.cell.idx <- c()
    temp.cells.with.guide <- input.pool.counts[temp.guide.idx]
    
    #   for each regulatory region overlapped
    #     find the overlaps with the bp.prop and use the max()
    #     sample cell indexes
    #     (store the cell indexes in a list?)
    #     add to already recorded indexes and keep track of the new indexes
    #     compute the dirichlet probabilities
    #     compute the counts given the probabilities and number of cells affected
    for(o in 1:nrow(temp.overl.fs.ordered)){
      
      temp.fs <- temp.overl.fs.ordered[o,,drop = F]
      temp.fs.ranges <- GRanges(seqnames = temp.fs$chrom,
                                ranges = IRanges(temp.fs$start, temp.fs$end))
      
      bp.overlaps <- as.data.frame(findOverlaps(temp.bp.prop.ranges, temp.fs.ranges,type = 'any'))
      
      temp.max.prop <- max(temp.bp.prop$bpProb[unique(bp.overlaps$queryHits)])
      
      temp.functional.cells <- round(temp.cells.with.guide * temp.max.prop)
      temp.fs.idx <- sample(c(1:temp.cells.with.guide), size = temp.functional.cells, replace = F)
      temp.nr.cells.affected <- length(which(! temp.fs.idx %in% temp.unique.cell.idx))
      temp.unique.cell.idx <- unique(c(temp.unique.cell.idx, temp.fs.idx))
      
      temp.functional.sorting.prob <- c()
      if(input.frame$screenType == 'selectionScreen'){
        # temp.functional.sorting.prob <- t(apply(t(effect.diff %*% t(fs.df$sortFactor[temp.overlaps$subjectHits])), 1, function(x){neg.sort.freq -  x}))
        temp.functional.sorting.prob <- (input.frame$negProp - (effect.diff * temp.fs$sortFactor)) * guide.disp[temp.guide.idx]
        # temp.functional.sorting.prob <- input.frame$negSortingFrequency - (effect.diff * temp.fs$sortFactor)
      } else {
        # temp.functional.sorting.prob <- t(apply(t(effect.diff %*% t(fs.df$sortFactor[temp.overlaps$subjectHits])), 1, function(x){neg.sort.freq +  x}))
        temp.functional.sorting.prob <- (input.frame$negProp + (effect.diff * temp.fs$sortFactor)) * guide.disp[temp.guide.idx]
        # temp.functional.sorting.prob <- input.frame$negSortingFrequency + (effect.diff * temp.fs$sortFactor)
      }
      
      temp.functional.sorting.prob.rdir <- generate_dirichlet_probs(temp.functional.sorting.prob)
      temp.functional.counts <- t(rmultinom(1, size = temp.functional.cells,
                                            prob = temp.functional.sorting.prob.rdir ))
      temp.fs.counts <- rbind(temp.fs.counts, temp.functional.counts)
    }
    
    #   compute the number of cells left with background sorting
    #   compute background counts
    #   combine background and FS counts
    #   update the repl.counts data frame
    nonFunctional.cells <- temp.cells.with.guide - length(temp.unique.cell.idx)
    nonFunctional.sorting.prob.rdir <- generate_dirichlet_probs(adj.bkg.freq[temp.guide.idx,])
    temp.bkg.counts <- t(rmultinom(1, size = nonFunctional.cells,
                               prob = nonFunctional.sorting.prob.rdir ))
    
    temp.all.counts <- temp.bkg.counts + colSums(temp.fs.counts)
    
    out.repl.counts[temp.overl$queryHits[1], ] <- temp.all.counts
    
  }
  
  return(out.repl.counts)
}


#' @title Comput the effect on sorting proportions for each base pairs overlapped by a guide(pair)
#' @param input.frame: list containing all the necessary parameters: $posSortingFrequency, $negSortingFrequency, $areaOfEffect_type
#' @param guide.info: guide data frame (1 row), $chrom, $start, $end, $sgTarget_1, $sgTarget_2
#' @param g1.efficiency: guide 1 efficiency
#' @param g2.efficiency: guide 2 efficiency
#' @return data frame: $bpProb
#' @export compute_paired_bp_sort_prop()

compute_paired_bp_sort_prop <- function(guide.info, input.frame, g1.efficiency, g2.efficiency){
  
  bp.coordinates <- c(guide.info$start:guide.info$end)
  out.df <- data.frame(chrom = rep(guide.info$chrom[1], length(bp.coordinates)),
                       start = bp.coordinates, end = bp.coordinates)
  
  if(input.frame$areaOfEffect_type == 'uniform'){
    
    combined.gEff <- g1.efficiency * g2.efficiency * input.frame$deletionProb
    out.df$bpProb <- combined.gEff
    return(out.df)
    
  } else if(input.frame$areaOfEffect_type == 'normal'){
    
    normal.sd <- input.frame$normal_areaOfEffect_sd
    
    # guide 1 effect scaled by distance
    g1.dist.scaling <- dnorm(bp.coordinates, mean = guide.info$sgTarget_1, sd = normal.sd) / dnorm(guide.info$sgTarget_1, mean = guide.info$sgTarget_1, sd = normal.sd)
    g1.probs <- g1.dist.scaling * g1.efficiency
    
    # guide 2 effect scaled by distance
    g2.dist.scaling <- dnorm(bp.coordinates, mean = guide.info$sgTarget_2, sd = normal.sd) / dnorm(guide.info$sgTarget_2, mean = guide.info$sgTarget_2, sd = normal.sd)
    g2.probs <- g2.dist.scaling * g2.efficiency
    
    # deletion probability, set the tails to 0 as they are not part of the deletion
    del.prob <- g1.efficiency * g2.efficiency * input.frame$deletionProb
    bp.del.prob <- rep(del.prob, length(bp.coordinates))
    bp.del.prob[c(1:input.frame$crisprEffectRange, (length(bp.del.prob) - input.frame$crisprEffectRange):(length(bp.del.prob)))] <- 0
    
    all.probs <- cbind(g1.probs, g2.probs, bp.del.prob)
    max.bp.prob <- unlist(apply(all.probs, 1, max))
    
    out.df$bpProb <- max.bp.prob
    return(out.df)
  }
  
}


#' @title adjust the counts for guides overlapping exons, use uniform model
#' @param input.frame: list containing all the necessary parameters: $posSortingFrequency, $negSortingFrequency, if dualCRISPR: $deletionSize
#' @param effect.diff: difference in dirichlet probabilities between null and alternative
#' @param repl.counts: data frame wih the null counts of the replicate
#' @param guide.ranges: guide info as GRanges object
#' @param all.guide.efficiencies: guide efficiencies
#' @param input.pool.counts: counts from the input pool
#' @param guide.disp: vector, contains the dispersion for each guide
#' @param bkg.freq: dirichlet background parameters (including dispersion)
#' @return data frame: repl.counts
#' @export compute_uniform_exon_overlap()
 
compute_uniform_exon_overlap <- function(input.frame, effect.diff, 
                                         repl.counts, guide.ranges, all.guide.efficiencies,
                                         input.pool.counts, guide.disp, bkg.freq){
  
  for(j in 1:nrow(input.frame$exon)){
  
    temp.exon.ranges <- GRanges(seqnames = input.frame$exon$chrom[j],
                                ranges = IRanges(input.frame$exon$start[j], input.frame$exon$end[j]))
    
    temp.overlaps <- as.data.frame(findOverlaps(guide.ranges, temp.exon.ranges, type = 'any'))
    
    unique.guide.overlaps <- unique(temp.overlaps$queryHits)
    overlap.guide.efficiency <- all.guide.efficiencies[unique.guide.overlaps]
    
    temp.exon.sort.prob <- c()
    
    for(g in 1:length(unique.guide.overlaps)){
      temp.cells.with.guide <- input.pool.counts[unique.guide.overlaps[g]]
      temp.guide.efficiency <- overlap.guide.efficiency[g]
      temp.guide.disp <- guide.disp[unique.guide.overlaps[g]]
      temp.bkg.freq <- bkg.freq[unique.guide.overlaps[g],]
      
      temp.functional.sorting.prob <- c()
      if(input.frame$screenType == 'selectionScreen'){ 
        temp.functional.sorting.prob <- t(input.frame$negProp - matrix(effect.diff)) * temp.guide.disp
      } else {
        temp.functional.sorting.prob <- t(input.frame$negProp + matrix(effect.diff)) * temp.guide.disp
      }

      temp.guide.counts <- generate_guide_counts(temp.cells.with.guide, temp.guide.efficiency, 
                                                 temp.functional.sorting.prob, temp.bkg.freq) #temp.functional.counts + temp.nonFunctional.counts
      repl.counts[unique.guide.overlaps[g],] <- temp.guide.counts
    }
  }
  
  return(repl.counts)
  
}


#' @title adjust the counts for guides overlapping exons, use uniform model
#' @param input.frame: list containing all the necessary parameters: $posSortingFrequency, $negSortingFrequency, if dualCRISPR: $deletionSize
#' @param effect.diff: difference in dirichlet probabilities between null and alternative
#' @param repl.counts: data frame wih the null counts of the replicate
#' @param guide.ranges: guide info as GRanges object
#' @param all.guide.efficiencies: guide efficiencies
#' @param input.pool.counts: counts from the input pool
#' @param sort.factor: strength of each enhancer
#' @param guide.disp: vector, contains the dispersion for each guide
#' @param bkg.freq: dirichlet background parameters (including dispersion)
#' @return data frame: repl.counts
#' @export compute_uniform_enhancer_overlap()

compute_uniform_enhancer_overlap <- function(input.frame, effect.diff, 
                                         repl.counts, guide.ranges, all.guide.efficiencies,
                                         input.pool.counts, sort.factor, guide.disp, bkg.freq){
  
  # for every enhancer, obtain the enhancer strength from the beta distribution
  # for every guide, get the guide efficiency
  # guides get sorted either according to enhancer strength or as negative sorting
  for(k in 1:nrow(input.frame$enhancer)){
    temp.enhancer.ranges <- GRanges(seqnames = input.frame$enhancer$chrom[k],
                                    ranges = IRanges(input.frame$enhancer$start[k], input.frame$enhancer$end[k]))
    
    temp.overlaps <- as.data.frame(findOverlaps(guide.ranges, temp.enhancer.ranges, type = 'any'))
    
    if(nrow(temp.overlaps) > 0){
      unique.overlaps <- unique(temp.overlaps$queryHits)
      
      guide.efficiency <- all.guide.efficiencies[unique.overlaps]
      filt.guide.disp <- guide.disp[unique.overlaps]
      filt.bkg.freq <- bkg.freq[unique.overlaps,]
      
      temp.enhancer.sort.prob <- c()
      for(e in 1:length(guide.efficiency)){
        temp.cells.with.guide <- input.pool.counts[unique.overlaps[e]]
        temp.guide.efficiency <- guide.efficiency[e]
        temp.guide.disp <- filt.guide.disp[e]
        temp.bkg.freq <- filt.bkg.freq[e,]
        
        temp.functional.sorting.prob <- c()
        if(input.frame$screenType == 'selectionScreen'){
          temp.functional.sorting.prob <- t(input.frame$negProp - matrix(effect.diff) * sort.factor[k]) * temp.guide.disp
        } else {
          temp.functional.sorting.prob <- t(sort.factor[k] * matrix(effect.diff) + input.frame$negProp) * temp.guide.disp
        }
        
        temp.guide.counts <- generate_guide_counts(temp.cells.with.guide, temp.guide.efficiency, 
                                                   temp.functional.sorting.prob, temp.bkg.freq) #temp.functional.counts + temp.nonFunctional.counts
        repl.counts[unique.overlaps[e],] <- temp.guide.counts
      }
    }
  }
  
  return(repl.counts)
  
}


#' @title give different genomic properties different sorting probabilities
#' @param cells.with.guide: integer, number of cells containing a guide
#' @param prct.cells.affected: between 0 and 1, efficiency of guide
#' @param temp.functional.sorting.prob: adjusted sorting probabilites
#' @param neg.sort.freq: background sorting frequency
#' @return vector: random sample of dirichlet probabilities
#' @export generate_guide_counts()

generate_guide_counts <- function(cells.with.guide, prct.cells.affected, 
                                  temp.functional.sorting.prob, neg.sort.freq, fs.strength = NULL){
  
  if(! is.null(fs.strength)){
    prct.cells.affected <- prct.cells.affected[order(fs.strength, decreasing = F)]
    temp.functional.sorting.prob <- temp.functional.sorting.prob[order(fs.strength, decreasing = F),,drop = F]
  }
  
  # identify the number of cells that are affected by each FS
  fs.idx <- c()
  nr.cells.affected <- c()
  for(i in 1:length(prct.cells.affected)){
    temp.functional.cells <- round(cells.with.guide * prct.cells.affected[i])
    temp.fs.idx <- sample(c(1:cells.with.guide), size = temp.functional.cells, replace = F)
    temp.nr.cells.affected <- length(which(! temp.fs.idx %in% fs.idx))
    nr.cells.affected <- c(nr.cells.affected, temp.nr.cells.affected)
    fs.idx <- unique(c(fs.idx, temp.fs.idx))
  }
  
  nonFunctional.cells <- cells.with.guide - sum(nr.cells.affected)
  nonFunctional.sorting.prob.rdir <- generate_dirichlet_probs(neg.sort.freq)
  pool.counts <- t(rmultinom(1, size = nonFunctional.cells,
                                           prob = nonFunctional.sorting.prob.rdir ))
  
  for(i  in 1:length(nr.cells.affected)){
    temp.functional.sorting.prob.rdir <- generate_dirichlet_probs(temp.functional.sorting.prob[i,])
    temp.functional.counts <- t(rmultinom(1, size = temp.functional.cells,
                                          prob = temp.functional.sorting.prob.rdir ))
    pool.counts <- pool.counts + temp.functional.counts
  }
  
  return(pool.counts)
  
}


#' @title give different genomic properties different sorting probabilities
#' @param input.sorting.prob: dirichlet sorting probabilitites from which to take a random sample
#' @return vector: random sample of dirichlet probabilities
#' @export generate_dirichlet_probs()

generate_dirichlet_probs <- function(input.sorting.prob){
  
  temp.sorting.prob.rdir <- rdirichlet(1, input.sorting.prob)
  
  if(length(which(is.na(temp.sorting.prob.rdir))) > 0){
    temp.nan.rows <- which(is.na(temp.sorting.prob.rdir[,1]))
    temp.fill.values <- rep(1 / ncol(temp.sorting.prob.rdir),
                            ncol(temp.sorting.prob.rdir))
    temp.sorting.prob.rdir[which(is.na(temp.sorting.prob.rdir[,1])),] <- temp.fill.values
  }
  
  return(temp.sorting.prob.rdir)
}


#' @title compute the shift in guide counts resulting from directly overlapping a functional sequence
#' @param repl.counts: data.frame with the counts for current replicate
#' @param fs.df: functional sequence data.frame: $chrom, $start, $end, $sortFactor
#' @param all.guide.efficiencies: guide efficiencies
#' @param screen.type: either 'selectionScreen' or 'FACS'
#' @param direct.fs.overlaps.list: list, each element is a df, overlaps of guide (queryHits) with fs (subjectHits)
#' @param effect.diff: difference in dirichlet probabilities between null and alternative
#' @param input.counts: counts per cell in the input pool
#' @param guide.disp: vector, contains the dispersion for each guide
#' @param adj.bkg.freq: dirichlet background parameters (including dispersion)
#' @param adj.fs.freq: dirichlet FS parameters (including dispersion)
#' @return data.frame with the updated counts
#' @export compute_normal_counts_directOverl()

compute_normal_counts_directOverl <- function(repl.counts, fs.df, all.guide.efficiencies, input.frame,
                                  direct.fs.overlaps.list, effect.diff, input.counts,
                                  adj.bkg.freq, adj.fs.freq, guide.disp){
  
  for(i in 1:length(direct.fs.overlaps.list)){
    
    temp.overlaps <- direct.fs.overlaps.list[[i]]
    guide.idx <- temp.overlaps$queryHits[1]
    
    temp.bkg.freq <- adj.bkg.freq[guide.idx, ]
    temp.fs.freq <- adj.fs.freq[guide.idx, ]
    temp.guide.disp <- guide.disp[guide.idx]
    
    # establish the dirichlet probabilities
    temp.functional.sorting.prob <- c()
    if(input.frame$screenType == 'selectionScreen'){
      temp.scaled.effect.diff <- t(effect.diff %*% t(fs.df$sortFactor[temp.overlaps$subjectHits]))
      temp.freq.diff <- t(apply(temp.scaled.effect.diff, 1, function(x){input.frame$negProp - x}))
      temp.functional.sorting.prob <- temp.freq.diff * temp.guide.disp
    } else {
      temp.scaled.effect.diff <- t(effect.diff %*% t(fs.df$sortFactor[temp.overlaps$subjectHits]))
      temp.freq.diff <- t(apply(temp.scaled.effect.diff, 1, function(x){input.frame$negProp + x}))
      temp.functional.sorting.prob <- temp.freq.diff * temp.guide.disp
    }
    
    guide.idx <- temp.overlaps$queryHits[1]
    
    temp.cells.with.guide <- input.counts[guide.idx]
    temp.guide.efficiency <- all.guide.efficiencies[guide.idx]
    
    # modified to account for dispersion adjusted background frequencies
    temp.guide.counts <- generate_guide_counts(temp.cells.with.guide, temp.guide.efficiency,
                                               temp.functional.sorting.prob, temp.bkg.freq, 
                                               fs.strength = fs.df$sortFactor[temp.overlaps$subjectHits])
    repl.counts[guide.idx,] <- temp.guide.counts
    
  }
  
  return(repl.counts)
}

#' @title compute the shift in guide counts resulting from indirectly overlapping a functional sequence
#' @param repl.counts: data.frame with the counts for current replicate
#' @param fs.df: functional sequence data.frame: $chrom, $start, $end, $sortFactor
#' @param all.guide.efficiencies: guide efficiencies
#' @param screen.type: either 'selectionScreen' or 'FACS'
#' @param indirect.fs.overlaps.list: list, each element is a df, overlaps of guide (queryHits) with fs (subjectHits)
#' @param effect.diff: difference in dirichlet probabilities between null and alternative
#' @param neg.sort.freq: sorting probabilities of the background
#' @param input.counts: counts per cell in the input pool
#' @param input.info: info data frame with the sgRNA target coordinates
#' @param normal.sd: standard deviation to use for the normal distr.
#' @param guide.disp: vector, contains the dispersion for each guide
#' @param adj.bkg.freq: dirichlet background parameters (including dispersion)
#' @param adj.fs.freq: dirichlet FS parameters (including dispersion)
#' @return data.frame with the updated counts
#' @export compute_normal_counts_indirectOverl()

compute_normal_counts_indirectOverl <- function(repl.counts, fs.df, all.guide.efficiencies, screen.type,
                                              indirect.fs.overlaps.list, effect.diff, neg.sort.prop, input.counts,
                                              input.info, normal.sd, adj.bkg.freq, adj.fs, guide.disp){
  
  for(i in 1:length(indirect.fs.overlaps.list)){
    
    temp.overlaps <- indirect.fs.overlaps.list[[i]]
    
    guide.idx <- temp.overlaps$queryHits[1]
    
    # for each overlap, it has to be determined what the distance to the closest start or end-site is from the sgTarget
    # use distance to compute the scaling of the effect.diff
    dist.scaling <- c()
    temp.functional.sorting.prob <- c()
    for(o in 1:nrow(temp.overlaps)){
      temp.min.dist <- min(c(abs(input.info$sgTarget[temp.overlaps$queryHits[o]] - fs.df$start[temp.overlaps$subjectHits[o]]), 
                             abs(input.info$sgTarget[temp.overlaps$queryHits[o]] - fs.df$end[temp.overlaps$subjectHits[o]])))
      
      # normalized scaling by distance
      temp.dist.scaling <- dnorm(temp.min.dist, mean = 0, sd = normal.sd) / dnorm(0, mean = 0, sd = normal.sd)
      dist.scaling <- c(dist.scaling, temp.dist.scaling)
      
      temp.enh.strength <- fs.df$sortFactor[temp.overlaps$subjectHits[o]]
      if(screen.type == 'selectionScreen'){
        temp.functional.sorting.prob <- rbind(temp.functional.sorting.prob, t(neg.sort.prop - effect.diff * temp.enh.strength) * guide.disp[guide.idx])
      } else {
        temp.functional.sorting.prob <- rbind(temp.functional.sorting.prob, t(temp.enh.strength * effect.diff + neg.sort.prop) * guide.disp[guide.idx])
      }
      
    }

    temp.cells.with.guide <- input.counts[guide.idx]
    temp.guide.efficiency <- all.guide.efficiencies[guide.idx]
    temp.cells.affected <- temp.guide.efficiency * dist.scaling
    temp.guide.counts <- generate_guide_counts(temp.cells.with.guide, temp.cells.affected,
                                               temp.functional.sorting.prob, adj.bkg.freq[guide.idx, ],
                                               fs.df$sortFactor[temp.overlaps$subjectHits])
    repl.counts[guide.idx,] <- temp.guide.counts
    
  }
  
  return(repl.counts)
}

  
#' @title Compute the area of effect using a normal distribution for single guide simulations
#' @param input.frame: list containing all the necessary parameters: $posSortingFrequency, $negSortingFrequency, if dualCRISPR: $deletionSize
#' @param effect.diff: difference in dirichlet probabilities between null and alternative
#' @param repl.counts: data frame wih the null counts of the replicate
#' @param guide.ranges: guide info as GRanges object
#' @param all.guide.efficiencies: guide efficiencies
#' @param input.info: info data frame with the sgRNA target coordinates
#' @param sort.factor: beta distributed fractions of signal strength relative to exons
#' @param input.pool.counts: counts per cell in the input pool
#' @param guide.disp: vector, contains the dispersion for each guide
#' @param adj.bkg.freq: dirichlet background parameters (including dispersion)
#' @param adj.fs.freq: dirichlet FS parameters (including dispersion)
#' @return data frame: repl.counts
#' @export compute_normal_areaOfEffect()

compute_normal_areaOfEffect <- function(input.frame, effect.diff, repl.counts, 
                                        guide.ranges, all.guide.efficiencies,
                                        input.info, sort.factor, input.pool.counts,
                                        adj.bkg.freq, adj.fs.freq, guide.disp){
  
  # combine exons and enhancers into one data frame
  all.enhancers <- input.frame$enhancer
  all.enhancers$sortFactor <- sort.factor
  
  all.exons <- input.frame$exon
  all.exons$sortFactor <- 1
  
  functional.sequences <- rbind(all.exons, all.enhancers)
  functional.sequences.ranges <- GRanges(seqnames = functional.sequences$chrom,
                                         ranges = IRanges(functional.sequences$start, functional.sequences$end))
  
  # compute all direct overlaps
  guide.sgTarget.ranges <- GRanges(seqnames = input.info$chrom,
                                   ranges = IRanges(input.info$sgTarget, input.info$sgTarget))
  direct.fs.overlaps <- as.data.frame(findOverlaps(guide.sgTarget.ranges, functional.sequences.ranges, type = 'within'))
  unique.direct.overlaps <- unique(direct.fs.overlaps$queryHits)
  
  # for all direct overlaps, split according to guide
  # if multiple overlaps; if exon present use exon,
  # else compute strongest signal based on distance and enhancer strength
  direct.fs.overlaps.list <- split(direct.fs.overlaps, direct.fs.overlaps$queryHits)
  
  repl.counts <- compute_normal_counts_directOverl(repl.counts, functional.sequences, all.guide.efficiencies, 
                                                   input.frame, direct.fs.overlaps.list, effect.diff, 
                                                   input.pool.counts, adj.bkg.freq, adj.fs.freq, guide.disp)
  
  # compute all overlaps
  # identify all guides only indirectly overlapping FS
  guide.area.ranges <- GRanges(seqnames = input.info$chrom,
                                   ranges = IRanges(input.info$start, input.info$end))
  all.fs.overlaps <- as.data.frame(findOverlaps(guide.area.ranges, functional.sequences.ranges, type = 'any'))
  indirect.fs.overlaps <- all.fs.overlaps[-match(unique.direct.overlaps, all.fs.overlaps$queryHits),]
  unique.indirect.overlaps.list <- split(indirect.fs.overlaps, indirect.fs.overlaps$queryHits)
  repl.counts <- compute_normal_counts_indirectOverl(repl.counts, functional.sequences, all.guide.efficiencies, 
                                                   input.frame$screenType, unique.indirect.overlaps.list, effect.diff, 
                                                   input.frame$negProp, input.pool.counts,
                                                   input.info, input.frame$normal_areaOfEffect_sd,
                                                   adj.bkg.freq, adj.fs.freq, guide.disp)
  
  return(repl.counts)

}

#' @title give different genomic properties different sorting probabilities, formerly 'full_replicate_simulation_sepDistrSampl_v2'
#' @param input.frame: list containing all the necessary parameters: $posSortingFrequency, $negSortingFrequency, if dualCRISPR: $deletionSize
#' @param input.info: info file
#' @param sim.nr: which simulation number
#' @return list: counts, enhancerStrength, guide_efficiencies
#' @export single_guide_replicate_simulation()

# give different genomic properties different sorting probabilities
# this method samples from either the null sorting distrubuton or samples
# from an enhancer-strength distribution
# x% of guides have an enhancer effect, remaining guides have no effect
# this implementation also records the guide efficiency
single_guide_replicate_simulation <- function(input.frame, input.info, sim.nr){
  
  out.sim.data <- c()
  out.repl.sort.fractions.filt <- c()
  
  guide.ranges <- GRanges(seqnames = input.info$chrom,
                          ranges = IRanges(input.info$start, input.info$end))
  
  
  effect.diff <- c()
  
  # Divide Negative Sorting Frequency to Obtain Proportions
  input.frame$negProp <- input.frame$negSortingFrequency / sum(input.frame$negSortingFrequency)
  input.frame$posProp <- input.frame$posSortingFrequency / sum(input.frame$posSortingFrequency)
  
  # does this lead to vectors like (4, 0)?
  if(input.frame$screenType == 'selectionScreen'){
    effect.diff <- input.frame$negProp - input.frame$posProp
  } else {
    effect.diff <- input.frame$posProp - input.frame$negProp
  }
  
  sort.factor <- rbeta(nrow(input.frame$enhancer), shape1 = input.frame$enhancerShape1, shape2 = input.frame$enhancerShape2)
  all.guide.efficiencies <- rbeta(nrow(input.info), shape1 = input.frame$guideShape1, shape2 = input.frame$guideShape2)
  
  # for the number of replicates
  for(i in 1:length(input.frame$inputPools)){
    print(paste0('Generating replicate ', i))
    
    # Scale the Negative Sorting Frequencies according to the Radical Fit
    # Dispersion Estimation Coefficients for MYC radical r2: -36.21266742  -0.03280368   3.28222794
    input.distr <- input.frame$inputPools[[i]]
    
    # Obtain Dispersion for Each Guide and Save in Vector
    print('Computing dispersions')
    guide.disp <- c()
    if(input.frame$dispersionType == 'independent'){
      guide.disp <- rep(sum(input.frame$negSortingFrequency), length(input.distr))
    } else if(input.frame$dispersionType == 'exponential'){
      
      # radical fit
      # guide_disp <- -36.21266742 + -0.03280368 * input.distr + 3.28222794 * sqrt(input.distr)
      
      # exponential fit with linear term
      # guide.disp <- -1.092370e+02 + -4.259984e-03 * input.distr + 2.139774e+01 * log(input.distr)
      
      # exponential fit without linear term
      guide.disp <- -89.17911 + 17.78130 * log(input.distr)
      
      # set floor value of 3 for dispersion to avoid negative probabilities
      guide.disp[guide.disp < 3] <- 3

    } else if(input.frame$dispersionType == 'spline'){
      disp.model <- readRDS(input.frame$splineModelLoc)
      total.counts.df <- data.frame(counts = input.distr, stringsAsFactors = F)
      spline.predict <- suppressWarnings(predict(disp.model, total.counts.df))
      guide.disp <- spline.predict
      guide.disp[guide.disp < 1] <- 1
    }

    # compute the per-guide dirichlet parameters by including the dispersion
    temp.adj.bkg.freq <- t(t(guide.disp)) %*% input.frame$negProp 
    temp.adj.fs.freq <- t(t(guide.disp)) %*% input.frame$posProp
    
    print('Generating sorting probabilities')
    temp.sort.prob <- matrix(0, nrow = length(input.distr), ncol = length(input.frame$negSortingFrequency))
    for (j in 1:length(input.distr)) {
       temp.sort.prob[j, ] <- rdirichlet(1, temp.adj.bkg.freq[j, ])
    }

    if(length(which(is.na(temp.sort.prob))) > 0){
      print('Adjust edge-case probabilities')
      temp.nan.rows <- which(is.na(temp.sort.prob[,1]))
      temp.fill.values <- rep(1 / ncol(temp.sort.prob), ncol(temp.sort.prob))
      temp.sort.prob[temp.nan.rows,] <- temp.fill.values
    }
    
    print('Generate intitial counts')
    repl.counts <- counts_from_probs(temp.sort.prob, input.frame$inputPools[[i]])
    

    print('Compute area of effect')
    # adjust counts of guides overlapping functional sequences
    if(input.frame$areaOfEffect_type == 'uniform'){
      repl.counts <- compute_uniform_exon_overlap(input.frame, effect.diff, 
                                                  repl.counts, guide.ranges, all.guide.efficiencies,
                                                  input.frame$inputPools[[i]], guide.disp, temp.adj.bkg.freq)
      
      repl.counts <- compute_uniform_enhancer_overlap(input.frame, effect.diff, 
                                                   repl.counts, guide.ranges, all.guide.efficiencies,
                                                   input.frame$inputPools[[i]], sort.factor, guide.disp, temp.adj.bkg.freq)
      
    } else if(input.frame$areaOfEffect_type == 'normal'){
      # test this first
      repl.counts <- compute_normal_areaOfEffect(input.frame, effect.diff, 
                                                 repl.counts, guide.ranges, all.guide.efficiencies,
                                                 input.info, sort.factor, input.frame$inputPools[[i]],
                                                 temp.adj.bkg.freq, temp.adj.fs.freq, guide.disp)
    }

    print('Perform sequencing')
    repl.sequenced <- experiment_sequencing(cbind(input.frame$inputPools[[i]], repl.counts),
                                            input.frame$pcrDupl, input.frame$seqDepth[[i]])
    repl.sort.fractions <- colSums(cbind(input.frame$inputPools[[i]], repl.counts)) / sum(colSums(cbind(input.frame$inputPools[[i]], repl.counts)))
    
    repl.sequenced.filtered <- repl.sequenced
    repl.sort.fractions.filt <- repl.sort.fractions
    if(input.frame$screenType == 'selectionScreen'){ 
      repl.sequenced.filtered <- repl.sequenced[,c(1,2)]
      repl.sort.fractions.filt <- repl.sort.fractions[c(1,2)]
    }
    
    colnames(repl.sequenced.filtered) <- paste(paste0('sim', sim.nr, '_repl', i), input.frame$poolNames, sep = '_')
    out.sim.data <- cbind(out.sim.data, repl.sequenced.filtered)
    
    out.repl.sort.fractions.filt <- c(out.repl.sort.fractions.filt, repl.sort.fractions.filt)
  }
  return(list(counts = out.sim.data, 
              enhancerStrength = sort.factor, 
              guide_efficiencies = all.guide.efficiencies,
              pool_proportions = out.repl.sort.fractions.filt))
}


# this version was used for when the area of effect for guides was assumed to be uniform only
simulate_data_old <- function(input.list){

  input.list <- set_default_flags(input.list)

  updated.info <- input.list$guides  #input.info # need to update

  if(input.list$crisprSystem %in% c('CRISPRi', 'CRISPRa', 'CRISPRcas9')){
    updated.info$start <- updated.info$start - round(input.list$crisprEffectRange / 2)
    updated.info$end <- updated.info$end + round(input.list$crisprEffectRange / 2)
  }

  dir.create(file.path(paste0(input.list$outDir, input.list$simName)))

  for(sim in 1:input.list$nrSims){
    print(paste0('Starting simulation ', sim, ' out of ', input.list$nrSims))
    input.list$inputPools <- list()
    input.list$enhancer <- generate_enhancers(input.list, updated.info, input.list$exon)

    print('Generating initial pool')
    for(repl in 1:length(input.list$inputCellNr)){
      if(input.list$randomizeRepl == 'yes'){
        randomized.counts <- sample(input.list$inputGuideDistr[, repl])
        input.list$inputPools[[paste0('repl',repl)]] <- rmultinom(1,
          input.list$inputCellNr[repl], randomized.counts)
      } else {
        input.list$inputPools[[paste0('repl',repl)]] <- rmultinom(1,
          input.list$inputCellNr[repl], input.list$inputGuideDistr[, repl])
      }
    }

    temp.sim <- full_replicate_simulation_sepDistrSampl(input.list, updated.info, sim)
    
    write.csv(temp.sim$counts, file = paste0(input.list$outDir, input.list$simName,'/',
      input.list$simName, '_sim', sim, '_counts.csv'), row.names = F)

    input.list$enhancer$enhancerStrength <- temp.sim$enhancerStrength
    ordered.enhancers <- input.list$enhancer[order(input.list$enhancer$start),]

    write.csv(ordered.enhancers, file = paste0(input.list$outDir, input.list$simName,'/',
      input.list$simName, '_sim', sim, '_enhancers.csv'), row.names = F)

    init.labels <- rep('chr', nrow(updated.info))
    final.info <- input.list$guides
    final.info$label <- init.labels
    exon.label.pos <- sim_enhancer_overlaps(updated.info, input.list$exon) # also works for exons
    final.info$label[exon.label.pos] <- 'exon'
    enh.label.pos <- sim_enhancer_overlaps(updated.info, ordered.enhancers)
    final.info$label[enh.label.pos] <- 'pos'
    random.neg <- sample(which(final.info$label == 'chr'), input.list$nrNeg)
    for(i in 1:length(random.neg)){
      final.info$chrom[random.neg[i]] <- NA
      final.info$start[random.neg[i]] <- NA
      final.info$end[random.neg[i]] <- NA
      final.info$label[random.neg[i]] <- 'neg'
    }
    
    final.info$guideEfficiency <- temp.sim$guide_efficiencies

    
    # final.info[random.neg, c(1:4)] <- data.frame(NA, NA, NA, 'neg', stringsAsFactors = F)
    write.csv(final.info, file = paste0(input.list$outDir, input.list$simName,'/',
      input.list$simName, '_sim', sim, '_info.csv'), row.names = F)
  }
}

# below was discontinued on 19/05/31 and re-written to simplify the process of data simulation
# this does change a few things, especially with the input and therefore is not fully backward compatible
simulate_data_discontinued <- function(input.list, two.distr = NULL, using.facs = 'yes'){

  input.list <- set_default_flags(input.list)

  updated.info <- input.list$guides  #input.info # need to update

  if(input.list$crisprSystem %in% c('CRISPRi', 'CRISPRa', 'CRISPRcas9')){
    updated.info$start <- updated.info$start - round(input.list$crisprEffectRange / 2)
    updated.info$end <- updated.info$end + round(input.list$crisprEffectRange / 2)
  }

  dir.create(file.path(paste0(input.list$outDir, input.list$simName)))

  for(sim in 1:input.list$nrSims){
    input.list$inputPools <- list()
    input.list$enhancer <- generate_enhancers(input.list, updated.info, input.list$exon)

    for(repl in 1:length(input.list$inputCellNr)){
      if(input.list$randomizeRepl == 'yes'){
        randomized.counts <- sample(input.list$inputGuideDistr[, repl])
        input.list$inputPools[[paste0('repl',repl)]] <- rmultinom(1,
          input.list$inputCellNr[repl], randomized.counts)
      } else {
        input.list$inputPools[[paste0('repl',repl)]] <- rmultinom(1,
          input.list$inputCellNr[repl], input.list$inputGuideDistr[, repl])
      }
    }

    temp.sim <- c()
    if(using.facs == 'yes'){
      if(is.null(two.distr)){
        temp.sim <- full_replicate_simulation(input.list, updated.info, sim)
      } else {
        temp.sim <- full_replicate_simulation_sepDistrSampl(input.list, updated.info, sim)
      }
    } else {
      print('not yet implemented')
      break()
      if(is.null(two.distr)){
        print('not yet implemented')
        break()
      } else {
        # ToDo
      }
    }
    write.csv(temp.sim$counts, file = paste0(input.list$outDir, input.list$simName,'/',
      input.list$simName, '_sim', sim, '_counts.csv'), row.names = F)

    input.list$enhancer$enhancerStrength <- temp.sim$enhancerStrength
    ordered.enhancers <- input.list$enhancer[order(input.list$enhancer$start),]

    write.csv(ordered.enhancers, file = paste0(input.list$outDir, input.list$simName,'/',
      input.list$simName, '_sim', sim, '_enhancers.csv'), row.names = F)

    init.labels <- rep('chr', nrow(updated.info))
    final.info <- input.list$guides
    final.info$label <- init.labels
    exon.label.pos <- sim_enhancer_overlaps(updated.info, input.list$exon) # also works for exons
    final.info$label[exon.label.pos] <- 'exon'
    enh.label.pos <- sim_enhancer_overlaps(updated.info, ordered.enhancers)
    final.info$label[enh.label.pos] <- 'pos'
    random.neg <- sample(which(final.info$label == 'chr'), input.list$nrNeg)
    final.info[random.neg, c(1:4)] <- data.frame(NA, NA, NA, 'neg', stringsAsFactors = F)
    write.csv(final.info, file = paste0(input.list$outDir, input.list$simName,'/',
      input.list$simName, '_sim', sim, '_info.csv'), row.names = F)
  }
}

# give different genomic properties different sorting probabilities
# this method samples from either the null sorting distrubuton or samples
# from an enhancer-strength distribution
# x% of guides have an enhancer effect, remaining guides have no effect
# this implementation also records the guide efficiency
# this version was used for when the area of effect for guides was assumed to be uniform only
full_replicate_simulation_sepDistrSampl <- function(input.frame, input.info, sim.nr){

  out.sim.data <- c()

  guide.ranges <- GRanges(seqnames = input.info$chrom,
    ranges = IRanges(input.info$start, input.info$end))


  effect.diff <- c()
  if(input.frame$screenType == 'selectionScreen'){ #'selectionScreen' %in% names(input.frame)){
    effect.diff <- input.frame$negSortingFrequency - input.frame$posSortingFrequency
  } else {
    effect.diff <- input.frame$posSortingFrequency - input.frame$negSortingFrequency
  }
  # exon.neg.diff <- input.frame$posSortingFrequency - input.frame$negSortingFrequency
  sort.factor <- rbeta(nrow(input.frame$enhancer), shape1 = input.frame$enhancerShape1, shape2 = input.frame$enhancerShape2)
  all.guide.efficiencies <- rbeta(nrow(input.info), shape1 = input.frame$guideShape1, shape2 = input.frame$guideShape2)

  # for the number of replicates
  for(i in 1:length(input.frame$inputPools)){
    print(paste0('Generating replicate ', i))
    # add counts, assuming negative sorting probabilities
    temp.sort.prob <- rdirichlet(length(input.frame$inputPools[[i]]),
      input.frame$negSortingFrequency)
    if(length(which(is.na(temp.sort.prob))) > 0){
      temp.nan.rows <- which(is.na(temp.sort.prob[,1]))
      temp.fill.values <- rep(1 / ncol(temp.sort.prob), ncol(temp.sort.prob))
      temp.sort.prob[which(is.na(temp.sort.prob[,1])),] <- temp.fill.values
    }
    repl.counts <- counts_from_probs(temp.sort.prob, input.frame$inputPools[[i]])

    for(j in 1:nrow(input.frame$exon)){
      temp.exon.ranges <- GRanges(seqnames = input.frame$exon$chrom[j],
        ranges = IRanges(input.frame$exon$start[j], input.frame$exon$end[j]))

      temp.overlaps <- as.data.frame(findOverlaps(guide.ranges, temp.exon.ranges, type = 'any'))

      unique.guide.overlaps <- unique(temp.overlaps$queryHits)
      overlap.guide.efficiency <- all.guide.efficiencies[unique.guide.overlaps] #rbeta(length(unique.guide.overlaps), shape1 = input.frame$guideShape1, shape2 = input.frame$guideShape2)

      temp.exon.sort.prob <- c()

      # since dealing with exons, there is no enhancer efficiency
      # for each guide, the efficient guide gets sorted according to exon efficiency
      # the remaining guides get sorted according to negative sorting frequencies
      for(g in 1:length(unique.guide.overlaps)){
        temp.cells.with.guide <- input.frame$inputPools[[i]][unique.guide.overlaps[g]]
        temp.guide.efficiency <- overlap.guide.efficiency[g]

        temp.functional.cells <- round(temp.cells.with.guide * temp.guide.efficiency)
        temp.nonFunctional.cells <- temp.cells.with.guide - temp.functional.cells

        if(input.frame$screenType == 'selectionScreen'){ #'selectionScreen' %in% names(input.frame))
          temp.functional.sorting.prob <- input.frame$negSortingFrequency - effect.diff # + input.frame$posSortingFrequency
        } else {
          temp.functional.sorting.prob <- effect.diff + input.frame$negSortingFrequency
        }
        # temp.functional.sorting.prob <- exon.neg.diff + input.frame$negSortingFrequency

        temp.functional.sorting.prob.rdir <- rdirichlet(1, temp.functional.sorting.prob)
        temp.nonFunctional.sorting.prob.rdir <- rdirichlet(1, input.frame$negSortingFrequency)

        if(length(which(is.na(temp.functional.sorting.prob.rdir))) > 0){
          temp.nan.rows <- which(is.na(temp.functional.sorting.prob.rdir[,1]))
          temp.fill.values <- rep(1 / ncol(temp.functional.sorting.prob.rdir),
            ncol(temp.functional.sorting.prob.rdir))
          temp.functional.sorting.prob.rdir[which(is.na(temp.functional.sorting.prob.rdir[,1])),] <- temp.fill.values
        }

        if(length(which(is.na(temp.nonFunctional.sorting.prob.rdir))) > 0){
          temp.nan.rows <- which(is.na(temp.nonFunctional.sorting.prob.rdir[,1]))
          temp.fill.values <- rep(1 / ncol(temp.nonFunctional.sorting.prob.rdir),
            ncol(temp.nonFunctional.sorting.prob.rdir))
          temp.nonFunctional.sorting.prob.rdir[which(is.na(temp.nonFunctional.sorting.prob.rdir[,1])),] <- temp.fill.values
        }

        temp.functional.counts <- t(rmultinom(1, size = temp.functional.cells,
          prob = temp.functional.sorting.prob.rdir ))
        temp.nonFunctional.counts <- t(rmultinom(1, size = temp.nonFunctional.cells,
          prob = temp.nonFunctional.sorting.prob.rdir ))

        temp.guide.counts <- temp.functional.counts + temp.nonFunctional.counts
        repl.counts[unique.guide.overlaps[g],] <- temp.guide.counts
      }
    }

    # for every enhancer, obtain the enhancer strength from the beta distribution
    # for every guide, get the guide efficiency
    # guides get sorted either according to enhancer strength or as negative sorting
    for(k in 1:nrow(input.frame$enhancer)){
      temp.enhancer.ranges <- GRanges(seqnames = input.frame$enhancer$chrom[k],
        ranges = IRanges(input.frame$enhancer$start[k], input.frame$enhancer$end[k]))

      temp.overlaps <- as.data.frame(findOverlaps(guide.ranges, temp.enhancer.ranges, type = 'any'))

      if(nrow(temp.overlaps) > 0){
        unique.overlaps <- unique(temp.overlaps$queryHits)

        guide.efficiency <- all.guide.efficiencies[unique.overlaps] #rbeta(length(unique.overlaps), shape1 = input.frame$guideShape1, shape2 = input.frame$guideShape2)
        temp.enhancer.sort.prob <- c()
        for(e in 1:length(guide.efficiency)){
          temp.cells.with.guide <- input.frame$inputPools[[i]][unique.overlaps[e]]
          temp.guide.efficiency <- guide.efficiency[e]

          temp.functional.cells <- round(temp.cells.with.guide * temp.guide.efficiency)
          temp.nonFunctional.cells <- temp.cells.with.guide - temp.functional.cells

          if(input.frame$screenType == 'selectionScreen'){ #'selectionScreen' %in% names(input.frame))
            temp.functional.sorting.prob <- input.frame$negSortingFrequency - effect.diff * sort.factor[k] # + input.frame$posSortingFrequency
            # temp.functional.sorting.prob <- sort.factor[k] * effect.diff + input.frame$posSortingFrequency
          } else {
            temp.functional.sorting.prob <- sort.factor[k] * effect.diff + input.frame$negSortingFrequency
          }
          # temp.functional.sorting.prob <- sort.factor[k] * exon.neg.diff + input.frame$negSortingFrequency

          temp.functional.sorting.prob.rdir <- rdirichlet(1, temp.functional.sorting.prob)
          temp.nonFunctional.sorting.prob.rdir <- rdirichlet(1, input.frame$negSortingFrequency)

          if(length(which(is.na(temp.functional.sorting.prob.rdir))) > 0){
            temp.nan.rows <- which(is.na(temp.functional.sorting.prob.rdir[,1]))
            temp.fill.values <- rep(1 / ncol(temp.functional.sorting.prob.rdir),
              ncol(temp.functional.sorting.prob.rdir))
            temp.functional.sorting.prob.rdir[which(is.na(temp.functional.sorting.prob.rdir[,1])),] <- temp.fill.values
          }

          if(length(which(is.na(temp.nonFunctional.sorting.prob.rdir))) > 0){
            temp.nan.rows <- which(is.na(temp.nonFunctional.sorting.prob.rdir[,1]))
            temp.fill.values <- rep(1 / ncol(temp.nonFunctional.sorting.prob.rdir),
              ncol(temp.nonFunctional.sorting.prob.rdir))
            temp.nonFunctional.sorting.prob.rdir[which(is.na(temp.nonFunctional.sorting.prob.rdir[,1])),] <- temp.fill.values
          }

          temp.functional.counts <- t(rmultinom(1, size = temp.functional.cells,
            prob = temp.functional.sorting.prob.rdir ))
          temp.nonFunctional.counts <- t(rmultinom(1, size = temp.nonFunctional.cells,
            prob = temp.nonFunctional.sorting.prob.rdir ))

          temp.guide.counts <- temp.functional.counts + temp.nonFunctional.counts
          repl.counts[unique.overlaps[e],] <- temp.guide.counts
        }
      }
    }
    repl.sequenced <- experiment_sequencing(cbind(input.frame$inputPools[[i]], repl.counts),
      input.frame$pcrDupl, input.frame$seqDepth[[i]])

    repl.sequenced.filtered <- repl.sequenced
    if(input.frame$screenType == 'selectionScreen'){ #'selectionScreen' %in% names(input.frame))
      repl.sequenced.filtered <- repl.sequenced[,c(1,2)]
    }

    colnames(repl.sequenced.filtered) <- paste(paste0('sim', sim.nr, '_repl', i), input.frame$poolNames, sep = '_')
    out.sim.data <- cbind(out.sim.data, repl.sequenced.filtered)
  }
  return(list(counts = out.sim.data, enhancerStrength = sort.factor, guide_efficiencies = all.guide.efficiencies))
}

# deprecated: did not record guide efficiency
# give different genomic properties different sorting probabilities
# this method samples from either the null sorting distrubuton or samples
# from an enhancer-strength distribution
# x% of guides have an enhancer effect, remaining guides have no effect
full_replicate_simulation_sepDistrSampl_v1 <- function(input.frame, input.info, sim.nr){

  out.sim.data <- c()

  guide.ranges <- GRanges(seqnames = input.info$chrom,
    ranges = IRanges(input.info$start, input.info$end))


  effect.diff <- c()
  if(input.frame$screenType == 'selectionScreen'){ # 'selectionScreen' %in% names(input.frame)){
    effect.diff <- input.frame$negSortingFrequency - input.frame$posSortingFrequency
  } else {
    effect.diff <- input.frame$posSortingFrequency - input.frame$negSortingFrequency
  }
  # exon.neg.diff <- input.frame$posSortingFrequency - input.frame$negSortingFrequency
  sort.factor <- rbeta(nrow(input.frame$enhancer), shape1 = input.frame$enhancerShape1, shape2 = input.frame$enhancerShape2)

  # for the number of replicates
  for(i in 1:length(input.frame$inputPools)){
    print(paste0('Generating replicate ', i))
    # add counts, assuming negative sorting probabilities
    temp.sort.prob <- rdirichlet(length(input.frame$inputPools[[i]]),
      input.frame$negSortingFrequency)
    if(length(which(is.na(temp.sort.prob))) > 0){
      temp.nan.rows <- which(is.na(temp.sort.prob[,1]))
      temp.fill.values <- rep(1 / ncol(temp.sort.prob), ncol(temp.sort.prob))
      temp.sort.prob[which(is.na(temp.sort.prob[,1])),] <- temp.fill.values
    }
    repl.counts <- counts_from_probs(temp.sort.prob, input.frame$inputPools[[i]])

    for(j in 1:nrow(input.frame$exon)){
      temp.exon.ranges <- GRanges(seqnames = input.frame$exon$chrom[j],
        ranges = IRanges(input.frame$exon$start[j], input.frame$exon$end[j]))

      temp.overlaps <- as.data.frame(findOverlaps(guide.ranges, temp.exon.ranges, type = 'any'))

      unique.guide.overlaps <- unique(temp.overlaps$queryHits)
      overlap.guide.efficiency <- rbeta(length(unique.guide.overlaps), shape1 = input.frame$guideShape1, shape2 = input.frame$guideShape2)

      temp.exon.sort.prob <- c()

      # since dealing with exons, there is no enhancer efficiency
      # for each guide, the efficient guide gets sorted according to exon efficiency
      # the remaining guides get sorted according to negative sorting frequencies
      for(g in 1:length(unique.guide.overlaps)){
        temp.cells.with.guide <- input.frame$inputPools[[i]][unique.guide.overlaps[g]]
        temp.guide.efficiency <- overlap.guide.efficiency[g]

        temp.functional.cells <- round(temp.cells.with.guide * temp.guide.efficiency)
        temp.nonFunctional.cells <- temp.cells.with.guide - temp.functional.cells

        if(input.frame$screenType == 'selectionScreen'){ #'selectionScreen' %in% names(input.frame)){
          temp.functional.sorting.prob <- input.frame$negSortingFrequency - effect.diff # + input.frame$posSortingFrequency
        } else {
          temp.functional.sorting.prob <- effect.diff + input.frame$negSortingFrequency
        }
        # temp.functional.sorting.prob <- exon.neg.diff + input.frame$negSortingFrequency

        temp.functional.sorting.prob.rdir <- rdirichlet(1, temp.functional.sorting.prob)
        temp.nonFunctional.sorting.prob.rdir <- rdirichlet(1, input.frame$negSortingFrequency)

        if(length(which(is.na(temp.functional.sorting.prob.rdir))) > 0){
          temp.nan.rows <- which(is.na(temp.functional.sorting.prob.rdir[,1]))
          temp.fill.values <- rep(1 / ncol(temp.functional.sorting.prob.rdir),
            ncol(temp.functional.sorting.prob.rdir))
          temp.functional.sorting.prob.rdir[which(is.na(temp.functional.sorting.prob.rdir[,1])),] <- temp.fill.values
        }

        if(length(which(is.na(temp.nonFunctional.sorting.prob.rdir))) > 0){
          temp.nan.rows <- which(is.na(temp.nonFunctional.sorting.prob.rdir[,1]))
          temp.fill.values <- rep(1 / ncol(temp.nonFunctional.sorting.prob.rdir),
            ncol(temp.nonFunctional.sorting.prob.rdir))
          temp.nonFunctional.sorting.prob.rdir[which(is.na(temp.nonFunctional.sorting.prob.rdir[,1])),] <- temp.fill.values
        }

        temp.functional.counts <- t(rmultinom(1, size = temp.functional.cells,
          prob = temp.functional.sorting.prob.rdir ))
        temp.nonFunctional.counts <- t(rmultinom(1, size = temp.nonFunctional.cells,
          prob = temp.nonFunctional.sorting.prob.rdir ))

        temp.guide.counts <- temp.functional.counts + temp.nonFunctional.counts
        repl.counts[unique.guide.overlaps[g],] <- temp.guide.counts
      }
    }

    # for every enhancer, obtain the enhancer strength from the beta distribution
    # for every guide, get the guide efficiency
    # guides get sorted either according to enhancer strength or as negative sorting
    for(k in 1:nrow(input.frame$enhancer)){
      temp.enhancer.ranges <- GRanges(seqnames = input.frame$enhancer$chrom[k],
        ranges = IRanges(input.frame$enhancer$start[k], input.frame$enhancer$end[k]))

      temp.overlaps <- as.data.frame(findOverlaps(guide.ranges, temp.enhancer.ranges, type = 'any'))

      if(nrow(temp.overlaps) > 0){
        unique.overlaps <- unique(temp.overlaps$queryHits)

        guide.efficiency <- rbeta(length(unique.overlaps), shape1 = input.frame$guideShape1, shape2 = input.frame$guideShape2)
        temp.enhancer.sort.prob <- c()
        for(e in 1:length(guide.efficiency)){
          temp.cells.with.guide <- input.frame$inputPools[[i]][unique.overlaps[e]]
          temp.guide.efficiency <- guide.efficiency[e]

          temp.functional.cells <- round(temp.cells.with.guide * temp.guide.efficiency)
          temp.nonFunctional.cells <- temp.cells.with.guide - temp.functional.cells

          if(input.frame$screenType == 'selectionScreen'){ #'selectionScreen' %in% names(input.frame)){
            temp.functional.sorting.prob <- input.frame$negSortingFrequency - effect.diff * sort.factor[k] # + input.frame$posSortingFrequency
            # temp.functional.sorting.prob <- sort.factor[k] * effect.diff + input.frame$posSortingFrequency
          } else {
            temp.functional.sorting.prob <- sort.factor[k] * effect.diff + input.frame$negSortingFrequency
          }
          # temp.functional.sorting.prob <- sort.factor[k] * exon.neg.diff + input.frame$negSortingFrequency

          temp.functional.sorting.prob.rdir <- rdirichlet(1, temp.functional.sorting.prob)
          temp.nonFunctional.sorting.prob.rdir <- rdirichlet(1, input.frame$negSortingFrequency)

          if(length(which(is.na(temp.functional.sorting.prob.rdir))) > 0){
            temp.nan.rows <- which(is.na(temp.functional.sorting.prob.rdir[,1]))
            temp.fill.values <- rep(1 / ncol(temp.functional.sorting.prob.rdir),
              ncol(temp.functional.sorting.prob.rdir))
            temp.functional.sorting.prob.rdir[which(is.na(temp.functional.sorting.prob.rdir[,1])),] <- temp.fill.values
          }

          if(length(which(is.na(temp.nonFunctional.sorting.prob.rdir))) > 0){
            temp.nan.rows <- which(is.na(temp.nonFunctional.sorting.prob.rdir[,1]))
            temp.fill.values <- rep(1 / ncol(temp.nonFunctional.sorting.prob.rdir),
              ncol(temp.nonFunctional.sorting.prob.rdir))
            temp.nonFunctional.sorting.prob.rdir[which(is.na(temp.nonFunctional.sorting.prob.rdir[,1])),] <- temp.fill.values
          }

          temp.functional.counts <- t(rmultinom(1, size = temp.functional.cells,
            prob = temp.functional.sorting.prob.rdir ))
          temp.nonFunctional.counts <- t(rmultinom(1, size = temp.nonFunctional.cells,
            prob = temp.nonFunctional.sorting.prob.rdir ))

          temp.guide.counts <- temp.functional.counts + temp.nonFunctional.counts
          repl.counts[unique.overlaps[e],] <- temp.guide.counts
        }
      }
    }
    repl.sequenced <- experiment_sequencing(cbind(input.frame$inputPools[[i]], repl.counts),
      input.frame$pcrDupl, input.frame$seqDepth[[i]])

    repl.sequenced.filtered <- repl.sequenced
    if(input.frame$screenType == 'selectionScreen'){ #'selectionScreen' %in% names(input.frame)){
      repl.sequenced.filtered <- repl.sequenced[,c(1,2)]
    }

    colnames(repl.sequenced.filtered) <- paste(paste0('sim', sim.nr, '_repl', i), input.frame$poolNames, sep = '_')
    out.sim.data <- cbind(out.sim.data, repl.sequenced.filtered)
  }
  return(list(counts = out.sim.data, enhancerStrength = sort.factor))
}

# give different genomic properties different sorting probabilities
# this method creates a compound distribution for obtinain the sorting probabilities
# deprecated
full_replicate_simulation <- function(input.frame, input.info, sim.nr){

  out.sim.data <- c()

  guide.ranges <- GRanges(seqnames = input.info$chrom,
    ranges = IRanges(input.info$start, input.info$end))

  exon.neg.diff <- input.frame$posSortingFrequency - input.frame$negSortingFrequency
  sort.factor <- rbeta(nrow(input.frame$enhancer), shape1 = input.frame$enhancerShape1, shape2 = input.frame$enhancerShape2)

  # for the number of replicates
  for(i in 1:length(input.frame$inputPools)){
    temp.sort.prob <- rdirichlet(length(input.frame$inputPools[[i]]),
      input.frame$negSortingFrequency)

    for(j in 1:nrow(input.frame$exon)){
      temp.exon.ranges <- GRanges(seqnames = input.frame$exon$chrom[j],
        ranges = IRanges(input.frame$exon$start[j], input.frame$exon$end[j]))

      temp.overlaps <- as.data.frame(findOverlaps(guide.ranges, temp.exon.ranges, type = 'any'))

      unique.guide.overlaps <- unique(temp.overlaps$queryHits)
      overlap.guide.efficiency <- rbeta(length(unique.guide.overlaps), shape1 = input.frame$guideShape1, shape2 = input.frame$guideShape2)

      temp.exon.sort.prob <- c()

      for(g in 1:length(unique.guide.overlaps)){
        temp.sorting.prob.guide <- overlap.guide.efficiency[g] * exon.neg.diff
        temp.sorting.prob <- temp.sorting.prob.guide + input.frame$negSortingFrequency
        guide.exon.sorting.prob <- rdirichlet(1,
          temp.sorting.prob)
        temp.exon.sort.prob <- rbind(temp.exon.sort.prob, guide.exon.sorting.prob)
      }

      temp.sort.prob[unique.guide.overlaps, ] <- temp.exon.sort.prob
    }

    # for every enhancer, obtain the forting frequency from the beta distribution
    # and reduce the initial sorting frequency by designated factor.
    # Note: with example frequencies of 39 and 13 for exon and neg, it does not
    # make sense to use the 39 as starting poin but instead the difference between pos and neg
    # therefore beta.out * (39 - 13)

    for(k in 1:nrow(input.frame$enhancer)){
      temp.enhancer.ranges <- GRanges(seqnames = input.frame$enhancer$chrom[k],
        ranges = IRanges(input.frame$enhancer$start[k], input.frame$enhancer$end[k]))

      temp.overlaps <- as.data.frame(findOverlaps(guide.ranges, temp.enhancer.ranges, type = 'any'))
      unique.overlaps <- unique(temp.overlaps$queryHits)

      guide.efficiency <- rbeta(length(unique.overlaps), shape1 = input.frame$guideShape1, shape2 = input.frame$guideShape2)
      temp.enhancer.sort.prob <- c()
      for(e in 1:length(guide.efficiency)){
        enhancer.sort.freq <- sort.factor[k] * guide.efficiency[e] * exon.neg.diff + input.frame$negSortingFrequency
        temp.enhancer.sort.prob.guide <- rdirichlet(1, enhancer.sort.freq)
        temp.enhancer.sort.prob <- rbind(temp.enhancer.sort.prob, temp.enhancer.sort.prob.guide)
      }

      temp.sort.prob[unique.overlaps, ] <- temp.enhancer.sort.prob
    }
    repl.counts <- counts_from_probs(temp.sort.prob, input.frame$inputPools[[i]])
    repl.sequenced <- experiment_sequencing(cbind(input.frame$inputPools[[i]], repl.counts),
      input.frame$pcrDupl, input.frame$seqDepth[[i]])
    colnames(repl.sequenced) <- paste(paste0('sim', sim.nr, '_repl', i), input.frame$poolNames, sep = '_')
    out.sim.data <- cbind(out.sim.data, repl.sequenced)
  }
  return(list(counts = out.sim.data, enhancerStrength = sort.factor))
}

generate_enhancers <- function(input.list, input.info, input.exon){
  min.start <- min(input.info$start)
  
  # Are enhancers allowed outside of this max end?
  max.end <- max(input.info$end) - 1.5 * input.list$enhancerSize

  potential.start.pos <- c(min.start:max.end)

  overlapping.exon <- which(potential.start.pos > (min(input.exon$start)- 1.5 * input.list$enhancerSize) &
    potential.start.pos < (max(input.exon$end) + 1.5 * input.list$enhancerSize) )

  # is this dropping potential enhancers that overlap with the positive controls?
  final.potential.start.pos <- potential.start.pos[ -overlapping.exon]

  start.pos <- sample(final.potential.start.pos, input.list$nrEnhancers)
  end.pos <- start.pos + input.list$enhancerSize

  # what does the input.info$chrom[1] do?
  out.enhancers <- data.frame(chrom = rep(input.info$chrom[1], input.list$nrEnhancers),
    start = start.pos, end = end.pos, stringsAsFactors = F)

  return(out.enhancers)
}

sim_enhancer_overlaps <- function(input.info, input.enh){

  info.ranges <- GRanges(seqnames = input.info$chrom,
    ranges = IRanges(input.info$start, input.info$end))

  enh.ranges <- GRanges(seqnames = input.enh$chrom,
    ranges = IRanges(input.enh$start, input.enh$end))

  overlaps <- as.data.frame(findOverlaps(info.ranges, enh.ranges, type = 'any'))

  return(unique(overlaps$queryHits))
}

create_ZINB_shape <- function(nrSim,in_eta,in_rate,in_dispersion){
  #given an eta, how many zeros will there be
  nr_fromZero <- rbinom(1, nrSim, in_eta)
  sample_zeros <- rep(0,nr_fromZero)

  nb_values <- rnbinom((nrSim-nr_fromZero),size = in_dispersion,mu = in_rate)

  out_zinb <- c(sample_zeros,nb_values)

  shuffled_out <- sample(out_zinb)

  return(shuffled_out)
}

estimate_ZINB_shape <- function(in_data,par){
  d <- par[1]
  rate <- par[2]
  eta <- par[3]

	in_row_counts <- unlist(in_data)

  non_zero_obs <- in_row_counts[which(in_row_counts > 0)]
  nr_zeros <- rep(0,(length(in_row_counts)-length(non_zero_obs)))

  zero_lik <- log(eta + (1-eta)*dnbinom(nr_zeros, size = d, mu = rate))
  non_zero_lik <- log(1-eta) + dnbinom(non_zero_obs, size = d, mu = rate, log = TRUE)

  lik_vec <- c(zero_lik,non_zero_lik)

  if(length(which(lik_vec == Inf)) > 0){
    print(paste('there are: ',length(which(lik_vec == Inf)),' Inf in lik_vec',par,sep = ''))
    print(paste('par/mu: ',par,sep = ''))
  }
  if(length(which(lik_vec == -Inf)) > 0){
    print(paste('there are: ',length(which(lik_vec == -Inf)),' -Inf in lik_vec',par,sep = ''))
    print(paste('par/mu: ',par,sep = ''))
  }

  result = -sum(lik_vec)
  result
}

obtain_ZINB_pars <- function(input.data){
  ZINB_optim_par_out <- optim(par=c(0.5,50,0.5),estimate_ZINB_shape,in_data = input.data,method = 'L-BFGS-B',lower = c(0.0001,0.0001,0.0001))
  return(list(eta = ZINB_optim_par_out$par[3], rate = ZINB_optim_par_out$par[2],
    dispersion = ZINB_optim_par_out$par[1]))
}

generate_tiling_deletion_positions <- function(input.list){

  start.positions <- seq(from = input.list$tilingStart, to = input.list$tilingEnd,
    by = input.list$stepSize)
  end.positions <- start.positions + input.list$tilingSize

  out.chrom <- rep(input.list$screenChrom, length(start.positions))

  out.df <- data.frame(chrom = out.chrom, start = start.positions,
    end = end.positions, stringsAsFactors = F)

  return(out.df)
}
