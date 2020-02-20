
suppressPackageStartupMessages(require(MCMCpack))
suppressPackageStartupMessages(require(transport))
suppressPackageStartupMessages(require(IRanges))
suppressPackageStartupMessages(require(GenomicRanges))

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
  if(nrow(out.list$inputGuideDistr) != nrow(out.list$guides)){
    if(nrow(out.list$inputGuideDistr) >= nrow(out.list$guides)){
      out.list$inputGuideDistr <- out.list$inputGuideDistr[sample(1:nrow(out.list$inputGuideDistr), size = nrow(out.list$guides)),]
    } else {
      out.list$inputGuideDistr <- out.list$inputGuideDistr[sample(1:nrow(out.list$inputGuideDistr), size = nrow(out.list$guides), replace = T),]
    }
  }
  if(! 'crisprEffectRange' %in% input.list.names){
    if(out.list$crisprSystem == 'CRISPRi' | out.list$crisprSystem == 'CRISPRa'){
      out.list$crisprEffectRange <- 1000
    } else if(out.list$crisprSystem == 'Cas9'){
      out.list$crisprEffectRange <- 20
    } else if(out.list$crisprSystem == 'dualCRISPR'){
      out.list$crisprEffectRange <- 0
    } else {
      print('Error: please specify a valid CRISPR system (Cas9, CRISPRi, CRISPRa, dualCRISPR)')
      break()
    }
  }
  if(! 'nrNeg' %in% input.list.names){
    out.list$nrNeg <- 300
  }
  if(! 'nrSims' %in% input.list.names){
    out.list$nrSims <- 1
  }
  if(! 'inputCellNr' %in% input.list.names){
    out.list$inputCellNr <- rep(50e6, ncol(out.list$inputGuideDistr))
  }
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
    } else{
      print('Error: please specify valid enhancerStrength (high, medium, low)')
      break()
    }
  }
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
    } else{
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
    temp.seq.depth <- out.list$seqDepth

    adj.seq.depth <- lapply(temp.seq.depth, function(x){
      return(c(x, x[length(x)]))
    })

    out.list$seqDepth <- adj.seq.depth
    if('selectionStrength' %in% input.list.names){
      if(out.list$selectionStrength == 'high'){
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
        if(out.list$selectionStrength == 'high'){
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

  return(out.list)
}

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

simulate_data <- function(input.list){

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
    
    browser()

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
  if(input.frame$screenType == 'selectionScreen'){ #'selectionScreen' %in% names(input.frame)){
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
  max.end <- max(input.info$end) - 1.5 * input.list$enhancerSize

  potential.start.pos <- c(min.start:max.end)

  overlapping.exon <- which(potential.start.pos > (min(input.exon$start)- 1.5 * input.list$enhancerSize) &
    potential.start.pos < (max(input.exon$end) + 1.5 * input.list$enhancerSize) )

  final.potential.start.pos <- potential.start.pos[ -overlapping.exon]

  start.pos <- sample(final.potential.start.pos, input.list$nrEnhancers)
  end.pos <- start.pos + input.list$enhancerSize

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
