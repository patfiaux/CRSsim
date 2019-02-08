
print('Loading packages...')
suppressPackageStartupMessages(require(IRanges))
suppressPackageStartupMessages(require(GenomicRanges))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(glmmTMB))

suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(pROC))
suppressPackageStartupMessages(require(MESS))
suppressPackageStartupMessages(require(edgeR))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(extraDistr))

#given step size, standard deviation around target site, number of targets, number of chromosomes
  #
  # generate a testing dataset consiting of both tiles and enhancers (two seperate files)
  # files contain format: chromosome  start end
#in_targets vector should be of equal length to in_chroms vector to specify the number of targets per chromosome
#in_enhancers: should be of equal length to in_chroms vector to specify the number of enhancers per chromosome
#targetName: file name for guide targets
#enhName: file name for enahcer file
#exampe:
  #simulate_sgPosition_data(15,10,c(2000,2000,2000),c(0,2,0),100,30,'testTargetPos_V2','testEnhancerPos_V2')
simulate_sgPosition_data <- function(in_stepSize,in_targetSd,in_targets,in_enhancers,in_enhancerLength,
  in_enhancerSD){ #,targetName,enhName){

  out_targetPos <- c()
  out_enhancerPos <- c()
  #for each chromosome
  for(j in 1:length(in_targets)){
    chrom_tileNr <- in_targets[j]
    chrom_enhancNr <- in_enhancers[j]
    chrom_out <- c()
    target_out <- c()
    end_out <- c()
    for(i in 1:chrom_tileNr){
      chrom_out <- c(chrom_out,paste('chr',j,sep = ''))
      target_out <- c(target_out,round(abs(rnorm(1,mean = i*in_stepSize,sd = in_targetSd))))
    }
    if(chrom_enhancNr > 0){
      enhancerLoc_options <- c(1:max(target_out))
      enhancer_start <- sample(enhancerLoc_options,chrom_enhancNr,replace = FALSE)
      enhancer_end <- enhancer_start + round(rnorm(length(enhancer_start),mean = in_enhancerLength, sd = in_enhancerSD))
      enh_chroms <- rep(paste('chr',j,sep = '') , chrom_enhancNr)
      temp_enhancer_out <- cbind.data.frame(enh_chroms,enhancer_start,enhancer_end,stringsAsFactors = FALSE)
      sorted_temp_enhancer_out <- temp_enhancer_out[order(temp_enhancer_out[,2]),]
      out_enhancerPos <- rbind.data.frame(out_enhancerPos,sorted_temp_enhancer_out)
    }
    #sorting
    temp_tilePos <- cbind.data.frame(chrom_out,target_out,stringsAsFactors = FALSE)
    sorted_temp_tilePos <- temp_tilePos[order(temp_tilePos[,2]),]
    out_targetPos <- rbind.data.frame(out_targetPos,sorted_temp_tilePos)
  }

  colnames(out_targetPos) <- c('chrom', 'start')
  out_targetPos$end  <- out_targetPos$start + 20
  colnames(out_enhancerPos) <- c('chrom', 'start', 'end')

  return(list(guide_targets = out_targetPos, enhancer_pos = out_enhancerPos))
  # save_a_df_as_CSV(targetName,out_targetPos)
  # save_a_df_as_CSV(enhName,out_enhancerPos)
}

#given tile deletion length, standard deviation around tile position, step size, number of tiles, number of chromosomes
  #
  # generate a testing dataset consiting of both tiles and enhancers (two seperate files)
  # files contain format: chromosome  start end
#in_tiles vector should be of equal length to in_chroms vector to specify the number of tiles per chromosome
#in_enhancers: should be of equal length to in_chroms vector to specify the number of enhancers per chromosome
#targetName: file name for guide targets
#enhName: file name for enahcer file
#exampe:
  #simulate_tilePosition_data(100,20,40,c(20,20,20),c(0,2,0),40,5,'testTilePos','testEnhancerPos')
  #simulate_tilePosition_data(1000,75,50,c(6000),c(5),500,300,'TilePos','EnhancerPos')
simulate_tilePosition_data <- function(in_tileLength,in_tileSd,in_tileStepSize,in_tiles,in_enhancers,in_enhancerLength,
  in_enhancerSD){ #,targetName,enhName){

  out_targetPos <- c()
  out_enhancerPos <- c()
  #for each chromosome
  for(j in 1:length(in_tiles)){
    chrom_tileNr <- in_tiles[j]
    chrom_enhancNr <- in_enhancers[j]
    chrom_out <- c()
    start_out <- c()
    end_out <- c()
    for(i in 1:chrom_tileNr){
      chrom_out <- c(chrom_out,paste('chr',j,sep = ''))
      start_out <- c(start_out,round(abs(rnorm(1,mean = i*in_tileStepSize,sd = in_tileSd))))
      end_out <- c(end_out, round(abs(rnorm(1,mean = ((i*in_tileStepSize)+in_tileLength),sd = in_tileSd))))
    }
    if(chrom_enhancNr > 0){
      enhancerLoc_options <- c(1:max(end_out))
      enhancer_start <- sample(enhancerLoc_options,chrom_enhancNr,replace = FALSE)
      enhancer_end <- enhancer_start + round(rnorm(length(enhancer_start),mean = in_enhancerLength, sd = in_enhancerSD))
      enh_chroms <- rep(paste('chr',j,sep = '') , chrom_enhancNr)
      temp_enhancer_out <- cbind.data.frame(enh_chroms,enhancer_start,enhancer_end,stringsAsFactors = FALSE)
      sorted_temp_enhancer_out <- temp_enhancer_out[order(temp_enhancer_out[,2]),]
      out_enhancerPos <- rbind.data.frame(out_enhancerPos,sorted_temp_enhancer_out)
    }
    #sorting
    temp_tilePos <- cbind.data.frame(chrom_out,start_out,end_out,stringsAsFactors = FALSE)
    sorted_temp_tilePos <- temp_tilePos[order(temp_tilePos[,2]),]
    out_targetPos <- rbind.data.frame(out_targetPos,sorted_temp_tilePos)
  }

  colnames(out_targetPos) <- c('chrom', 'start', 'end')
  colnames(out_enhancerPos) <- c('chrom', 'start', 'end')
  return(list(guide_targets = out_targetPos, enhancer_pos = out_enhancerPos))
  # save_a_df_as_CSV(targetName,out_targetPos)
  # save_a_df_as_CSV(enhName,out_enhancerPos)
}

# output the analysis evaluation results in csv format for comparison of different data treatments
analysis_eval_out <- function(spec.file, score.file = NULL, label.file = NULL, data.dir = NULL){

    # running regular analysis:
    # read in specification file
    if(is.null(spec.file)) stop('Error in "analyze_data()": specification file has to be provided')
    analysis.specs <- read_analysis_specs(spec.file, data.dir)
    analysis.specs.names <- names(analysis.specs)

    # if no score file is provided, scores have to be calculated
    final.score.list <- list()  #instantiate the final score list (specifiy format?)

    if(is.null(score.file)){
      final.score.list <- calculate_scores(analysis.specs, label.file = label.file)  # needs to be implemented
    } else {
      final.score.list <- obtain_computed_scores(score.file)  # needs to be implemented
    }

    if('evaluatePerformance' %in% analysis.specs.names){
      performanceEvaluation_recording(final.score.list, analysis.specs)  # needs to be implemented
    }

    score_plotting(final.score.list, analysis.specs)  # needs to be implemented

}

# plotting of both AUC and prAUC of the scores based on labels
# input:
#   input.score.list: list of scores in score output format (rawScores, formatScores, log2_rate_ratio, chromosome, labels, start, end)
#   input.specs: list of specifications
# output: AUC and prAUC plot of scores
#   writing of evaluation df.
#     pre_bin_type, pre_bin_size, pre_bin_step, pre_bin_methd,
#     post_bin_methd, post_bin_size, Method, AUC, prAUC
performanceEvaluation_recording <- function(input.score.list, input.specs){
  pre_bin_type <- c()
  pre_bin_size <- c()
  pre_bin_step <- c()
  pre_bin_methd <- c()
  post_bin_methd <- c()
  post_bin_size	<- c()
  Method <- c()
  AUC <- c()
  prAUC <- c()

  methd.auc <- list()
  methd.prAuc <- list()

  score.names <- names(input.score.list)
  for(i in 1:length(score.names)){
    temp.score.df <- input.score.list[[i]]
    temp.score.name <- score.names[i]  #NB, viterbi;
    print(paste('Evaluate per guide performance for:', temp.score.name, sep = ' '))
    temp.pos.labels <- input.specs$positiveLabels
    temp.neg.labels <- input.specs$negativeLabels
    temp.labels <- temp.score.df$label
    temp.scores <- temp.score.df$formatScores
    if(temp.score.name == 'viterbi'){
      temp.viterbi.df <- data.frame(formatScores = temp.scores,
        label = temp.score.df$origLabel, viterbi_label = temp.score.df$viterbiLabels, stringsAsFactors = 'F')

      # assumption: enhancers are encoded as '1', non-enhancer as '2'
      ordered.temp.viterbi.df <- temp.viterbi.df[with(temp.viterbi.df, order(viterbi_label, -formatScores)),]
      ordered.temp.viterbi.df$newScores <- c(nrow(ordered.temp.viterbi.df):1)

      temp.scores <- ordered.temp.viterbi.df$newScores
      temp.labels <- ordered.temp.viterbi.df$label
    }
    if('simulated_data' %in% names(input.specs)){
      all.labels <- unique(temp.labels)
      temp.neg.labels <- all.labels[-which(all.labels %in% c(temp.pos.labels, 'exon') )]
    }
    auc.list <- AUC_x_y(temp.scores, temp.labels, temp.pos.labels, temp.neg.labels)
    prAUC.list <- precisionRecall_AUC(temp.scores, temp.labels, temp.pos.labels, temp.neg.labels)
    methd.auc[[temp.score.name]] <- auc.list
    methd.prAuc[[temp.score.name]] <- prAUC.list

    method.name <- temp.score.name  #strsplit(temp.score.name,'_')[[1]][1]
    Method <- c(Method, method.name)

    if(input.specs$dataArranging == 'none'){
      pre_bin_type <- c(pre_bin_type, 'none')
      pre_bin_size <- c(pre_bin_size, NA)
      pre_bin_step <- c(pre_bin_step, NA)
      pre_bin_methd <- c(pre_bin_methd, NA)
    } else {
      if(input.specs$binBy == 'mean'){
        pre_bin_methd <- c(pre_bin_methd, 'mean')
      } else {
        pre_bin_methd <- c(pre_bin_methd, 'sum')
      }
      if('preBinning' %in% names(input.specs)){
        pre_bin_type <- c(pre_bin_type,'nonOverl')
        pre_bin_size <- c(pre_bin_size, input.specs$preBinning)
        pre_bin_step <- c(pre_bin_step, NA)
      } else {
        pre_bin_type <- c(pre_bin_type,'overl')
        pre_bin_size <- c(pre_bin_size, input.specs$preOverlapBinning[2])
        pre_bin_step <- c(pre_bin_step, input.specs$preOverlapBinning[1])
      }
    }

    if('postScoringAnalysis' %in% names(input.specs)){
      temp.method.name <- strsplit(temp.score.name,'_')[[1]]
      if(length(temp.method.name) > 1){
        post.score.method <- temp.method.name[2]
        if(post.score.method == 'alphaRRA'){
          post_bin_methd <- c(post_bin_methd, 'alphaRRA')
          post_bin_size	<- c(post_bin_size, input.specs$binSize)
        } else if(post.score.method == 'CRESTrra'){
          post_bin_methd <- c(post_bin_methd, 'CRESTrra')
          post_bin_size	<- c(post_bin_size, input.specs$binSize)
        } else if(post.score.method == 'guideWindowRRA'){
          post_bin_methd <- c(post_bin_methd, 'guideWindowRRA')
          post_bin_size	<- c(post_bin_size, input.specs$guidePerFixedWindow)
        } else if(post.score.method == 'slidingWindow'){
          post_bin_methd <- c(post_bin_methd, 'slidingWindow')
          post_bin_size	<- c(post_bin_size, input.specs$guidePerSlidingWindow)
        } else if(post.score.method == 'genomeScores'){
          post_bin_methd <- c(post_bin_methd, 'genomeScores')
          post_bin_size	<- c(post_bin_size, round(mean(temp.score.df$end - temp.score.df$start)))
        } else {
          print('WARNING! No correct post scoring method found!')
        }
      } else {
        post_bin_methd <- c(post_bin_methd, 'none')
        post_bin_size	<- c(post_bin_size, NA)
      }
    } else {
      post_bin_methd <- c(post_bin_methd, 'none')
      post_bin_size	<- c(post_bin_size, NA)
    }
    AUC <- c(AUC, auc.list$num_AUC)
    prAUC <- c(prAUC, round(prAUC.list$num_AUC, 3))
  }
  plot_methodEval(methd.auc,'AUCeval',input.specs$dataName,'bottomright','Fraction True negative included','Fraction True positives included')
  plot_methodEval(methd.prAuc,'prAUCeval',input.specs$dataName,'bottomleft','Recall','Precision')

  out.df <- cbind.data.frame(pre_bin_type, pre_bin_size, pre_bin_step, pre_bin_methd,
    post_bin_methd, post_bin_size, Method, AUC, prAUC, stringsAsFactors = FALSE)
  write.table(out.df, file = 'perGuide_method_eval.csv', sep = ',', row.names = FALSE, col.names = FALSE)

  return(cbind.data.frame(Method, prAUC, stringsAsFactors = F))

}

# given a csv file containing the positive regions, calculate the top score
# for each region, corresponding to a specific bin size
# out: region_scores , region_labels, region_start, region_end
pos_region_scores <- function(input.score.df, input.specs){

  pos.regions <- read.csv(input.specs$pos_regions, header = T, stringsAsFactors = F)

  # get max for each enhancer
  pos.range <- GRanges(seqnames = pos.regions$chrom,
        ranges = IRanges(pos.regions$start, pos.regions$end))
  guide.range <- GRanges(seqnames = input.score.df$chrom,
        ranges = IRanges(input.score.df$start, input.score.df$end))

  enhancer.overlaps <- as.data.frame(findOverlaps(pos.range, guide.range, type = 'any'))
  enhancer.score.list <- split(enhancer.overlaps, enhancer.overlaps$queryHits)

  enhancer.scores <- unlist(lapply(enhancer.score.list, function(x){
      temp.score <- max(input.score.df$formatScores[x$subjectHits])
      temp.score
      }))

  enhancer.label <- rep('pos', length(enhancer.scores))
  enhancer.info <- pos.regions[unique(enhancer.overlaps$queryHits),]

  # remove tiles overlaping enhancers
  reduced.score.df <- input.score.df[-enhancer.overlaps$subjectHits,]

  total.scores <- c(enhancer.scores, reduced.score.df$formatScores)
  total.labels <- c(enhancer.label, reduced.score.df$label)
  final.start <- c(enhancer.info$start, reduced.score.df$start)

  out.df <- data.frame(region_scores = total.scores, region_labels = total.labels,
    region_start = final.start, stringsAsFactors = F)

  return(out.df)
}

# plotting of both AUC and prAUC of the scores based on labels
# input:
#   input.score.list: list of scores in score output format (rawScores, formatScores, log2_rate_ratio, chromosome, labels, start, end)
#   input.specs: list of specifications
# output: AUC and prAUC plot of scores
#   writing of evaluation df.
#     pre_bin_type, pre_bin_size, pre_bin_step, pre_bin_methd,
#     post_bin_methd, post_bin_size, Method, AUC, prAUC
performanceEvaluation_perRegion_recording <- function(input.score.list, input.specs){
  pre_bin_type <- c()
  pre_bin_size <- c()
  pre_bin_step <- c()
  pre_bin_methd <- c()
  post_bin_methd <- c()
  post_bin_size	<- c()
  Method <- c()
  AUC <- c()
  prAUC <- c()

  methd.auc <- list()
  methd.prAuc <- list()

  score.names <- names(input.score.list)
  for(i in 1:length(score.names)){
    init.score.df <- input.score.list[[i]]
    # ensure that only real genomic regions are being used
    temp.score.df <- c()
    if(length(which(is.na(init.score.df$start))) > 0){
      temp.score.df <- init.score.df[-which(is.na(init.score.df$start)),]
    } else {
      temp.score.df <- init.score.df
    }
    temp.score.name <- score.names[i]  #NB, viterbi;
    print(paste('Evaluate per region performance for:', temp.score.name, sep = ' '))
    temp.pos.labels <- input.specs$positiveLabels
    temp.neg.labels <- input.specs$negativeLabels

    if(temp.score.name == 'viterbi'){
      temp.viterbi.df <- data.frame(origformatScores = temp.score.df$formatScores,
        label = temp.score.df$origLabel, viterbi_label = temp.score.df$viterbiLabels,
        start = temp.score.df$start, end = temp.score.df$end, chrom = temp.score.df$chrom, stringsAsFactors = 'F')

      # assumption: enhancers are encoded as '1', non-enhancer as '2'
      ordered.temp.viterbi.df <- temp.viterbi.df[with(temp.viterbi.df, order(viterbi_label, -origformatScores)),]
      ordered.temp.viterbi.df$formatScores <- c(nrow(ordered.temp.viterbi.df):1)

      # temp.scores <- ordered.temp.viterbi.df$newScores
      # temp.labels <- ordered.temp.viterbi.df$label

      temp.score.df <- ordered.temp.viterbi.df
    }

    if('simulated_data' %in% names(input.specs)){
      temp.score.df <- temp.score.df[which(temp.score.df$label %in% c(temp.pos.labels, temp.neg.labels) ),]
    }

    if(length(unique(temp.score.df$chrom)) > 1){
      print('Warning: per region AUC can currently only be calculated for one chromosome')
      break
    }

    region.scores <- c()
    sub.score.names <- strsplit(temp.score.name, '_')[[1]]
    if(temp.score.name %in% c('nbGlmmUntr', 'dirichletMultinom')){
      region.scores <- data.frame(region_scores = temp.score.df$formatScores,
        region_labels = temp.score.df$label, region_start = temp.score.df$start,
        stringsAsFactors = F)
    } else if(length(sub.score.names) > 1 & sub.score.names[2] == 'genomeScores'){
      region.scores <- data.frame(region_scores = temp.score.df$formatScores,
        region_labels = temp.score.df$label, region_start = temp.score.df$start,
        stringsAsFactors = F)
    } else {
      # take max score for each region and determine how well method did in detecting all positive regions
      region.scores <- pos_region_scores(temp.score.df, input.specs)
    }


    auc.list <- AUC_x_y_perRegion(region.scores$region_scores, region.scores$region_labels,
      temp.pos.labels, temp.neg.labels)
    prAUC.list <- precisionRecall_AUC_perRegion(region.scores$region_scores, region.scores$region_labels,
      temp.pos.labels, temp.neg.labels)
    methd.auc[[temp.score.name]] <- auc.list
    methd.prAuc[[temp.score.name]] <- prAUC.list

    method.name <- temp.score.name  #strsplit(temp.score.name,'_')[[1]][1]
    Method <- c(Method, method.name)

    if(input.specs$dataArranging == 'none'){
      pre_bin_type <- c(pre_bin_type, 'none')
      pre_bin_size <- c(pre_bin_size, NA)
      pre_bin_step <- c(pre_bin_step, NA)
      pre_bin_methd <- c(pre_bin_methd, NA)
    } else {
      if(input.specs$binBy == 'mean'){
        pre_bin_methd <- c(pre_bin_methd, 'mean')
      } else {
        pre_bin_methd <- c(pre_bin_methd, 'sum')
      }
      if('preBinning' %in% names(input.specs)){
        pre_bin_type <- c(pre_bin_type,'nonOverl')
        pre_bin_size <- c(pre_bin_size, input.specs$preBinning)
        pre_bin_step <- c(pre_bin_step, NA)
      } else {
        pre_bin_type <- c(pre_bin_type,'overl')
        pre_bin_size <- c(pre_bin_size, input.specs$preOverlapBinning[2])
        pre_bin_step <- c(pre_bin_step, input.specs$preOverlapBinning[1])
      }
    }

    if('postScoringAnalysis' %in% names(input.specs)){
      temp.method.name <- strsplit(temp.score.name,'_')[[1]]
      if(length(temp.method.name) > 1){
        post.score.method <- temp.method.name[2]
        if(post.score.method == 'alphaRRA'){
          post_bin_methd <- c(post_bin_methd, 'alphaRRA')
          post_bin_size	<- c(post_bin_size, input.specs$binSize)
        } else if(post.score.method == 'CRESTrra'){
          post_bin_methd <- c(post_bin_methd, 'CRESTrra')
          post_bin_size	<- c(post_bin_size, input.specs$binSize)
        } else if(post.score.method == 'guideWindowRRA'){
          post_bin_methd <- c(post_bin_methd, 'guideWindowRRA')
          post_bin_size	<- c(post_bin_size, input.specs$guidePerFixedWindow)
        } else if(post.score.method == 'slidingWindow'){
          post_bin_methd <- c(post_bin_methd, 'slidingWindow')
          post_bin_size	<- c(post_bin_size, input.specs$guidePerSlidingWindow)
        } else if(post.score.method == 'genomeScores'){
          post_bin_methd <- c(post_bin_methd, 'genomeScores')
          post_bin_size	<- c(post_bin_size, round(mean(temp.score.df$end - temp.score.df$start)))
        } else {
          print('WARNING! No correct post scoring method found!')
        }
      } else {
        post_bin_methd <- c(post_bin_methd, 'none')
        post_bin_size	<- c(post_bin_size, NA)
      }
    } else {
      post_bin_methd <- c(post_bin_methd, 'none')
      post_bin_size	<- c(post_bin_size, NA)
    }
    AUC <- c(AUC, auc.list$num_AUC)
    prAUC <- c(prAUC, round(prAUC.list$num_AUC, 3))
  }
  plot_methodEval(methd.auc,'perRegion_AUCeval',input.specs$dataName,'bottomright','Fraction True negative included','Fraction True positives included')
  plot_methodEval(methd.prAuc,'perRegion_prAUCeval',input.specs$dataName,'bottomleft','Recall','Precision')

  out.df <- cbind.data.frame(pre_bin_type, pre_bin_size, pre_bin_step, pre_bin_methd,
    post_bin_methd, post_bin_size, Method, AUC, prAUC, stringsAsFactors = FALSE)
  write.table(out.df, file = 'perRegion_method_eval.csv', sep = ',', row.names = FALSE, col.names = FALSE)

  return(cbind.data.frame(Method, prAUC, stringsAsFactors = F))

}

# given a csv file containing the positive regions, calculate the top score
# for each region, corresponding to a specific bin size
# out: region_scores , region_labels, region_start, region_end
pos_element_scores <- function(input.score.df, input.specs){

  pos.regions <- read.csv(input.specs$pos_regions, header = T, stringsAsFactors = F)

  guide.effect.size <- c()
  if(input.specs$crisprSystem == 'dualCRISPR'){
    guide.effect.size <- input.specs$deletionSize
  } else {
    guide.effect.size <- input.specs$crisprEffectRange
  }

  # get max for each enhancer
  pos.range <- GRanges(seqnames = pos.regions$chrom,
        ranges = IRanges(pos.regions$start - guide.effect.size, pos.regions$end + guide.effect.size))
  guide.range <- GRanges(seqnames = input.score.df$chrom,
        ranges = IRanges(input.score.df$start, input.score.df$end))

  enhancer.overlaps <- as.data.frame(findOverlaps(pos.range, guide.range, type = 'any'))
  enhancer.score.list <- split(enhancer.overlaps, enhancer.overlaps$queryHits)

  enhancer.scores <- unlist(lapply(enhancer.score.list, function(x){
      temp.score <- max(input.score.df$formatScores[x$subjectHits])
      temp.score
      }))

  enhancer.label <- rep('pos', length(enhancer.scores))
  enhancer.info <- pos.regions[unique(enhancer.overlaps$queryHits),]

  # remove tiles overlaping enhancers
  reduced.score.df <- input.score.df[-enhancer.overlaps$subjectHits,]
  if('pos' %in% unique(reduced.score.df$label)){
    reduced.score.df <- reduced.score.df[-which(reduced.score.df$label == 'pos'),]
  }

  reduced.start.pos <- seq(from = min(reduced.score.df$start), to = max(reduced.score.df$end),
    by = guide.effect.size)
  reduced.end.pos <- reduced.start.pos + guide.effect.size - 1

  reduced.window.df <- data.frame(chrom = rep(reduced.score.df$chrom[1], length(reduced.start.pos)),
    start = reduced.start.pos, end = reduced.end.pos, stringsAsFactors = F)

  reduced.window.guides <- GRanges(seqnames = reduced.window.df$chrom,
    ranges = IRanges(reduced.window.df$start, reduced.window.df$end))
  reduced.score.ranges <- GRanges(seqnames = reduced.score.df$chrom,
    ranges = IRanges(reduced.score.df$start, reduced.score.df$end))

  reduced.overlaps <- as.data.frame(findOverlaps(reduced.window.guides,
    reduced.score.ranges, type = 'any'))

  reduced.score.list <- split(reduced.overlaps, reduced.overlaps$queryHits)

  reduced.window.scores <- unlist(lapply(reduced.score.list, function(x){
    temp.max.score <- max(reduced.score.df$formatScores[x$subjectHits])
    # temp.labels.present <- which(input.specs$labelHierarchy %in% temp.label)
    # out.label <- input.specs$labelHierarchy[max(temp.labels.present)]
    return(temp.max.score)
    }))

  final.reduced.window.df <- reduced.window.df[unique(reduced.overlaps$queryHits), ]
  final.reduced.window.df$formatScores <- reduced.window.scores
  final.reduced.window.df$label <- rep(reduced.score.df$label[1], nrow(final.reduced.window.df))

  total.scores <- c(enhancer.scores, final.reduced.window.df$formatScores)
  total.labels <- c(enhancer.label, final.reduced.window.df$label)
  final.start <- c(enhancer.info$start, final.reduced.window.df$start)

  out.df <- data.frame(region_scores = total.scores, region_labels = total.labels,
    region_start = final.start, stringsAsFactors = F)

  return(out.df)
}

# plotting of both AUC and prAUC of the scores based on labels
# input:
#   input.score.list: list of scores in score output format (rawScores, formatScores, log2_rate_ratio, chromosome, labels, start, end)
#   input.specs: list of specifications
# output: AUC and prAUC plot of scores
#   writing of evaluation df.
#     pre_bin_type, pre_bin_size, pre_bin_step, pre_bin_methd,
#     post_bin_methd, post_bin_size, Method, AUC, prAUC
performanceEvaluation_perElement_recording <- function(input.score.list, input.specs){
  pre_bin_type <- c()
  pre_bin_size <- c()
  pre_bin_step <- c()
  pre_bin_methd <- c()
  post_bin_methd <- c()
  post_bin_size	<- c()
  Method <- c()
  AUC <- c()
  prAUC <- c()

  methd.auc <- list()
  methd.prAuc <- list()

  score.names <- names(input.score.list)
  for(i in 1:length(score.names)){
    init.score.df <- input.score.list[[i]]
    # ensure that only real genomic regions are being used
    temp.score.df <- c()
    if(length(which(is.na(init.score.df$start))) > 0){
      temp.score.df <- init.score.df[-which(is.na(init.score.df$start)),]
    } else {
      temp.score.df <- init.score.df
    }
    temp.score.name <- score.names[i]  #NB, viterbi;
    print(paste('Evaluate per element performance for:', temp.score.name, sep = ' '))
    temp.pos.labels <- input.specs$positiveLabels
    temp.neg.labels <- input.specs$negativeLabels

    if(temp.score.name == 'viterbi'){
      temp.viterbi.df <- data.frame(origformatScores = temp.score.df$formatScores,
        label = temp.score.df$origLabel, viterbi_label = temp.score.df$viterbiLabels,
        start = temp.score.df$start, end = temp.score.df$end, chrom = temp.score.df$chrom, stringsAsFactors = 'F')

      # assumption: enhancers are encoded as '1', non-enhancer as '2'
      ordered.temp.viterbi.df <- temp.viterbi.df[with(temp.viterbi.df, order(viterbi_label, -origformatScores)),]
      ordered.temp.viterbi.df$formatScores <- c(nrow(ordered.temp.viterbi.df):1)

      temp.score.df <- ordered.temp.viterbi.df
    }

    if('simulated_data' %in% names(input.specs)){
      temp.score.df <- temp.score.df[which(temp.score.df$label %in% c(temp.pos.labels, temp.neg.labels) ),]
    }

    if(length(unique(temp.score.df$chrom)) > 1){
      print('Warning: per element AUC can currently only be calculated for one chromosome')
      break
    }

    # region.scores <- c()
    # sub.score.names <- strsplit(temp.score.name, '_')[[1]]
    # if(temp.score.name %in% c('nbGlmmUntr', 'dirichletMultinom')){
    #   region.scores <- data.frame(region_scores = temp.score.df$formatScores,
    #     region_labels = temp.score.df$label, region_start = temp.score.df$start,
    #     stringsAsFactors = F)
    # } else if(length(sub.score.names) > 1 & sub.score.names[2] == 'genomeScores'){
    #   region.scores <- data.frame(region_scores = temp.score.df$formatScores,
    #     region_labels = temp.score.df$label, region_start = temp.score.df$start,
    #     stringsAsFactors = F)
    # } else {
    #   # take max score for each region and determine how well method did in detecting all positive regions
    #   region.scores <- pos_element_scores(temp.score.df, input.specs)
    # }
    region.scores <- pos_element_scores(temp.score.df, input.specs)

    auc.list <- AUC_x_y_perRegion(region.scores$region_scores, region.scores$region_labels,
      temp.pos.labels, temp.neg.labels)
    prAUC.list <- precisionRecall_AUC_perRegion(region.scores$region_scores, region.scores$region_labels,
      temp.pos.labels, temp.neg.labels)
    methd.auc[[temp.score.name]] <- auc.list
    methd.prAuc[[temp.score.name]] <- prAUC.list

    method.name <- temp.score.name  #strsplit(temp.score.name,'_')[[1]][1]
    Method <- c(Method, method.name)

    if(input.specs$dataArranging == 'none'){
      pre_bin_type <- c(pre_bin_type, 'none')
      pre_bin_size <- c(pre_bin_size, NA)
      pre_bin_step <- c(pre_bin_step, NA)
      pre_bin_methd <- c(pre_bin_methd, NA)
    } else {
      if(input.specs$binBy == 'mean'){
        pre_bin_methd <- c(pre_bin_methd, 'mean')
      } else {
        pre_bin_methd <- c(pre_bin_methd, 'sum')
      }
      if('preBinning' %in% names(input.specs)){
        pre_bin_type <- c(pre_bin_type,'nonOverl')
        pre_bin_size <- c(pre_bin_size, input.specs$preBinning)
        pre_bin_step <- c(pre_bin_step, NA)
      } else {
        pre_bin_type <- c(pre_bin_type,'overl')
        pre_bin_size <- c(pre_bin_size, input.specs$preOverlapBinning[2])
        pre_bin_step <- c(pre_bin_step, input.specs$preOverlapBinning[1])
      }
    }

    if('postScoringAnalysis' %in% names(input.specs)){
      temp.method.name <- strsplit(temp.score.name,'_')[[1]]
      if(length(temp.method.name) > 1){
        post.score.method <- temp.method.name[2]
        if(post.score.method == 'alphaRRA'){
          post_bin_methd <- c(post_bin_methd, 'alphaRRA')
          post_bin_size	<- c(post_bin_size, input.specs$binSize)
        } else if(post.score.method == 'CRESTrra'){
          post_bin_methd <- c(post_bin_methd, 'CRESTrra')
          post_bin_size	<- c(post_bin_size, input.specs$binSize)
        } else if(post.score.method == 'guideWindowRRA'){
          post_bin_methd <- c(post_bin_methd, 'guideWindowRRA')
          post_bin_size	<- c(post_bin_size, input.specs$guidePerFixedWindow)
        } else if(post.score.method == 'slidingWindow'){
          post_bin_methd <- c(post_bin_methd, 'slidingWindow')
          post_bin_size	<- c(post_bin_size, input.specs$guidePerSlidingWindow)
        } else if(post.score.method == 'genomeScores'){
          post_bin_methd <- c(post_bin_methd, 'genomeScores')
          post_bin_size	<- c(post_bin_size, round(mean(temp.score.df$end - temp.score.df$start)))
        } else {
          print('WARNING! No correct post scoring method found!')
        }
      } else {
        post_bin_methd <- c(post_bin_methd, 'none')
        post_bin_size	<- c(post_bin_size, NA)
      }
    } else {
      post_bin_methd <- c(post_bin_methd, 'none')
      post_bin_size	<- c(post_bin_size, NA)
    }
    AUC <- c(AUC, auc.list$num_AUC)
    prAUC <- c(prAUC, round(prAUC.list$num_AUC, 3))
  }
  plot_methodEval(methd.auc,'perElement_AUCeval',input.specs$dataName,'bottomright','Fraction True negative included','Fraction True positives included')
  plot_methodEval(methd.prAuc,'perElement_prAUCeval',input.specs$dataName,'bottomleft','Recall','Precision')

  out.df <- cbind.data.frame(pre_bin_type, pre_bin_size, pre_bin_step, pre_bin_methd,
    post_bin_methd, post_bin_size, Method, AUC, prAUC, stringsAsFactors = FALSE)
  write.table(out.df, file = 'perElement_method_eval.csv', sep = ',', row.names = FALSE, col.names = FALSE)

  return(cbind.data.frame(Method, prAUC, stringsAsFactors = F))

}

# Old version: performed binning of negative region. However, This does now accurately reflect
# the performance of the methods in my opinion. New method is looking at per guide score
# in negative regions
# given a csv file containing the positive regions, calculate the top score
# for each region, corresponding to a specific bin size
# out: region_scores , region_labels, region_start, region_end
pos_region_scores_old <- function(input.score.df, input.specs){

  pos.regions <- read.csv(input.specs$pos_regions, header = T, stringsAsFactors = F)
  # pos.regions <- input.specs$pos_regions

  total.region.start <- min(input.score.df$start)
  total.region.end <- max(input.score.df$end)

  # get max for each enhancer
  pos.range <- GRanges(seqnames = pos.regions$chrom,
        ranges = IRanges(pos.regions$start, pos.regions$end))
  guide.range <- GRanges(seqnames = input.score.df$chrom,
        ranges = IRanges(input.score.df$start, input.score.df$end))

  enhancer.overlaps <- as.data.frame(findOverlaps(pos.range, guide.range, type = 'any'))
  enhancer.score.list <- split(enhancer.overlaps, enhancer.overlaps$queryHits)

  enhancer.scores <- unlist(lapply(enhancer.score.list, function(x){
      temp.score <- max(input.score.df$formatScores[x$subjectHits])
      temp.score
      }))

  enhancer.label <- rep('pos', length(enhancer.scores))
  enhancer.info <- pos.regions[unique(enhancer.overlaps$queryHits),]

  # remove tiles overlaping enhancers
  reduced.score.df <- input.score.df[-enhancer.overlaps$subjectHits,]
  reduced.guide.range <- GRanges(seqnames = reduced.score.df$chrom,
        ranges = IRanges(reduced.score.df$start, reduced.score.df$end))

  binned.region <- seq(total.region.start, total.region.end, by = input.specs$pos_region_size)
  binned.range <- GRanges(seqnames = rep(pos.regions$chrom[1], length(binned.region)),
      ranges = IRanges(binned.region, binned.region + input.specs$pos_region_size - 1))

  bin.enhancer.overlaps <- as.data.frame(findOverlaps(binned.range, pos.range, type = 'any'))

  reduced.binned.region <- binned.region[-unique(bin.enhancer.overlaps$queryHits)]
  reduced.binned.range <- GRanges(seqnames = rep(pos.regions$chrom[1], length(reduced.binned.region)),
      ranges = IRanges(reduced.binned.region, reduced.binned.region + input.specs$pos_region_size - 1))

  guide.overlaps <- as.data.frame(findOverlaps(reduced.binned.range, reduced.guide.range, type = 'any'))

  bin.score.list <- split(guide.overlaps, guide.overlaps$queryHits)

  bin.scores <- unlist(lapply(bin.score.list, function(x){
    temp.score <- max(reduced.score.df$formatScores[x$subjectHits])
    temp.score
    }))

  binned.region.start <- reduced.binned.region[unique(guide.overlaps$queryHits)]

  bin.labels <- rep('neg', length(bin.scores))

  total.scores <- c(enhancer.scores, bin.scores)
  total.labels <- c(enhancer.label, bin.labels)
  final.start <- c(enhancer.info$start, binned.region.start)

  # out.df <- data.frame(region_scores = bin.scores, region_labels = region.labels,
  #   region_start = bin.coverage.start,
  #   stringsAsFactors = F)

  out.df <- data.frame(region_scores = total.scores, region_labels = total.labels,
    region_start = final.start, stringsAsFactors = F)

  return(out.df)
}

# Older version: in case of simulations is removing exons as opposed to keeping only
# regions which are of interest (positives and negatives)
# plotting of both AUC and prAUC of the scores based on labels
# input:
#   input.score.list: list of scores in score output format (rawScores, formatScores, log2_rate_ratio, chromosome, labels, start, end)
#   input.specs: list of specifications
# output: AUC and prAUC plot of scores
#   writing of evaluation df.
#     pre_bin_type, pre_bin_size, pre_bin_step, pre_bin_methd,
#     post_bin_methd, post_bin_size, Method, AUC, prAUC
performanceEvaluation_perRegion_recording_old <- function(input.score.list, input.specs){
  pre_bin_type <- c()
  pre_bin_size <- c()
  pre_bin_step <- c()
  pre_bin_methd <- c()
  post_bin_methd <- c()
  post_bin_size	<- c()
  Method <- c()
  AUC <- c()
  prAUC <- c()

  methd.auc <- list()
  methd.prAuc <- list()

  score.names <- names(input.score.list)
  for(i in 1:length(score.names)){
    init.score.df <- input.score.list[[i]]
    # ensure that only real genomic regions are being used
    temp.score.df <- c()
    if(length(which(is.na(init.score.df$start))) > 0){
      temp.score.df <- init.score.df[-which(is.na(init.score.df$start)),]
    } else {
      temp.score.df <- init.score.df
    }
    temp.score.name <- score.names[i]  #NB, viterbi;
    print(paste('Evaluate per region performance for:', temp.score.name, sep = ' '))
    temp.pos.labels <- input.specs$positiveLabels
    temp.neg.labels <- input.specs$negativeLabels

    if(temp.score.name == 'viterbi'){
      temp.viterbi.df <- data.frame(origformatScores = temp.score.df$formatScores,
        label = temp.score.df$origLabel, viterbi_label = temp.score.df$viterbiLabels,
        start = temp.score.df$start, end = temp.score.df$end, chrom = temp.score.df$chrom, stringsAsFactors = 'F')

      # assumption: enhancers are encoded as '1', non-enhancer as '2'
      ordered.temp.viterbi.df <- temp.viterbi.df[with(temp.viterbi.df, order(viterbi_label, -origformatScores)),]
      ordered.temp.viterbi.df$formatScores <- c(nrow(ordered.temp.viterbi.df):1)

      # temp.scores <- ordered.temp.viterbi.df$newScores
      # temp.labels <- ordered.temp.viterbi.df$label

      temp.score.df <- ordered.temp.viterbi.df
    }

    if('simulated_data' %in% names(input.specs)){
      if('exon' %in% unique(temp.score.df$label)){
        temp.score.df <- temp.score.df[which(temp.score.df$label %in% c(temp.pos.labels, temp.neg.labels) ),]
      }
    }

    if(length(unique(temp.score.df$chrom)) > 1){
      print('Warning: per region AUC can currently only be calculated for one chromosome')
      break
    }
    # take max score for each region and determine how well method did in detecting all positive regions
    region.scores <- pos_region_scores(temp.score.df, input.specs)

    auc.list <- AUC_x_y_perRegion(region.scores$region_scores, region.scores$region_labels,
      temp.pos.labels, temp.neg.labels)
    prAUC.list <- precisionRecall_AUC_perRegion(region.scores$region_scores, region.scores$region_labels,
      temp.pos.labels, temp.neg.labels)
    methd.auc[[temp.score.name]] <- auc.list
    methd.prAuc[[temp.score.name]] <- prAUC.list

    method.name <- temp.score.name  #strsplit(temp.score.name,'_')[[1]][1]
    Method <- c(Method, method.name)

    if(input.specs$dataArranging == 'none'){
      pre_bin_type <- c(pre_bin_type, 'none')
      pre_bin_size <- c(pre_bin_size, NA)
      pre_bin_step <- c(pre_bin_step, NA)
      pre_bin_methd <- c(pre_bin_methd, NA)
    } else {
      if(input.specs$binBy == 'mean'){
        pre_bin_methd <- c(pre_bin_methd, 'mean')
      } else {
        pre_bin_methd <- c(pre_bin_methd, 'sum')
      }
      if('preBinning' %in% names(input.specs)){
        pre_bin_type <- c(pre_bin_type,'nonOverl')
        pre_bin_size <- c(pre_bin_size, input.specs$preBinning)
        pre_bin_step <- c(pre_bin_step, NA)
      } else {
        pre_bin_type <- c(pre_bin_type,'overl')
        pre_bin_size <- c(pre_bin_size, input.specs$preOverlapBinning[2])
        pre_bin_step <- c(pre_bin_step, input.specs$preOverlapBinning[1])
      }
    }

    if('postScoringAnalysis' %in% names(input.specs)){
      temp.method.name <- strsplit(temp.score.name,'_')[[1]]
      if(length(temp.method.name) > 1){
        post.score.method <- temp.method.name[2]
        if(post.score.method == 'alphaRRA'){
          post_bin_methd <- c(post_bin_methd, 'alphaRRA')
          post_bin_size	<- c(post_bin_size, input.specs$binSize)
        } else if(post.score.method == 'CRESTrra'){
          post_bin_methd <- c(post_bin_methd, 'CRESTrra')
          post_bin_size	<- c(post_bin_size, input.specs$binSize)
        } else if(post.score.method == 'guideWindowRRA'){
          post_bin_methd <- c(post_bin_methd, 'guideWindowRRA')
          post_bin_size	<- c(post_bin_size, input.specs$guidePerFixedWindow)
        } else if(post.score.method == 'slidingWindow'){
          post_bin_methd <- c(post_bin_methd, 'slidingWindow')
          post_bin_size	<- c(post_bin_size, input.specs$guidePerSlidingWindow)
        } else {
          print('WARNING! No correct post scoring method found!')
        }
      } else {
        post_bin_methd <- c(post_bin_methd, 'none')
        post_bin_size	<- c(post_bin_size, NA)
      }
    } else {
      post_bin_methd <- c(post_bin_methd, 'none')
      post_bin_size	<- c(post_bin_size, NA)
    }
    AUC <- c(AUC, auc.list$num_AUC)
    prAUC <- c(prAUC, round(prAUC.list$num_AUC, 3))
  }
  plot_methodEval(methd.auc,'perRegion_AUCeval',input.specs$dataName,'bottomright','Fraction True negative included','Fraction True positives included')
  plot_methodEval(methd.prAuc,'perRegion_prAUCeval',input.specs$dataName,'bottomleft','Recall','Precision')

  out.df <- cbind.data.frame(pre_bin_type, pre_bin_size, pre_bin_step, pre_bin_methd,
    post_bin_methd, post_bin_size, Method, AUC, prAUC, stringsAsFactors = FALSE)
  write.table(out.df, file = 'perRegion_method_eval.csv', sep = ',', row.names = FALSE, col.names = FALSE)

  return(cbind.data.frame(Method, prAUC, stringsAsFactors = F))

}

# given tile deletion length, standard deviation around tile position, step size, number of tiles, number of chromosomes
#  generate a testing dataset consiting of both tiles and enhancers (two seperate files)
# input:
#   in_tileLength: length of tiles
#   in_tileSd: standard deviation used for varying tile length
#   in_tileStepSize: step size between tile start positions
#   in_tiles: of equal length to in_chroms vector to specify the number of tiles per chromosome
#   in_enhancers: of equal length to in_chroms vector to specify the number of enhancers per chromosome
#   in_enhancerLength: length of enhancers
#   in_enhancerSD: standard deviation used for varying enhancer length
#   targetName: name of file contining locations of gual guide targets
#   enhName: name of file contining locations of enhancers
# output: data frame, colnames: # files contain format: chromosome  start end
#exampe:
  #simulate_tilePosition_data(1000,75,50,c(6000),c(5),500,300,'TilePos','EnhancerPos')
simulate_tilePosition_data <- function(in_tileLength,in_tileSd,in_tileStepSize,in_tiles,in_enhancers,in_enhancerLength,
  in_enhancerSD,targetName,enhName){

  out_targetPos <- c()
  out_enhancerPos <- c()
  #for each chromosome
  for(j in 1:length(in_tiles)){
    chrom_tileNr <- in_tiles[j]
    chrom_enhancNr <- in_enhancers[j]
    chrom_out <- c()
    start_out <- c()
    end_out <- c()
    for(i in 1:chrom_tileNr){
      chrom_out <- c(chrom_out,paste('chr',j,sep = ''))
      start_out <- c(start_out,round(abs(rnorm(1,mean = i*in_tileStepSize,sd = in_tileSd))))
      end_out <- c(end_out, round(abs(rnorm(1,mean = ((i*in_tileStepSize)+in_tileLength),sd = in_tileSd))))
    }
    if(chrom_enhancNr > 0){
      enhancerLoc_options <- c(1:max(end_out))
      enhancer_start <- sample(enhancerLoc_options,chrom_enhancNr,replace = FALSE)
      enhancer_end <- enhancer_start + round(rnorm(length(enhancer_start),mean = in_enhancerLength, sd = in_enhancerSD))
      enh_chroms <- rep(paste('chr',j,sep = '') , chrom_enhancNr)
      temp_enhancer_out <- cbind.data.frame(enh_chroms,enhancer_start,enhancer_end,stringsAsFactors = FALSE)
      sorted_temp_enhancer_out <- temp_enhancer_out[order(temp_enhancer_out[,2]),]
      out_enhancerPos <- rbind.data.frame(out_enhancerPos,sorted_temp_enhancer_out)
    }
    #sorting
    temp_tilePos <- cbind.data.frame(chrom_out,start_out,end_out,stringsAsFactors = FALSE)
    sorted_temp_tilePos <- temp_tilePos[order(temp_tilePos[,2]),]
    out_targetPos <- rbind.data.frame(out_targetPos,sorted_temp_tilePos)
  }
  write.csv(x= out_targetPos, file = paste(targetName,'.csv', sep = ''))
  write.csv(x= out_enhancerPos, file = paste(enhName,'.csv', sep = ''))
}

evaluate_multi_simulations <- function(spec.file = NULL, score.file = NULL, label.file = NULL, data.dir = NULL){

  main.analysis.dir <- getwd()
  analysis.specs <- read_analysis_specs(spec.file, data.dir)

  if('noPostScoreAnalysis' %in% names(analysis.specs)){
    analysis.specs$postScoringAnalysis <- NULL
  }
  analysis.specs.names <- names(analysis.specs)

  combined.guide.prAUC <- c()
  combined.region.prAUC <- c()
  combined.element.prAUC <- c()

  for(sim in 1:analysis.specs$nrSims){
    sim.dir <- paste0(main.analysis.dir, '/', analysis.specs$simName, '_sim', sim)
    dir.create(file.path(sim.dir))
    setwd(file.path(sim.dir))
    print(paste0('analyzing simulation:', sim))
    analysis.specs$CountFileLoc <- paste0(analysis.specs$dataDir, analysis.specs$simName, '/',
      analysis.specs$simName, '_sim', sim, '_counts.csv')
    analysis.specs$sgRNAInfoFileLoc <- paste0(analysis.specs$dataDir, analysis.specs$simName, '/',
      analysis.specs$simName, '_sim', sim, '_info.csv')
    analysis.specs$pos_regions <- paste0(analysis.specs$dataDir, analysis.specs$simName, '/',
      analysis.specs$simName, '_sim', sim, '_enhancers.csv')

    # if no score file is provided, scores have to be calculated
    final.score.list <- list()  #instantiate the final score list (specifiy format?)

    if(is.null(score.file)){
      final.score.list <- calculate_scores(analysis.specs, label.file = label.file)  # needs to be implemented
    } else {
      final.score.list <- obtain_computed_scores(score.file)  # needs to be implemented
    }

    print('Genome scores obtained ...')

    if('evaluate_perGuide_Performance' %in% analysis.specs.names){
      # print('Evaluate per guide performance ...')
      # performanceEvaluation(final.score.list, analysis.specs)
      guide.prAUC <- performanceEvaluation_recording(final.score.list, analysis.specs)
      guide.prAUC.sim <- cbind(guide.prAUC, rep(paste0('sim',sim), nrow(guide.prAUC)))
      combined.guide.prAUC <- rbind(combined.guide.prAUC, guide.prAUC.sim)
    }
    if('evaluate_perRegion_Performance' %in% analysis.specs.names){
      # print('Evaluate per region performance ...')
      region.prAUC <- performanceEvaluation_perRegion_recording(final.score.list, analysis.specs)
      region.prAUC.sim <- cbind(region.prAUC, rep(paste0('sim',sim), nrow(region.prAUC)))
      combined.region.prAUC <- rbind(combined.region.prAUC, region.prAUC.sim)
    }
    if('evaluate_perElement_Performance' %in% analysis.specs.names){
      # print('Evaluate per region performance ...')
      element.prAUC <- performanceEvaluation_perElement_recording(final.score.list, analysis.specs)
      element.prAUC.sim <- cbind(element.prAUC, rep(paste0('sim',sim), nrow(element.prAUC)))
      combined.element.prAUC <- rbind(combined.element.prAUC, element.prAUC.sim)
    }

    if('postScoringAnalysis' %in% analysis.specs.names){
      print('Saving genome scores ...')
      save_all_scores(final.score.list, analysis.specs)
    }

    if('createBedGraph' %in% analysis.specs.names){
      print('Creating Bedgraphs ...')
      create_bedgraphs(final.score.list, analysis.specs)
    }

    if('scorePlotting' %in% analysis.specs.names){
      print('Score plotting ...')
      score_plotting(final.score.list, analysis.specs)
    }

    setwd('../')
  }

  colnames(combined.guide.prAUC) <- c('method', 'prAUC', 'sim')
  combined.guide.prAUC <- as.data.frame(combined.guide.prAUC)
  combined.guide.prAUC$prAUC <- as.numeric(combined.guide.prAUC$prAUC)
  total.guide.prAUC.plot <- ggplot(data = combined.guide.prAUC) + geom_point(aes(x = method, y = prAUC, colour = as.factor(sim))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  pdf(paste0(analysis.specs$simName,'_guide_prAUC.pdf'), width = 14, height = 7)
  print(total.guide.prAUC.plot)
  dev.off()
  write.csv(combined.guide.prAUC, file = paste0(analysis.specs$simName,'_guide_prAUC.csv'), row.names = F)

  colnames(combined.region.prAUC) <- c('method', 'prAUC', 'sim')
  combined.region.prAUC <- as.data.frame(combined.region.prAUC, stringsAsFactors = F)
  combined.region.prAUC$prAUC <- as.numeric(combined.region.prAUC$prAUC)
  total.region.prAUC.plot <- ggplot(data = combined.region.prAUC) + geom_point(aes(x = method, y = prAUC, colour = as.factor(sim))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  pdf(paste0(analysis.specs$simName,'_region_prAUC.pdf'), width = 14, height = 7)
  print(total.region.prAUC.plot)
  dev.off()
  write.csv(combined.region.prAUC, file = paste0(analysis.specs$simName,'_region_prAUC.csv'), row.names = F)

  colnames(combined.element.prAUC) <- c('method', 'prAUC', 'sim')
  combined.element.prAUC <- as.data.frame(combined.element.prAUC, stringsAsFactors = F)
  combined.element.prAUC$prAUC <- as.numeric(combined.element.prAUC$prAUC)
  total.element.prAUC.plot <- ggplot(data = combined.element.prAUC) + geom_point(aes(x = method, y = prAUC, colour = as.factor(sim))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  pdf(paste0(analysis.specs$simName,'_element_prAUC.pdf'), width = 14, height = 7)
  print(total.element.prAUC.plot)
  dev.off()
  write.csv(combined.element.prAUC, file = paste0(analysis.specs$simName,'_element_prAUC.csv'), row.names = F)

}

# master analysis function
# input: spec.file
#     specification file contains instructions about analysis types to use, plotting and evaluations
#   optional argument: scoreFile, labelFile
#     score.file: contains the location of the file containing already computed scores (avoids re-calculation)
#     label.file: contains the genomic coordinates for specific labels
#     data.dir: contains the location of the data directory. Used in case data only referenced by name, not loaction)
# output: none
#   There is no direct output to this function
#   Analyzing, file saving and plotting are all done during the function call
# per analysis run, only one type of data re-arrangement can be made
# per analysis run, only one type of filtering can be made
analyze_data <- function(spec.file = NULL, label.file = NULL, data.dir = NULL){

  # read in specification file
  if(is.null(spec.file)) stop('Error in "analyze_data()": specification file has to be provided')
  analysis.specs <- read_analysis_specs(spec.file, data.dir)

  if('noPostScoreAnalysis' %in% names(analysis.specs)){
    analysis.specs$postScoringAnalysis <- NULL
  }
  analysis.specs.names <- names(analysis.specs)

  # if no score file is provided, scores have to be calculated
  final.score.list <- list()  #instantiate the final score list (specifiy format?)

  final.score.list <- calculate_scores(analysis.specs, label.file = label.file)

  print('Genome scores obtained ...')
  if('evaluate_perGuide_Performance' %in% analysis.specs.names){
    # print('Evaluate per guide performance ...')
    # performanceEvaluation(final.score.list, analysis.specs)
    guide.prAUC <- performanceEvaluation_recording(final.score.list, analysis.specs)
  }

  if('evaluate_perRegion_Performance' %in% analysis.specs.names){
    # print('Evaluate per region performance ...')
    region.prAUC <- performanceEvaluation_perRegion_recording(final.score.list, analysis.specs)
  }

  if('postScoringAnalysis' %in% analysis.specs.names){
    print('Saving genome scores ...')
    save_all_scores(final.score.list, analysis.specs)
  }

  if('createBedGraph' %in% analysis.specs.names){
    print('Creating Bedgraphs ...')
    create_bedgraphs(final.score.list, analysis.specs)
  }

  if('scorePlotting' %in% analysis.specs.names){
    print('Score plotting ...')
    score_plotting(final.score.list, analysis.specs)
  }

}

# can only be done for targeting guides. Non-targeting fuides will be removed prior to
#  peak calling
# input: data frame
#   *Scores, chrom, label, start, end
# * scores will be adjusted based on input to 'genomeScores' or 'guideScores'
# output: list
#   peakData: peakScores, chrom, label, start, end, peakNr
#   peakRange: peakScores, chrom, label, start, end
identify_peaks <- function(input.df, input.name, input.peak.range, peak.minimum, out.name){
  input.df <- read.csv(input.df, header = T, stringsAsFactors = F)
  input.df <- input.df[-which(is.na(input.df$start)),]
  if('genomeScore' %in% colnames(input.df)){
    input.df <- input.df[which(input.df$genomeScore > peak.minimum),]
    input.df$scores <- input.df$genomeScore
  } else if('guideScore' %in% colnames(input.df)){
    input.df <- input.df[which(input.df$guideScore > peak.minimum),]
    input.df$scores <- input.df$guideScore
  } else {
    print('unknown file format!')
    print('need columns: chrom, start, end and either guideScore or genomeScore')
    break
  }
  output.df <- c()
  output.peak.ranges <- c()
  row.index <- 1
  peak.index <- 1  # keeps track of the rows that are in the currently established peak
  peak.nr <- 1  # kees track of how many peaks have been recorded so far
  current.peak <- GRanges(seqnames = input.df$chrom[peak.index],
    ranges = IRanges(input.df$start[peak.index] - input.peak.range, input.df$end[peak.index] + input.peak.range))

  while(row.index < nrow(input.df)){
    row.index <- row.index + 1
    cur.del <- GRanges(seqnames = input.df$chrom[row.index],
      ranges = IRanges(input.df$start[row.index] - input.peak.range, input.df$end[row.index] + input.peak.range))  # current deletion range

    cur.peak.overl <- as.data.frame(findOverlaps(current.peak, cur.del, type = 'any'))

    # if there is no overlap, add the processed current peak as its own peak
    # else, add the deletion to the set which makes up the current peak
    if(nrow(cur.peak.overl) == 0){
      final.peak <- input.df[peak.index,]
      final.peak$peakNr <- rep(paste0(input.name, '_peak_',peak.nr), nrow(final.peak))
      output.df <- rbind(output.df, final.peak)
      peak.range <- cbind.data.frame(peakScores = max(final.peak$scores),
        chrom = final.peak$chrom[1], label = paste0(input.name, '_peak_',peak.nr),
        start = min(final.peak$start) - input.peak.range, end = max(final.peak$end) + input.peak.range)
      output.peak.ranges <- rbind(output.peak.ranges, peak.range)
      peak.nr <- peak.nr + 1
      peak.index <- row.index
      current.peak <- GRanges(seqnames = input.df$chrom[peak.index],
        ranges = IRanges(input.df$start[peak.index] - input.peak.range, input.df$end[peak.index] + input.peak.range))
    } else {
      peak.index <- c(peak.index, row.index)
      current.peak <- GRanges(seqnames = input.df$chrom[peak.index],
        ranges = IRanges(input.df$start[peak.index] - input.peak.range, input.df$end[peak.index] + input.peak.range))
    }
  }
  final.peak <- input.df[peak.index,]
  final.peak$peakNr <- rep(paste0(input.name, '_peak_',peak.nr), nrow(final.peak))
  output.df <- rbind(output.df, final.peak)
  peak.range <- cbind.data.frame(peakScores = mean(final.peak$scores),
    chrom = final.peak$chrom[1], label = paste0(input.name, '_peak_',peak.nr),
    start = min(final.peak$start) - input.peak.range, end = max(final.peak$end) + input.peak.range)
  output.peak.ranges <- rbind(output.peak.ranges, peak.range)

  write.csv(output.df, file = paste0(out.name, '_peakData.csv'), row.names = F)
  write.csv(output.peak.ranges, file = paste0(out.name, '_peakRanges.csv'), row.names = F)
  #
  # out.list <- list(peakData = output.df, peakRange = output.peak.ranges)
  # return(out.list)

}

convert_to_bedgraph <- function(input.file, input.name, out.name){
  input.df <- read.csv(input.file, header = T, stringsAsFactors = F)
  bedgraph.colors <- col2rgb(rainbow(1, s = 1, v = 1, start = 0, end = 1, alpha = 1))
  score.df <- input.df[which(! is.na(input.df$start)),]

  if('genomeScore' %in% colnames(score.df)){
    score.df$scores <- score.df$genomeScore
  } else if('guideScore' %in% colnames(score.df)){
    score.df$scores <- score.df$guideScore
  } else if('peakScores' %in% colnames(score.df)){
    score.df$scores <- score.df$peakScores
  } else {
    print('unknown file format!')
    print('need columns: chrom, start, end and either guideScore or genomeScore')
    break
  }

  bg.color <- paste(bedgraph.colors[1,1], bedgraph.colors[2,1], bedgraph.colors[3,1], sep = ',')
  first.chrom <- score.df[which(input.df$chrom == input.df$chrom[1]),]
  header1 <- paste0('browser position ', first.chrom$chrom[1], ':', min(first.chrom$start, na.rm = T),
    '-', max(first.chrom$end, na.rm = T), '\n')
  header2 <- paste0("track type=bedGraph name='", out.name, '_', input.name, "_bg' description='",
    out.name, '_', input.name, "_bg' visibility=full color=", bg.color, '\n')

  score.gr <- GRanges(seqnames = score.df$chrom, ranges = IRanges(score.df$start, score.df$end),
    score = score.df$scores)
  score.gr <- sortSeqlevels(score.gr)
  score.gr <- sort(score.gr)
  df.score <- as.data.frame(score.gr)
  df.score.final <- df.score[,c(1,2,3,6)]

  out.file <- paste0(out.name, '_', input.name, ".bedgraph")
  cat(header1, file = out.file)
  cat(header2, file = out.file, append = TRUE)

  write.table(df.score.final, file = out.file, append = TRUE, row.names = F, col.names = F, quote = F, sep = '\t')

}

# given a list of scores, create bedgraphs for each one of them
create_bedgraphs <- function(input.score.list, input.specs){
  score.names <- names(input.score.list)
  nr.bg <- length(score.names)  # create unique colors for each bedgraph to plot
  bedgraph.colors <- col2rgb(rainbow(nr.bg,s = 1, v = 1, start = 0, end =  max(1, nr.bg - 1)/nr.bg, alpha = 1))
  for(i in 1:length(score.names)){
    intit.temp.score.df <- input.score.list[[i]]
    temp.score.df <- intit.temp.score.df[which(! is.na(intit.temp.score.df$start)),]
    temp.score.name <- score.names[i]
    print(paste0('Creating Bedgraph for: ', temp.score.name))
    temp.bg.color <- paste(bedgraph.colors[1,i], bedgraph.colors[2,i], bedgraph.colors[3,i], sep = ',')
    first.chrom <- temp.score.df[which(temp.score.df$chrom == temp.score.df$chrom[1]),]
    temp.header1 <- paste0('browser position ', first.chrom$chrom[1], ':', min(first.chrom$start, na.rm = T),
      '-', max(first.chrom$end, na.rm = T), '\n')
    temp.header2 <- paste0("track type=bedGraph name='", input.specs$dataName, temp.score.name, "_bg' description='",
      input.specs$dataName, temp.score.name, "_bg' visibility=full color=", temp.bg.color, '\n')

    temp.score.gr <- GRanges(seqnames = temp.score.df$chrom, ranges = IRanges(temp.score.df$start, temp.score.df$end),
      score = temp.score.df$formatScores)
    temp.score.gr <- sortSeqlevels(temp.score.gr)
    temp.score.gr <- sort(temp.score.gr)
    df.temp.score <- as.data.frame(temp.score.gr)
    df.temp.score.final <- df.temp.score[,c(1,2,3,6)]

    tmp.file <- paste0(input.specs$dataName, temp.score.name, ".bedgraph")
    cat(temp.header1, file = tmp.file)
    cat(temp.header2, file = tmp.file, append = TRUE)

    write.table(df.temp.score.final, file = tmp.file, append = TRUE, row.names = F, col.names = F, quote = F, sep = '\t')

    temp.method.name <- strsplit(score.names[i],'_')[[1]]
    if('viterbi' %in% temp.method.name | 'frwdBkwd' %in% temp.method.name | 'nbGlmm' %in% temp.method.name | 'nbGlmmDev' %in% temp.method.name){
      tmp.peaks <- df.temp.score.final[which(df.temp.score.final[,4] > 0),]

      tmp.file.peaks <- paste0(input.specs$dataName, temp.score.name, "_peaks.bedgraph")
      temp.header2.peaks <- paste0("track type=bedGraph name='", temp.score.name, "_peaks_bg' description='",
        temp.score.name, "_peaks_bg' visibility=full color=", temp.bg.color, '\n')

      cat(temp.header1, file = tmp.file.peaks)
      cat(temp.header2.peaks, file = tmp.file.peaks, append = TRUE)
      write.table(tmp.peaks, file = tmp.file.peaks, append = TRUE, row.names = F, col.names = F, quote = F, sep = '\t')
    }

    # viterbi needs explicit labelling of the best hidden path
    if('viterbi' %in% temp.method.name){
      viterbi.file <- paste0(input.specs$dataName, temp.score.name, "_hiddenPath.bedgraph")
      temp.header2.hidden <- paste0("track type=bedGraph name='", temp.score.name, "_hiddenPath_bg' description='",
        temp.score.name, "_hiddenPath_bg' visibility=full color=", temp.bg.color, '\n')

      cat(temp.header1, file = viterbi.file)
      cat(temp.header2.hidden, file = viterbi.file, append = TRUE)
      hidden.path <- rep(0, nrow(temp.score.df))
      hidden.path[which(temp.score.df$viterbiLabels == 1)] <- 1
      temp.viterbi.gr <- GRanges(seqnames = temp.score.df$chrom, ranges = IRanges(temp.score.df$start, temp.score.df$end),
        score = hidden.path)
      temp.viterbi.gr <- sortSeqlevels(temp.viterbi.gr)
      temp.viterbi.gr <- sort(temp.viterbi.gr)
      df.temp.viterbi <- as.data.frame(temp.viterbi.gr)
      df.temp.viterbi.final <- df.temp.viterbi[,c(1,2,3,6)]

      write.table(df.temp.viterbi.final, file = viterbi.file, append = TRUE, row.names = F, col.names = F, quote = F, sep = '\t')
    }

  }
}

score_plotting <- function(input.scores, input.specs){

  scores.to.plot <- input.scores

  # below is depricated: FDR scores are currently not calculated (conflict with HMM mehtods)
  # # for alphaRRA we want both p-values and FDR andjusted p-values
  # # need to duplicate results from alphaRRA and convert formatScores to -log10 of rawScores
  # if('postScoringAnalysis' %in% names(input.specs)){
  #   scores.to.duplicate <- which(unlist(lapply(strsplit(names(scores.to.plot), '_'), function(x){x[2] == 'alphaRRA'})) == TRUE)
  #   if(length(scores.to.duplicate) >= 1){
  #     score.names <- names(scores.to.plot)
  #     for(i in 1:length(scores.to.duplicate)){
  #       temp.name <- paste0(score.names[scores.to.duplicate[i]], '_FDR')
  #       temp.scores <- scores.to.plot[[scores.to.duplicate[i]]]
  #       temp.scores$formatScores <- -log10(temp.scores$rawScores)
  #       scores.to.plot[[temp.name]]<- temp.scores
  #     }
  #   }
  # }
  score.names <- names(scores.to.plot)

  # need to account for viterbi labelling
  orig.label.hierarchy <- input.specs$labelHierarchy

  for(i in 1:length(score.names)){
    raw.methd.scores <- scores.to.plot[[i]]
    na.rows <- which(is.na(raw.methd.scores$rawScores))

    # if not enough reads present, DESeq2 will return NA
    # remove these rows before plotting
    methd.scores <- c()
    if(length(na.rows) > 0){
      methd.scores <- raw.methd.scores[-na.rows,]
    } else {
      methd.scores <- raw.methd.scores
    }

    methd.name <- score.names[i]
    if('viterbiLabels' %in% colnames(methd.scores)){
      input.specs$labelHierarchy <- c(paste(orig.label.hierarchy, 'nonEnh', sep = '_'), paste(orig.label.hierarchy, 'enh', sep = '_'))
    } else {
      input.specs$labelHierarchy <- orig.label.hierarchy
    }
    if('plotAllChrom' %in% names(input.specs)){
      all.chroms <- unique(methd.scores$chrom)
      plotMultiChrom(methd.scores, input.specs, methd.name, all.chroms, 'plotAllChrom')
    }
    if('plotChroms' %in% names(input.specs)){
      chroms.to.plot <- input.specs$plotChroms
      plotMultiChrom(methd.scores, input.specs, methd.name, chroms.to.plot, 'plotChroms')
    }
    if('plotSeparateChroms' %in% names(input.specs)){
      chroms.to.plot <- unique(methd.scores$chrom)

      # if there is a 'NA' chromosome, remove it since there is no point in
      # plotting non-targeting guides by themselves
      if(NA %in% chroms.to.plot){
        chroms.to.plot <- chroms.to.plot[-which(is.na(chroms.to.plot))]
      }
      if('NA' %in% chroms.to.plot){
        chroms.to.plot <- chroms.to.plot[-which(chroms.to.plot == 'NA')]
      }
      plotSingleChrom(methd.scores, input.specs, methd.name, chroms.to.plot, 'plotSeparateChroms')
    }
    if('zoomChromosomePlot' %in% names(input.specs)){
      zoom.file <- read.csv(input.specs$zoomRange, as.is = TRUE, stringsAsFactors = FALSE)
      plotZoomChrom(methd.scores, input.specs, methd.name, chroms.to.plot, 'plotZoomChroms', zoom.file)
    }
    if('tileZoomChromosomePlot' %in% names(input.specs)){
      zoom.file <- read.csv(input.specs$tileZoomRange, as.is = TRUE, stringsAsFactors = FALSE)
      plotTileZoomChrom(methd.scores, input.specs, methd.name, chroms.to.plot, 'plotTileZoomChroms', zoom.file)
    }
  }
}

# plot each zoom regions individually
# input:
#   input.scores: data frame containing scores (rawScores, formatScores, log2_rate_ratio, chromosome, label, start, end)
#   input.specs: specification list
#   input.name: name of scores
#   input.chroms: chromosomes to be plotted next to each other
#   plot.type: specifying whether it's plotting all chromosomes or a subset of them
#   zoom.regions: regions which shuld be zoomed into (cols; chrom, start, end)
# output:
#   pdf containing scores per chromosome
plotTileZoomChrom <- function(input.scores, input.specs, input.name, input.chroms, plot.type, zoom.regions){

  chrom.df <- filter(input.scores, grepl('chr',chrom))
  na.df <- c()

  if(nrow(chrom.df) < nrow(input.scores)){
    na.df <- filter(input.scores, grepl('NA',chrom))
  }

  for(i in 1:nrow(zoom.regions)){
    zoom.df <- input.scores[input.scores$chrom == zoom.regions$chrom[i] & input.scores$start >= zoom.regions$start[i] & input.scores$end <= zoom.regions$end[i], ]
    zoom.x.start <- zoom.df$start
    zoom.x.end <- zoom.df$end
    zoom.y.pos <- zoom.df$formatScores
    zoom.labels <- zoom.df$label
    zoom.inter.chr.space <- round(1/8* (max(zoom.x.end) - min(zoom.x.start)))

    all.gene.pos <- input.specs$plotGenesData
    zoom.gene.df <- all.gene.pos[all.gene.pos$chrom == zoom.regions$chrom[i] & all.gene.pos$start >= zoom.regions$start[i] & all.gene.pos$end <= zoom.regions$end[i], ]
    gene.pos.df <- c()
    if(nrow(zoom.gene.df) > 0){
      gene.pos.df <- zoom.gene.df
    }

    if(nrow(chrom.df) < nrow(input.scores)){
      all.na.labels <- unique(na.df$label)
      for(lab in all.na.labels){
        na.lab.df <- na.df[na.df$label == lab,]
        na.lab.start <- max(zoom.x.end) + zoom.inter.chr.space
        na.lab.stepSize <- median(sort(zoom.x.start)[2:length(zoom.x.start)] - sort(zoom.x.start)[1:(length(zoom.x.start) - 1)])  # adjust step size to median step size of targeting guides
        na.lab.x.pos <- na.lab.start + na.lab.stepSize*c(1:nrow(na.lab.df))

        zoom.x.start <- c(zoom.x.start, na.lab.x.pos)
        zoom.x.end <- c(zoom.x.end, na.lab.x.pos+3)
        zoom.y.pos <- c(zoom.y.pos, na.lab.df$formatScores)
        zoom.labels <- c(zoom.labels, na.lab.df$label)
      }
    }
    chrom.range <- NULL
    pdf(paste(input.specs$dataName, input.name, plot.type, 'tileZoom', i, '.pdf',sep = '_'),width = 14, height = 7)
    create_tile_chromosome_scorePlot(zoom.x.start, zoom.x.end, zoom.y.pos, zoom.labels, gene.pos.df, chrom.range, input.specs)
    dev.off()
  }
}

# plot each zoom regions individually
# input:
#   input.scores: data frame containing scores (rawScores, formatScores, log2_rate_ratio, chromosome, label, start, end)
#   input.specs: specification list
#   input.name: name of scores
#   input.chroms: chromosomes to be plotted next to each other
#   plot.type: specifying whether it's plotting all chromosomes or a subset of them
#   zoom.regions: regions which shuld be zoomed into (cols; chrom, start, end)
# output:
#   pdf containing scores per chromosome
plotZoomChrom <- function(input.scores, input.specs, input.name, input.chroms, plot.type, zoom.regions){

  chrom.df <- filter(input.scores, grepl('chr',chrom))
  na.df <- c()

  if(nrow(chrom.df) < nrow(input.scores)){
    na.df <- filter(input.scores, grepl('NA',chrom))
  }

  for(i in 1:nrow(zoom.regions)){
    zoom.df <- input.scores[input.scores$chrom == zoom.regions$chrom[i] & input.scores$start >= zoom.regions$start[i] & input.scores$end <= zoom.regions$end[i], ]
    zoom.x.pos <- (zoom.df$start + zoom.df$end) / 2
    zoom.y.pos <- zoom.df$formatScores
    zoom.labels <- zoom.df$label
    zoom.inter.chr.space <- round(1/8* (max(zoom.x.pos) - min(zoom.x.pos)))

    gene.pos.df <- c()
    if(! is.null(input.specs$plotGenesData)){
      all.gene.pos <- input.specs$plotGenesData
      zoom.gene.df <- all.gene.pos[all.gene.pos$chrom == zoom.regions$chrom[i] & all.gene.pos$start >= zoom.regions$start[i] & all.gene.pos$end <= zoom.regions$end[i], ]
      if(nrow(zoom.gene.df) > 0){
        gene.pos.df <- zoom.gene.df
      }
    }

    if(nrow(chrom.df) < nrow(input.scores)){
      all.na.labels <- unique(na.df$label)
      for(lab in all.na.labels){
        na.lab.df <- na.df[na.df$label == lab,]
        na.lab.start <- max(zoom.x.pos) + zoom.inter.chr.space
        na.lab.stepSize <- median(sort(zoom.x.pos)[2:length(zoom.x.pos)] - sort(zoom.x.pos)[1:(length(zoom.x.pos) - 1)])  # adjust step size to median step size of targeting guides
        na.lab.x.pos <- na.lab.start + na.lab.stepSize*c(1:nrow(na.lab.df))

        zoom.x.pos <- c(zoom.x.pos, na.lab.x.pos)
        zoom.y.pos <- c(zoom.y.pos, na.lab.df$formatScores)
        zoom.labels <- c(zoom.labels, na.lab.df$label)
      }
    }
    chrom.range <- NULL
    pdf(paste(input.specs$dataName, input.name, plot.type, 'zoom', i, '.pdf',sep = '_'),width = 14, height = 7)
    create_chromosome_scorePlot(zoom.x.pos, zoom.y.pos, zoom.labels, gene.pos.df, chrom.range, input.specs)
    dev.off()
  }
}

# plot chromosomes individually, option to zoom in into specific region
# input:
#   input.scores: data frame containing scores (rawScores, formatScores, log2_rate_ratio, chromosome, label, start, end)
#   input.specs: specification list
#   input.name: name of scores
#   input.chroms: chromosomes to be plotted next to each other
#   plot.type: specifying whether it's plotting all chromosomes or a subset of them
# output:
#   pdf containing scores per chromosome
plotSingleChrom <- function(input.scores, input.specs, input.name, input.chroms, plot.type){

  chrom.df <- filter(input.scores, grepl('chr',chrom))
  chrom.midpoints <- (chrom.df$start + chrom.df$end) / 2
  all.chroms <- input.chroms[grep('chr', input.chroms, perl = TRUE)]
  na.df <- c()

  if(nrow(chrom.df) < nrow(input.scores)){
    na.df <- filter(input.scores, grepl('NA',chrom))
  }

  for(chrom in all.chroms){
    temp.chrom.df <- chrom.df[chrom.df$chrom == chrom,]
    temp.x.pos <- chrom.midpoints[chrom.df$chrom == chrom]
    temp.y.pos <- temp.chrom.df$formatScores
    temp.labels <- temp.chrom.df$label
    temp.inter.chr.space <- round(1/8* (max(temp.x.pos) - min(temp.x.pos))) + 1

    gene.pos.df <- c()
    gene.pos.df <- get_gene_Position(input.specs, chrom, temp.x.pos, gene.pos.df, 0, 0)

    if(! is.null(na.df) ){
      all.na.labels <- unique(na.df$label)
      for(lab in all.na.labels){
        na.lab.df <- na.df[na.df$label == lab,]
        na.lab.start <- max(temp.x.pos) + temp.inter.chr.space
        na.lab.stepSize <- median(sort(temp.x.pos)[2:length(temp.x.pos)] - sort(temp.x.pos)[1:(length(temp.x.pos) - 1)])  # adjust step size to median step size of targeting guides
        if(is.na(na.lab.stepSize)){
          na.lab.stepSize <- 1
        }
        na.lab.x.pos <- na.lab.start + na.lab.stepSize*c(1:nrow(na.lab.df))

        temp.x.pos <- c(temp.x.pos, na.lab.x.pos)
        temp.y.pos <- c(temp.y.pos, na.lab.df$formatScores)
        temp.labels <- c(temp.labels, na.lab.df$label)
      }
    }
    chrom.range <- NULL

    if('showFDR' %in% names(input.specs)){
      fdrs <- c()
      for(i in 1:length(input.specs$showFDR)){
        # amongs the rows which are below a certain FDR in the rawScore, return the minimum raw
        temp.fdr <- min(input.scores$formatScores[which(input.scores$rawScores <= input.specs$showFDR[i])])
        fdrs <- c(fdrs, temp.fdr)
      }
      input.specs$FDRcutoffs <- fdrs
    }

    pdf(paste(input.specs$dataName, input.name, plot.type, chrom, '.pdf',sep = '_'),width = 14, height = 7)
    create_chromosome_scorePlot(temp.x.pos, temp.y.pos, temp.labels, gene.pos.df, chrom.range, input.specs)
    dev.off()
  }
}

# plot multiple chromosome scores next to each other
# chromosome coordinates will be adjusted
# input:
#   input.scores: data frame containing scores (rawScores, formatScores, log2_rate_ratio, chromosome, label, start, end)
#   input.specs: specification list
#   input.name: name of scores
#   input.chroms: chromosomes to be plotted next to each other
#   plot.type: specifying whether it's plotting all chromosomes or a subset of them
# output:
#   pdf containing scores per chromosome
plotMultiChrom <- function(input.scores, input.specs, input.name, input.chroms, plot.type){

  chrom.df <- filter(input.scores, grepl('chr',chrom))
  all.chroms <- input.chroms[grep('chr', input.chroms, perl = TRUE)]
  na.df <- c()

  if(nrow(chrom.df) < nrow(input.scores)){
    na.df <- filter(input.scores, grepl('NA',chrom))
  }

  x.pos <- c()
  y.pos <- c()
  final.labels <- c()
  chrom.range <- c()
  inter.chr.space <- c()
  gene.pos.df <- c()  # data frame containing the coordinates for genes which are to be plotted

  for(chrom in all.chroms){
    temp.chrom.df <- chrom.df[chrom.df$chrom == chrom,]
    sorted.temp.chrom.df <- temp.chrom.df[order(temp.chrom.df$start),]
    temp.x.pos <- (sorted.temp.chrom.df$start + sorted.temp.chrom.df$end) / 2
    temp.y.pos <- sorted.temp.chrom.df$formatScores
    temp.labels <- sorted.temp.chrom.df$label

    if(length(x.pos) == 0){
      #have x-position be 1, and keep spacing between guides
      #define spacing between chromosomes
      shift.temp.x.pos <- temp.x.pos - temp.x.pos[1] + 1
      # inter.chr.space <- round(0.25*mean(shift.temp.x.pos))
      inter.chr.space <- round(1/8* (max(temp.x.pos) - min(temp.x.pos)))
      #if specific gene are to be plotted, check if for current chromosome a
      # gene is to be plotted, adjust gene position relative to the shift made
      gene.pos.df <- get_gene_Position(input.specs, chrom, shift.temp.x.pos, gene.pos.df, temp.x.pos[1], 0)

    } else {
      #have x-position be 1, and keep spacing between guides
      # start where the last chromosome left off with plus the interchromosome space to visually keep them apart
      coord.shift <- inter.chr.space + max(x.pos)
      shift.temp.x.pos <- temp.x.pos - temp.x.pos[1] + 1 + coord.shift
      gene.pos.df <- get_gene_Position(input.specs, chrom, shift.temp.x.pos, gene.pos.df, temp.x.pos[1], coord.shift)
    }
    x.pos <- c(x.pos, shift.temp.x.pos)
    y.pos <- c(y.pos, temp.y.pos)
    final.labels <- c(final.labels, temp.labels)

    chrom.range <- get_chrom_range(chrom, shift.temp.x.pos, gene.pos.df, chrom.range)
  }

  if(nrow(chrom.df) < nrow(input.scores)){
    all.na.labels <- unique(na.df$label)
    for(lab in all.na.labels){
      na.lab.df <- na.df[na.df$label == lab,]
      na.lab.start <- max(c(chrom.range[,3],x.pos)) + inter.chr.space
      na.lab.stepSize <- median(sort(x.pos)[2:length(x.pos)] - sort(x.pos)[1:(length(x.pos) - 1)])  # adjust step size to median step size of targeting guides
      na.lab.x.pos <- na.lab.start + na.lab.stepSize*c(1:nrow(na.lab.df))

      x.pos <- c(x.pos, na.lab.x.pos)
      y.pos <- c(y.pos, na.lab.df$formatScores)
      final.labels <- c(final.labels, na.lab.df$label)
    }
  }

  if('showFDR' %in% names(input.specs)){
    fdrs <- c()
    for(i in 1:length(input.specs$showFDR)){
      # amongs the rows which are below a certain FDR in the rawScore, return the minimum raw
      temp.fdr <- min(input.scores$formatScores[which(input.scores$rawScores <= input.specs$showFDR[i])])
      fdrs <- c(fdrs, temp.fdr)
    }
    input.specs$FDRcutoffs <- fdrs
  }

  pdf(paste(input.specs$dataName, input.name, plot.type, '.pdf',sep = '_'),width = 14, height = 7)
  create_chromosome_scorePlot(x.pos, y.pos, final.labels, gene.pos.df, chrom.range, input.specs)
  dev.off()
}

#plot the p-values along the chromosome of interest
create_tile_chromosome_scorePlot <- function(input.tile.start, input.tile.end, input.y.values,
  input.labels, input.gene.pos, input.chrom.range, input.specs){

  chrom.info <- cbind.data.frame(coords = input.tile.start, coordsEnd = input.tile.end, scores = input.y.values, labels = input.labels)

  out.chrom.plot <- ggplot()+
    geom_segment(data = chrom.info[chrom.info$label == input.specs$labelHierarchy[1], ], aes(x = coords, xend = coordsEnd,
      y = scores, yend = scores, colour = labels), size = 1)

  for(i in 2:length(input.specs$labelHierarchy)){
    if(input.specs$labelHierarchy[i] %in% unique(input.labels)){
      out.chrom.plot <- out.chrom.plot + geom_segment(data = chrom.info[chrom.info$label == input.specs$labelHierarchy[i], ],
        aes(x = coords, xend = coordsEnd, y = scores, yend = scores, colour = labels), size = 1.5)
    }
  }

  # keep chromosome range in a list to avoid overrriding of level specific variables
  if(! is.null(input.chrom.range)){
    chrom.range.list <- list()
    for(i in 1:nrow(input.chrom.range)){
      chrom.range <- cbind.data.frame(pos = c(input.chrom.range[i,2],input.chrom.range[i,3]),
      val = c((min(input.y.values)-2),(min(input.y.values)-2)))
      chrom.range.list[[i]] <- as.data.frame(chrom.range)
      out.chrom.plot <- out.chrom.plot + geom_line(data = chrom.range.list[[i]], aes(x = pos,y = val, size = 3),colour = 'red')
    }
  }

  # keep gene positions in a list to avoid overrriding of level specific variables
  if(! is.null(input.gene.pos)){
  gene.pos.list <- list()
    for(i in 1:nrow(input.gene.pos)){
      gene.pos <- cbind.data.frame(pos = c(input.gene.pos[i,2],input.gene.pos[i,3]), val = c(0,0))
      gene.pos.list[[i]] <- as.data.frame(gene.pos)
      out.chrom.plot <- out.chrom.plot + geom_line(data = gene.pos.list[[i]], aes(x = pos,y = val, size = 4),colour = 'orange')
    }
  }
  out.chrom.plot <- out.chrom.plot + theme_bw() + scale_colour_discrete(name = 'Scores')
  print(out.chrom.plot)
}

#plot the p-values along the chromosome of interest
create_chromosome_scorePlot <- function(input.x.values, input.y.values, input.labels, input.gene.pos, input.chrom.range, input.specs){
  chrom.info <- cbind.data.frame(coords = input.x.values, scores = input.y.values, labels = input.labels)

  out.chrom.plot <- ggplot()+
    geom_point(data = chrom.info[chrom.info$label == input.specs$labelHierarchy[1], ], aes(x = coords, y = scores,
      colour = labels), size = 1)

  for(i in 2:length(input.specs$labelHierarchy)){
    if(input.specs$labelHierarchy[i] %in% unique(input.labels)){
      out.chrom.plot <- out.chrom.plot + geom_point(data = chrom.info[chrom.info$label == input.specs$labelHierarchy[i], ],
        aes(x = coords, y = scores, colour = labels), size = 1.5)
    }
  }

  # keep chromosome range in a list to avoid overrriding of level specific variables
  if(! is.null(input.chrom.range)){
    chrom.range.list <- list()
    for(i in 1:nrow(input.chrom.range)){
      chrom.range <- cbind.data.frame(pos = c(input.chrom.range[i,2],input.chrom.range[i,3]),
      val = c((min(input.y.values)-2),(min(input.y.values)-2)))
      chrom.range.list[[i]] <- as.data.frame(chrom.range)
      out.chrom.plot <- out.chrom.plot + geom_line(data = chrom.range.list[[i]], aes(x = pos,y = val, size = 3),colour = 'red')
    }
  }

  # keep gene positions in a list to avoid overrriding of level specific variables
  if(! is.null(input.gene.pos)){
  gene.pos.list <- list()
    for(i in 1:nrow(input.gene.pos)){
      gene.pos <- cbind.data.frame(pos = c(input.gene.pos[i,2],input.gene.pos[i,3]), val = c(0,0))
      gene.pos.list[[i]] <- as.data.frame(gene.pos)
      out.chrom.plot <- out.chrom.plot + geom_line(data = gene.pos.list[[i]], aes(x = pos,y = val, size = 4),colour = 'orange')
    }
  }

  if('showFDR' %in% names(input.specs)){
    fdr.cutoff.list <- list()
    for(i in 1:length(input.specs$showFDR)){
      fdr.line <- cbind.data.frame(pos = c(min(chrom.info$coords),max(chrom.info$coords)),
        FDRcutoff = c(input.specs$FDRcutoffs[i], input.specs$FDRcutoffs[i]))
      fdr.cutoff.list[[i]] <- fdr.line
      out.chrom.plot <- out.chrom.plot + geom_line(data = fdr.cutoff.list[[i]],
        aes(x = pos, y = FDRcutoff, size = 0.8), colour = 'blue')
    }
  }

  out.chrom.plot <- out.chrom.plot + theme_bw() + scale_colour_discrete(name = 'Scores')
  print(out.chrom.plot)
}

#plot the p-values along the chromosome of interest
create_chromosome_pValuePlot_old <- function(in_x_values,in_y_values,in_colors,in_labels,in_colorType,input.gene.pos,in_ChrSectrLabel){
  chrom_p_loc <- cbind.data.frame(mPoints = in_x_values, scores = in_y_values,in_colors)
  colnames(chrom_p_loc) <- c('mPoints','pVal','color')

  out_chrom_P <- ggplot(chrom_p_loc, aes(x=mPoints, y=scores, colour=color))+
    geom_point(data=chrom_p_loc[which(chrom_p_loc[,3] == in_labels[1]), ], alpha=0.7)

  for(i in 2:length(in_labels)){
    out_chrom_P <- out_chrom_P + geom_point(data=chrom_p_loc[which(chrom_p_loc[,3] == in_labels[i]), ], colour = in_colorType[i], size=1.8)
  }

  if(! is.null(in_ChrSectrLabel)){
    sectr_pos_list <- list()
    for(i in 1:nrow(in_ChrSectrLabel)){
      temp_sectr_pos <- cbind(c(in_ChrSectrLabel[i,2],in_ChrSectrLabel[i,3]),c(-2,-2))
      colnames(temp_sectr_pos) <-c('pos','val')
      sectr_pos_list[[i]] <- as.data.frame(temp_sectr_pos)
      # print(sectr_pos_list[[i]])
      out_chrom_P <- out_chrom_P + geom_line(data = sectr_pos_list[[i]], aes(x = pos,y = val, size = 3),colour = 'red')
    }
  }

  gene_pos_list <- list()
  for(i in 1:nrow(in_genePosDf)){
    temp_gene_pos <- cbind(c(in_genePosDf[i,2],in_genePosDf[i,3]),c(0,0))
    colnames(temp_gene_pos) <-c('pos','val')
    gene_pos_list[[i]] <- as.data.frame(temp_gene_pos)
    out_chrom_P <- out_chrom_P + geom_line(data = gene_pos_list[[i]], aes(x = pos,y = val, size = 4),colour = 'orange')
  }

  out_chrom_P <- out_chrom_P + theme_bw()
  print(out_chrom_P)
}

# for a given chromosoem, create the sector delineating start and end (for plotting)
# input:
#   input.chrom: chromsome for which sector to be obtained (ex: 'chr6')
#   input.coords: chromosome coordinates of all scores (used to adjust gene coordinates)
#   input.gene.pos: exisiting gene coordinates for other chromsomes, already adjusted (cols: chrom, start, end)
#   input.range: existing sectors from other chromosomes
# output: data frame
#   cols: chromRange, start, end
get_chrom_range <- function(input.chrom, input.coords, input.gene.pos, input.range){
  chrom.gene.pos <- input.gene.pos[input.gene.pos$chrom == input.chrom,]
  chrom.range <- cbind.data.frame(chromRange = input.chrom,
    start = min(min(input.coords), min(chrom.gene.pos[,2])),
    end =  max(max(input.coords), max(chrom.gene.pos[,3])), stringsAsFactors = FALSE)
  out.chr.range <- rbind(input.range, chrom.range)
  return(out.chr.range)
}

# given a chromosome to plot, obtaing the gene coordinates, if available, and adjust accordingly
# input:
#   input.specs: specification list
#   input.chrom: chromsome for which gene coordinates are to be obtained (ex: 'chr6')
#   input.coords: chromosome coordinates of all scores (used to adjust gene coordinates)
#   input.gene.pos: exisiting gene coordinates for other chromsomes, already adjusted (cols: chrom, start, end)
#   input.initial.shift: if multiple chromosomes are plotted next to each other , this shoul dbe the minimum x-value, to start the coordinates at 1
#   input.additional.coord.shift: number by which genes should be shifted after setting coordinates using initial.shift, used if multiple chromosomes are plotted at the same time
# output: data dataframe
#   cols. chrom, start, end
get_gene_Position <- function(input.specs, input.chrom, input.coords, input.gene.pos, input.initial.shift, input.additional.coord.shift){
  sorted.input.coords <- sort(input.coords)
  if('plotGenes' %in% names(input.specs)){
    chrom.gene.rows <- which(input.specs$plotGenesData[,1] == input.chrom)
    if(length(chrom.gene.rows) > 0){
      temp.gene.pos <- input.specs$plotGenesData[chrom.gene.rows,]
      new.gene.pos <- cbind.data.frame(chrom = temp.gene.pos[,1],
        start = (temp.gene.pos[,2] - input.initial.shift + 1 + input.additional.coord.shift),
        end = (temp.gene.pos[,3] - input.initial.shift + 1 + input.additional.coord.shift), stringsAsFactors = FALSE)
      gene.pos.df <- rbind(input.gene.pos, new.gene.pos)
    } else {  # if no gene is specified for said chromosome simply give a minimal range
      new.gene.pos <- cbind.data.frame(chrom = input.chrom,
        start = sorted.input.coords[1], end = sorted.input.coords[1] + 1, stringsAsFactors = FALSE)
      gene.pos.df <- rbind(input.gene.pos, new.gene.pos)
    }
  } else {
    new.gene.pos <- cbind.data.frame(chrom = input.chrom,
      start = sorted.input.coords[1], end = sorted.input.coords[1] + 1, stringsAsFactors = FALSE)
    gene.pos.df <- rbind(input.gene.pos, new.gene.pos)
  }
  return(gene.pos.df)
}

# save te final scores used for plotting, AUC etc.
# input: list with scores per method, input specifications
# output: .csv for each methods scores
save_all_scores <- function(input.score.list, input.specs){
  score.names <- names(input.score.list)
  for(i in 1:length(score.names)){
    temp.score.df <- input.score.list[[i]]
    temp.score.name <- score.names[i]
    method.type <- strsplit(temp.score.name, '_')[[1]]
    method.identifier <- method.type[1]
    if(length(method.type) > 1){
      if(method.type[2] == 'genomeScores'){
        out.csv <- data.frame(genomeScore = temp.score.df$formatScores, chrom = temp.score.df$chrom,
          start = temp.score.df$start, end = temp.score.df$end, label = temp.score.df$label,
          log2_FC = temp.score.df$log2_rate_ratio, nrSupportGuides = temp.score.df$nrSupportGuides,
          stringsAsFactors = F)
      } else {
        out.csv <- data.frame(genomeScore = temp.score.df$formatScores, chrom = temp.score.df$chrom,
          start = temp.score.df$start, end = temp.score.df$end, label = temp.score.df$label,
          log2_FC = temp.score.df$log2_rate_ratio, stringsAsFactors = F)
      }
      write.csv(out.csv, file = paste(input.specs$dataName, method.identifier, method.type[2], '.csv', sep = '_'), row.names = FALSE)
    }
  }
}

# save te final scores used for plotting, AUC etc.
# input: list with scores per method, input specifications
# output: .csv for each methods scores
save_all_guide_scores <- function(input.score.list, input.specs){
  score.names <- names(input.score.list)
  for(i in 1:length(score.names)){
    temp.score.df <- input.score.list[[i]]
    temp.score.name <- score.names[i]
    method.identifier <- strsplit(temp.score.name, '_')[[1]][1]
    out.csv <- temp.score.df
    out.names <- names(out.csv)
    out.names[which(out.names == 'formatScores')] <- 'guideScore'
    names(out.csv) <- out.names

    # out.csv <- data.frame(guideScore = temp.score.df$formatScores, chrom = temp.score.df$chrom,
    #   start = temp.score.df$start, end = temp.score.df$end, label = temp.score.df$label, stringsAsFactors = F)
    write.csv(out.csv, file = paste(input.specs$dataName, method.identifier, 'guideScores.csv', sep = '_'), row.names = FALSE)
  }
}

# depricated: no longer contains all features necessary (can't correctly analyze simulated data)
# plotting of both AUC and prAUC of the scores based on labels
# input:
#   input.score.list: list of scores in score output format (rawScores, formatScores, log2_rate_ratio, chromosome, labels, start, end)
#   input.specs: list of specifications
# output: AUC and prAUC plot of scores
performanceEvaluation <- function(input.score.list, input.specs){

  methd.auc <- list()
  methd.prAuc <- list()

  score.names <- names(input.score.list)
  for(i in 1:length(score.names)){
    temp.score.df <- input.score.list[[i]]
    temp.score.name <- score.names[i]
    methd.auc[[temp.score.name]] <- AUC_x_y(temp.score.df$formatScores, temp.score.df$label, input.specs$positiveLabels, input.specs$negativeLabels)
    methd.prAuc[[temp.score.name]] <- precisionRecall_AUC(temp.score.df$formatScores, temp.score.df$label, input.specs$positiveLabels, input.specs$negativeLabels)
  }
  plot_methodEval(methd.auc,'AUCeval',input.specs$dataName,'bottomright','Fraction True negative included','Fraction True positives included')
  plot_methodEval(methd.prAuc,'prAUCeval',input.specs$dataName,'bottomleft','Recall','Precision')

}

# plot either AUC or prAUCs for multiple method scores
# input:
#   in_list: list containing:xAxis, yAxis , num_AUC (numerical AUC value)
#   in_evalType: either AUC or prAUC, used for labelling plot
#   in_datName: unique identifier of the data set
#   lgndPos: placement of the legend on the plot
#   xLabl: X-axis label
#   yLabl: y-axis label
# output: plot or either an AUC or prAUC of all scores
plot_methodEval <- function(in_list,in_evalType,in_datName,lgndPos,xLabl,yLabl){
  nr.aucs <- length(in_list)
  lgnd_colors <- rainbow(nr.aucs,s = 1, v = 1, start = 0, end =  max(1, nr.aucs - 1)/nr.aucs, alpha = 1)
  lgnd_txt <- c()
  methd.names <- names(in_list)
  pdf(paste(in_datName,'_',in_evalType,'.pdf',sep = ''))
  plot(x = in_list[[1]]$xAxis, y = in_list[[1]]$yAxis, type = 'l', col = lgnd_colors[1],main = in_evalType,xlab = xLabl, ylab = yLabl)
  lgnd_txt <- paste(methd.names[1],': ', in_list[[1]]$num_AUC,sep = '')
  if(length(in_list) > 1)
  for(i in 2:length(in_list)){
    lines(x = in_list[[i]]$xAxis, y = in_list[[i]]$yAxis, type = 'l', col = lgnd_colors[i])
    lgnd_txt <- c(lgnd_txt,paste(methd.names[i],': ', in_list[[i]]$num_AUC,sep = ''))
  }
  legend(lgndPos,legend=lgnd_txt, text.col=lgnd_colors[1:length(in_list)])
  dev.off()
}

# calculating area under the curve using regular AUC, looking at %negative control vs %positive control included
# input:
#   formatted.pvals: formatted p-values as input. Largest one is treated as most significant one
#   in_label: vector of equal length as formatted.pvals, assigns leable with each p-value
#   pos_labels: label identifier for positive labels (ex.: 'pos', 'exon', ..)
#   neg_labels: label identifier for negative labels (ex.: 'neg', 'nonTargeting', ..)
#output: list containing: xAxis, yAxis, AUC (already numeric)
AUC_x_y <- function(formatted.pvals,in_label,pos_labels,neg_labels){

  comb_df <- cbind.data.frame(-formatted.pvals,in_label,stringsAsFactors = FALSE)
  sorted_comb_df <- comb_df[order(comb_df[,1]),]

  #x-value represents the number of false positives incorporated
  x_axis <- 0
  #y-value represents the % of total true positives found
  y_axis <- 0

  total_nr_pos <- length(which(in_label %in% pos_labels))  #count_elemnt_occurance(pos_labels,in_label)
  total_nr_neg <- length(which(in_label %in% neg_labels))  #$count_elemnt_occurance(neg_labels,in_label)

  #not all positives have been found
  all_pos_NOT_found <- TRUE
  posit_tracker <- 1
  nr_pos_found <- 0
  nr_neg_found <- 0

  #category is either true positive (1) or false positive (0)
  auc_category <- c()

  while(all_pos_NOT_found){

    temp_row <- sorted_comb_df[posit_tracker,]

    if(temp_row[2] %in% pos_labels){
      nr_pos_found <- nr_pos_found + 1
      x_axis <- c(x_axis,nr_neg_found)
      y_axis <- c(y_axis,nr_pos_found)
      auc_category <- c(auc_category, 1)
      if(nr_pos_found == total_nr_pos){
        all_pos_NOT_found <- FALSE
        break
      }
    }
    else if(temp_row[2] %in% neg_labels){
      nr_neg_found <- nr_neg_found + 1
      x_axis <- c(x_axis,nr_neg_found)
      y_axis <- c(y_axis,nr_pos_found)
      auc_category <- c(auc_category, 0)
    }
    posit_tracker <- posit_tracker + 1
  }

  if(nr_neg_found < total_nr_neg){
    x_axis <- c(x_axis,c((nr_neg_found+1):total_nr_neg))
    y_axis <- c(y_axis,rep(total_nr_pos,length(c((nr_neg_found+1):total_nr_neg))))
  }

  #following example on r bloggers
  add_negatives <- rep(0,total_nr_neg-nr_neg_found)
  auc_category <- c(auc_category,add_negatives)
  auc_prediction <- rev(seq_along(auc_category))
  roc_obj <- roc(auc_category, auc_prediction)

  return(list(xAxis = (x_axis/total_nr_neg), yAxis = (y_axis/total_nr_pos), num_AUC = round(pROC::auc(roc_obj),3)))
}

#calculating precision recall AUC. recall (x) vs precision (y)
#using linear function to calculate area under the curve
#input: p-values (not negative ones, i.e. -log10), labels, label numbers
#output: list with 3 elements: xAxis, yAxis, num_AUC (already numeric)
precisionRecall_AUC <- function(formatted.pvals,in_label,pos_labels,neg_labels){

  comb_df <- cbind.data.frame(-formatted.pvals,in_label,stringsAsFactors = FALSE)
  sorted_comb_df <- comb_df[order(comb_df[,1]),]

  #x-value represents the number of false positives incorporated
  x_axis <- 0
  #y-value represents the % of total true positives found
  y_axis <- 1

  #not all positives have been found
  nr_pos_found <- 0
  total_pos <- length(which(in_label %in% pos_labels))
  total_found <- 0

  for(i in 1:length(formatted.pvals)){

    temp_row <- sorted_comb_df[i,]

    if(temp_row[2] %in% pos_labels){
      total_found <- total_found + 1
      nr_pos_found <- nr_pos_found + 1
      x_axis <- c(x_axis,nr_pos_found/total_pos)
      y_axis <- c(y_axis,nr_pos_found/total_found)
      if(nr_pos_found == total_pos){
        break
      }
    }
    else if(temp_row[2] %in% neg_labels){
      total_found <- total_found + 1
      x_axis <- c(x_axis,nr_pos_found/total_pos)
      y_axis <- c(y_axis,nr_pos_found/total_found)
    }
  }

  prAUC = round(MESS::auc(x= x_axis, y = y_axis,type='linear'),3)

  x_axis <- c(x_axis,x_axis[length(x_axis)])
  y_axis <- c(y_axis,0)

  return(list(xAxis = x_axis, yAxis = y_axis, num_AUC = prAUC))
}

# calculating area under the curve using regular AUC, looking at %negative control vs %positive control included
# input:
#   formatted.pvals: formatted p-values as input. Largest one is treated as most significant one
#   in_label: vector of equal length as formatted.pvals, assigns leable with each p-value
#   pos_labels: label identifier for positive labels (ex.: 'pos', 'exon', ..)
#   neg_labels: label identifier for negative labels (ex.: 'neg', 'nonTargeting', ..)
#output: list containing: xAxis, yAxis, AUC (already numeric)
AUC_x_y_perRegion <- function(formatted.pvals,in_label,pos_labels,neg_labels){

  comb_df <- cbind.data.frame(formatted.pvals,in_label,stringsAsFactors = FALSE)
  sorted_comb_df <- comb_df[order(comb_df[,1], comb_df[,2], decreasing = TRUE),]

  #x-value represents the number of false positives incorporated
  x_axis <- 0
  #y-value represents the % of total true positives found
  y_axis <- 0

  total_nr_pos <- length(which(in_label %in% pos_labels))  #count_elemnt_occurance(pos_labels,in_label)
  total_nr_neg <- length(which(in_label %in% neg_labels))  #$count_elemnt_occurance(neg_labels,in_label)

  #not all positives have been found
  all_pos_NOT_found <- TRUE
  posit_tracker <- 1
  nr_pos_found <- 0
  nr_neg_found <- 0

  #category is either true positive (1) or false positive (0)
  auc_category <- c()

  while(all_pos_NOT_found){

    temp_row <- sorted_comb_df[posit_tracker,]

    if(temp_row[2] %in% pos_labels){
      nr_pos_found <- nr_pos_found + 1
      x_axis <- c(x_axis,nr_neg_found)
      y_axis <- c(y_axis,nr_pos_found)
      auc_category <- c(auc_category, 1)
      if(nr_pos_found == total_nr_pos){
        all_pos_NOT_found <- FALSE
        break
      }
    }
    else if(temp_row[2] %in% neg_labels){
      nr_neg_found <- nr_neg_found + 1
      x_axis <- c(x_axis,nr_neg_found)
      y_axis <- c(y_axis,nr_pos_found)
      auc_category <- c(auc_category, 0)
    }
    posit_tracker <- posit_tracker + 1
  }

  if(nr_neg_found < total_nr_neg){
    x_axis <- c(x_axis,c((nr_neg_found+1):total_nr_neg))
    y_axis <- c(y_axis,rep(total_nr_pos,length(c((nr_neg_found+1):total_nr_neg))))
  }

  #following example on r bloggers
  add_negatives <- rep(0,total_nr_neg-nr_neg_found)
  auc_category <- c(auc_category,add_negatives)
  auc_prediction <- rev(seq_along(auc_category))
  roc_obj <- roc(auc_category, auc_prediction)

  return(list(xAxis = (x_axis/total_nr_neg), yAxis = (y_axis/total_nr_pos), num_AUC = round(pROC::auc(roc_obj),3)))
}

#calculating precision recall AUC. recall (x) vs precision (y)
#using linear function to calculate area under the curve
#input: p-values (not negative ones, i.e. -log10), labels, label numbers
#output: list with 3 elements: xAxis, yAxis, num_AUC (already numeric)
precisionRecall_AUC_perRegion <- function(formatted.pvals,in_label,pos_labels,neg_labels){

  comb_df <- cbind.data.frame(formatted.pvals,in_label,stringsAsFactors = FALSE)
  sorted_comb_df <- comb_df[order(comb_df[,1], comb_df[,2], decreasing = TRUE),]

  #x-value represents the number of false positives incorporated
  x_axis <- 0
  #y-value represents the % of total true positives found
  y_axis <- 1

  #not all positives have been found
  nr_pos_found <- 0
  total_pos <- length(which(in_label %in% pos_labels))
  total_found <- 0

  for(i in 1:length(formatted.pvals)){

    temp_row <- sorted_comb_df[i,]

    if(temp_row[2] %in% pos_labels){
      total_found <- total_found + 1
      nr_pos_found <- nr_pos_found + 1
      x_axis <- c(x_axis,nr_pos_found/total_pos)
      y_axis <- c(y_axis,nr_pos_found/total_found)
      if(nr_pos_found == total_pos){
        break
      }
    }
    else if(temp_row[2] %in% neg_labels){
      total_found <- total_found + 1
      x_axis <- c(x_axis,nr_pos_found/total_pos)
      y_axis <- c(y_axis,nr_pos_found/total_found)
    }
  }

  prAUC = round(MESS::auc(x= x_axis, y = y_axis,type='linear'),3)

  x_axis <- c(x_axis,x_axis[length(x_axis)])
  y_axis <- c(y_axis,0)

  return(list(xAxis = x_axis, yAxis = y_axis, num_AUC = prAUC))
}

# given the specifications, calculate the scores for the specified methods
# input: analysis specification file, label file (optional)
#   for further details on input files, see analyze_data() documentation
# output: list containing dataframes for each analysis method, each dataframe having the following columns:
#   $method:
#     rawScores: calculated probabilities of p-values
#     formatScores: adjusted scores such that they are signed by ratio and can be ranked from most to least significant
#     log2_rate_ratio: log2 ratio used for sign-adjusting scores, either based on count ratio or calculated rates
#     chromosome, labels, sgRNA targets (can be 1 or 2 columns) (initial values from the info file)
calculate_scores <- function(analysis.specs, label.file = NULL){

  analysis.specs.names <- names(analysis.specs)
  count.data <- read.csv(analysis.specs$CountFileLoc,as.is = TRUE)
  data.info <- read.csv(analysis.specs$sgRNAInfoFileLoc,as.is = TRUE, stringsAsFactors = FALSE)

  #convert NA in chromosome column to string 'NA'
  if(length(which(is.na(data.info[,1]))) > 0){
    data.info[which(is.na(data.info[,1])),1] <- 'NA'
  }
  #if a label file is supplied, add corresponding note to the specification file
  if(! is.null(label.file)){
    analysis.specs$labelFile <- label.file
  }

  #extract samples of interest
  # output: list
  #   $counts, contains counts of relevant samples
  #   $info, contains chromsome and positional information
  #   $specs, contains specs and added $subGroup1 and $subGroup2
  sample.list <- extract_samples(count.data,data.info,analysis.specs)
  sample.counts <- sample.list$counts
  sample.info <- sample.list$info
  sample.specs <- sample.list$specs

  # re-arrange data if specified
  # output: list
  #   $counts, contains re-arranged counts
  #   $info, contains chromsome and positional information about re-arranged data
  #   $specs, contains specs about re-arranged data
  arranged.data.list <- reArrange_data(sample.counts, sample.info, sample.specs)
  arranged.data.counts <- arranged.data.list$counts
  arranged.data.info <- arranged.data.list$info
  arranged.data.specs <- arranged.data.list$specs

  #filter data if specified
  filtered.data.list <- filtering_counts(arranged.data.counts, arranged.data.info, arranged.data.specs)
  filtered.data.counts <- filtered.data.list$counts
  filtered.data.info <- filtered.data.list$info
  filtered.data.specs <- filtered.data.list$specs

  #added this line
  filtered.data.specs$new_grp1 <- filtered.data.specs$subGroup1
  filtered.data.specs$new_grp2 <- filtered.data.specs$subGroup2

  # if specified; counts, info and rows removed will be saved for reproducability
  save_analysisData(filtered.data.counts, filtered.data.info, filtered.data.specs)

  #output scores from specified methods
  method.score.list <- score_calculation(filtered.data.counts, filtered.data.info, filtered.data.specs)

  print('Guide scores obtained ...')
  if('saveGuideScores' %in% names(filtered.data.specs)){
    print('Saving guide scores ...')
    save_all_guide_scores(method.score.list, filtered.data.specs)
  }
  #apply additional methods on the output scores for high-level scores
  post.method.score.list <- post_score_calculation(method.score.list, filtered.data.specs)

  return(post.method.score.list)
}

#for all scores, apply additional post-scoring methods
# obtain the scores for all specified analysis methods
# input: list, containing the following:
#   rawScores, formatScores, log2_rate_ratio, chromosome, label, sgRNA targets (can be 1 or 2 columns)
# output: list, each method adds the following to the
#   rawScores, formatScores, log2_rate_ratio, chromosome, label, sgRNA targets (can be 1 or 2 columns)
#   for more description see: score_calculation
post_score_calculation <- function(input.list, input.specs){
  if('postScoringAnalysis' %in% names(input.specs)){
    out.list <- list()
    if('postScoringMAGeCK' %in% names(input.specs)){
      pSM.path <- input.specs$MAGeCKrra
      out.list <- c(out.list, RRA_wrapper(input.list, input.specs,pSM.path, 'MAGeCKrra'))
    }
    if('postScoringCREST' %in% names(input.specs)){
      pSC.path <- input.specs$CRESTrra
      out.list <- c(out.list, RRA_wrapper(input.list, input.specs, pSC.path, 'CRESTrra'))
    }
    if('postScoreAlphaRRA' %in% names(input.specs)){
      out.list <- c(out.list, alphaRRA_wrapper(input.list, input.specs))
    }
    if('postScoreGuideWindowRRA' %in% names(input.specs)){
      out.list <- c(out.list, guideWindowRRA_wrapper(input.list, input.specs))
    }
    if('postScoreSlidingWindow' %in% names(input.specs)){
      out.list <- c(out.list, slidingWindow_wrapper(input.list, input.specs))
    }
    if('RELICS_genomeScoring' %in% names(input.specs)){
      out.list <- c(out.list, RELICS_genomeScoring_wrapper(input.list, input.specs))
    }
    return(out.list)
  } else {
    return(input.list)
  }
}

#apply sliding window to all score datasets
#input:
#   input.list: all method-specifi results to be analyzed using a specific RRA
#   input.specs: specification list
#output: list
#   rawScores, formatScores, log2_rate_ratio, chromosome, label, sgRNA targets (can be 1 or 2 columns)
RELICS_genomeScoring_wrapper <- function(input.list, input.specs){
  out.list <- list()
  list.names <- names(input.list)

  for(i in 1:length(input.list)){
    if(! list.names[i] %in% c('viterbi', 'frwdBkwd', 'dirichletMultinom', 'nbGlmmUntr')){
      temp.list.name <- paste(list.names[i], 'genomeScores', sep = '_')
      out.list[[temp.list.name]] <- RELICS_genomeScoring(input.list[[i]], input.specs, list.names[i])
    } else {
      temp.list.name <- list.names[i]
      out.list[[temp.list.name]] <- input.list[[i]]
    }
  }

  return(out.list)
}

# slide window along scores and claculate sum or average, depending on type of analysis
# input.window, input.ll.df, input.label.hierarchy, input.gene, input.pdf.name
RELICS_genomeScoring <- function(input.df, input.specs, analysis.name){

  # generate a window of effect for each guide
  if(input.specs$crisprSystem %in% c('CRISPRcas9', 'CRISPRi', 'CRISPRa')){
    input.df$start <- input.df$start - round(input.specs$crisprEffectRange / 2)
    input.df$end <- input.df$end + round(input.specs$crisprEffectRange / 2)
  }

  # need to sort the data according to start of guide
  targeting.window.scores <- c()

  #separate out chromsome from non-chromosome guides
  targeting.df <- c()  # data frame containing only valid chromosomes
  avg.overlaps <- c()
  non.targeting.df <- c()
  if('NA' %in% input.df$chrom){
    targeting.df <- filter(input.df, grepl('chr',chrom))
    non.targeting.df <- filter(input.df, grepl('NA',chrom))
  } else {
    targeting.df <- input.df
  }

  scoring.type <- c()

  if(analysis.name %in% c('RELICS','nbGlmm','nbGlmmDev', 'allDataNbGlmmCategoryParallel',
    'allDataNbGlmm', 'numIntPoisGlmm', 'uniPoisGlmm', 'multiPoisGlmm', 'allDataMultiPoisGlmmPerRepl',
    'allDataMultiPoisGlmmAcrossRepl', 'reducedDataMultiPoisGlmmPerRepl', 'reducedDataMultiPoisGlmmAcrossRepl')){
      scoring.type <- 1
    } else if(analysis.name %in% 'FoldChange'){
      scoring.type <- 2
    } else { # DESeq2, edgeR, ... p-value based
      scoring.type <- 3
    }

  # score.regions <- analysis.name %in% c('RELICS','nbGlmm','nbGlmmDev', 'allDataNbGlmmCategoryParallel',
  #   'allDataNbGlmm', 'numIntPoisGlmm', 'uniPoisGlmm', 'multiPoisGlmm', 'allDataMultiPoisGlmmPerRepl',
  #   'allDataMultiPoisGlmmAcrossRepl', 'reducedDataMultiPoisGlmmPerRepl', 'reducedDataMultiPoisGlmmAcrossRepl')

  if(length(unique(targeting.df$chrom)) == 1){
    sorted.targeting.chrom <- targeting.df[order(targeting.df$start),]
    # temp.score.list <- score_regions(sorted.df, input.specs$labelHierarchy, score.regions)
    # targeting.window.scores <- temp.score.list$out_df
    # targeting.window.scores <- score_regions(sorted.targeting.chrom,
    #   input.specs$labelHierarchy, score.regions)
    targeting.window.scores <- score_regions(sorted.targeting.chrom,
      input.specs$labelHierarchy, scoring.type)
    # avg.overlaps <- temp.score.list$avg_overlaps
  } else {
    all.chroms <- unique(targeting.df$chrom)
    for(i in 1:length(all.chroms)){
      targeting.chrom <- targeting.df[(which(targeting.df$chrom == all.chroms[i])),]
      sorted.targeting.chrom <- targeting.chrom[order(targeting.chrom$start),]
      # temp.score.list <- score_regions(sorted.targeting.chrom,
      #   input.specs$labelHierarchy, score.regions)
      # temp.targeting.window.scores <- temp.score.list$out_df
      # targeting.window.scores <- score_regions(sorted.targeting.chrom,
      #   input.specs$labelHierarchy, score.regions)
      targeting.window.scores <- score_regions(sorted.targeting.chrom,
        input.specs$labelHierarchy, scoring.type)
      # avg.overlaps <- c(avg.overlaps, temp.score.list$avg_overlaps)
      targeting.window.scores <- rbind(targeting.window.scores, temp.targeting.window.scores)
    }
  }
  if(nrow(targeting.df) != nrow(input.df)){
    avg.overlap <- round(median(targeting.window.scores$nrSupportGuides))
    # avg.overlap <- round(mean(avg.overlaps))
    non.targeting.window.scores <- score_nonTarget_regions(scoring.type,
      non.targeting.df, avg.overlap, non.targeting.df$label[1])
    targeting.window.scores <- rbind(targeting.window.scores, non.targeting.window.scores)
  }

  return(targeting.window.scores)
}

score_nonTarget_regions <- function(scoringType, input.df, avg.overlap, input.label){

  row.seq <- seq(1, nrow(input.df), by = avg.overlap)
  out.raw.scores <- vector('numeric', length = length(row.seq))
  out.format.scores <- vector('numeric', length = length(row.seq))
  out.log2.fc <- vector('numeric', length = length(row.seq))

  overlapping <- c(0:(avg.overlap - 1))
  for(i in 1:(length(row.seq) - 1)){
    temp.rows <- row.seq[i] + overlapping
    if(scoringType == 1){ # llr based
      out.raw.scores[i] <- sum(input.df$formatScores[temp.rows])
      out.format.scores[i] <- sum(input.df$formatScores[temp.rows])
      out.log2.fc[i] <- 1
    } else if(scoringType == 2){
      out.raw.scores[i] <- mean(input.df$formatScores[temp.rows])
      out.format.scores[i] <- mean(input.df$formatScores[temp.rows])
      out.log2.fc[i] <- mean(input.df$formatScores[temp.rows])
    } else {
      temp.df <- input.df[temp.rows,]
      avg.sign <- mean(temp.df$log2_rate_ratio)
      temp.pval <- pchisq(-2 * sum(log(temp.df$rawScores)), df = 2 * nrow(temp.df), lower.tail=F)

      out.raw.scores[i] <- temp.pval
      out.format.scores[i] <- -log10(temp.pval) * sign(avg.sign)
      out.log2.fc[i] <- avg.sign
    }
  }
  last.index <- length(row.seq)
  temp.rows <- c(row.seq[last.index]:nrow(input.df))
  if(scoringType == 1){ # llr based
    out.raw.scores[last.index] <- sum(input.df$formatScores[temp.rows])
    out.format.scores[last.index] <- sum(input.df$formatScores[temp.rows])
    out.log2.fc[last.index] <- 1
  } else if(scoringType == 2){
    out.raw.scores[last.index] <- mean(input.df$formatScores[temp.rows])
    out.format.scores[last.index] <- mean(input.df$formatScores[temp.rows])
    out.log2.fc[last.index] <- mean(input.df$formatScores[temp.rows])
  } else {
    temp.df <- input.df[temp.rows,]
    avg.sign <- mean(temp.df$log2_rate_ratio)
    temp.pval <- pchisq(-2 * sum(log(temp.df$rawScores)), df = 2 * nrow(temp.df), lower.tail=F)

    out.raw.scores[last.index] <- temp.pval
    out.format.scores[last.index] <- -log10(temp.pval) * sign(avg.sign)
    out.log2.fc[last.index] <- avg.sign
  }

  out.df <- data.frame(rawScores = out.raw.scores, formatScores = out.format.scores,
    log2_rate_ratio = rep(1, length(out.raw.scores)), chrom = rep('NA', length(out.raw.scores)),
    label = rep(input.label, length(out.raw.scores)), start = rep(NA, length(out.raw.scores)),
    end = rep(NA, length(out.raw.scores)), nrSupportGuides = c(rep(avg.overlap,
      (length(out.raw.scores)-1)), length(temp.rows)), stringsAsFactors = F)

  return(out.df)
}

score_regions <- function(input.df, input.label.hierarchy, scoringType){

  all.breaks <- unique(c(input.df$start, input.df$end))
  sorted.breaks <- sort(all.breaks)
  all.region.breaks <- sorted.breaks
  start.regions <- all.region.breaks[c(1:(length(all.region.breaks) - 1))]
  end.regions <- all.region.breaks[c(2:length(all.region.breaks))] - 1

  region.ranges <- GRanges(seqnames = rep(input.df$chrom[1], length(start.regions)),
    ranges = IRanges(start.regions, end.regions))

  # adjust score ranges such that they fall within, not onto, boarders of genomic ranges
  score.ranges <- GRanges(seqnames = input.df$chrom,
    ranges = IRanges(input.df$start+1, input.df$end-1))

  score.overlaps <- as.data.frame(findOverlaps(region.ranges, score.ranges, type = 'any'))

  region.overlap.list <- split(score.overlaps, score.overlaps$queryHits)

  # all.region.scores <- vector('numeric', length = length(start.regions))
  all.region.rawScores <- vector('numeric', length = length(start.regions))
  all.region.formatScores <- vector('numeric', length = length(start.regions))
  all.region.foldChange <- vector('numeric', length = length(start.regions))

  all.region.guide.supports <- vector('numeric', length = length(start.regions))
  unique.regions <- unique(score.overlaps$queryHits)

  region.labels <- unlist(lapply(region.overlap.list, function(x){
    temp.label <- input.df$label[x$subjectHits]
    temp.labels.present <- which(input.label.hierarchy %in% temp.label)
    out.label <- input.label.hierarchy[max(temp.labels.present)]
    return(out.label)
    }))

  avg.region.overlaps <- unlist(lapply(region.overlap.list, function(x){
    temp.overlaps <- length(input.df$formatScores[x$subjectHits])
    return(temp.overlaps)
    }))

  # if the input is based on bayes factor, add the scores, else take the mean
  region.scores <- c()
  # if(scoringType == 1) {
  #   region.scores <- unlist(lapply(region.overlap.list, function(x){
  #     temp.scores <- input.df$formatScores[x$subjectHits]
  #     out.score <- sum(temp.scores)
  #     return(out.score)
  #     }))
  # } else if(scoringType == 2){ # fold change
  #   region.scores <- unlist(lapply(region.overlap.list, function(x){
  #     temp.df <- input.df[x$subjectHits,]
  #     out.score <- mean(temp.df$formatScores)
  #     return(out.score)
  #     }))
  # } else { # p-value based result
  #   region.scores <- unlist(lapply(region.overlap.list, function(x){
  #     temp.df <- input.df[x$subjectHits,]
  #     avg.sign <- mean(temp.df$log2_rate_ratio)
  #     temp.pval <- pchisq(-2 * sum(log(temp.df$rawScores)), df = 2 * nrow(temp.df), lower.tail=F)
  #     out.pval <- -log10(temp.pval) * sign(avg.sign)
  #     return(out.pval)
  #     }))
  # }

  if(scoringType == 1) {
    region.scores <- do.call(rbind, lapply(region.overlap.list, function(x){
      temp.scores <- input.df$formatScores[x$subjectHits]
      out.score <- sum(temp.scores)
      out.values <- c(out.score, out.score, rep(1, length(out.score)) )
      return(out.values)
      }))
  } else if(scoringType == 2){ # fold change
    region.scores <- do.call(rbind, lapply(region.overlap.list, function(x){
      temp.df <- input.df[x$subjectHits,]
      out.score <- mean(temp.df$formatScores)
      out.values <- c(out.score, out.score, out.score)
      return(out.values)
      }))
  } else { # p-value based result, use Fisher's method to combine
    region.scores <- do.call(rbind, lapply(region.overlap.list, function(x){
      temp.df <- input.df[x$subjectHits,]
      avg.sign <- mean(temp.df$log2_rate_ratio)
      temp.pval <- pchisq(-2 * sum(log(temp.df$rawScores)), df = 2 * nrow(temp.df), lower.tail=F)
      out.pval <- -log10(temp.pval) * sign(avg.sign)
      out.values <- c(temp.pval, out.pval, avg.sign)
      return(out.values)
      }))
  }

  # all.region.scores[unique.regions] <- region.scores

  all.region.rawScores[unique.regions] <- region.scores[,1]
  all.region.formatScores[unique.regions] <- region.scores[,2]
  all.region.foldChange[unique.regions] <- region.scores[,3]

  all.region.labels <- rep('chr', length(all.region.rawScores))
  all.region.labels[unique.regions] <- region.labels
  all.region.guide.supports[unique.regions] <- avg.region.overlaps

  # out.df <- data.frame(rawScores = all.region.scores, formatScores = all.region.scores,
  #   log2_rate_ratio = rep(1, length(all.region.scores)),
  #   chrom = rep(input.df$chrom[1], length(all.region.scores)),
  #   label = all.region.labels, start = start.regions, end = end.regions,
  #   nrSupportGuides = all.region.guide.supports, stringsAsFactors = F)

  out.df <- data.frame(rawScores = all.region.rawScores, formatScores = all.region.formatScores,
    log2_rate_ratio = all.region.foldChange,
    chrom = rep(input.df$chrom[1], length(all.region.rawScores)),
    label = all.region.labels, start = start.regions, end = end.regions,
    nrSupportGuides = all.region.guide.supports, stringsAsFactors = F)

  return(out.df)
  # return(list(out_df = out.df, avg_overlaps = avg.region.overlaps))
}

RELICS_genomeScoring_weighted <- function(input.df, input.specs, analysis.name){

  # generate a window of effect for each guide
  if(input.specs$crisprSystem %in% c('CRISPRcas9', 'CRISPRi', 'CRISPRa')){
    input.df$start <- input.df$start - round(input.specs$crisprEffectRange / 2)
    input.df$end <- input.df$end + round(input.specs$crisprEffectRange / 2)
  }

  # need to sort the data according to start of guide
  targeting.window.scores <- c()

  #separate out chromsome from non-chromosome guides
  targeting.df <- c()  # data frame containing only valid chromosomes
  avg.overlaps <- c()
  non.targeting.df <- c()
  if('NA' %in% input.df$chrom){
    targeting.df <- filter(input.df, grepl('chr',chrom))
    non.targeting.df <- filter(input.df, grepl('NA',chrom))
  } else {
    targeting.df <- input.df
  }

  scoring.type <- c()

  if(analysis.name %in% c('RELICS','nbGlmm','nbGlmmDev', 'allDataNbGlmmCategoryParallel',
    'allDataNbGlmm', 'numIntPoisGlmm', 'uniPoisGlmm', 'multiPoisGlmm', 'allDataMultiPoisGlmmPerRepl',
    'allDataMultiPoisGlmmAcrossRepl', 'reducedDataMultiPoisGlmmPerRepl', 'reducedDataMultiPoisGlmmAcrossRepl')){
      scoring.type <- 1
    } else if(analysis.name %in% 'FoldChange'){
      scoring.type <- 2
    } else { # DESeq2, edgeR, ... p-value based
      scoring.type <- 3
    }

  # score.regions <- analysis.name %in% c('RELICS','nbGlmm','nbGlmmDev', 'allDataNbGlmmCategoryParallel',
  #   'allDataNbGlmm', 'numIntPoisGlmm', 'uniPoisGlmm', 'multiPoisGlmm', 'allDataMultiPoisGlmmPerRepl',
  #   'allDataMultiPoisGlmmAcrossRepl', 'reducedDataMultiPoisGlmmPerRepl', 'reducedDataMultiPoisGlmmAcrossRepl')

  if(length(unique(targeting.df$chrom)) == 1){
    sorted.targeting.chrom <- targeting.df[order(targeting.df$start),]
    # temp.score.list <- score_regions(sorted.df, input.specs$labelHierarchy, score.regions)
    # targeting.window.scores <- temp.score.list$out_df
    # targeting.window.scores <- score_regions(sorted.targeting.chrom,
    #   input.specs$labelHierarchy, score.regions)
    targeting.window.scores <- score_regions_weighted(sorted.targeting.chrom,
      input.specs$labelHierarchy, scoring.type)
    # avg.overlaps <- temp.score.list$avg_overlaps
  } else {
    all.chroms <- unique(targeting.df$chrom)
    for(i in 1:length(all.chroms)){
      targeting.chrom <- targeting.df[(which(targeting.df$chrom == all.chroms[i])),]
      sorted.targeting.chrom <- targeting.chrom[order(targeting.chrom$start),]
      # temp.score.list <- score_regions(sorted.targeting.chrom,
      #   input.specs$labelHierarchy, score.regions)
      # temp.targeting.window.scores <- temp.score.list$out_df
      # targeting.window.scores <- score_regions(sorted.targeting.chrom,
      #   input.specs$labelHierarchy, score.regions)
      targeting.window.scores <- score_regions_weighted(sorted.targeting.chrom,
        input.specs$labelHierarchy, scoring.type)
      # avg.overlaps <- c(avg.overlaps, temp.score.list$avg_overlaps)
      targeting.window.scores <- rbind(targeting.window.scores, temp.targeting.window.scores)
    }
  }
  if(nrow(targeting.df) != nrow(input.df)){
    avg.overlap <- round(median(targeting.window.scores$nrSupportGuides))
    # avg.overlap <- round(mean(avg.overlaps))
    non.targeting.window.scores <- score_nonTarget_regions_weighted(scoring.type,
      non.targeting.df, avg.overlap, non.targeting.df$label[1])
    targeting.window.scores <- rbind(targeting.window.scores, non.targeting.window.scores)
  }

  return(targeting.window.scores)
}

score_nonTarget_regions_weighted <- function(scoringType, input.df, avg.overlap, input.label){

  row.seq <- seq(1, nrow(input.df), by = avg.overlap)
  out.raw.scores <- vector('numeric', length = length(row.seq))
  out.format.scores <- vector('numeric', length = length(row.seq))
  out.log2.fc <- vector('numeric', length = length(row.seq))

  overlapping <- c(0:(avg.overlap - 1))
  for(i in 1:(length(row.seq) - 1)){
    temp.rows <- row.seq[i] + overlapping
    if(scoringType == 1){ # llr based
      out.raw.scores[i] <- sum(input.df$formatScores[temp.rows])
      out.format.scores[i] <- sum(input.df$formatScores[temp.rows])
      out.log2.fc[i] <- 1
    } else if(scoringType == 2){
      out.raw.scores[i] <- mean(input.df$formatScores[temp.rows])
      out.format.scores[i] <- mean(input.df$formatScores[temp.rows])
      out.log2.fc[i] <- mean(input.df$formatScores[temp.rows])
    } else {
      temp.df <- input.df[temp.rows,]
      avg.sign <- mean(temp.df$log2_rate_ratio * -log10(temp.df$rawScores))
      temp.pval <- pchisq(-2 * sum(log(temp.df$rawScores)), df = 2 * nrow(temp.df), lower.tail=F)

      out.raw.scores[i] <- temp.pval
      out.format.scores[i] <- -log10(temp.pval) * sign(avg.sign)
      out.log2.fc[i] <- avg.sign
    }
  }
  last.index <- length(row.seq)
  temp.rows <- c(row.seq[last.index]:nrow(input.df))
  if(scoringType == 1){ # llr based
    out.raw.scores[last.index] <- sum(input.df$formatScores[temp.rows])
    out.format.scores[last.index] <- sum(input.df$formatScores[temp.rows])
    out.log2.fc[last.index] <- 1
  } else if(scoringType == 2){
    out.raw.scores[last.index] <- mean(input.df$formatScores[temp.rows])
    out.format.scores[last.index] <- mean(input.df$formatScores[temp.rows])
    out.log2.fc[last.index] <- mean(input.df$formatScores[temp.rows])
  } else {
    temp.df <- input.df[temp.rows,]
    avg.sign <- mean(temp.df$log2_rate_ratio * -log10(temp.df$rawScores))
    temp.pval <- pchisq(-2 * sum(log(temp.df$rawScores)), df = 2 * nrow(temp.df), lower.tail=F)

    out.raw.scores[last.index] <- temp.pval
    out.format.scores[last.index] <- -log10(temp.pval) * sign(avg.sign)
    out.log2.fc[last.index] <- avg.sign
  }

  out.df <- data.frame(rawScores = out.raw.scores, formatScores = out.format.scores,
    log2_rate_ratio = rep(1, length(out.raw.scores)), chrom = rep('NA', length(out.raw.scores)),
    label = rep(input.label, length(out.raw.scores)), start = rep(NA, length(out.raw.scores)),
    end = rep(NA, length(out.raw.scores)), nrSupportGuides = c(rep(avg.overlap,
      (length(out.raw.scores)-1)), length(temp.rows)), stringsAsFactors = F)

  return(out.df)
}

score_regions_weighted <- function(input.df, input.label.hierarchy, scoringType){

  all.breaks <- unique(c(input.df$start, input.df$end))
  sorted.breaks <- sort(all.breaks)
  all.region.breaks <- sorted.breaks
  start.regions <- all.region.breaks[c(1:(length(all.region.breaks) - 1))]
  end.regions <- all.region.breaks[c(2:length(all.region.breaks))] - 1

  region.ranges <- GRanges(seqnames = rep(input.df$chrom[1], length(start.regions)),
    ranges = IRanges(start.regions, end.regions))

  # adjust score ranges such that they fall within, not onto, boarders of genomic ranges
  score.ranges <- GRanges(seqnames = input.df$chrom,
    ranges = IRanges(input.df$start+1, input.df$end-1))

  score.overlaps <- as.data.frame(findOverlaps(region.ranges, score.ranges, type = 'any'))

  region.overlap.list <- split(score.overlaps, score.overlaps$queryHits)

  # all.region.scores <- vector('numeric', length = length(start.regions))
  all.region.rawScores <- vector('numeric', length = length(start.regions))
  all.region.formatScores <- vector('numeric', length = length(start.regions))
  all.region.foldChange <- vector('numeric', length = length(start.regions))

  all.region.guide.supports <- vector('numeric', length = length(start.regions))
  unique.regions <- unique(score.overlaps$queryHits)

  region.labels <- unlist(lapply(region.overlap.list, function(x){
    temp.label <- input.df$label[x$subjectHits]
    temp.labels.present <- which(input.label.hierarchy %in% temp.label)
    out.label <- input.label.hierarchy[max(temp.labels.present)]
    return(out.label)
    }))

  avg.region.overlaps <- unlist(lapply(region.overlap.list, function(x){
    temp.overlaps <- length(input.df$formatScores[x$subjectHits])
    return(temp.overlaps)
    }))

  # if the input is based on bayes factor, add the scores, else take the mean
  region.scores <- c()

  if(scoringType == 1) {
    region.scores <- do.call(rbind, lapply(region.overlap.list, function(x){
      temp.scores <- input.df$formatScores[x$subjectHits]
      out.score <- sum(temp.scores)
      out.values <- c(out.score, out.score, rep(1, length(out.score)) )
      return(out.values)
      }))
  } else if(scoringType == 2){ # fold change
    region.scores <- do.call(rbind, lapply(region.overlap.list, function(x){
      temp.df <- input.df[x$subjectHits,]
      out.score <- mean(temp.df$formatScores)
      out.values <- c(out.score, out.score, out.score)
      return(out.values)
      }))
  } else { # p-value based result, use Fisher's method to combine
    region.scores <- do.call(rbind, lapply(region.overlap.list, function(x){
      temp.df <- input.df[x$subjectHits,]
      avg.sign <- mean(temp.df$log2_rate_ratio * -log10(temp.df$rawScores))
      temp.pval <- pchisq(-2 * sum(log(temp.df$rawScores)), df = 2 * nrow(temp.df), lower.tail=F)
      out.pval <- -log10(temp.pval) * sign(avg.sign)
      out.values <- c(temp.pval, out.pval, avg.sign)
      return(out.values)
      }))
  }

  # all.region.scores[unique.regions] <- region.scores

  all.region.rawScores[unique.regions] <- region.scores[,1]
  all.region.formatScores[unique.regions] <- region.scores[,2]
  all.region.foldChange[unique.regions] <- region.scores[,3]

  all.region.labels <- rep('chr', length(all.region.rawScores))
  all.region.labels[unique.regions] <- region.labels
  all.region.guide.supports[unique.regions] <- avg.region.overlaps

  # out.df <- data.frame(rawScores = all.region.scores, formatScores = all.region.scores,
  #   log2_rate_ratio = rep(1, length(all.region.scores)),
  #   chrom = rep(input.df$chrom[1], length(all.region.scores)),
  #   label = all.region.labels, start = start.regions, end = end.regions,
  #   nrSupportGuides = all.region.guide.supports, stringsAsFactors = F)

  out.df <- data.frame(rawScores = all.region.rawScores, formatScores = all.region.formatScores,
    log2_rate_ratio = all.region.foldChange,
    chrom = rep(input.df$chrom[1], length(all.region.rawScores)),
    label = all.region.labels, start = start.regions, end = end.regions,
    nrSupportGuides = all.region.guide.supports, stringsAsFactors = F)

  return(out.df)
  # return(list(out_df = out.df, avg_overlaps = avg.region.overlaps))
}

#apply sliding window to all score datasets
#input:
#   input.list: all method-specifi results to be analyzed using a specific RRA
#   input.specs: specification list
#output: list
#   rawScores, formatScores, log2_rate_ratio, chromosome, label, sgRNA targets (can be 1 or 2 columns)
slidingWindow_wrapper <- function(input.list, input.specs){
  out.list <- list()
  list.names <- names(input.list)

  for(i in 1:length(input.list)){
    if(! list.names[i] %in% c('viterbi', 'frwdBkwd')){
      temp.list.name <- paste(list.names[i], 'slidingWindow', sep = '_')
      out.list[[temp.list.name]] <- calculate_sliding_window(input.list[[i]], input.specs, list.names[i])
    } else {
      temp.list.name <- list.names[i]
      out.list[[temp.list.name]] <- input.list[[i]]
    }
  }

  return(out.list)
}

# slide window along scores and claculate sum or average, depending on type of analysis
# input.window, input.ll.df, input.label.hierarchy, input.gene, input.pdf.name
calculate_sliding_window <- function(input.df, input.specs, analysis.name){

  # need to sort the data according to start of guide
  targeting.window.scores <- c()

  #separate out chromsome from non-chromosome guides
  targeting.df <- c()  # data frame containing only valid chromosomes
  non.targeting.df <- c()
  if('NA' %in% input.df$chrom){
    targeting.df <- filter(input.df, grepl('chr',chrom))
    non.targeting.df <- filter(input.df, grepl('NA',chrom))
  } else {
    targeting.df <- input.df
  }

  if(length(unique(targeting.df$chrom)) == 1){
    sorted.df <- targeting.df[order(targeting.df$start),]
    targeting.window.scores <- c()
    if(analysis.name %in% c('nbGlmm')){
      targeting.window.scores <- slide_ll_window(input.specs$guidePerSlidingWindow, sorted.df,
        input.specs$labelHierarchy, input.specs$maxWindowSize)
    } else {
      targeting.window.scores <- slide_score_window(input.specs$guidePerSlidingWindow, sorted.df,
        input.specs$labelHierarchy, input.specs$maxWindowSize)
    }
  } else {
    all.chroms <- unique(targeting.df$chrom)
    for(i in 1:length(all.chroms)){
      targeting.chrom <- targeting.df[(which(targeting.df$chrom == all.chroms[i])),]
      sorted.targeting.chrom <- targeting.chrom[order(targeting.chrom$start),]
      if(analysis.name %in% c('nbGlmm')){
        temp.targeting.window.scores <- slide_ll_window(input.specs$guidePerSlidingWindow,
          sorted.targeting.chrom, input.specs$labelHierarchy, input.specs$maxWindowSize)
        targeting.window.scores <- rbind(targeting.window.scores, temp.targeting.window.scores)
      } else {
        temp.targeting.window.scores <- slide_score_window(input.specs$guidePerSlidingWindow,
          sorted.targeting.chrom, input.specs$labelHierarchy, input.specs$maxWindowSize)
        targeting.window.scores <- rbind(targeting.window.scores, temp.targeting.window.scores)
      }
    }
  }

  if(nrow(targeting.df) != nrow(input.df)){
    # non.targeting.df <- input.df[which(input.df$label %in% c('neg')),]
    non.targeting.sorted.df <- non.targeting.df[order(non.targeting.df$start),]

    non.targeting.window.scores <- c()
    if(analysis.name %in% c('nbGlmm')){
      non.targeting.window.scores <- slide_ll_window_nonTargeting(input.specs$guidePerSlidingWindow,
        non.targeting.sorted.df, input.specs$labelHierarchy, input.specs$maxWindowSize)
    } else {
      non.targeting.window.scores <- slide_score_window_nonTargeting(input.specs$guidePerSlidingWindow,
        non.targeting.sorted.df, input.specs$labelHierarchy, input.specs$maxWindowSize)
    }

    targeting.window.scores <- rbind(targeting.window.scores, non.targeting.window.scores)
  }

  # df.to.plot <- cbind.data.frame(chrom = rep(sorted.df$chrom[1], nrow(targeting.window.scores)),
  #   start = targeting.window.scores$start, end = targeting.window.scores$end,
  #   label = targeting.window.scores$window_label, formatScores = targeting.window.scores$window_score,
  #   x_pos = targeting.window.scores$window_pos, y_pos = targeting.window.scores$window_score, stringsAsFactors = F)
  #
  # out.score.list <- list()
  # out.score.list$windowDF <- df.to.plot

  # plot_scores_distr(df.to.plot, input.label.hierarchy, input.gene = input.gene)
  #
  # pdf(paste0('window_',input.pdf.name), width = 14, height = 7)
  # plot_scores_distr(df.to.plot, input.label.hierarchy, input.gene = input.gene)
  # dev.off()

  return(targeting.window.scores)
}

reduce_ll_window_size <- function(input.window.df, input.label.hierarchy, input.window.max,
  input.row.start, input.row.max){


  out.list <- list()
  out.list$endSliding <- FALSE

  # if the first region is already too large ...
  if(input.window.df$end[1] - input.window.df$start[1] > input.window.max){
    new.window.start <- input.row.start
    reduced.window.end <- nrow(input.window.df) # know already the first one didn not work
    for(i in 2:reduced.window.end){
      new.window.start <- new.window.start + 1
      # if the region (already too large) has changed
      if(input.window.df$end[1] != input.window.df$end[i]){
        out.list$end <- max(input.window.df$end[1:(i - 1)])
        out.list$label <- input.label.hierarchy[max(which(input.label.hierarchy %in% input.window.df$label[1:(i - 1)]))]
        out.list$score <- sum(input.window.df$formatScores[1:(i - 1)])
        out.list$newWindowStart <- new.window.start
        return(out.list)
      }
    }
    out.list$end <- max(input.window.df$end[1:(i - 1)])
    out.list$label <- input.label.hierarchy[max(which(input.label.hierarchy %in% input.window.df$label[1:(i - 1)]))]
    out.list$score <- sum(input.window.df$formatScores[1:(i - 1)])
    out.list$newWindowStart <- input.row.start + 1
    out.list$endSliding <- TRUE # if the end of the file is reached and the same region
                                # continues, combine into one window
  } else {
    new.window.start <- input.row.max
    reduced.window.end <- nrow(input.window.df) - 1 # know already the first one didn not work
    for(i in reduced.window.end:1){
      if(max(input.window.df$end[1:reduced.window.end]) - input.window.df$start[1] > input.window.max){
        new.window.start <- new.window.start - 1
        reduced.window.end <- reduced.window.end - 1
      } else {
        out.list$end <- max(input.window.df$end[1:reduced.window.end])
        out.list$label <- input.label.hierarchy[max(which(input.label.hierarchy %in% input.window.df$label[1:reduced.window.end]))]
        out.list$score <- sum(input.window.df$formatScores[1:reduced.window.end])
        out.list$newWindowStart <- new.window.start
        return(out.list)
      }
    }
  }

  return(out.list)
}

slide_ll_window <- function(input.window, input.df, input.label.hierarchy, input.window.max){
  window.size <- input.window - 1
  window.pos <- vector('numeric', length = nrow(input.df))
  window.label <- vector('numeric', length = nrow(input.df))
  window.start <- vector('numeric', length = nrow(input.df))
  window.end <- vector('numeric', length = nrow(input.df))

  window.ll <- vector(mode = 'numeric', length = nrow(input.df))
  row.position <- 1
  continue <- TRUE
  nr.rows <- nrow(input.df)
  final.row <- nr.rows
  nr.windows.generated <- 0
  #initialize the first instance

  while(continue){  # row.position + window.size <= nrow(input.df)
    nr.windows.generated <- nr.windows.generated + 1
    curr.last.row <- row.position + window.size  # current last low to be included in window
    if(curr.last.row > nr.rows){  #reaching end of file
      # if there are more, smaller windows possible toward the end
      if(max(input.df$end[row.position:final.row]) - input.df$start[row.position] > input.window.max){
        window.specs <- reduce_ll_window_size(input.df[row.position:final.row, ],
          input.label.hierarchy, input.window.max, row.position, final.row)
        window.start[nr.windows.generated] <- input.df$start[row.position]
        window.end[nr.windows.generated] <- window.specs$end
        window.label[nr.windows.generated] <- window.specs$label
        window.ll[nr.windows.generated] <- window.specs$score
        if(window.specs$endSliding){
          continue <- FALSE
        } else {
          row.position <- window.specs$newWindowStart
        }
      } else {
        window.start[nr.windows.generated] <- input.df$start[row.position]
        window.end[nr.windows.generated] <- max(input.df$end[row.position:final.row])
        window.label[nr.windows.generated] <- input.label.hierarchy[max(which(input.label.hierarchy %in% input.df$label[row.position:final.row]))]
        window.ll[nr.windows.generated] <- sum(input.df$formatScores[row.position:final.row])
        continue <- FALSE
      }
    } else {
      if(max(input.df$end[row.position:curr.last.row]) - input.df$start[row.position] > input.window.max){
        window.specs <- reduce_ll_window_size(input.df[row.position:curr.last.row, ],
          input.label.hierarchy, input.window.max, row.position, curr.last.row)
        window.start[nr.windows.generated] <- input.df$start[row.position]
        window.end[nr.windows.generated] <- window.specs$end
        window.label[nr.windows.generated] <- window.specs$label
        window.ll[nr.windows.generated] <- window.specs$score
        row.position <- window.specs$newWindowStart
      } else {
        window.start[nr.windows.generated] <- input.df$start[row.position]
        window.end[nr.windows.generated] <- max(input.df$end[row.position:curr.last.row])
        window.label[nr.windows.generated] <- input.label.hierarchy[max(which(input.label.hierarchy %in% input.df$label[row.position:curr.last.row]))]
        window.ll[nr.windows.generated] <- sum(input.df$formatScores[row.position:curr.last.row])
        row.position <- row.position + 1
      }
    }
  }

  out.df <- data.frame(rawScores = window.ll[1:nr.windows.generated],
    formatScores = window.ll[1:nr.windows.generated], log2_rate_ratio = rep(1, nr.windows.generated),
    chrom = rep(input.df$chrom[1], nr.windows.generated),
    label = window.label[1:nr.windows.generated], start = window.start[1:nr.windows.generated],
    end = window.end[1:nr.windows.generated], stringsAsFactors = F)

  return(out.df)
}

slide_ll_window_nonTargeting <- function(input.window, input.df, input.label.hierarchy, input.window.max){
  window.size <- window.size <- min(input.window - 1, nrow(input.df) - 1)
  window.pos <- vector('numeric', length = (nrow(input.df) - window.size) )
  window.label <- vector('numeric', length = (nrow(input.df) - window.size) )
  window.start <- vector('numeric', length = (nrow(input.df) - window.size) )
  window.end <- vector('numeric', length = (nrow(input.df) - window.size) )

  window.ll <- vector(mode = 'numeric', length = (nrow(input.df) - window.size))
  row.position <- 1
  while(row.position + window.size <= nrow(input.df)){
    # temp.window.min <- min(input.df$start[row.position:(row.position+window.size)])
    # temp.window.max <- max(input.df$end[row.position:(row.position+window.size)])
    window.start[row.position] <- NA #temp.window.min
    window.end[row.position] <- NA #temp.window.max
    window.label[row.position] <- input.label.hierarchy[max(which(input.label.hierarchy %in% input.df$label[row.position:(row.position+window.size)]))]
    window.ll[row.position] <- sum(input.df$formatScores[row.position:(row.position+window.size)])
    row.position <- row.position + 1
  }

  out.df <- data.frame(rawScores = window.ll, formatScores = window.ll,
    log2_rate_ratio = rep(1, length(window.ll)), chrom = rep(input.df$chrom[1], length(window.ll)),
    label = window.label, start = window.start, end = window.end, stringsAsFactors = F)

  return(out.df)
}

reduce_score_window_size <- function(input.window.df, input.label.hierarchy, input.window.max,
  input.row.start, input.row.max){


  out.list <- list()
  out.list$endSliding <- FALSE

  # if the first region is already too large ...
  if(input.window.df$end[1] - input.window.df$start[1] > input.window.max){
    new.window.start <- input.row.start
    reduced.window.end <- nrow(input.window.df) # know already the first one didn not work
    for(i in 2:reduced.window.end){
      new.window.start <- new.window.start + 1
      # if the region (already too large) has changed
      if(input.window.df$end[1] != input.window.df$end[i]){
        out.list$end <- max(input.window.df$end[1:(i - 1)])
        out.list$label <- input.label.hierarchy[max(which(input.label.hierarchy %in% input.window.df$label[1:(i - 1)]))]
        out.list$score <- mean(input.window.df$formatScores[1:(i - 1)])
        out.list$newWindowStart <- new.window.start
        return(out.list)
      }
    }
    out.list$end <- max(input.window.df$end[1:(i - 1)])
    out.list$label <- input.label.hierarchy[max(which(input.label.hierarchy %in% input.window.df$label[1:(i - 1)]))]
    out.list$score <- mean(input.window.df$formatScores[1:(i - 1)])
    out.list$newWindowStart <- input.row.start + 1
    out.list$endSliding <- TRUE # if the end of the file is reached and the same region
                                # continues, combine into one window
  } else {
    new.window.start <- input.row.max
    reduced.window.end <- nrow(input.window.df) - 1 # know already the first one didn not work
    for(i in reduced.window.end:1){
      if(max(input.window.df$end[1:reduced.window.end]) - input.window.df$start[1] > input.window.max){
        new.window.start <- new.window.start - 1
        reduced.window.end <- reduced.window.end - 1
      } else {
        out.list$end <- max(input.window.df$end[1:reduced.window.end])
        out.list$label <- input.label.hierarchy[max(which(input.label.hierarchy %in% input.window.df$label[1:reduced.window.end]))]
        out.list$score <- mean(input.window.df$formatScores[1:reduced.window.end])
        out.list$newWindowStart <- new.window.start
        return(out.list)
      }
    }
  }

  return(out.list)
}

slide_score_window <- function(input.window, input.df, input.label.hierarchy, input.window.max){
  window.size <- input.window - 1
  window.pos <- vector('numeric', length = nrow(input.df))
  window.label <- vector('numeric', length = nrow(input.df))
  window.start <- vector('numeric', length = nrow(input.df))
  window.end <- vector('numeric', length = nrow(input.df))

  window.ll <- vector(mode = 'numeric', length = nrow(input.df))
  row.position <- 1
  continue <- TRUE
  nr.rows <- nrow(input.df)
  final.row <- nr.rows
  nr.windows.generated <- 0
  #initialize the first instance

  while(continue){  # row.position + window.size <= nrow(input.df)
    nr.windows.generated <- nr.windows.generated + 1
    curr.last.row <- row.position + window.size  # current last low to be included in window
    if(curr.last.row > nr.rows){  #reaching end of file

      # if there are more, smaller windows possible toward the end
      if(max(input.df$end[row.position:final.row]) - input.df$start[row.position] > input.window.max){
        window.specs <- reduce_score_window_size(input.df[row.position:final.row, ],
          input.label.hierarchy, input.window.max, row.position, final.row)
        window.start[nr.windows.generated] <- input.df$start[row.position]
        window.end[nr.windows.generated] <- window.specs$end
        window.label[nr.windows.generated] <- window.specs$label
        window.ll[nr.windows.generated] <- window.specs$score
        if(window.specs$endSliding){
          continue <- FALSE
        } else {
          row.position <- window.specs$newWindowStart
        }
      } else {
        window.start[nr.windows.generated] <- input.df$start[row.position]
        window.end[nr.windows.generated] <- max(input.df$end[row.position:final.row])
        window.label[nr.windows.generated] <- input.label.hierarchy[max(which(input.label.hierarchy %in% input.df$label[row.position:final.row]))]
        window.ll[nr.windows.generated] <- mean(input.df$formatScores[row.position:final.row])
        continue <- FALSE
      }
    } else {
      if(max(input.df$end[row.position:curr.last.row]) - input.df$start[row.position] > input.window.max){
        window.specs <- reduce_score_window_size(input.df[row.position:curr.last.row, ],
          input.label.hierarchy, input.window.max, row.position, curr.last.row)
        window.start[nr.windows.generated] <- input.df$start[row.position]
        window.end[nr.windows.generated] <- window.specs$end
        window.label[nr.windows.generated] <- window.specs$label
        window.ll[nr.windows.generated] <- window.specs$score
        row.position <- window.specs$newWindowStart
      } else {
        window.start[nr.windows.generated] <- input.df$start[row.position]
        window.end[nr.windows.generated] <- max(input.df$end[row.position:curr.last.row])
        window.label[nr.windows.generated] <- input.label.hierarchy[max(which(input.label.hierarchy %in% input.df$label[row.position:curr.last.row]))]
        window.ll[nr.windows.generated] <- mean(input.df$formatScores[row.position:curr.last.row])
        row.position <- row.position + 1
      }
    }
  }

  out.df <- data.frame(rawScores = window.ll[1:nr.windows.generated],
    formatScores = window.ll[1:nr.windows.generated], log2_rate_ratio = rep(1, nr.windows.generated),
    chrom = rep(input.df$chrom[1], nr.windows.generated),
    label = window.label[1:nr.windows.generated], start = window.start[1:nr.windows.generated],
    end = window.end[1:nr.windows.generated], stringsAsFactors = F)

  return(out.df)
}

slide_score_window_nonTargeting <- function(input.window, input.df, input.label.hierarchy, input.window.max){
  window.size <- window.size <- min(input.window - 1, nrow(input.df) - 1)
  window.pos <- vector('numeric', length = (nrow(input.df) - window.size) )
  window.label <- vector('numeric', length = (nrow(input.df) - window.size) )
  window.start <- vector('numeric', length = (nrow(input.df) - window.size) )
  window.end <- vector('numeric', length = (nrow(input.df) - window.size) )

  window.ll <- vector(mode = 'numeric', length = (nrow(input.df) - window.size))
  row.position <- 1
  while(row.position + window.size <= nrow(input.df)){
    # temp.window.min <- min(input.df$start[row.position:(row.position+window.size)])
    # temp.window.max <- max(input.df$end[row.position:(row.position+window.size)])
    window.start[row.position] <- NA #temp.window.min
    window.end[row.position] <- NA #temp.window.max
    window.label[row.position] <- input.label.hierarchy[max(which(input.label.hierarchy %in% input.df$label[row.position:(row.position+window.size)]))]
    window.ll[row.position] <- mean(input.df$formatScores[row.position:(row.position+window.size)])
    row.position <- row.position + 1
  }

  out.df <- data.frame(rawScores = window.ll, formatScores = window.ll,
    log2_rate_ratio = rep(1, length(window.ll)), chrom = rep(input.df$chrom[1], length(window.ll)),
    label = window.label, start = window.start, end = window.end, stringsAsFactors = F)

  return(out.df)
}

# incorrect version:
# can't handle spread out, spares regions at the end on the data frame
slide_ll_window_incorrect <- function(input.window, input.df, input.label.hierarchy, input.window.max){
  window.size <- min(input.window - 1, nrow(input.df) - 1)
  window.pos <- vector('numeric', length = (nrow(input.df) - window.size) )
  window.label <- vector('numeric', length = (nrow(input.df) - window.size) )
  window.start <- vector('numeric', length = (nrow(input.df) - window.size) )
  window.end <- vector('numeric', length = (nrow(input.df) - window.size) )

  window.ll <- vector(mode = 'numeric', length = (nrow(input.df) - window.size))
  row.position <- 1
  while(row.position + window.size <= nrow(input.df)){
    # check if window size within acceptable boundary
    temp.window.min <- min(input.df$start[row.position:(row.position+window.size)])
    temp.window.max <- max(input.df$end[row.position:(row.position+window.size)])
    window.start[row.position] <- temp.window.min
    window.end[row.position] <- temp.window.max
    window.label[row.position] <- input.label.hierarchy[max(which(input.label.hierarchy %in% input.df$label[row.position:(row.position+window.size)]))]
    window.ll[row.position] <- sum(input.df$formatScores[row.position:(row.position+window.size)])
    row.position <- row.position + 1
  }

  initial.out.df <- data.frame(rawScores = window.ll, formatScores = window.ll,
    log2_rate_ratio = rep(1, length(window.ll)), chrom = rep(input.df$chrom[1], length(window.ll)),
    label = window.label, start = window.start, end = window.end, stringsAsFactors = F)

  out.df <- c()
  if(length(which(initial.out.df$end - initial.out.df$start > input.window.max))){
    out.df <- initial.out.df[-which(initial.out.df$end - initial.out.df$start > input.window.max),]
  } else {
    out.df <- initial.out.df
  }

  return(out.df)
}

# incorrect version:
# can't handle spread out, spares regions at the end on the data frame
slide_score_window_incorrect <- function(input.window, input.df, input.label.hierarchy, input.window.max){
  window.size <- window.size <- min(input.window - 1, nrow(input.df) - 1)
  window.pos <- vector('numeric', length = (nrow(input.df) - window.size) )
  window.label <- vector('numeric', length = (nrow(input.df) - window.size) )
  window.start <- vector('numeric', length = (nrow(input.df) - window.size) )
  window.end <- vector('numeric', length = (nrow(input.df) - window.size) )

  window.ll <- vector(mode = 'numeric', length = (nrow(input.df) - window.size))
  row.position <- 1
  while(row.position + window.size <= nrow(input.df)){
    temp.window.min <- min(input.df$start[row.position:(row.position+window.size)])
    temp.window.max <- max(input.df$end[row.position:(row.position+window.size)])
    window.start[row.position] <- temp.window.min
    window.end[row.position] <- temp.window.max
    window.label[row.position] <- input.label.hierarchy[max(which(input.label.hierarchy %in% input.df$label[row.position:(row.position+window.size)]))]
    window.ll[row.position] <- mean(input.df$formatScores[row.position:(row.position+window.size)])
    row.position <- row.position + 1
  }

  initial.out.df <- data.frame(rawScores = window.ll, formatScores = window.ll,
    log2_rate_ratio = rep(1, length(window.ll)), chrom = rep(input.df$chrom[1], length(window.ll)),
    label = window.label, start = window.start, end = window.end, stringsAsFactors = F)

  out.df <- c()
  if(length(which(initial.out.df$end - initial.out.df$start > input.window.max))){
    out.df <- initial.out.df[-which(initial.out.df$end - initial.out.df$start > input.window.max),]
  } else {
    out.df <- initial.out.df
  }

  return(out.df)
}

# group the bins as fixed guide blocks, non-overlaping
# this way no score gets counted twice
# apply own RRA to all score datasets
#input:
#   input.list: all method-specifi results to be analyzed using a specific RRA
#   input.specs: specification list
#output: list
#   rawScores, formatScores, log2_rate_ratio, chromosome, label, sgRNA targets (can be 1 or 2 columns)
guideWindowRRA_wrapper <- function(input.list, input.specs){
  out.list <- list()
  list.names <- names(input.list)

  for(i in 1:length(input.list)){
    if(! list.names[i] %in% c('viterbi', 'frwdBkwd')){
      temp.list.name <- paste(list.names[i], 'guideWindowRRA', sep = '_')
      print(paste('Running guide-window RRA for:', list.names[i], sep = ' '))
      out.list[[temp.list.name]] <- guideWindowRRA(input.list[[i]], input.specs)
    } else {
      temp.list.name <- list.names[i]
      out.list[[temp.list.name]] <- input.list[[i]]
    }
  }
  print('Finished RRA permutations!')
  return(out.list)
}

# this version returns 1 as value for the bin if no siginifcant guide is present
# input:
#   input.df: contains resutls on which rra should be perfomred: rawScores, formatScores, log2_rate_ratio, chromosome, labels, sgRNA targets (can be 1 or 2 columns)
#   input.specs: specifcation list
# output: data frame; rawScores, formatScores, log2_rate_ratio, chromosome, label, sgRNA targets (can be 1 or 2 columns)
#   rawScores: permutation p-Values
#   formatScores: Benjamini-Hochberg adjusted p-values, with -log10 formating
guideWindowRRA <- function(input.df, input.specs){

  # rank according to the formatted scores, and normalize by total number of guides
  rank.df <- input.df
  rank.df$rank <- rank(-rank.df$formatScores)
  rank.df$normRank <- rank.df$rank/nrow(rank.df)

  #separate out chromsome from non-chromosome guides
  na.df <- c()  # data frame containing only guides not targeting chromosome
  all.chrom.df <- c()  # data frame containing only valid chromosomes
  if('NA' %in% rank.df$chrom){
    all.chrom.df <- filter(rank.df, grepl('chr',chrom))
    na.df <- filter(rank.df, grepl('NA',chrom))
  } else {
    all.chrom.df <- rank.df
  }

  all.chrom.ranges <- c()
  all.chrom.labels <- c()
  all.chroms <- unique(all.chrom.df$chrom)
  for(chrom in all.chroms){
    chrom.df <- all.chrom.df[all.chrom.df$chrom == chrom,]
    sorted.chrom.df <- chrom.df[order(chrom.df$start),]

    sorted.chrom.df$guideWindow <- rep(paste0(chrom,'_', c(1:ceiling(nrow(sorted.chrom.df) / input.specs$guidePerFixedWindow))),
      each = input.specs$guidePerFixedWindow)[1:nrow(sorted.chrom.df)]  # adding the block label for splitting

    all.chrom.ranges <- rbind(all.chrom.ranges, sorted.chrom.df)

  }

  if(! is.null(na.df)){
    all.na.lables <- unique(na.df[,5])
    for(lab in all.na.lables){
      lab.na.df <- na.df[na.df$label == lab,]
      lab.na.df$guideWindow <- rep(paste0(lab, c(1:ceiling(nrow(lab.na.df) / input.specs$guidePerFixedWindow))),
        each = input.specs$guidePerFixedWindow)[1:nrow(lab.na.df)]
      all.chrom.ranges <- rbind(all.chrom.ranges, lab.na.df)
    }
  }

  #what is the highest/worst possible rank given threshold
  score.threshold.cutoff <- -log10(input.specs$scoreThresh)
  all.chrom.ranges.scores <- all.chrom.ranges$formatScores
  all.chrom.ranges.norm.ranks <- rank(-all.chrom.ranges.scores)/length(all.chrom.ranges.scores)
  all.chrom.ranges$norm_ranks <- all.chrom.ranges.norm.ranks
  all.chrom.ranges$score_ID <- c(1:nrow(all.chrom.ranges))
  rank.thresh <- max(all.chrom.ranges.norm.ranks[which(all.chrom.ranges.scores > score.threshold.cutoff)])

  # create a list where each element contains the indeces for the tiles which overlap the bin
  bin.tile.list <- split(all.chrom.ranges, all.chrom.ranges$guideWindow)

  # obtain the rho scores according to the cutoff
  bin.rho.scores <- unlist(lapply(bin.tile.list, function(x){
    # temp.ranks <- all.chrom.ranges.norm.ranks[x$subjectHits]
    temp.ranks <- x$norm_ranks
    signif.ranks <- temp.ranks[which(temp.ranks <= rank.thresh)]
    all.len <- length(temp.ranks)
    signif.len <- length(signif.ranks)
    if(signif.len > 0){
      return(min(pbeta(sort(signif.ranks), seq(1, signif.len, 1), seq(all.len, all.len - signif.len + 1, -1))))
    } else {
      return(1)
    }
    }))

  bin.info.list <- lapply(bin.tile.list, function(x){
    # temp.ranks <- all.chrom.ranges.norm.ranks[x$subjectHits]
    bin.start <- min(x$start)
    bin.end <- max(x$end)
    bin.chr <- x$chrom[1]
    bin.label <- input.specs$labelHierarchy[max(which(input.specs$labelHierarchy %in% x$label))]
    return(data.frame(chrom = bin.chr, start = bin.start, end = bin.end, label = bin.label, stringsAsFactors = F))
    })

  bin.info <- do.call(rbind.data.frame, bin.info.list)
  rownames(bin.info) <- c()

  permut.df <- data.frame(matrix(NA, nrow = length(bin.tile.list), ncol = input.specs$rraPermutNr))
  ptm <- proc.time()
  for(i in seq_len(input.specs$rraPermutNr)){
    rand.ranks <- sample(all.chrom.ranges.norm.ranks, size = length(all.chrom.ranges.norm.ranks), replace = FALSE)
    permut.df[[i]] <- unlist(lapply(bin.tile.list, function(x){
      temp.ranks <- rand.ranks[x$score_ID]
      signif.ranks <- temp.ranks[which(temp.ranks <= rank.thresh)]
      all.len <- length(temp.ranks)
      signif.len <- length(signif.ranks)
      if(signif.len > 0){
        return(min(pbeta(sort(signif.ranks), seq(1, signif.len, 1), seq(all.len, all.len - signif.len + 1, -1))))
      } else{
        return(1)
      }
      }))
  }

  permut.time <- proc.time() - ptm

  null.pvals <- unlist(permut.df)
  sorted.null.pvals <- sort(null.pvals)

  permut.pvals <- (findInterval(bin.rho.scores, sorted.null.pvals) + 1) / length(sorted.null.pvals)

  adj.pvals <- p.adjust(permut.pvals, method = 'BH')

  rra.out.df <- cbind.data.frame(rawScores = adj.pvals, formatScores = -log10(permut.pvals),
    log2_rate_ratio = rep(1, length(adj.pvals)),
    chrom = bin.info$chrom, label = bin.info$label, start = bin.info$start,
    end = bin.info$end, stringsAsFactors = F)

  return(rra.out.df)
}

#apply own RRA to all score datasets
#input:
#   input.list: all method-specifi results to be analyzed using a specific RRA
#   input.specs: specification list
#output: list
#   rawScores, formatScores, log2_rate_ratio, chromosome, label, sgRNA targets (can be 1 or 2 columns)
alphaRRA_wrapper <- function(input.list, input.specs){
  out.list <- list()
  list.names <- names(input.list)

  for(i in 1:length(input.list)){
    if(! list.names[i] %in% c('viterbi', 'frwdBkwd')){
      temp.list.name <- paste(list.names[i], 'alphaRRA', sep = '_')
      print(paste('Running alpha RRA for:', list.names[i], sep = ' '))
      out.list[[temp.list.name]] <- alphaRRA(input.list[[i]], input.specs)
    } else {
      temp.list.name <- list.names[i]
      out.list[[temp.list.name]] <- input.list[[i]]
    }
  }
  print('Finished RRA permutations!')
  return(out.list)
}

# this version returns 1 as value for the bin if no siginifcant guide is present
# input:
#   input.df: contains resutls on which rra should be perfomred: rawScores, formatScores, log2_rate_ratio, chromosome, labels, sgRNA targets (can be 1 or 2 columns)
#   input.specs: specifcation list
# output: data frame; rawScores, formatScores, log2_rate_ratio, chromosome, label, sgRNA targets (can be 1 or 2 columns)
#   rawScores: permutation p-Values
#   formatScores: Benjamini-Hochberg adjusted p-values, with -log10 formating
alphaRRA <- function(input.df, input.specs){

  # rank according to the formatted scores, and normalize by total number of guides
  rank.df <- input.df
  rank.df$rank <- rank(-rank.df$formatScores)
  rank.df$normRank <- rank.df$rank/nrow(rank.df)

  #separate out chromsome from non-chromosome guides
  na.df <- c()  # data frame containing only guides not targeting chromosome
  all.chrom.df <- c()  # data frame containing only valid chromosomes
  if('NA' %in% rank.df$chrom){
    all.chrom.df <- filter(rank.df, grepl('chr',chrom))
    na.df <- filter(rank.df, grepl('NA',chrom))
  } else {
    all.chrom.df <- rank.df
  }

  # adjust info, if single guide, add additional column with identical values
  all.chrom.info <- adjustInfo(all.chrom.df)

  all.chrom.ranges <- GRanges()
  all.bin.ranges <- GRanges()
  all.chrom.labels <- c()
  all.chroms <- unique(all.chrom.df$chrom)
  for(chrom in all.chroms){
    chrom.df <- all.chrom.df[all.chrom.df$chrom == chrom,]
    chrom.info <- all.chrom.info[all.chrom.info$chrom == chrom,]
    range.list <- obtain_data_GRanges(chrom.info, input.specs, -chrom.df$formatScores)  # keep negative sign!
    chrom.ranges <- range.list$gRanges
    all.chrom.ranges <- suppressWarnings(c(all.chrom.ranges, chrom.ranges)) # formerly 'append'
    bin.ranges <- obtain_bin_GRanges(range.list$rangeStartEnd, chrom, input.specs)
    all.bin.ranges <- suppressWarnings(c(all.bin.ranges, bin.ranges))
    chrom.labels <- obtain_bin_labels(bin.ranges, chrom.info, input.specs, chrom)
    all.chrom.labels <- c(all.chrom.labels,chrom.labels)
  }

  granges.bin.info <- as.data.frame(all.bin.ranges)
  all.bin.info <- cbind.data.frame(chrom = as.character(granges.bin.info$seqnames), label = all.chrom.labels,
    start = granges.bin.info$start, end = granges.bin.info$end, stringsAsFactors = F)

  if(! is.null(na.df)){
    all.na.lables <- unique(na.df[,5])
    for(lab in all.na.lables){
      lab.na.df <- na.df[na.df$label == lab,]
      na.bin.list <- rra_na_binGeneration(input.specs$rraNaBin, nrow(lab.na.df), input.specs)

      na.lab.info <- cbind.data.frame(chrom = lab.na.df$chrom, label = lab.na.df$label,
        start = na.bin.list$naRanges, end = na.bin.list$naRanges + input.specs$CRISPReffectRange,
        stringsAsFactors = F)

      na.lab.range <- GRanges(seqnames = paste(na.lab.info[,1],lab,sep = '_'),
          ranges = IRanges(na.lab.info[,3],na.lab.info[,4]), score = lab.na.df$formatScores)

      na.bins <- na.bin.list$naBins
      na.bin.ranges <- GRanges(seqnames = rep(paste('NA',lab,sep = '_'),nrow(na.bins)), ranges = IRanges(na.bins[,1], na.bins[,2]))
      na.bin.ranges$names <- paste(rep(paste('NA',lab,sep = '_'),nrow(na.bins)), na.bins[,1], na.bins[,2], sep = '_')
      all.chrom.ranges <- suppressWarnings(c(all.chrom.ranges, na.lab.range)) # formerly used append, caused warnings
      all.bin.ranges <- suppressWarnings(c(all.bin.ranges, na.bin.ranges))
      all.chrom.labels <- c(all.chrom.labels,rep(lab,nrow(na.bins)))
      temp.na.info <- cbind.data.frame(chrom = rep('NA', nrow(na.bins)), label = rep(lab,nrow(na.bins)),
        start = rep(NA, nrow(na.bins)), end = rep(NA, nrow(na.bins)), stringsAsFactors = F)
      all.bin.info <- rbind(all.bin.info,temp.na.info)
    }
  }

  bin.overlaps <- as.data.frame(findOverlaps(all.bin.ranges, all.chrom.ranges, type = 'any'))
  overlapping.bin.info <- all.bin.info[unique(bin.overlaps$queryHits),]

  #what is the highest/worst possible rank given threshold
  score.threshold.cutoff <- -log10(input.specs$scoreThresh)
  all.chrom.ranges.scores <- all.chrom.ranges$score
  all.chrom.ranges.norm.ranks <- rank(-all.chrom.ranges.scores)/length(all.chrom.ranges.scores)
  rank.thresh <- max(all.chrom.ranges.norm.ranks[which(all.chrom.ranges.scores > score.threshold.cutoff)])

  # create a list where each element contains the indeces for the tiles which overlap the bin
  bin.tile.list <- split(bin.overlaps, bin.overlaps$queryHits)

  # obtain the rho scores according to the cutoff
  bin.rho.scores <- unlist(lapply(bin.tile.list, function(x){
    temp.ranks <- all.chrom.ranges.norm.ranks[x$subjectHits]
    signif.ranks <- temp.ranks[which(temp.ranks <= rank.thresh)]
    all.len <- length(temp.ranks)
    signif.len <- length(signif.ranks)
    if(signif.len > 0){
      return(min(pbeta(sort(signif.ranks), seq(1, signif.len, 1), seq(all.len, all.len - signif.len + 1, -1))))
    } else {
      return(1)
    }
    }))

  print('Running RRA permutations....')
  permut.df <- data.frame(matrix(NA, nrow = length(bin.tile.list), ncol = input.specs$rraPermutNr))
  ptm <- proc.time()
  for(i in seq_len(input.specs$rraPermutNr)){
    rand.ranks <- sample(all.chrom.ranges.norm.ranks, size = length(all.chrom.ranges.norm.ranks), replace = FALSE)
    permut.df[[i]] <- unlist(lapply(bin.tile.list, function(x){
      temp.ranks <- rand.ranks[x$subjectHits]
      signif.ranks <- temp.ranks[which(temp.ranks <= rank.thresh)]
      all.len <- length(temp.ranks)
      signif.len <- length(signif.ranks)
      if(signif.len > 0){
        return(min(pbeta(sort(signif.ranks), seq(1, signif.len, 1), seq(all.len, all.len - signif.len + 1, -1))))
      } else{
        return(1)
      }
      }))
  }
  permut.time <- proc.time() - ptm
  print('Finished RRA permutations!')

  null.pvals <- unlist(permut.df)
  sorted.null.pvals <- sort(null.pvals)

  permut.pvals <- (findInterval(bin.rho.scores, sorted.null.pvals) + 1) / length(sorted.null.pvals)

  adj.pvals <- p.adjust(permut.pvals, method = 'BH')

  rra.out.df <- cbind.data.frame(rawScores = adj.pvals, formatScores = -log10(permut.pvals), log2_rate_ratio = rep(1, length(adj.pvals)),
    chrom = overlapping.bin.info$chrom, label = overlapping.bin.info$label, start = overlapping.bin.info$start,
    end = overlapping.bin.info$end, stringsAsFactors = F)

  return(rra.out.df)
}

#apply RRA to all score datasets
#input:
#   input.list: all method-specifi results to be analyzed using a specific RRA
#   input.specs: specification list
#output: list
#   rawScores, formatScores, log2_rate_ratio, chromosome, label, sgRNA targets (can be 1 or 2 columns)
RRA_wrapper <- function(input.list, input.specs, input.rra.path, input.name){

  out.list <- list()
  list.names <- names(input.list)

  for(i in 1:length(input.list)){
    if(! list.names[i] %in% c('viterbi', 'frwdBkwd')){
      temp.list.name <- paste(list.names[i], input.name, sep = '_')
      print(paste("Running MAGeCK's RRA for:", list.names[i], sep = ' '))
      out.list[[temp.list.name]] <- RRA(input.list[[i]], input.specs, list.names[i], input.rra.path, input.name)
    } else {
      temp.list.name <- list.names[i]
      out.list[[temp.list.name]] <- input.list[[i]]
    }
  }

  return(out.list)
}

# apply specified RRA to data
# input:
#   input.df: contains resutls on which rra should be perfomred: rawScores, formatScores, log2_rate_ratio, chromosome, labels, sgRNA targets (can be 1 or 2 columns)
#   input.specs: specifcation list
#   input.df.name: name of method on which RRA is performed (to save RRA files as method specific)
#   input.rra.path: RRA specific path
# output: data frame; rawScores, formatScores, log2_rate_ratio, chromosome, label, sgRNA targets (can be 1 or 2 columns)
# re-label bins based on either pre-exisiting input labels and their position, or from input label files
RRA <- function(input.df, input.specs, input.df.name, input.rra.path, input.RRA.name){

  #separate out chromsome from non-chromosome guides
  na.df <- c()  # data frame containing only guides not targeting chromosome
  all.chrom.df <- c()  # data frame containing only valid chromosomes
  if('NA' %in% input.df$chrom){
    all.chrom.df <- filter(input.df, grepl('chr',chrom))
    na.df <- filter(input.df, grepl('NA',chrom))
  } else {
    all.chrom.df <- input.df
  }

  # adjust info, if single guide, add additional column with identical values
  all.chrom.info <- adjustInfo(all.chrom.df)

  final.bin.labels <- c()  # labels associated with each bin
  final.bin.names <- c()  # name of each bin, used for matching with labels
  final.chrom.bin.means <- c()
  all.chrom.bins.overl <- c()  #overlap of all chromosome bins with guides

  # obtain binning info from chromsome targeting guides
  # set up RRA analysis file
  mean.bin.overlaps <- c()  # average number of overlaps per bin
  subjHitTracker <- 0  # ordering of guides across chromosomes seems to matter. This makes sure that ordering keeps chromosomes separate

  all.chroms <- unique(all.chrom.df$chrom)
  for(chrom in all.chroms){
    chrom.df <- all.chrom.df[all.chrom.df$chrom == chrom,]
    chrom.info <- all.chrom.info[all.chrom.info$chrom == chrom,]
    range.list <- obtain_data_GRanges(chrom.info, input.specs, chrom.df$formatScores)
    chrom.ranges <- range.list$gRanges
    bin.ranges <- obtain_bin_GRanges(range.list$rangeStartEnd, chrom, input.specs)
    chrom.labels <- obtain_bin_labels(bin.ranges, chrom.info, input.specs, chrom)
    chrom.overlaps <- as.data.frame(findOverlaps(bin.ranges,chrom.ranges,type = 'any'))
    # extract labels from chromsome bins which overlap, and the label names used for RRA
    final.bin.labels <- c(final.bin.labels, chrom.labels[chrom.overlaps$queryHits])
    final.bin.names <- c(final.bin.names, bin.ranges$names[chrom.overlaps$queryHits])
    #obtain average number of overlaps per bin (used to obtain average per na-bin)
    chom.bin.mean <- mean(table(chrom.overlaps$queryHits))
    final.chrom.bin.means <- c(final.chrom.bin.means, chom.bin.mean)

    chrom.bins.overl <- data.frame(bin = bin.ranges$names[chrom.overlaps$queryHits],
      score = chrom.ranges$score[chrom.overlaps$subjectHits],
      subjectHits = (chrom.overlaps$subjectHits + subjHitTracker))
    subjHitTracker <- subjHitTracker + max(chrom.overlaps$subjectHits)
    all.chrom.bins.overl <- rbind(all.chrom.bins.overl, chrom.bins.overl)
    # chrom.rra.format.df <- RRAformat_overlaps(bin.ranges, chrom.ranges, chrom.overlaps)
    # final.rra.format.df <- rbind(final.rra.format.df, chrom.rra.format.df)
  }

  #add NA part back into this
  if(! is.null(na.df)){
    na.chrom.bins.overl.list <- na_df_binning(final.chrom.bin.means, na.df, input.specs, subjHitTracker)
    all.chrom.bins.overl <- rbind(all.chrom.bins.overl, na.chrom.bins.overl.list$naBinOverl)
    final.bin.labels <- c(final.bin.labels, na.chrom.bins.overl.list$naLabels)
    final.bin.names <- c(final.bin.names, na.chrom.bins.overl.list$naBinNames)
  }

  final.rra.format.df <- RRAformat_overlaps(all.chrom.bins.overl)

  rra.format.filename <- paste(input.specs$dataName,input.df.name,input.RRA.name,'rankedBins.bed',sep = '_')
  write.table(final.rra.format.df, file = rra.format.filename, append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

  rra.max.percentile <- input.specs$RRAmaxPercent
  rra.command <- paste(input.rra.path,'/RRA -i ', rra.format.filename, ' -o ', input.specs$dataName,
    input.df.name, input.RRA.name, '_binPvalue.bed -p ', rra.max.percentile, sep = '')
  system(rra.command)

  rra.out.name <- paste(input.specs$dataName,input.df.name, input.RRA.name,'_binPvalue.bed',sep = '')
  #example output file format:
  #   group_id items_in_group   lo_value          p     FDR goodsgrna
  # 1 chr6_31139084_31139134             17 1.3031e-13 1.2338e-07 6.7e-05        14
  # 2 chr6_31138934_31138984             18 5.4450e-13 1.2338e-07 6.7e-05        14
  rra.out.file <- read.table(rra.out.name, head=TRUE, stringsAsFactors = FALSE)

  #use final.rra.format.df$bin to access names which are to be matched.
  #re-format rra output, add correct labels
  rra.out.df <- format_rra_output(rra.out.file,final.bin.labels,final.bin.names)

  return(rra.out.df)
}

# binning of data frame containing non-targeting guides
# non-targeting guides can be either negative controls, or guides which could not be
#   assigned to a specific location due to ambigous hits
# input:
#   input.mean.bin.guides: vector containing the mean number of guides a bin contains across all chromosomes
#   input.na.df: data frame contianint the non-targeting guide information (rawScores, formatScores, log2_rate_ratio, chromosome, label, sgRNA targets (can be 1 or 2 columns))
#   input.specs: specification list
#   input.tracker: tracker which ensures that each non-targeting label category stays together for RRA analysis
# output: list
#   naBinOverl: data frame containing the columns necessary for creating the RRA input file (bin, score, subjectHits)
#   naLabels: labels associated ith each non-targeting category
#   naBinNames: name of each bin, used for mathcing against RRA names and re-assigning the labels
na_df_binning <- function(input.mean.bin.guides, input.na.df, input.specs, input.tracker){

  output.na.bins.overl <- c()
  output.na.bin.names <- c()
  output.na.bin.labels <- c()

  bin.mean <- round(mean(input.mean.bin.guides))
  na.tracker <- input.tracker
  #treat each label within na separately
  #account for that when reading back into in 'format_rra_output'
  all.na.lables <- unique(input.na.df[,5])
  for(lab in all.na.lables){
    lab.na.df <- input.na.df[input.na.df$label == lab,]
    na.bin.list <- rra_na_binGeneration(bin.mean, nrow(lab.na.df), input.specs)

    na.lab.info <- cbind.data.frame(chrom = lab.na.df$chrom, label = lab.na.df$label,
      start = na.bin.list$naRanges, end = na.bin.list$naRanges + input.specs$CRISPReffectRange)

    na.lab.range <- GRanges(seqnames = paste(na.lab.info[,1],lab,sep = '_'),
        ranges = IRanges(na.lab.info[,3],na.lab.info[,4]), score = lab.na.df$formatScores)

    na.bins <- na.bin.list$naBins
    na.bin.ranges <- GRanges(seqnames = rep(paste('NA',lab,sep = '_'),nrow(na.bins)), ranges = IRanges(na.bins[,1], na.bins[,2]))
    na.bin.ranges$names <- paste(rep(paste('NA',lab,sep = '_'),nrow(na.bins)), na.bins[,1], na.bins[,2], sep = '_')

    na.lab.overlaps <- as.data.frame(findOverlaps(na.bin.ranges,na.lab.range,type = 'any'))

    output.na.bin.labels <- c(output.na.bin.labels, rep(lab, nrow(na.lab.overlaps)))
    output.na.bin.names <- c(output.na.bin.names, na.bin.ranges$names[na.lab.overlaps$queryHits])

    #format: ?
    na.lab.bins.overl <- data.frame(bin = na.bin.ranges$names[na.lab.overlaps$queryHits],
      score = na.lab.range$score[na.lab.overlaps$subjectHits],
      subjectHits = (na.lab.overlaps$subjectHits + na.tracker))
    na.tracker <- na.tracker + max(na.lab.overlaps$subjectHits)  # need to add tracker to the start
    output.na.bins.overl <- rbind(output.na.bins.overl, na.lab.bins.overl)

  }
  return(list(naBinOverl = output.na.bins.overl, naLabels =  output.na.bin.labels, naBinNames = output.na.bin.names))
}

#given, mean number of guides per bin, the length of na controls and specs
# implementation idea from: https://stackoverflow.com/questions/31120552/r-generate-a-sequence-of-numbers
# output: list:
  # naRanges: dataframe for ranges for each control
  # naBins: dataframe containing all the bins
  #   $naRanges
  #  [1]  1  1  1  1  1  1  1  1  1  1 22 22 22 22 22 22 22 22 22 22 43 43 43 43 43 43 43 43 43 43 64 64 64 64 64
  #
  # $naBins
  #   bin.range.setup[, 1] bin.range.setup[, 1] + input.specs$CRISPReffectRange
  # 1                    1                                                   21
  # 2                   22                                                   42
  # 3                   43                                                   63
  # 4                   64                                                   84
rra_na_binGeneration <- function(input.mean, input.na.length, input.specs){
  nr.na.bins <- ceiling(input.na.length/input.mean)  #since taking the ceiling, at most this many
  bin.starts <- rep(1,nr.na.bins) * c(0:(nr.na.bins - 1)) * (input.specs$CRISPReffectRange+1) + 1
  bin.range.setup <- cbind.data.frame(bin.starts, rep(input.mean,length(bin.starts))) #expand each column one item column 2 times
  full.na.bin.ranges <- rep(bin.range.setup[,1], bin.range.setup[,2])
  out.na.ranges <- full.na.bin.ranges[1:input.na.length]
  out.na.bins <- cbind.data.frame(bin.range.setup[,1], bin.range.setup[,1] + input.specs$CRISPReffectRange)
  out.list <- list()
  out.list$naRanges <- out.na.ranges
  out.list$naBins <- out.na.bins
  return(out.list)
}

# given the output from rra, adjust the outputformat and add the correct labels
# input:
#   input.rra.df: output file from running rra (group_id, items_in_group, lo_value, p, FDR, goodsgrna)
#   input.labels: labels per bin
#   input.label.bins: bins used in rra, row-positions maching up to the labels ('chr6_31138934_31138984')
# output: data frame
#   columns: rawScores, formatScores (-log10 p), log2_rate_ratio, chromosome, label, sgRNA targets (can be 1 or 2 columns)
format_rra_output <- function(input.rra.df, input.labels, input.label.bins){

  chrom.df <- filter(input.rra.df, grepl('chr',group_id))
  chrom.rra.matches <- match(chrom.df[,1], input.label.bins)
  chrom.labels <- input.labels[chrom.rra.matches]
  chrom.coords <- do.call(rbind, strsplit(chrom.df[,1],'_'))

  # rra outputs are given a ratio of 1, since the ratio was already used for the ranking
  out.df <- cbind.data.frame(rawScores = chrom.df$FDR, formatScores = -log10(chrom.df$p),
    log2_rate_ratio = rep(1,nrow(chrom.df)), chrom = chrom.coords[,1],
    label = chrom.labels, start = as.numeric(chrom.coords[,2]), end = as.numeric(chrom.coords[,3]), stringsAsFactors = FALSE)

  if(nrow(input.rra.df) > nrow(chrom.df)){
    na.df <- filter(input.rra.df, grepl('NA',group_id))
    na.rra.matches <- match(na.df[,1], input.label.bins)
    na.labels <- input.labels[na.rra.matches]
    na.coords <- as.data.frame(do.call(rbind, strsplit(na.df[,1],'_')), stringsAsFactors = FALSE)
    label.order <- order(na.labels)
    ordered.na.df <- na.df[label.order,]
    ordered.na.coords <- na.coords[label.order,]

    # rra outputs are given a ratio of 1, since the ratio was already used for the ranking
    na.df <- cbind.data.frame(rawScores = ordered.na.df$p, formatScores = -log10(ordered.na.df$p),
      log2_rate_ratio = rep(1,nrow(na.df)), chrom = ordered.na.coords[,1],
      label = ordered.na.coords[,2], start = rep(NA,nrow(na.df)),
      end = rep(NA,nrow(na.df)), stringsAsFactors = FALSE)

    out.df <- rbind(out.df, na.df)
  }
  return(out.df)
}

#input: data frame of bins overlapped by guides across all chromsome
#   format: bin, score, subjectHits
#output: RRA formatted data frame:
  # sgRNA: not sure if necessary, might be just for labelling
  # symbol: what bin is overlapped
  # pool: 'list' is only thing specified here, not sure how used
  # p.high: value equal to the scroes from the input.sgRanges
  # prob: not sure why all is set to 1
  # chose: seems to delineat the first time a given sgRNA shows up
RRAformat_overlaps <- function(input.bins.overl){

  #order according to the position of the deletion start
  temp.bins.overl.ordered <- input.bins.overl[order(input.bins.overl$subjectHits),]
  temp.bins.overl.ordered$rank <- c(1:nrow(temp.bins.overl.ordered))
  #not sure what below does, seems like first time a particular guide appears
  index <- do.call(rbind, lapply(split(temp.bins.overl.ordered[,c(3,4)],
    temp.bins.overl.ordered$subjectHits), function(x) x[1,]))
  temp.bins.overl.ordered$chosen = 0
  temp.bins.overl.ordered$chosen[index$rank] = 1

  temp.out.df <- data.frame(sgrna=1:nrow(temp.bins.overl.ordered),
    symbol=temp.bins.overl.ordered$bin,
    pool="list", p.high=temp.bins.overl.ordered$score, prob=1,
    chosen=temp.bins.overl.ordered$chosen)

  return(temp.out.df)
}

# given a set of bins, re-assign original labels to bin-specific labels
# input:
#   input.bins:  bins of the chromosome of interest. Are in GRanges format
#   input.info: chromsome specific info
#   input.specs: specification list
#   input.chrom: chromsome currently labelled
# output: vectro containing the labels, vector format (c('pos','neg',...))
obtain_bin_labels <- function(input.bins, input.info, input.specs, input.chrom){

  out.labels <- c()
  #if a label file is provided, use it to re-assign bin labels, else use info labels
  if('labelFile' %in% names(input.specs)){
    label.file <- input.specs$labelFile
    check_label_match(label.file,input.info, input.specs$labelHierarchy[1])  # need to check!
    label.df <- obtain_label_df(label.file, input.specs)  # need to check!
    out.labels <- assignLabels(input.bins, label.df, input.specs, input.chrom)  # need to implement
  } else {
    label.df <- info_to_label(input.info)
    out.labels <- assignLabels(input.bins, label.df, input.specs, input.chrom)
  }
  return(out.labels)
}

# given bins for a specific chromosome and labels, assign each bin a label
# input:
#   input.bins: bins of the chromosome of interest. Are in GRanges format
#   input.labels: data frame containing information about the labels
#     colnames: label, chrom, start, end
#   input.specs: specification file, used for obtianing label hierarchy
#   input.chrom: chromosome currently being re-labelled
# output: vector continaing label names (c('pos','pos','neg'...))
assignLabels <- function(input.bins, input.labels, input.specs, input.chrom){

  bins.df <- as.data.frame(input.bins)
  bin.labels <- rep(input.specs$labelHierarchy[1], nrow(bins.df))  # initiate as default bin label

  #if there are labels for the given chromosome, assing labels per bin, else return the default label
  if(input.chrom %in% input.labels[,2]){
    chrom.labels <- input.labels[input.labels[,2] == input.chrom,]
    granges.label <- GRanges(seqnames = chrom.labels[,2], ranges = IRanges(chrom.labels[,3],chrom.labels[,4]), label = chrom.labels[,1])
    bin.labels <- per_label_overlap(input.bins, granges.label, input.specs, bin.labels)
  } else {
    return(bin.labels)
  }
}

# for a given chromosome, find the chromosome specific labels and overlap them with the chromomsome bins
# input:
#   input.bins: bins of the chromosome of interest. Are in GRanges format
#   input.label.ranges:chromosome specific labels, GRanges format
#   input.specs: specification file, used for obtianing label hierarchy
#   input.bin.labels: vector of labels, is continously updataed for each label processed
# output: vector continaing label names (c('pos','pos','neg'...))
per_label_overlap <- function(input.bins, input.label.ranges, input.specs, input.bin.labels){

  bin.labels <- input.bin.labels

  for(labs in input.specs$labelHierarchy){
    if(labs %in% input.label.ranges$label){
      lab.input.label.ranges <- input.label.ranges[input.label.ranges$label == labs,]
      lab.overlaps <- as.data.frame(findOverlaps(input.bins, lab.input.label.ranges,type = 'any'))
      bin.labels[lab.overlaps$queryHits] <- labs
    }
  }
  return(bin.labels)
}

# converts regular input file to label.df format.
# input.labels: data frame containing information about the labels (info was already adjusted)
#   columns: chrom, label, start, end
# output: data frame with following columns: label, chrom, start, end
info_to_label <- function(input.info){
  out.df <- cbind.data.frame(label = input.info[,2], chrom = input.info[,1], start = input.info[,3], end = input.info[,4], stringsAsFactors = FALSE)
  return(out.df)
}

# given an input label file, extract the label coordinates
# input:
#   input.label.file: file containing the location of the label files (pos:../EnhancerPos.csv)
#   input.specs: specification list
# output: label data frame,
#   columns: label, chrom, start, end
obtain_label_df <- function(input.label.file, input.specs){
  read.labels <- scan(input.label.file,what='character')

  label.coords <- c()
  for(lab in read.labels){
    temp.lab <- strsplit(lab, ':')[[1]][1]
    temp.file <- strsplit(lab, ':')[[1]][2]
    if(temp.lab != input.specs$labelHierarchy[1]){
      temp.lab.coord <- read.csv(temp.file, as.is = TRUE)
      temp.out.lab.coord <- cbind.data.frame(label = rep(temp.lab,nrow(temp.lab.coord)), chrom = temp.lab.coord[,1],
        start = temp.lab.coord[,2], end = temp.lab.coord[,3])
      label.coords <- rbind(label.coords, temp.out.lab.coord)
    }
  }
  return(label.coords)
}

# check that the number of labels in the label file matches up with the ones from the data
# input:
#   input.label.file: location and filename containing label-specific paths to csv files contianing label location
#   input.info: info file (cols: chrom, label, start, end)
#   default.label: first label of the $labelHierarchy, the default label
# output: none if tests are passed
check_label_match <- function(input.label.file, input.info, default.label){
  read.labels <- scan(input.label.file,what='character')
  info.labels <- unique(input.info[,2])

  labels <- c()
  label.coords <- c()
  for(lab in read.labels){
    temp.lab.name <- strsplit(lab, ':')[[1]][1]
    labels <- c(labels, temp.lab.name)
    if(lab == default.label){
      stopifnot(strsplit(lab, ':')[[1]][2] == '-')  # if the lowest hierarchy label is not the default label
    }
    # temp.lab.coord <- read.csv(strsplit(lab, ':')[[1]][2], as.is = TRUE)
  }
  stopifnot(length(labels) >= length(info.labels))  # 'Error in "check_label_match": provided label file specifiers must be same as in data')
  stopifnot(all((info.labels %in% labels) == TRUE))  # 'Error in "check_label_match": provided labels must match (upper and lower case is active) labels of data')

}

# obtain the bins for the genomic region of interest
# input:
#   input.range: range for the bins (vector of format: c(start, end))
#   input.chrom: chromosome for which the bins are created
#   input.specs: specifications file for bin size
# output: GRanges object:
#   example: 'temp.bin.ranges':
#   GRanges object with 5 ranges and 1 metadata column:
#     seqnames               ranges strand |                  names
#        <Rle>            <IRanges>  <Rle> |            <character>
#   [1]     chr6 [30132134, 30132184]      * | chr6_30132134_30132184
#   [2]     chr6 [30132184, 30132234]      * | chr6_30132184_30132234
obtain_bin_GRanges <- function(input.range, input.chrom, input.specs){
  temp.range.bin.start <- seq(input.range[1], input.range[2], by = input.specs$binSize)
  temp.range.bin.end <- temp.range.bin.start + input.specs$binSize

  temp.bin.ranges <- GRanges(seqnames = rep(input.chrom,length(temp.range.bin.start)),
    ranges = IRanges(temp.range.bin.start, temp.range.bin.end))
  temp.bin.ranges$names <- paste(input.chrom, paste(temp.range.bin.start, temp.range.bin.end,
    sep = '_'), sep = '_')
  return(temp.bin.ranges)
}

# given information of a chromosome, and position-specific scores, convert to GRanges
#   and add scores as additional information
#input:
#   input.chrom.info: chromosome specific info; 'chrom', 'label', and either one or two columns
#   input.specs: data specifications,
#   input.scores: chromosome specific scores, note that RRA goes after negatively enriched scores
#output: list:
  # gRanges: contains the GRanges object of the guide ranges
  # rangeStartEnd: contains max and min position of targeting sgRNAs
obtain_data_GRanges <- function(input.chrom.info, input.specs, input.scores){
  out.sgRange <-c()
  out.sgRange <- GRanges(seqnames = input.chrom.info$chrom,
    ranges = IRanges(input.chrom.info$start, input.chrom.info$end),
    score = -input.scores)
  temp.range.start <- min(input.chrom.info$start)
  temp.range.end <- max(input.chrom.info$end)

  out.list <- list()
  out.list$gRanges <- out.sgRange
  out.list$rangeStartEnd <- c(temp.range.start, temp.range.end)
  return(out.list)
}

# extract the chromsome information from the score output format
# if dealing with single guide data, add a column with identical target site to simulate dual guides
# for downstream analysis consistency
# input: data frame; rawScores, formatScores, log2_rate_ratio, chromosome, labels, sgRNA targets (can be 1 or 2 columns)
# output: data frame; chromosome, labels, sgRNA target1 ($start), sgRNA target2 ($end)
adjustInfo <- function(input.df){
  output.df <- c()
  #if data is from dual guides (7 columns, extract chromosome information), else add column
  if('end' %in% colnames(input.df)){
    output.df <- input.df[,c(4:ncol(input.df))]
  } else {
    output.df <- input.df[,c(4:ncol(input.df))]
    output.df$end <- output.df$start + 1
  }
  return(output.df)
}

# obtain the scores for all specified analysis methods
# input: count df, info df, specification list
# output: list, each method adds the following to the
#   $method, where each $mthod contains the following rows:
#     rawScores: calculated probabilities of p-values
#     formatScores: adjusted scores such that they are signed by ratio and can be ranked from most to least significant
#     log2_rate_ratio: log2 ratio used for sign-adjusting scores, either based on count ratio or calculated rates
#     chromosome, labels, sgRNA targets (can be 1 or 2 columns) (all values from the info file, can be altered in methods like CREST analysis)
score_calculation <- function(input.counts, input.info, input.specs){
  out.list <- list()
  if('NB' %in% input.specs$Method){
    print('#####################################')
    print('Running NB analysis')
    out.list$NB <- NB_analysis(input.counts, input.info, input.specs)
    print('Finished NB analysis')
  }
  if('ZINB' %in% input.specs$Method){
    print('#####################################')
    print('Running ZINB analysis')
    out.list$ZINB <- ZINB_analysis(input.counts, input.info, input.specs)
    print('Finished ZINB analysis')
  }
  if('FoldChange' %in% input.specs$Method){
    print('#####################################')
    print('Running FoldChange analysis')
    out.list$FoldChange <- fold_change_analysis(input.counts, input.info, input.specs)
    print('Finished FoldChange analysis')
  }
  if('edgeR' %in% input.specs$Method){
    print('#####################################')
    print('Running edgeR analysis')
    out.list$edgeR <- edgeR_analysis(input.counts, input.info, input.specs)
    print('Finished edgeR analysis')
  }
  if('DESeq2' %in% input.specs$Method){
    print('#####################################')
    print('Running DESeq2 analysis')
    out.list$DESeq2 <- DESeq2_analysis(input.counts, input.info, input.specs)
    print('Finished DESeq2 analysis')
  }
  if('NB-GLMM' %in% input.specs$Method){
    print('#####################################')
    print('Running NB-GLMM analysis')
    out.list$nbGlmm <- nb_tmb_glmm_analysis(input.counts, input.info, input.specs)
    print('Finished NB-GLMM analysis')
  }
  if('RELICS-search' %in% input.specs$Method){
    print('#####################################')
    print('Running RELICS-search ')
    out.list$RELICS <- RELICS_search(input.counts, input.info, input.specs)
    print('Finished RELICS-search')
  }
  if('NB-GLMM-Untrained' %in% input.specs$Method){
    print('#####################################')
    print('Running NB-GLMM-Untrained')
    out.list$nbGlmmUntr <- nbGlmm_untrained(input.counts, input.info, input.specs)
    print('NB-GLMM-Untrained')
  }
  dirichlet_multinomial_analysis
  if('Dirichlet-Multinom' %in% input.specs$Method){
    print('#####################################')
    print('Running Dirichlet Multinomial Analysis')
    out.list$dirichletMultinom <- dirichlet_multinomial_analysis(input.counts, input.info, input.specs)
    print('Finished Dirichlet Multinomial Analysis')
  }
  if('NB-GLMM-Dev' %in% input.specs$Method){
    print('#####################################')
    print('Running NB-GLMM-Dev analysis')
    out.list$nbGlmmDev <- nb_deviance_glmm_analysis(input.counts, input.info, input.specs)
    print('Finished NB-GLMM-Dev analysis')
  }
  if('allData-NB-GLMM-category-parallel' %in% input.specs$Method){
    print('#####################################')
    print('Running allData-NB-GLMM-category-parallel analysis')
    out.list$allDataNbGlmmCategoryParallel <- allData_nb_tmb_glmm_categoryAnalysis_parallel(input.counts, input.info, input.specs)
    print('Finished allData-NB-GLMM-category-parallel analysis')
  }
  if('allData-NB-GLMM' %in% input.specs$Method){
    print('#####################################')
    print('Running allData-NB-GLMM analysis')
    out.list$allDataNbGlmm <- allData_nb_tmb_glmm_analysis(input.counts, input.info, input.specs)
    print('Finished allData-NB-GLMM analysis')
  }
  if('NumInt-Poisson-GLMM' %in% input.specs$Method){
    print('#####################################')
    print('Running NumInt-Poisson-GLMM analysis')
    out.list$numIntPoisGlmm <- poisson_univariate_integration_tmb_glmm_analysis(input.counts, input.info, input.specs)
    print('Finished NumInt-Poisson-GLMM analysis')
  }
  if('Uni-Poisson-GLMM' %in% input.specs$Method){
    print('#####################################')
    print('Running Uni-Poisson-GLMM analysis')
    out.list$uniPoisGlmm <- poisson_1D_laplace_tmb_glmm_analysis(input.counts, input.info, input.specs)
    print('Finished Uni-Poisson-GLMM analysis')
  }
  if('Multi-Poisson-GLMM' %in% input.specs$Method){
    print('#####################################')
    print('Running Multi-Poisson-GLMM analysis')
    out.list$multiPoisGlmm <- poisson_multiDim_laplace_tmb_glmm_analysis(input.counts, input.info, input.specs)
    print('Finished Multi-Poisson-GLMM analysis')
  }
  if('allData-Multi-Poisson-GLMM-perRepl' %in% input.specs$Method){
    print('#####################################')
    print('allData-Multi-Poisson-GLMM-perRepl')
    out.list$allDataMultiPoisGlmmPerRepl <- allData_poisson_multiDim_laplace_tmb_glmm_analysis_perRepl(input.counts, input.info, input.specs)
    print('Finished allData-Multi-Poisson-GLMM-perRepl')
  }
  if('allData-Multi-Poisson-GLMM-acrossRepl' %in% input.specs$Method){
    print('#####################################')
    print('allData-Multi-Poisson-GLMM-acrossRepl')
    out.list$allDataMultiPoisGlmmAcrossRepl <- allData_poisson_multiDim_laplace_tmb_glmm_analysis_acrossRepl(input.counts, input.info, input.specs)
    print('Finished allData-Multi-Poisson-GLMM-acrossRepl')
  }
  if('reducedData-Multi-Poisson-GLMM-perRepl' %in% input.specs$Method){
    print('#####################################')
    print('reducedData-Multi-Poisson-GLMM-perRepl')
    out.list$reducedDataMultiPoisGlmmPerRepl <- reducedData_poisson_multiDim_laplace_tmb_glmm_analysis_perRepl(input.counts, input.info, input.specs)
    print('Finished reducedData-Multi-Poisson-GLMM-perRepl')
  }
  if('reducedData-Multi-Poisson-GLMM-acrossRepl' %in% input.specs$Method){
    print('#####################################')
    print('reducedData-Multi-Poisson-GLMM-acrossRepl')
    out.list$reducedDataMultiPoisGlmmAcrossRepl <- reducedData_poisson_multiDim_laplace_tmb_glmm_analysis_acrossRepl(input.counts, input.info, input.specs)
    print('Finished reducedData-Multi-Poisson-GLMM-acrossRepl')
  }
  if('Viterbi' %in% input.specs$Method){
    print('#####################################')
    print('Running Viterbi analysis')
    out.list$viterbi <- viterbi_hmm_tmb_glmm_analysis(input.counts, input.info, input.specs)
    print('Finished Viterbi analysis')
  }
  if('Frwd-Bkwd' %in% input.specs$Method){
    print('#####################################')
    print('Running Frwd-Bkwd analysis')
    out.list$frwdBkwd <- frwdBkwd_hmm_tmb_glmm_analysis(input.counts, input.info, input.specs)
    print('Finished Frwd-Bkwd analysis')
  }
  # if('CrestAnalysis' %in% input.specs$Method){
  #   out.list$CREST <- crest_seq_analysis(input.counts, input.info, input.specs)
  # }
  return(out.list)
}

glmm_tmb_parEst <- function(input.counts, input.specs){

  pool.IDs <- colnames(input.counts)
  guide.nrs <- c(1:nrow(input.counts))
  guide.names <- paste0('guide', guide.nrs)
  pool.nrs <- c(1:ncol(input.counts))
  pool.names <- paste0('pool', pool.nrs)
  # positive.types <- c(1:(ncol(input.counts) - 1))
  #
  #change to : intercept, pos1, pos
  # pos.type.label <- c('interc', paste0('pos', positive.types))
  # nr.positives <- length(which(input.type == 'pos'))

  neg.type.label <- c('interc', paste0('neg', pool.nrs[-length(pool.nrs)]))
  # nr.negatives <- length(which(input.type == 'neg'))

  guide.labels <- rep(guide.names, each = length(pool.names))
  pool.labels <- rep(pool.names, length(guide.names))
  type.labels <- rep(neg.type.label, length(guide.nrs))
  # type.labels <- c(rep(pos.type.label, nr.positives), rep(neg.type.label, nr.negatives))

  long.df <- cbind.data.frame(observed = c(t(input.counts)), guide_labels = guide.labels,
    pool_labels = pool.labels, type_labels = type.labels, stringsAsFactors = T)

  tmb.nb2.model <- glmmTMB(observed ~ type_labels + (1|guide_labels),
    data = long.df, family=nbinom2)

  model.dispersion <- sigma(tmb.nb2.model) # residual variance, theta corresponding to the NB variance
  model.tau <- exp(getME(tmb.nb2.model, name = 'theta'))  # SD of the random effect
  #  no idea why this works, no idea how else to get random effect SD

  nr.beta.neg <- ncol(input.counts) - 1
  tmb.betas.orig <- tmb.nb2.model$fit$par[2:(nr.beta.neg + 1)]
  beta.neg <- tmb.betas.orig[1:nr.beta.neg]

  # nr.beta.neg <- ncol(input.counts) - 1
  # nr.beta.pos <- ncol(input.counts) - 1
  # nr.total.beta <- nr.beta.neg + nr.beta.pos
  # tmb.betas.orig <- tmb.nb2.model$fit$par[2:(nr.total.beta + 1)]
  # beta.neg <- tmb.betas.orig[1:nr.beta.neg] #+ c(0, rep(tmb.betas.orig[1], nr.beta.neg - 1))
  # beta.pos <- tmb.betas.orig[(nr.beta.neg + 1):nr.total.beta] #beta.neg + c(0, tmb.betas.orig[(nr.beta.neg + 1):nr.total.beta])

  if('searchSilencers' %in% names(input.specs)){
    beta.pos <- beta.neg - (beta.pos - beta.neg)  # inverse the effect of the positives
  }

  tmb.neg_logLik <- logLik(tmb.nb2.model)[1]

  out.parameters <- c(tmb.nb2.model$fit$par[1], beta.neg, model.tau, model.dispersion, tmb.neg_logLik)
  out.df <- t(as.data.frame(out.parameters))
  colnames(out.df) <- c('beta_intercept', paste('beta_neg_pred', c(1:length(beta.neg)), sep = '_'),
    'tau', 'disp', 'neg_ll')

  return(out.df)

  # out.list <- list(betas = tmb.betas.orig[2:(nr.beta.neg + 1)], beta_pos_pred = beta.pos,
  #   beta_neg_pred = beta.neg, model_tau = model.tau, model_dispersion = model.dispersion,
  #   neg_logLik = tmb.neg_logLik, pool_IDs = pool.IDs, beta_intercept = tmb.nb2.model$fit$par[1])
  #
  # final.list <- list(out_df = out.df, out_list = out.list)

}

glmm_null_pars <- function(input.counts, input.specs){

  filtered.counts <- c()
  if(length(which(rowSums(input.counts) == 0)) > 0){
    filtered.counts <- input.counts[which(rowSums(input.counts) > 0),]
  } else {
    filtered.counts <- input.counts
  }

  null.pars <- glmm_tmb_parEst(filtered.counts, input.specs)

  return(null.pars)

  # total.counts <- rowSums(filtered.counts)
  #
  # out.optim <- optim(input.pars, dirichlet_null_logLik, method = 'Nelder-Mead',
  #   input.counts = filtered.counts, input.size = total.counts)
  #
  # out.ll <- out.optim$value
  # out.disp <- out.optim$par[1]
  # out.prop <- c(out.optim$par[2:length(out.optim$par)], 1 - sum(out.optim$par[2:length(out.optim$par)]))
  #
  # out.vector <- c(out.ll, c(out.disp, out.prop))
  # return(out.vector)
}

glmm_likelihood_calculation <- function(input.counts, input.betas, input.intercept,
  input.tau, input.disp){

  guide.ll <- vector('numeric', length = nrow(input.counts))

  for(i in 1:nrow(input.counts)){
    temp.guide.counts <- unlist(input.counts[i,])
    # add the 0 because first pool is intercept
    temp_integration_function <- RELICS_function_generator(temp.guide.counts,
      c(0, input.betas), input.tau, input.disp, input.intercept)
    temp.integr <- log(integrate(temp_integration_function,
      lower = (-5*input.tau + input.intercept),
      upper = (5*input.tau + input.intercept) )$value[1])
    if(temp.integr == -Inf){
      guide.ll[i] <- -745
    } else {
      guide.ll[i] <- temp.integr
    }
  }

  return(-sum(guide.ll))
}

glmm_null_ll <- function(input.counts, input.pars){

  nr.pars <- length(input.pars)
  input.disp <- input.pars[nr.pars - 1]
  input.tau <- input.pars[nr.pars - 2]
  input.intercept <- input.pars[1]
  input.betas <- input.pars[-c(1, (nr.pars - 2), (nr.pars - 1), nr.pars)]

  filtered.counts <- c()
  if(length(which(rowSums(input.counts) == 0)) > 0){
    filtered.counts <- input.counts[which(rowSums(input.counts) > 0),]
    if(nrow(filtered.counts) == 0){
      return(0)
    }
  } else {
    filtered.counts <- input.counts
  }

  neg.ll <- glmm_likelihood_calculation(filtered.counts, input.betas, input.intercept,
    input.tau, input.disp)

  return(neg.ll)
}

optim_glmm_ll <- function(par, input.intercept, input.counts, input.tau, input.disp){

  betas <- par

  out.neg.ll <- glmm_likelihood_calculation(input.counts, betas, input.intercept,
    input.tau, input.disp)

  return(out.neg.ll)
}

glmm_alt_ll <- function(input.counts, input.pars){

  nr.pars <- length(input.pars)
  input.disp <- input.pars[nr.pars - 1]
  input.tau <- input.pars[nr.pars - 2]
  input.intercept <- input.pars[1]
  input.betas <- input.pars[-c(1, (nr.pars - 2), (nr.pars - 1), nr.pars)]

  filtered.counts <- c()
  if(length(which(rowSums(input.counts) == 0)) > 0){
    filtered.counts <- input.counts[which(rowSums(input.counts) > 0),]
    if(nrow(filtered.counts) == 0){
      return(c(0, input.betas))
    }
  } else {
    filtered.counts <- input.counts
  }

  optim.glmm.pars <- optim(input.betas, optim_glmm_ll, input.intercept = input.intercept,
    method = 'Nelder-Mead', input.counts = filtered.counts, input.tau = input.tau,
    input.disp = input.disp)

  out.neg.ll <- optim.glmm.pars$value
  out.pars <- optim.glmm.pars$par

  return(c(out.neg.ll, out.pars))
}

nbGlmm_untrained <- function(input.counts, input.info, input.specs){

  neg.counts <- input.counts[which(input.info$label %in% input.specs$glmm_untrainedNull), ]

  if(input.specs$crisprSystem %in% c('CRISPRcas9', 'CRISPRi', 'CRISPRa')){
    input.info$start <- input.info$start - round(input.specs$crisprEffectRange / 2)
    input.info$end <- input.info$end + round(input.specs$crisprEffectRange / 2)
  }

  genome.regions <- extract_all_genome_regions(input.counts, input.info, input.specs)

  final.null.ll <- c()
  final.alt.ll <- c()
  final.alt.pars <- list()
  final.llr <- c()
  final.pVal <- c()
  total.df <- 0  # total degrees of freedom

  final.null.par.list <- list()
  # for each replicate
  for(rep in 1:length(input.specs$repl_groups)){

    print(paste0('replicate: ', rep))
    ptm <- proc.time()
    temp.repl.cols <- input.specs$repl_groups[[rep]]
    temp.repl.neg <- neg.counts[, temp.repl.cols]

    # temp.null.init.alpha <- colSums(temp.repl.neg) / sum(temp.repl.neg)
    # temp.null.init.disp <- 1

    repl.glmm.null <- glmm_null_pars(temp.repl.neg, input.specs)
    final.null.par.list[[rep]] <- repl.glmm.null

    # repl.dirichlet.null.ll <- repl.dirichlet.null[1]
    # repl.dirichlet.null.disp <- repl.dirichlet.null[2]
    # repl.dirichlet.null.pars <- repl.dirichlet.null[3:length(repl.dirichlet.null)]

    repl.null.ll <- c()
    repl.alt.ll <- c()
    repl.alt.pars <- c()
    repl.info <- c()

    for(chrom in 1:length(genome.regions)){
      temp.genome.region <- genome.regions[[chrom]]

      print(paste0('Section: ', chrom))

      temp.glmm.null <- c()
      temp.glmm.alt <- c()

      if('nrFreeCores' %in% names(input.specs)){

        nr.cores <- max((detectCores() - input.specs$nrFreeCores), 1)
        cl <- c()
        if('parallelOutfile' %in% names(input.specs)){
          cl <- makeCluster(nr.cores, outfile = input.specs$parallelOutfile)
          setDefaultCluster(cl=cl)
        } else {
          cl <- makeCluster(nr.cores)
          setDefaultCluster(cl=cl)
        }
        # cl <- makeCluster(nr.cores, outfile = input.specs$parallelOutfile)
        # setDefaultCluster(cl=cl)

        clusterExport(cl, c('glmm_null_ll', 'glmm_alt_ll', 'glmm_likelihood_calculation',
          'RELICS_function_generator', 'optim_glmm_ll'), envir=environment())

        print('Computing Null')
        temp.glmm.null <- do.call(c, parLapply(cl, temp.genome.region$region_counts, function(x){
          repl.counts <- x[, temp.repl.cols]
          temp.chrom.glmm.out <- glmm_null_ll(repl.counts, repl.glmm.null)
          print(paste(temp.chrom.glmm.out, collapse = '_'))
          temp.chrom.glmm.out
        }))

        print('Computing alt.')
        temp.glmm.alt <- do.call(rbind, parLapply(cl, temp.genome.region$region_counts, function(x){
          repl.counts <- x[, temp.repl.cols]
          temp.chrom.glmm.out <- glmm_alt_ll(repl.counts, repl.glmm.null)
          print(paste(temp.chrom.glmm.out, collapse = '_'))
          temp.chrom.glmm.out
        }))

        stopCluster(cl)

      } else {
        print('Computing Null')
        temp.glmm.null <- do.call(c, lapply(temp.genome.region$region_counts, function(x){
          repl.counts <- x[, temp.repl.cols]
          temp.chrom.glmm.out <- glmm_null_ll(repl.counts, repl.glmm.null)
          temp.chrom.glmm.out
          }))

        print('Computing alt.')
        temp.glmm.alt <- do.call(rbind, lapply(temp.genome.region$region_counts, function(x){
          repl.counts <- x[, temp.repl.cols]
          temp.chrom.glmm.out <- glmm_alt_ll(repl.counts, repl.glmm.null)
          temp.chrom.glmm.out
          }))

      }

      repl.null.ll <- c(repl.null.ll, temp.glmm.null)
      repl.alt.ll <- c(repl.alt.ll, temp.glmm.alt[,1])
      repl.alt.pars <- rbind(repl.alt.pars, temp.glmm.alt[,c(2:ncol(temp.glmm.alt))])
      repl.info <- rbind(repl.info, temp.genome.region$region_info_df)
    }


    final.null.ll <- cbind(final.null.ll, repl.null.ll)
    final.alt.ll <- cbind(final.alt.ll, repl.alt.ll)

    repl.llr <- repl.null.ll - repl.alt.ll
    final.llr <- cbind(final.llr, repl.llr)

    total.df <- total.df + (length(temp.repl.cols) - 1)

    parallel.runtime <- proc.time() - ptm
    print(paste0('GLMM Untrained Replicate runtime: ', round(parallel.runtime[3], 3) ))

    if('save_glmmUntrained_info' %in% names(input.specs)){
      save_glmmUntrained_repl_alt_pars(repl.alt.pars, repl.alt.ll, rep, input.specs)
    }
  }

  final.pVal <- pchisq(2 * rowSums(final.llr), df = total.df, lower.tail=F)

  out.df <- data.frame(rawScores = rowSums(final.llr), formatScores = -log10(final.pVal),
  log2_rate_ratio = rep(1, length(final.pVal)), stringsAsFactors = F)

  out.df <- cbind(out.df, repl.info)

  if('save_glmmUntrained_info' %in% names(input.specs)){
    save_glmmUntrained_info(final.llr, final.null.par.list, final.null.ll, input.specs)
  }

  return(out.df)
}

save_glmmUntrained_repl_alt_pars <- function(input.repl.pars, input.repl.ll, repl.nr, input.specs){

  out.df <- cbind(input.repl.ll, input.repl.pars)
  colnames(out.df) <- c('neg_ll', paste0('beta_', c(1:ncol(input.repl.pars) ) ))
  write.csv(out.df, file = paste(input.specs$dataName, 'repl', repl.nr, 'alt_glmmUntrPars.csv', sep = '_'), row.names = FALSE)

}

save_glmmUntrained_info <- function(input.repl.llr, input.null.pars, input.null.ll, input.specs){

  colnames(input.repl.llr) <- paste0('llr_repl', c(1:ncol(input.repl.llr)))
  write.csv(input.repl.llr, file = paste(input.specs$dataName, 'repl_glmmUntr_llr.csv', sep = '_'), row.names = FALSE)

  colnames(input.null.ll) <- paste0('null_neg_ll_repl', c(1:ncol(input.null.ll)))
  write.csv(input.null.ll, file = paste(input.specs$dataName, 'repl_null_glmmUntr_neg_ll.csv', sep = '_'), row.names = FALSE)

  for(i in 1:length(input.null.pars)){

    temp.df <- as.data.frame(t(input.null.pars[[i]]), stringsAsFactors = F)
    write.csv(temp.df, file = paste(input.specs$dataName, 'repl', i, '_null_glmmUntr_neg_pars.csv', sep = '_'), row.names = FALSE)
  }

}

# incomplete
dirichlet_dispersion_behavior <- function(input.counts){

  ordered.counts <- input.counts[order(rowSums(input.counts)),]
}

# for the null distribution:
# calculate the dirichlet -log likelihood given a matrix of counts , a
# dispersion and an alpha for each column
dirichlet_null_logLik <- function(par, input.counts, input.size){
  disp <- par[1]
  proportions <- par[2:length(par)]

  total.prop <- abs(c(proportions, 1 - sum(proportions)))

  alpha <- total.prop * disp

  dir.log.lik <- ddirmnom(input.counts, size = input.size, alpha, log=T)

  if(length(which(dir.log.lik == -Inf)) > 0){
    dir.log.lik[which(dir.log.lik == -Inf)] <- -745
  }
  -sum(dir.log.lik)
}

optim_dirichlet_null_pars <- function(input.counts, input.pars){

  filtered.counts <- c()
  if(length(which(rowSums(input.counts) == 0)) > 0){
    filtered.counts <- input.counts[which(rowSums(input.counts) > 0),]
  } else {
    filtered.counts <- input.counts
  }

  total.counts <- rowSums(filtered.counts)

  out.optim <- optim(input.pars, dirichlet_null_logLik, method = 'Nelder-Mead',
    input.counts = filtered.counts, input.size = total.counts)

  out.ll <- out.optim$value
  out.disp <- out.optim$par[1]
  out.prop <- c(out.optim$par[2:length(out.optim$par)], 1 - sum(out.optim$par[2:length(out.optim$par)]))

  out.vector <- c(out.ll, c(out.disp, out.prop))
  return(out.vector)
}

dirichlet_null_ll <- function(input.counts, input.pars){

  filtered.counts <- c()
  if(length(which(rowSums(input.counts) == 0)) > 0){
    filtered.counts <- input.counts[which(rowSums(input.counts) > 0),]
    if(nrow(filtered.counts) == 0){
      return(0)
    }
  } else {
    filtered.counts <- input.counts
  }

  total.counts <- rowSums(filtered.counts)

  out.ll <- dirichlet_null_logLik(par = input.pars, input.counts = filtered.counts,
    input.size = total.counts)

  return(out.ll)
}

# for the null distribution:
# calculate the dirichlet -log likelihood given a matrix of counts , a
# dispersion and an alpha for each column
dirichlet_alt_logLik <- function(par, input.counts, input.size, input.disp){
  disp <- input.disp
  proportions <- par

  total.prop <- abs(c(proportions, 1 - sum(proportions)))

  alpha <- total.prop * disp

  dir.log.lik <- ddirmnom(input.counts, size = input.size, alpha, log=T)

  if(length(which(dir.log.lik == -Inf)) > 0){
    dir.log.lik[which(dir.log.lik == -Inf)] <- -745
  }

  -sum(dir.log.lik)
}

optim_dirichlet_alt_pars <- function(input.counts, input.pars, input.disp){

  filtered.counts <- c()
  if(length(which(rowSums(input.counts) == 0)) > 0){
    filtered.counts <- input.counts[which(rowSums(input.counts) > 0),]
    if(nrow(filtered.counts) == 0){
      return(c(0, input.pars, 1 - sum(input.pars)))
    }
  } else {
    filtered.counts <- input.counts
  }

  total.counts <- rowSums(filtered.counts)

  out.optim <- optim(input.pars, dirichlet_alt_logLik, method = 'Nelder-Mead',
    input.counts = filtered.counts, input.size = total.counts, input.disp = input.disp)

  out.ll <- out.optim$value
  out.disp <- out.optim$par[1]
  out.prop <- c(out.optim$par[2:length(out.optim$par)], 1 - sum(out.optim$par[2:length(out.optim$par)]))

  out.vector <- c(out.ll, c(out.disp, out.prop))
  return(out.vector)
}

# can only be computed for individual chromosomes
extract_genome_regions <- function(input.counts, input.info, input.specs){

  combined.df <- cbind(input.counts, input.info)

  all.breaks <- unique(c(combined.df$start, combined.df$end))
  sorted.breaks <- sort(all.breaks)
  all.region.breaks <- sorted.breaks
  start.regions <- all.region.breaks[c(1:(length(all.region.breaks) - 1))]
  end.regions <- all.region.breaks[c(2:length(all.region.breaks))] - 1

  region.ranges <- GRanges(seqnames = rep(combined.df$chrom[1], length(start.regions)),
    ranges = IRanges(start.regions, end.regions))

  # adjust data ranges such that they fall within, not onto, boarders of genomic ranges
  data.ranges <- GRanges(seqnames = combined.df$chrom,
    ranges = IRanges(combined.df$start+1, combined.df$end-1))

  data.overlaps <- as.data.frame(findOverlaps(region.ranges, data.ranges, type = 'any'))

  overlapping.start <- start.regions[unique(data.overlaps$queryHits)]
  overlapping.end <- end.regions[unique(data.overlaps$queryHits)]

  region.overlap.list <- split(data.overlaps, data.overlaps$queryHits)

  region.counts <- lapply(region.overlap.list, function(x){
    temp.label <- input.counts[x$subjectHits,]
    return(temp.label)
    })

  region.labels <- unlist(lapply(region.overlap.list, function(x){
    temp.label <- combined.df$label[x$subjectHits]
    temp.labels.present <- which(input.specs$labelHierarchy %in% temp.label)
    out.label <- input.specs$labelHierarchy[max(temp.labels.present)]
    return(out.label)
    }))

  region.overlaps <- unlist(lapply(region.overlap.list, function(x){
    temp.overlaps <- length(combined.df$label[x$subjectHits])
    return(temp.overlaps)
    }))

  region.info.df <- data.frame(chrom = rep(combined.df$chrom[1], length(overlapping.start)),
    label = region.labels, start = overlapping.start, end = overlapping.end,
    nrSupportGuides = region.overlaps, stringsAsFactors = F)

  out.list <- list(region_counts = region.counts, region_info_df = region.info.df)

  return(out.list)
}

extract_nonTargeting_genome_regions <- function(input.counts, avg.overlap,
  input.specs, input.label){

  row.seq <- seq(1, nrow(input.counts), by = avg.overlap)
  nr.nonTargeting <- length(row.seq)
  out.scores <- vector('numeric', length = length(row.seq))
  overlapping <- c(0:(avg.overlap - 1))

  count.list <- list()
  info.df <- data.frame(chrom = rep('NA', nr.nonTargeting),
    label = rep(input.label, nr.nonTargeting),
    start = rep(NA, nr.nonTargeting), end = rep(NA, nr.nonTargeting),
    nrSupportGuides = rep(avg.overlap, nr.nonTargeting), stringsAsFactors = F)

  for(i in 1:(length(row.seq) - 1)){
    temp.rows <- row.seq[i] + overlapping
    temp.counts <- input.counts[temp.rows, ]
    count.list[[i]] <- temp.counts
  }

  last.index <- length(row.seq)
  temp.rows <- c(row.seq[last.index]:nrow(info.df))
  temp.counts <- input.counts[temp.rows, ]
  count.list[[last.index]] <- temp.counts
  info.df$nrSupportGuides[last.index] <- length(temp.rows)

  out.list <- list(region_counts = count.list, region_info_df = info.df)

  return(out.list)

}

# wrapper function to extract per-chromosome genome regions and add in non-targeting negative controls
extract_all_genome_regions <- function(input.counts, input.info, input.specs){

  avg.overlaps <- c()
  targeting.info <- c()  # data frame containing only valid chromosomes
  non.targeting.info <- c()
  targeting.counts <- c()
  non.targeting.counts <- c()
  if('NA' %in% input.info$chrom){
    na.rows <- which(input.info$chrom == 'NA')
    targeting.info <- input.info[-na.rows,] #filter(input.df, grepl('chr',chrom))
    non.targeting.info <- input.info[na.rows,] #filter(input.df, grepl('NA',chrom))
    targeting.counts <- input.counts[-na.rows,]
    non.targeting.counts <- input.counts[na.rows,]

  } else {
    targeting.info <- input.info
    targeting.counts <- input.counts
  }

  chrom.list <- list()
  all.chroms <- unique(targeting.info$chrom)

  if(length(all.chroms) == 1){
    sorted.targeting.counts <- targeting.counts[order(targeting.info$start),]
    sorted.targeting.chrom <- targeting.info[order(targeting.info$start),]

    chrom.list[[all.chroms[1]]] <- extract_genome_regions(sorted.targeting.counts,
      sorted.targeting.chrom, input.specs)

    avg.overlaps <- c(avg.overlaps, chrom.list[[all.chroms[1]]]$region_info_df$nrSupportGuides)
  } else {
    for(i in 1:length(all.chroms)){
      temp.chrom.rows <- which(targeting.info$chrom == all.chroms[i])
      targeting.counts.chrom <- targeting.counts[temp.chrom.rows, ]
      targeting.info.chrom <- targeting.info[temp.chrom.rows, ]

      sorted.targeting.counts.chrom <- targeting.counts.chrom[order(targeting.info.chrom$start),]
      sorted.targeting.info.chrom <- targeting.info.chrom[order(targeting.info.chrom$start),]

      chrom.list[[all.chroms[i]]] <- extract_genome_regions(sorted.targeting.counts.chrom,
        sorted.targeting.info.chrom, input.specs)

      avg.overlaps <- c(avg.overlaps, chrom.list[[all.chroms[i]]]$region_info_df$nrSupportGuides)
    }
  }
  if(nrow(targeting.info) != nrow(input.info)){
    avg.overlap <- round(median(avg.overlaps))
    non.targeting.chrom.list <- extract_nonTargeting_genome_regions(non.targeting.counts,
      avg.overlap, input.specs, non.targeting.info$label[1])

    chrom.list[['na_chrom']] <- non.targeting.chrom.list
  }
  return(chrom.list)
}

dirichlet_multinomial_analysis <- function(input.counts, input.info, input.specs){

  neg.counts <- input.counts[which(input.info$label %in% input.specs$dirichlet_null), ]

  if(input.specs$crisprSystem %in% c('CRISPRcas9', 'CRISPRi', 'CRISPRa')){
    input.info$start <- input.info$start - round(input.specs$crisprEffectRange / 2)
    input.info$end <- input.info$end + round(input.specs$crisprEffectRange / 2)
  }

  genome.regions <- extract_all_genome_regions(input.counts, input.info, input.specs)

  final.null.ll <- c()
  final.alt.ll <- c()
  final.alt.pars <- list()
  final.llr <- c()
  final.pVal <- c()
  total.df <- 0  # total degrees of freedom

  final.null.par.list <- list()
  # for each replicate
  for(rep in 1:length(input.specs$repl_groups)){

    ptm <- proc.time()

    temp.repl.cols <- input.specs$repl_groups[[rep]]
    temp.repl.neg <- neg.counts[, temp.repl.cols]

    temp.null.init.alpha <- colSums(temp.repl.neg) / sum(temp.repl.neg)
    temp.null.init.disp <- 1

    repl.dirichlet.null <- optim_dirichlet_null_pars(temp.repl.neg,
      c(temp.null.init.disp, temp.null.init.alpha[1:(length(temp.null.init.alpha)-1)]))
    final.null.par.list[[rep]] <- repl.dirichlet.null

    repl.dirichlet.null.ll <- repl.dirichlet.null[1]
    repl.dirichlet.null.disp <- repl.dirichlet.null[2]
    repl.dirichlet.null.pars <- repl.dirichlet.null[3:length(repl.dirichlet.null)]

    if('plot_dirichlet_dispersion' %in% names(input.specs)){
      dirichlet_dispersion_behavior()
    }

    repl.null.ll <- c()
    repl.alt.ll <- c()
    repl.alt.pars <- c()
    repl.info <- c()

    for(chrom in 1:length(genome.regions)){
      temp.genome.region <- genome.regions[[chrom]]

      temp.dirichlet.null <- c()
      temp.dirichlet.alt <- c()

      if('nrFreeCores' %in% names(input.specs)){

        nr.cores <- max((detectCores() - input.specs$nrFreeCores), 1)
        cl <- makeCluster(nr.cores)
        setDefaultCluster(cl=cl)

        clusterExport(cl, c('dirichlet_null_ll', 'optim_dirichlet_alt_pars',
        'dirichlet_null_logLik', 'dirichlet_alt_logLik', "ddirmnom"), envir=environment())

        temp.dirichlet.null <- do.call(c, parLapply(cl, temp.genome.region$region_counts, function(x){
          repl.counts <- x[, temp.repl.cols]
          temp.chrom.dirichlet.out <- dirichlet_null_ll(repl.counts,
            c(repl.dirichlet.null.disp, repl.dirichlet.null.pars[1:(length(repl.dirichlet.null.pars) - 1)]))
          temp.chrom.dirichlet.out
          }))

        temp.dirichlet.alt <- do.call(rbind, parLapply(cl, temp.genome.region$region_counts, function(x){
          repl.counts <- x[, temp.repl.cols]
          temp.chrom.dirichlet.out <- optim_dirichlet_alt_pars(repl.counts,
            repl.dirichlet.null.pars[1:(length(repl.dirichlet.null.pars)-1)],
            repl.dirichlet.null.disp)
          temp.chrom.dirichlet.out
          }))

          stopCluster(cl)

        } else {
          temp.dirichlet.null <- do.call(c, lapply(temp.genome.region$region_counts, function(x){
            repl.counts <- x[, temp.repl.cols]
            temp.chrom.dirichlet.out <- dirichlet_null_ll(repl.counts,
              c(repl.dirichlet.null.disp, repl.dirichlet.null.pars[1:(length(repl.dirichlet.null.pars) - 1)]))
            temp.chrom.dirichlet.out
            }))

          temp.dirichlet.alt <- do.call(rbind, lapply(temp.genome.region$region_counts, function(x){
            repl.counts <- x[, temp.repl.cols]
            temp.chrom.dirichlet.out <- optim_dirichlet_alt_pars(repl.counts,
              repl.dirichlet.null.pars[1:(length(repl.dirichlet.null.pars)-1)],
              repl.dirichlet.null.disp)
            temp.chrom.dirichlet.out
            }))
        }

      repl.null.ll <- c(repl.null.ll, temp.dirichlet.null)
      repl.alt.ll <- c(repl.alt.ll, temp.dirichlet.alt[,1])
      repl.alt.pars <- rbind(repl.alt.pars, temp.dirichlet.alt[,c(2:ncol(temp.dirichlet.alt))])
      repl.info <- rbind(repl.info, temp.genome.region$region_info_df)
    }

    final.null.ll <- cbind(final.null.ll, repl.null.ll)
    final.alt.ll <- cbind(final.alt.ll, repl.alt.ll)

    repl.llr <- repl.null.ll - repl.alt.ll
    final.llr <- cbind(final.llr, repl.llr)

    total.df <- total.df + (length(temp.repl.cols) - 1)

    if('save_dirichlet_info' %in% names(input.specs)){
      save_dirichlet_repl_alt_pars(repl.alt.pars, repl.alt.ll, rep, input.specs)
    }

    parallel.runtime <- proc.time() - ptm
    print(paste0('GLMM Untrained Replicate runtime: ', round(parallel.runtime[3], 3) ))
  }

  final.pVal <- pchisq(2 * rowSums(final.llr), df = total.df, lower.tail=F)

  out.df <- data.frame(rawScores = rowSums(final.llr), formatScores = -log10(final.pVal),
  log2_rate_ratio = rep(1, length(final.pVal)), stringsAsFactors = F)

  out.df <- cbind(out.df, repl.info)

  if('save_dirichlet_info' %in% names(input.specs)){
    save_dirichlet_info(final.llr, final.null.par.list, final.null.ll, input.specs)
  }

  return(out.df)
}

save_dirichlet_repl_alt_pars <- function(input.repl.pars, input.repl.ll, repl.nr, input.specs){

  out.df <- cbind(input.repl.ll, input.repl.pars)
  colnames(out.df) <- c('neg_ll', paste0('alpha_', c(1:ncol(input.repl.pars)) ))
  write.csv(out.df, file = paste(input.specs$dataName, 'repl', repl.nr, 'alt_pars.csv', sep = '_'), row.names = FALSE)

}

save_dirichlet_info <- function(input.repl.alt.llr, input.null.pars, input.null.ll, input.specs){

  colnames(input.repl.alt.llr) <- paste0('llr_repl', c(1:ncol(input.repl.alt.llr)))
  write.csv(input.repl.alt.llr, file = paste(input.specs$dataName, 'repl_alt_llr.csv', sep = '_'), row.names = FALSE)

  colnames(input.null.ll) <- paste0('null_neg_ll_repl', c(1:ncol(input.null.ll)))
  write.csv(input.null.ll, file = paste(input.specs$dataName, 'repl_null_neg_ll.csv', sep = '_'), row.names = FALSE)

  for(i in 1:length(input.null.pars)){
    temp.df <- as.data.frame(t(input.null.pars[[i]]), stringsAsFactors = F)
    colnames(temp.df) <- c('neg_ll', 'dispersion', paste0('alpha_', c(1:(ncol(temp.df) - 2)) ))
    write.csv(temp.df, file = paste(input.specs$dataName, 'repl', i, '_null_pars.csv', sep = '_'), row.names = FALSE)
  }

}

# taken from: http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf
eexp <- function(x){
  if(x == 'LOGZERO'){
    0
  } else {
    exp(x)
  }
}

eln <- function(x){
  if(x == 0){
    'LOGZERO'
  } else {
    log(x)
  }
}

elnsum <- function(eln.x, eln.y){
  if(eln.x == 'LOGZERO' | eln.y == 'LOGZERO'){
    if(eln.x == 'LOGZERO'){
      eln.y
    } else {
      eln.x
    }
  } else {
    if(eln.x > eln.y){
      eln.x + eln(1 + exp(eln.y - eln.x))
    } else {
      eln.y + eln(1 + exp(eln.x - eln.y))
    }
  }
}

# during parallelization, , subfunction eln can't be called
# this function incorporates the if-else statement directly
elnsum_forParallel <- function(eln.x, eln.y){
  if(eln.x == 'LOGZERO' | eln.y == 'LOGZERO'){
    if(eln.x == 'LOGZERO'){
      eln.y
    } else {
      eln.x
    }
  } else {
    if(eln.x > eln.y){
      eln.x + (if((1 + exp(eln.y - eln.x)) == 0) 'LOGZERO' else (1 + exp(eln.y - eln.x))) #eln(1 + exp(eln.y - eln.x))
    } else {
      eln.y + (if((1 + exp(eln.x - eln.y)) == 0) 'LOGZERO' else (1 + exp(eln.x - eln.y))) # eln(1 + exp(eln.x - eln.y))
    }
  }
}

elnproduct <- function(eln.x, eln.y){
  if(eln.x == 'LOGZERO' | eln.y == 'LOGZERO'){
    'LOGZERO'
  } else {
    eln.x + eln.y
  }
}

emission_nonLog <- function(input.state, input.counts, input.par.list){

  if(input.state == 1){
    temp_pos_integration_function <- nb_integration_function_generator(input.counts,
      input.par.list$beta_pos_pred, input.par.list$model_tau, input.par.list$model_dispersion)
    return( integrate(temp_pos_integration_function, lower = -5*input.par.list$model_tau,
      upper = 5*input.par.list$model_tau)$value[1] )
  } else {
    temp_neg_integration_function <- nb_integration_function_generator(input.counts,
      input.par.list$beta_neg_pred, input.par.list$model_tau, input.par.list$model_dispersion)
    return( integrate(temp_neg_integration_function, lower = -5*input.par.list$model_tau,
      upper = 5*input.par.list$model_tau)$value[1] )
  }
}

frwdBkwd_hmm_tmb_glmm_analysis <- function(input.counts, input.info, input.specs){

  out.score.list <- list()

  sorted.counts <- input.counts[order(input.info$start),]
  sorted.info <- input.info[order(input.info$start),]

  glmm.counts <- sorted.counts[which(! is.na(sorted.info$start)), ]
  glmm.info <- sorted.info[which(! is.na(sorted.info$start)), ]

  rev.glmm.counts <- glmm.counts[seq(dim(glmm.counts)[1],1),]
  rev.glmm.info <- glmm.info[seq(dim(glmm.info)[1],1),]

  pos.counts <- sorted.counts[which(sorted.info$label %in% input.specs$glmm_positiveTraining), ]
  neg.counts <- sorted.counts[which(sorted.info$label %in% input.specs$glmm_negativeTraining), ]

  par.list <- list()

  # estimate nb_glmm_tmb parameters for each replicate
  for(i in 1:length(input.specs$repl_groups)){
    # estimate sample-specific, type specific paramters
    guide.type <- c(rep('pos', nrow(pos.counts)), rep('neg', nrow(neg.counts)))
    temp.label.counts <- rbind(pos.counts[,input.specs$repl_groups[[i]]], neg.counts[,input.specs$repl_groups[[i]]])
    par.list[[i]] <- estimate_glmm_tmb_nb_pars_jointly(temp.label.counts, guide.type, input.specs)
  }

  # max.probs[i, j]: contains the maximum obtainable probability of seeing the observed
  # state given the previous state and the transition and emission probabilites
  frwd.max.probs <- matrix(0, ncol = nrow(glmm.counts), nrow = 2)
  bkwd.max.probs <- matrix(0, ncol = nrow(glmm.counts), nrow = 2)

  transition.prob <- c()
  start.state.prob <- c()
  if(input.specs$transition_type == 'log'){
    transition.prob <- log( data.frame(input.specs$transition_prob, 1 - input.specs$transition_prob) )
    start.state.prob <- log(input.specs$hmmStartStateProbs)
  } else {
    transition.prob <- data.frame(input.specs$transition_prob, 1 - input.specs$transition_prob)
    start.state.prob <- input.specs$hmmStartStateProbs
  }

  # initialize hmm for forward and backward
  frwd.initial.prob <- vector(mode="numeric", length=length(input.specs$hmmStates))
  bkwd.initial.prob <- vector(mode="numeric", length=length(input.specs$hmmStates))
  for(state in 1:length(input.specs$hmmStates)){
    for(r in 1:length(input.specs$repl_groups)){
      frwd.temp.guide.counts <- unlist(glmm.counts[1, input.specs$repl_groups[[r]]])
      # technically: elnproduct(eln(frwd.initial.prob[state]), eln(emission_nonLog(state, frwd.temp.guide.counts, par.list[[r]])))
      # not sure how to initialize and add to inititalization, hence workaround below
      frwd.initial.prob[state] <- frwd.initial.prob[state] + eln(emission_nonLog(state, frwd.temp.guide.counts, par.list[[r]]))
      bkwd.temp.guide.counts <- unlist(rev.glmm.counts[1, input.specs$repl_groups[[r]]])
      bkwd.initial.prob[state] <- bkwd.initial.prob[state] + eln(emission_nonLog(state, bkwd.temp.guide.counts, par.list[[r]]))
    }
    frwd.max.probs[state, 1] <- elnproduct(start.state.prob[state], frwd.initial.prob[state])
    bkwd.max.probs[state, 1] <- elnproduct(start.state.prob[state], bkwd.initial.prob[state])
  }

  # for every guide(-pair) compute and record the highest probability of moving
  # from any previous state into a specific current state, given the observations
  for(g in 2:(nrow(glmm.counts) ) ){  # all possible guide(-pairs)
    for(state in 1:length(input.specs$hmmStates)){  # all the possible current states
      # frwd.temp.max.prob <- vector(mode="numeric", length=length(input.specs$hmmStates))  # contains the comulative probability of moving into state s across replicates
      # bkwd.temp.max.prob <- vector(mode="numeric", length=length(input.specs$hmmStates))
      frwd.temp.emiss.prob <- 0  # contains the comulative probability of moving into state s across replicates
      bkwd.temp.emiss.prob <- 0
      for(r in 1:length(input.specs$repl_groups)){  # all replicates
        frwd.temp.guide.counts <- unlist(glmm.counts[g, input.specs$repl_groups[[r]]])  # g-1 because initialized start, but also starting with counts
        frwd.temp.emiss.prob <- elnproduct(frwd.temp.emiss.prob, eln(emission_nonLog(state, frwd.temp.guide.counts, par.list[[r]])))
        bkwd.temp.guide.counts <- unlist(rev.glmm.counts[g, input.specs$repl_groups[[r]]])  # g-1 because initialized start, but also starting with counts
        bkwd.temp.emiss.prob <- elnproduct(bkwd.temp.emiss.prob, eln(emission_nonLog(state, bkwd.temp.guide.counts, par.list[[r]])))
      }
      frwd.temp.prev.trans.prob <- 'LOGZERO'  # contains the comulative probability of moving into state s across replicates
      bkwd.temp.prev.trans.prob <- 'LOGZERO'
      for(p in 1:length(input.specs$hmmStates)){  # all the possible previous states
        frwd.temp.prev.trans.prob <- elnsum(frwd.temp.prev.trans.prob, elnproduct(frwd.max.probs[p, (g - 1)], transition.prob[p, state]))
        bkwd.temp.prev.trans.prob <- elnsum(bkwd.temp.prev.trans.prob, elnproduct(bkwd.max.probs[p, (g - 1)], transition.prob[p, state]))
      }
      frwd.max.probs[state, g] <- elnproduct(frwd.temp.emiss.prob, frwd.temp.prev.trans.prob)
      # frwd.max.states[state, g-1] <- which(frwd.temp.max.prob == frwd.max.probs[state, g])[1]
      bkwd.max.probs[state, g] <- elnproduct(bkwd.temp.emiss.prob, bkwd.temp.prev.trans.prob)
      # bkwd.max.states[state, g-1] <- which(bkwd.temp.max.prob == bkwd.max.probs[state, g])[1]
    }
  }

  #########################################
  # check if sum(frwd.max.probs[, ncol(frwd.max.probs)])  == sum(bkwd.max.probs[, ncol(bkwd.max.probs)])
  #########################################
  frwd.prob <- sum(frwd.max.probs[, ncol(frwd.max.probs)])
  bkwd.prob <- sum(bkwd.max.probs[, ncol(bkwd.max.probs)])

  # flip the backward probs:
  flip.frwd.max.probs <- bkwd.max.probs[, seq(dim(bkwd.max.probs)[2],1)]

  posterior.probs <- frwd.max.probs + flip.frwd.max.probs - frwd.prob
  posterior.probs.ratio <- posterior.probs[1,] - posterior.probs[2,]

  # since there are no posterior probs for the non-targeting guides their
  # probabilities are calculated using the regular glmm log-likelihood ratio

  neg.counts <- sorted.counts[which( is.na(sorted.info$start)), ]
  neg.info <- sorted.info[which( is.na(sorted.info$start)), ]

  guide.pos.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(neg.counts)))
  guide.neg.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(neg.counts)))
  for(i in 1:nrow(neg.counts)){
    for(j in 1:length(input.specs$repl_groups)){
      temp.guide.counts <- unlist(neg.counts[i,input.specs$repl_groups[[j]]])
      temp_pos_integration_function <- nb_integration_function_generator(temp.guide.counts,
        par.list[[j]]$beta_pos_pred, par.list[[j]]$model_tau, par.list[[j]]$model_dispersion)
      temp_neg_integration_function <- nb_integration_function_generator(temp.guide.counts,
        par.list[[j]]$beta_neg_pred, par.list[[j]]$model_tau, par.list[[j]]$model_dispersion)
      guide.pos.ll[i,j] <- log(integrate(temp_pos_integration_function, lower = -5*par.list[[j]]$model_tau,
        upper = 5*par.list[[j]]$model_tau)$value[1])
      guide.neg.ll[i,j] <- log(integrate(temp_neg_integration_function, lower = -5*par.list[[j]]$model_tau,
        upper = 5*par.list[[j]]$model_tau)$value[1])
    }
  }

  adjusted.guide.pos.ll <- guide.pos.ll
  adjusted.guide.neg.ll <- guide.neg.ll

  # check below when debugging: how many fall in this category?
  for(i in 1:ncol(adjusted.guide.pos.ll)){
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == -Inf),i] <- min(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] > -Inf), i])
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == Inf),i] <- max(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] < Inf), i])

    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == -Inf),i] <- min(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] > -Inf), i])
    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == Inf),i] <- max(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] < Inf), i])
  }

  total.guide.pos.ll <- rowSums(adjusted.guide.pos.ll)
  total.guide.neg.ll <- rowSums(adjusted.guide.neg.ll)

  ll.guide.ratio <- total.guide.pos.ll - total.guide.neg.ll

  out.prob <- c(posterior.probs.ratio, ll.guide.ratio)

  par.summary <- summary_glmm_nb_pars(par.list)
  write.csv(par.summary, file = paste0(input.specs$dataName,'_glmm_nb_tmb_frwdBkwd_Summary.csv') , row.names = F)

  out.df <- cbind.data.frame(rawScores = out.prob, formatScores = out.prob,
    log2_rate_ratio = rep(1, length(out.prob)), stringsAsFactors = FALSE)
  #add the info part
  out.df <- cbind(out.df, sorted.info, stringsAsFactors = FALSE)

  return(out.df)

}

# states: 1 == enhancer, 2 == non-enhancer
# given a state, calculate the probability of emitting the observed counts
emission <- function(input.state, input.counts, input.par.list){

  if(input.state == 1){
    temp_pos_integration_function <- nb_integration_function_generator(input.counts,
      input.par.list$beta_pos_pred, input.par.list$model_tau, input.par.list$model_dispersion)
    return( log(integrate(temp_pos_integration_function, lower = -5*input.par.list$model_tau,
      upper = 5*input.par.list$model_tau)$value[1]) )
  } else {
    temp_neg_integration_function <- nb_integration_function_generator(input.counts,
      input.par.list$beta_neg_pred, input.par.list$model_tau, input.par.list$model_dispersion)
    return( log(integrate(temp_neg_integration_function, lower = -5*input.par.list$model_tau,
      upper = 5*input.par.list$model_tau)$value[1]) )
  }
}

viterbi_hmm_tmb_glmm_analysis <- function(input.counts, input.info, input.specs){

  out.score.list <- list()

  sorted.counts <- input.counts[order(input.info$start),]
  sorted.info <- input.info[order(input.info$start),]

  glmm.counts <- sorted.counts[which(! is.na(sorted.info$start)), ]
  glmm.info <- sorted.info[which(! is.na(sorted.info$start)), ]

  pos.counts <- sorted.counts[which(sorted.info$label %in% input.specs$glmm_positiveTraining), ]
  neg.counts <- sorted.counts[which(sorted.info$label %in% input.specs$glmm_negativeTraining), ]

  par.list <- list()

  # estimate nb_glmm_tmb parameters for each replicate
  for(i in 1:length(input.specs$repl_groups)){
    # estimate sample-specific, type specific paramters
    guide.type <- c(rep('pos', nrow(pos.counts)), rep('neg', nrow(neg.counts)))
    temp.label.counts <- rbind(pos.counts[,input.specs$repl_groups[[i]]], neg.counts[,input.specs$repl_groups[[i]]])
    par.list[[i]] <- estimate_glmm_tmb_nb_pars_jointly(temp.label.counts, guide.type, input.specs)
  }

  # max.states[i, j]: contains the numerical label of previous state, j-1, which maximized
  # the probability to move into this state j
  max.states <- matrix(0, ncol = nrow(glmm.counts) - 1, nrow = 2)  # account for initial state when hmm is started, hence g+1 length

  # max.probs[i, j]: contains the maximum obtainable probability of seeing the observed
  # state given the previous state and the transition and emission probabilites
  max.probs <- matrix(0, ncol = nrow(glmm.counts), nrow = 2)

  transition.prob <- c()
  start.state.prob <- c()
  if(input.specs$transition_type == 'log'){
    transition.prob <- log( data.frame(input.specs$transition_prob, 1 - input.specs$transition_prob) )
    start.state.prob <- log(input.specs$hmmStartStateProbs)
  } else {
    transition.prob <- data.frame(input.specs$transition_prob, 1 - input.specs$transition_prob)
    start.state.prob <- input.specs$hmmStartStateProbs
  }

  # initialize hmm
  initial.prob <- vector(mode="numeric", length=length(input.specs$hmmStates))
  for(state in 1:length(input.specs$hmmStates)){
    for(r in 1:length(input.specs$repl_groups)){
      temp.guide.counts <- unlist(glmm.counts[1, input.specs$repl_groups[[r]]])
      initial.prob[state] <- initial.prob[state] + emission(state, temp.guide.counts, par.list[[r]])
    }
    max.probs[state, 1] <- start.state.prob[state] + initial.prob[state]
  }

  # for every guide(-pair) compute and record the highest probability of moving
  # from any previous state into a specific current state, given the observations
  for(g in 2:(nrow(glmm.counts) ) ){ # all possible guide(-pairs)
    for(state in 1:length(input.specs$hmmStates)){ # all the possible current states
      temp.max.prob <- vector(mode="numeric", length=length(input.specs$hmmStates))  # contains the comulative probability of moving into state s across replicates
      for(p in 1:length(input.specs$hmmStates)){ # all the possible previous states
        for(r in 1:length(input.specs$repl_groups)){ # all replicates
          temp.guide.counts <- unlist(glmm.counts[g, input.specs$repl_groups[[r]]])  # g-1 because initialized start, but also starting with counts
          temp.max.prob[p] <- temp.max.prob[p] + emission(state, temp.guide.counts, par.list[[r]])
        }
        temp.max.prob[p] <- temp.max.prob[p] + max.probs[p, (g - 1)] + transition.prob[p, state]
      }
      max.probs[state, g] <- max(temp.max.prob)
      max.states[state, g-1] <- which(temp.max.prob == max.probs[state, g])[1]
    }
  }

  # no sink, the max(max.prob[,ncol(max.prob)]) represents the last state
  # use the above to then step through max.states and determine the sequence of
  # hidden states
  # because: max.states[state, g-1] the first hidden state is the first state of max.states
  out.states <- viterbi_backtrack(max.probs, max.states)

  guide.pos.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(sorted.counts)))
  guide.neg.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(sorted.counts)))
  for(i in 1:nrow(sorted.counts)){
    for(j in 1:length(input.specs$repl_groups)){
      temp.guide.counts <- unlist(sorted.counts[i,input.specs$repl_groups[[j]]])
      temp_pos_integration_function <- nb_integration_function_generator(temp.guide.counts,
        par.list[[j]]$beta_pos_pred, par.list[[j]]$model_tau, par.list[[j]]$model_dispersion)
      temp_neg_integration_function <- nb_integration_function_generator(temp.guide.counts,
        par.list[[j]]$beta_neg_pred, par.list[[j]]$model_tau, par.list[[j]]$model_dispersion)
      guide.pos.ll[i,j] <- log(integrate(temp_pos_integration_function, lower = -5*par.list[[j]]$model_tau,
        upper = 5*par.list[[j]]$model_tau)$value[1])
      guide.neg.ll[i,j] <- log(integrate(temp_neg_integration_function, lower = -5*par.list[[j]]$model_tau,
        upper = 5*par.list[[j]]$model_tau)$value[1])
    }
  }

  adjusted.guide.pos.ll <- guide.pos.ll
  adjusted.guide.neg.ll <- guide.neg.ll

  # check below when debugging: how many fall in this category?
  for(i in 1:ncol(adjusted.guide.pos.ll)){
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == -Inf),i] <- min(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] > -Inf), i])
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == Inf),i] <- max(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] < Inf), i])

    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == -Inf),i] <- min(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] > -Inf), i])
    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == Inf),i] <- max(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] < Inf), i])
  }

  total.guide.pos.ll <- rowSums(adjusted.guide.pos.ll)
  total.guide.neg.ll <- rowSums(adjusted.guide.neg.ll)

  ll.guide.ratio <- total.guide.pos.ll - total.guide.neg.ll

  par.summary <- summary_glmm_nb_pars(par.list)
  write.csv(par.summary, file = paste0(input.specs$dataName,'_glmm_nb_tmb_viterbi_Summary.csv') , row.names = F)

  out.df <- cbind.data.frame(rawScores = ll.guide.ratio, formatScores = ll.guide.ratio,
    log2_rate_ratio = rep(1, length(ll.guide.ratio)), stringsAsFactors = FALSE)
  #add the info part
  out.df <- cbind(out.df, sorted.info, stringsAsFactors = FALSE)

  out.df$viterbiLabels <- c(out.states, rep(2, length(which(is.na(sorted.info$start)))))
  out.df$origLabel <- out.df$label
  viterbi.state.labels <- rep('nonEnh', nrow(out.df))
  viterbi.state.labels[which(out.df$viterbiLabels == 1)] <- 'enh'
  out.df$label <- paste(out.df$origLabel, viterbi.state.labels, sep = '_')

  return(out.df)
}

viterbi_backtrack <- function(input.probs, input.states){

  best.state <- which(input.probs[, ncol(input.probs)] == max(input.probs[, ncol(input.probs)]))
  state.sequence <- vector('numeric', length = ncol(input.probs) - 1)
  state.index <- 1
  state.sequence[state.index] <- best.state
  state.index <- state.index + 1
  for(i in (ncol(input.states)):1){ # don't incorporate the start state
    state.sequence[state.index] <- input.states[best.state, i]
    best.state <- input.states[best.state, i]
    state.index <- state.index + 1
  }

  return(rev(state.sequence))
}

# this estimates the parameters for the positive and negative controls
# input.type describes if the guide is from a positive pool or not
# input.counts contain first all positive, then all negative guides!
estimate_glmm_tmb_nb_pars_allData <- function(input.counts, input.specs){

  pool.IDs <- colnames(input.counts)
  guide.nrs <- c(1:nrow(input.counts))
  guide.names <- paste0('guide', guide.nrs)
  pool.nrs <- c(1:ncol(input.counts))
  pool.names <- paste0('pool', pool.nrs)

  guide.labels <- rep(guide.names, each = length(pool.names))
  pool.labels <- rep(pool.names, length(guide.names))

  long.df <- cbind.data.frame(observed = c(t(input.counts)), guide_labels = guide.labels,
    pool_labels = pool.labels, stringsAsFactors = T)

  tmb.nb2.model <- glmmTMB(observed ~ pool_labels + (1|guide_labels),
    data = long.df, family=nbinom2)

  model.dispersion <- sigma(tmb.nb2.model) # residual variance, theta corresponding to the NB variance
  model.tau <- exp(getME(tmb.nb2.model, name = 'theta'))  # SD of the random effect
  #  no idea why this works, no idea how else to get random effect SD

  nr.beta.neg <- length(unique(long.df$pool_labels))
  nr.total.beta <- nr.beta.neg #+ nr.beta.pos
  tmb.betas.orig <- tmb.nb2.model$fit$par[1:nr.total.beta]
  beta.neg <- tmb.betas.orig[1:nr.beta.neg] + c(0, rep(tmb.betas.orig[1], nr.beta.neg - 1))

  tmb.neg_logLik <- logLik(tmb.nb2.model)[1]

  return(list(betas = tmb.betas.orig[1:nr.beta.neg],
    beta_pred = beta.neg, model_tau = model.tau, model_dispersion = model.dispersion,
    neg_logLik = tmb.neg_logLik, pool_IDs = pool.IDs))

}

category_optimization_function_parallel <- function(par, input.counts, input.default.ll,
  input.specs, input.models, integration_fun, log_sum_fun){ #, input.beta.list

  # print('category_optimization_function')
  nr.categories <- length(input.specs$categoryProbs)
  # nr.betas.per.repl <- length(input.specs$repl_groups[[1]]) - 1
  # compute the ll for each category, for each guide separately
  per.guide.lls <- data.frame(matrix(ncol = nr.categories, nrow = nrow(input.counts)))

  # need to normalize the category probabilities to sum to 1
  temp.cat.probs <- par[1:nr.categories]
  total.cat.probs <- sum(temp.cat.probs)
  cat.prob <- par[1:nr.categories]/total.cat.probs

  # extract the betas
  cat.beta.probs.orig <- par[(nr.categories + 1):length(par)]
  repl.beta.list <- list()

  cat.beta.probs <- cat.beta.probs.orig

  for(repl in 1:length(input.specs$repl_groups)){
    repl.beta.list[[repl]] <- list()
    # repl.beta.list[[repl]][[1]] <- par.list[[repl]]$beta_pred
    temp.nr.repl.beta <- length(input.specs$repl_groups[[repl]]) - 1 # input is identical, so not optimized over
    for(cat in 1:(nr.categories - 1) ){
      tmp.beta.pos <- c(1:temp.nr.repl.beta) #+ (cat - 2) * temp.nr.repl.beta #nr.par.betas.per.repl * (cat - 2) + (nr.par.betas.per.repl * (nr.categories - 1)) * (repl - 1)
      temp.betas <- cat.beta.probs[tmp.beta.pos]
      cat.beta.probs <- cat.beta.probs[-tmp.beta.pos] # remove extracted betas
      # add beta_0
      temp.betas <- c(input.models[[repl]]$beta_pred[1], temp.betas)
      repl.beta.list[[repl]][[cat]] <- temp.betas
    }
  }

  # for each guide
    # for each category
      # multiply the probability of being in the category with the probability of
      # the observed data given the category fixed effects, across all replicates
  for(i in 1:nrow(input.counts)){
    # the first category is already precomputed because it's the default one
    per.guide.lls[i, 1] <- log(cat.prob[1])
    per.guide.lls[i, 1] <- per.guide.lls[i, 1] + sum(input.default.ll[i, ])

    for(ct in 2:nr.categories){
      per.guide.lls[i, ct] <- log(cat.prob[ct])
      for(repl in 1:length(input.specs$repl_groups)){
        temp.guide.counts <- unlist(input.counts[i,input.specs$repl_groups[[repl]]])
        temp.category.betas <- repl.beta.list[[repl]][[ct-1]]
        temp_integration_function <- integration_fun(temp.guide.counts,
          temp.category.betas, input.models[[repl]]$model_tau, input.models[[repl]]$model_dispersion)
        temp.category.ll <- log(integrate(temp_integration_function, lower = -5*input.models[[repl]]$model_tau,
          upper = 5*input.models[[repl]]$model_tau)$value[1])
        per.guide.lls[i, ct] <- per.guide.lls[i, ct] + temp.category.ll
      }
    }
  }

  # sum across the categories
  out.cat.sums <- vector('numeric', length = nrow(per.guide.lls))
  for(row in 1:nrow(input.counts)){
    out.cat.sums[row] <- per.guide.lls[row, 1]
    for(cat in 2:nr.categories){
      out.cat.sums[row] <- log_sum_fun(out.cat.sums[row], per.guide.lls[row, cat])
    }
  }

  sum(out.cat.sums)
}

allData_nb_tmb_glmm_categoryAnalysis_parallel <- function(input.counts, input.info, input.specs){

  out.score.list <- list()

  glmm.info <- input.info
  glmm.counts <- input.counts

  par.list <- list()

  print('computing per-replicate parameters')
  # estimate the parameters for all the data
  for(i in 1:length(input.specs$repl_groups)){
    # estimate sample-specific, type specific paramters
    temp.label.counts <- glmm.counts[,input.specs$repl_groups[[i]]]
    par.list[[i]] <- estimate_glmm_tmb_nb_pars_allData(temp.label.counts, input.specs)
  }

  # optimize over the probabilities of being in a category (gamma)
  # optimize over the alpha parameters for each category for n-1 pools
  print('computing default category ll')
  default.category.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(glmm.counts)))
  for(i in 1:nrow(glmm.counts)){
    for(j in 1:length(input.specs$repl_groups)){
      temp.guide.counts <- unlist(glmm.counts[i,input.specs$repl_groups[[j]]])
      temp_neg_integration_function <- nb_integration_function_generator(temp.guide.counts,
        par.list[[j]]$beta_pred, par.list[[j]]$model_tau, par.list[[j]]$model_dispersion)
      default.category.ll[i,j] <- log(integrate(temp_neg_integration_function, lower = -5*par.list[[j]]$model_tau,
        upper = 5*par.list[[j]]$model_tau)$value[1])
    }
  }

  nr.categories <- length(input.specs$categoryProbs)

  # structure:
  #   repl1: 'category 2 beta', 'category 3 beta' ..
  #   repl2: 'category 2 beta', 'category 3 beta' ..
  optim.beta.init <- c()
  for(repl in 1:length(input.specs$repl_groups)){
    temp.repl.betas <- par.list[[repl]]$beta_pred[2:length(par.list[[repl]]$beta_pred)]
    optim.beta.init <- c(optim.beta.init, rep(temp.repl.betas, nr.categories - 1))
  }

  lower.catProbs <- rep(0.001, length(input.specs$categoryProbs))
  upper.catProbs <- rep(0.999, length(input.specs$categoryProbs))
  lower.betas <- rep(0, length(optim.beta.init))
  max.count <- max(glmm.counts)
  upper.betas <- rep(log(max.count), length(optim.beta.init))

  print('setting up parallelization')

  nr.cores <- max((detectCores() - input.specs$nrFreeCores), 1)
  cl <- makeCluster(nr.cores)
  setDefaultCluster(cl=cl)

  print('optimizing category pars..')

  ptm <- proc.time()
  optim.out <- optimParallel(par = c(input.specs$categoryProbs, optim.beta.init), category_optimization_function_parallel,
    input.counts = glmm.counts, input.default.ll = default.category.ll, input.models = par.list,
    integration_fun = nb_integration_function_generator, log_sum_fun = elnsum_forParallel,
    input.specs = input.specs, method = 'L-BFGS-B', lower = c(lower.catProbs, pmax(lower.betas, optim.beta.init - 2) ),
    upper = c(upper.catProbs, optim.beta.init + 2), control=list(fnscale=-1))
  parallelOptim.runtime <- proc.time() - ptm

  print(paste0('optimParallel runtime: ', parallelOptim.runtime[3]))
  stopCluster(cl)

  optim.out.pars <- optim.out$par
  optim.out.value <- optim.out$value

  # need to nromalize category probabilities to 1
  optim.cat.probs <- optim.out.pars[c(1:nr.categories)] / sum(optim.out.pars[c(1:nr.categories)])

  # extract betas
  # pars: gamma_1, gamma_2, .., beta_1_repl1_cat1, beta_2_repl1_cat1, beta_1_repl1_cat2, beta_2_repl1_cat2 |
  # beta_1_repl2_cat1, beta_2_repl2_cat1, beta_1_repl2_cat2, beta_2_repl2_cat2
  optim.cat.beta.list <- list()
  nr.par.betas.per.repl <- length(input.specs$repl_groups[[1]]) - 1 # nr of beta paremeter estimated per repl, per pool
  cat.beta.probs.orig <- optim.out.pars[(nr.categories + 1):length(optim.out.pars)]
  cat.beta.probs <- cat.beta.probs.orig

  #first category is default category
  for(repl in 1:length(input.specs$repl_groups)){
    optim.cat.beta.list[[repl]] <- list()
    optim.cat.beta.list[[repl]][[1]] <- par.list[[repl]]$beta_pred
    temp.nr.repl.beta <- length(par.list[[repl]]$beta_pred) - 1 # input is identical, so not optimized over
    for(cat in 2:nr.categories){
      tmp.beta.pos <- c(1:temp.nr.repl.beta)
      temp.betas <- cat.beta.probs[tmp.beta.pos]
      cat.beta.probs <- cat.beta.probs[-tmp.beta.pos] # remove extracted betas
      # add beta_0
      temp.betas.complete <- c(par.list[[repl]]$beta_pred[1], temp.betas)
      optim.cat.beta.list[[repl]][[cat]] <- temp.betas.complete
    }
  }

  print('computing per-category likelihood')
  guide.category.ll <- data.frame(matrix(ncol = nr.categories, nrow = nrow(glmm.counts), 0))

  for(i in 1:nrow(glmm.counts)){
    for(j in 1:length(input.specs$repl_groups)){
      temp.guide.counts <- unlist(glmm.counts[i,input.specs$repl_groups[[j]]])
      for(cat in 1:nr.categories){
        temp.betas <- optim.cat.beta.list[[j]][[cat]]
        temp_integration_function <- nb_integration_function_generator(temp.guide.counts,
          temp.betas, par.list[[j]]$model_tau, par.list[[j]]$model_dispersion)
        guide.cat.ll <- log(integrate(temp_integration_function, lower = -5*par.list[[j]]$model_tau,
          upper = 5*par.list[[j]]$model_tau)$value[1]) + log(optim.cat.probs[cat])

        guide.category.ll[i, cat] <- guide.category.ll[i, cat] + guide.cat.ll
      }
    }
  }

  adjusted.guide.category.ll <- guide.category.ll
  for(i in 1:ncol(adjusted.guide.category.ll)){
    adjusted.guide.category.ll[which(adjusted.guide.category.ll[,i] == -Inf),i] <- min(adjusted.guide.category.ll[which(adjusted.guide.category.ll[,i] > -Inf), i])
    adjusted.guide.category.ll[which(adjusted.guide.category.ll[,i] == Inf),i] <- max(adjusted.guide.category.ll[which(adjusted.guide.category.ll[,i] < Inf), i])
  }

  # as initial outpur, give category 2 vs default (category 1), add remaining ones at end of df
  colnames(adjusted.guide.category.ll) <- paste(rep('category', nr.categories),
    c(1:nr.categories), sep = '_')

  out.df <- cbind.data.frame(rawScores = adjusted.guide.category.ll$category_2 - adjusted.guide.category.ll$category_1,
    formatScores = adjusted.guide.category.ll$category_2 - adjusted.guide.category.ll$category_1,
    log2_rate_ratio = rep(1, nrow(adjusted.guide.category.ll)), stringsAsFactors = FALSE)
  #add the info part
  out.df <- cbind(out.df, glmm.info, stringsAsFactors = FALSE)

  out.df[names(adjusted.guide.category.ll)] <- adjusted.guide.category.ll[,c(1:nr.categories)]

  summary.df <- summary_glmm_nb_categoryPars(par.list, optim.cat.beta.list, optim.cat.probs, input.specs)
  write.csv(summary.df, file = paste0(input.specs$dataName,'_glmm_allData_nb_catParll_tmb_Summary.csv') , row.names = F)

  return(out.df)
}

summary_glmm_nb_categoryPars <- function(input.model.list, input.beta.list, input.cat.probs, input.specs){

  longest.betas <- max(unlist(lapply(input.beta.list, function(x){length(x[[1]])})))

  out.df <- c()
  for(repl in 1:length(input.specs$repl_groups)){
    temp.info <- c(round(input.model.list[[repl]]$model_tau, 3),
    round(input.model.list[[repl]]$model_dispersion, 3))
    for(cat in 1:length(input.specs$categoryProbs)){
      temp.cat.info <- c(temp.info, input.cat.probs[cat], round(input.beta.list[[repl]][[cat]], 3))
      length(temp.cat.info) <- longest.betas + 3
      out.df <- rbind(out.df, temp.cat.info)
    }
  }

  colnames(out.df) <- c('tau', 'dispersion', 'category_prob', 'input_beta', paste('beta', c(2:longest.betas), sep = '_'))
  return(out.df)
}

reducedData_poisson_multiDim_laplace_tmb_glmm_analysis_acrossRepl <- function(input.counts, input.info, input.specs){

  out.score.list <- list()

  glmm.info <- input.info
  glmm.counts <- input.counts

  start.time <- proc.time()

  pos.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_reducedData_positiveTraining), ]
  full.neg.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_reducedData_negativeTraining), ]
  rand.rows <- sample(c(1:nrow(full.neg.counts)), round(nrow(full.neg.counts) * input.specs$glmm_neg_trainingFraction))
  neg.counts <- full.neg.counts[rand.rows, ]

  par.list <- list()

  for(i in 1:length(input.specs$multivar_repl_groups)){
    # estimate sample-specific, type specific paramters
    guide.type <- c(rep('pos', nrow(pos.counts)), rep('neg', nrow(neg.counts)))
    temp.label.counts <- rbind(pos.counts[,input.specs$multivar_repl_groups[[i]]], neg.counts[,input.specs$multivar_repl_groups[[i]]])
    par.list[[i]] <- estimate_multivar_glmm_tmb_poisson_pars(temp.label.counts, guide.type, input.specs)
  }

  par.est.time <- proc.time()

  print('allData_poisson_multiDim_laplace_tmb_glmm_analysis_acrossRepl')
  par.diff <- par.est.time - start.time
  print(paste0('time to estimate pars: ', par.diff[3]))

  guide.pos.ll <- data.frame(matrix(ncol = length(input.specs$multivar_repl_groups), nrow = nrow(glmm.counts)))
  guide.neg.ll <- data.frame(matrix(ncol = length(input.specs$multivar_repl_groups), nrow = nrow(glmm.counts)))

  # browser()

  for(i in 1:nrow(glmm.counts)){
    # print(paste('processed row:', i, 'out of', nrow(glmm.counts), sep = ' '))
    for(j in 1:length(input.specs$multivar_repl_groups)){
      temp.guide.counts <- unlist(glmm.counts[i,input.specs$multivar_repl_groups[[j]]])
      temp.log.const.two <- sum(- lfactorial(temp.guide.counts))
      guide.pos.ll[i,j] <- laplace_multiDimensional(temp.guide.counts, par.list[[j]]$beta_pos_pred,
        par.list[[j]]$model_cov, par.list[[j]]$pool_sds, par.list[[j]]$inverted_sigma,
        par.list[[j]]$log_const_one, temp.log.const.two)
      guide.neg.ll[i,j] <- laplace_multiDimensional(temp.guide.counts, par.list[[j]]$beta_neg_pred,
        par.list[[j]]$model_cov, par.list[[j]]$pool_sds, par.list[[j]]$inverted_sigma,
        par.list[[j]]$log_const_one, temp.log.const.two)
    }
  }

  log.calc.time <- proc.time()
  calc.diff <- log.calc.time - par.est.time
  print(paste0('time to calculate loglik: ', calc.diff[3]))

  adjusted.guide.pos.ll <- guide.pos.ll
  adjusted.guide.neg.ll <- guide.neg.ll

  # browser()
  # check below when debugging: how many fall in this category?
  for(i in 1:ncol(adjusted.guide.pos.ll)){
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == -Inf),i] <- min(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] > -Inf), i])
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == Inf),i] <- max(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] < Inf), i])

    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == -Inf),i] <- min(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] > -Inf), i])
    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == Inf),i] <- max(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] < Inf), i])
  }

  total.guide.pos.ll <- rowSums(adjusted.guide.pos.ll)
  total.guide.neg.ll <- rowSums(adjusted.guide.neg.ll)

  ll.guide.ratio <- total.guide.pos.ll - total.guide.neg.ll

  par.summary <- summary_glmm_multiPoisson_pars(par.list)
  write.csv(par.summary, file = paste0(input.specs$dataName,'_glmm_reducedData_acrossRepl_multiPois_tmb_Summary.csv') , row.names = F)

  var.covar <- c()
  for(i in 1:length(par.list)){
    var.covar <- rbind(var.covar, par.list[[i]]$model_cov, rep(0, nrow(par.list[[i]]$model_cov)))
  }
  write.csv(var.covar, file = paste0(input.specs$dataName,'_glmm_reducedData_acrossRepl_multiPois_Covar_tmb_Summary.csv') , row.names = F)

  out.df <- cbind.data.frame(rawScores = ll.guide.ratio, formatScores = ll.guide.ratio,
    log2_rate_ratio = rep(1, length(ll.guide.ratio)), stringsAsFactors = FALSE)
  #add the info part
  out.df <- cbind(out.df, glmm.info, stringsAsFactors = FALSE)

  return(out.df)

}

reducedData_poisson_multiDim_laplace_tmb_glmm_analysis_perRepl <- function(input.counts, input.info, input.specs){

  out.score.list <- list()

  glmm.info <- input.info
  glmm.counts <- input.counts

  start.time <- proc.time()

  pos.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_reducedData_positiveTraining), ]
  full.neg.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_reducedData_negativeTraining), ]
  rand.rows <- sample(c(1:nrow(full.neg.counts)), round(nrow(full.neg.counts) * input.specs$glmm_neg_trainingFraction))
  neg.counts <- full.neg.counts[rand.rows, ]

  par.list <- list()

  for(i in 1:length(input.specs$repl_groups)){
    # estimate sample-specific, type specific paramters
    guide.type <- c(rep('pos', nrow(pos.counts)), rep('neg', nrow(neg.counts)))
    temp.label.counts <- rbind(pos.counts[,input.specs$repl_groups[[i]]], neg.counts[,input.specs$repl_groups[[i]]])
    par.list[[i]] <- estimate_multivar_glmm_tmb_poisson_pars(temp.label.counts, guide.type, input.specs)
  }

  par.est.time <- proc.time()

  print('allData_poisson_multiDim_laplace_tmb_glmm_analysis_perRepl')
  par.diff <- par.est.time - start.time
  print(paste0('time to estimate pars: ', par.diff[3]))

  guide.pos.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(glmm.counts)))
  guide.neg.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(glmm.counts)))

  # browser()

  for(i in 1:nrow(glmm.counts)){
    # print(paste('processed row:', i, 'out of', nrow(glmm.counts), sep = ' '))
    for(j in 1:length(input.specs$repl_groups)){
      temp.guide.counts <- unlist(glmm.counts[i,input.specs$repl_groups[[j]]])
      temp.log.const.two <- sum(- lfactorial(temp.guide.counts))
      guide.pos.ll[i,j] <- laplace_multiDimensional(temp.guide.counts, par.list[[j]]$beta_pos_pred,
        par.list[[j]]$model_cov, par.list[[j]]$pool_sds, par.list[[j]]$inverted_sigma,
        par.list[[j]]$log_const_one, temp.log.const.two)
      guide.neg.ll[i,j] <- laplace_multiDimensional(temp.guide.counts, par.list[[j]]$beta_neg_pred,
        par.list[[j]]$model_cov, par.list[[j]]$pool_sds, par.list[[j]]$inverted_sigma,
        par.list[[j]]$log_const_one, temp.log.const.two)
    }
  }

  log.calc.time <- proc.time()
  calc.diff <- log.calc.time - par.est.time
  print(paste0('time to calculate loglik: ', calc.diff[3]))

  adjusted.guide.pos.ll <- guide.pos.ll
  adjusted.guide.neg.ll <- guide.neg.ll

  # browser()
  # check below when debugging: how many fall in this category?
  for(i in 1:ncol(adjusted.guide.pos.ll)){
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == -Inf),i] <- min(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] > -Inf), i])
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == Inf),i] <- max(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] < Inf), i])

    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == -Inf),i] <- min(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] > -Inf), i])
    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == Inf),i] <- max(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] < Inf), i])
  }

  total.guide.pos.ll <- rowSums(adjusted.guide.pos.ll)
  total.guide.neg.ll <- rowSums(adjusted.guide.neg.ll)

  ll.guide.ratio <- total.guide.pos.ll - total.guide.neg.ll

  par.summary <- summary_glmm_multiPoisson_pars(par.list)
  write.csv(par.summary, file = paste0(input.specs$dataName,'_glmm_reducedData_perRepl_multiPois_tmb_Summary.csv') , row.names = F)

  var.covar <- c()
  for(i in 1:length(par.list)){
    var.covar <- rbind(var.covar, par.list[[i]]$model_cov, rep(0, nrow(par.list[[i]]$model_cov)))
  }
  write.csv(var.covar, file = paste0(input.specs$dataName,'_glmm_reducedData_perRepl_multiPois_Covar_tmb_Summary.csv') , row.names = F)

  out.df <- cbind.data.frame(rawScores = ll.guide.ratio, formatScores = ll.guide.ratio,
    log2_rate_ratio = rep(1, length(ll.guide.ratio)), stringsAsFactors = FALSE)
  #add the info part
  out.df <- cbind(out.df, glmm.info, stringsAsFactors = FALSE)

  return(out.df)

}

allData_poisson_multiDim_laplace_tmb_glmm_analysis_acrossRepl <- function(input.counts, input.info, input.specs){

  out.score.list <- list()

  glmm.info <- input.info
  glmm.counts <- input.counts

  start.time <- proc.time()

  pos.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_allData_positiveTraining), ]
  neg.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_allData_negativeTraining), ]

  par.list <- list()

  for(i in 1:length(input.specs$multivar_repl_groups)){
    # estimate sample-specific, type specific paramters
    guide.type <- c(rep('pos', nrow(pos.counts)), rep('neg', nrow(neg.counts)))
    temp.label.counts <- rbind(pos.counts[,input.specs$multivar_repl_groups[[i]]], neg.counts[,input.specs$multivar_repl_groups[[i]]])
    par.list[[i]] <- estimate_multivar_glmm_tmb_poisson_pars(temp.label.counts, guide.type, input.specs)
  }

  par.est.time <- proc.time()

  print('allData_poisson_multiDim_laplace_tmb_glmm_analysis_acrossRepl')
  par.diff <- par.est.time - start.time
  print(paste0('time to estimate pars: ', par.diff[3]))

  guide.pos.ll <- data.frame(matrix(ncol = length(input.specs$multivar_repl_groups), nrow = nrow(glmm.counts)))
  guide.neg.ll <- data.frame(matrix(ncol = length(input.specs$multivar_repl_groups), nrow = nrow(glmm.counts)))

  # browser()

  for(i in 1:nrow(glmm.counts)){
    # print(paste('processed row:', i, 'out of', nrow(glmm.counts), sep = ' '))
    for(j in 1:length(input.specs$multivar_repl_groups)){
      temp.guide.counts <- unlist(glmm.counts[i,input.specs$multivar_repl_groups[[j]]])
      temp.log.const.two <- sum(- lfactorial(temp.guide.counts))
      guide.pos.ll[i,j] <- laplace_multiDimensional(temp.guide.counts, par.list[[j]]$beta_pos_pred,
        par.list[[j]]$model_cov, par.list[[j]]$pool_sds, par.list[[j]]$inverted_sigma,
        par.list[[j]]$log_const_one, temp.log.const.two)
      guide.neg.ll[i,j] <- laplace_multiDimensional(temp.guide.counts, par.list[[j]]$beta_neg_pred,
        par.list[[j]]$model_cov, par.list[[j]]$pool_sds, par.list[[j]]$inverted_sigma,
        par.list[[j]]$log_const_one, temp.log.const.two)
    }
  }

  log.calc.time <- proc.time()
  calc.diff <- log.calc.time - par.est.time
  print(paste0('time to calculate loglik: ', calc.diff[3]))

  adjusted.guide.pos.ll <- guide.pos.ll
  adjusted.guide.neg.ll <- guide.neg.ll

  # browser()
  # check below when debugging: how many fall in this category?
  for(i in 1:ncol(adjusted.guide.pos.ll)){
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == -Inf),i] <- min(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] > -Inf), i])
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == Inf),i] <- max(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] < Inf), i])

    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == -Inf),i] <- min(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] > -Inf), i])
    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == Inf),i] <- max(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] < Inf), i])
  }

  total.guide.pos.ll <- rowSums(adjusted.guide.pos.ll)
  total.guide.neg.ll <- rowSums(adjusted.guide.neg.ll)

  ll.guide.ratio <- total.guide.pos.ll - total.guide.neg.ll

  par.summary <- summary_glmm_multiPoisson_pars(par.list)
  write.csv(par.summary, file = paste0(input.specs$dataName,'_glmm_allData_acrossRepl_multiPois_tmb_Summary.csv') , row.names = F)

  var.covar <- c()
  for(i in 1:length(par.list)){
    var.covar <- rbind(var.covar, par.list[[i]]$model_cov, rep(0, nrow(par.list[[i]]$model_cov)))
  }
  write.csv(var.covar, file = paste0(input.specs$dataName,'_glmm_allData_acrossRepl_multiPois_Covar_tmb_Summary.csv') , row.names = F)

  out.df <- cbind.data.frame(rawScores = ll.guide.ratio, formatScores = ll.guide.ratio,
    log2_rate_ratio = rep(1, length(ll.guide.ratio)), stringsAsFactors = FALSE)
  #add the info part
  out.df <- cbind(out.df, glmm.info, stringsAsFactors = FALSE)

  return(out.df)

}

allData_poisson_multiDim_laplace_tmb_glmm_analysis_perRepl <- function(input.counts, input.info, input.specs){

  out.score.list <- list()

  glmm.info <- input.info
  glmm.counts <- input.counts

  start.time <- proc.time()

  pos.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_allData_positiveTraining), ]
  neg.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_allData_negativeTraining), ]

  par.list <- list()

  for(i in 1:length(input.specs$repl_groups)){
    # estimate sample-specific, type specific paramters
    guide.type <- c(rep('pos', nrow(pos.counts)), rep('neg', nrow(neg.counts)))
    temp.label.counts <- rbind(pos.counts[,input.specs$repl_groups[[i]]], neg.counts[,input.specs$repl_groups[[i]]])
    par.list[[i]] <- estimate_multivar_glmm_tmb_poisson_pars(temp.label.counts, guide.type, input.specs)
  }

  par.est.time <- proc.time()

  print('allData_poisson_multiDim_laplace_tmb_glmm_analysis_perRepl')
  par.diff <- par.est.time - start.time
  print(paste0('time to estimate pars: ', par.diff[3]))

  guide.pos.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(glmm.counts)))
  guide.neg.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(glmm.counts)))

  # browser()

  for(i in 1:nrow(glmm.counts)){
    # print(paste('processed row:', i, 'out of', nrow(glmm.counts), sep = ' '))
    for(j in 1:length(input.specs$repl_groups)){
      temp.guide.counts <- unlist(glmm.counts[i,input.specs$repl_groups[[j]]])
      temp.log.const.two <- sum(- lfactorial(temp.guide.counts))
      guide.pos.ll[i,j] <- laplace_multiDimensional(temp.guide.counts, par.list[[j]]$beta_pos_pred,
        par.list[[j]]$model_cov, par.list[[j]]$pool_sds, par.list[[j]]$inverted_sigma,
        par.list[[j]]$log_const_one, temp.log.const.two)
      guide.neg.ll[i,j] <- laplace_multiDimensional(temp.guide.counts, par.list[[j]]$beta_neg_pred,
        par.list[[j]]$model_cov, par.list[[j]]$pool_sds, par.list[[j]]$inverted_sigma,
        par.list[[j]]$log_const_one, temp.log.const.two)
    }
  }

  log.calc.time <- proc.time()
  calc.diff <- log.calc.time - par.est.time
  print(paste0('time to calculate loglik: ', calc.diff[3]))

  adjusted.guide.pos.ll <- guide.pos.ll
  adjusted.guide.neg.ll <- guide.neg.ll

  # browser()
  # check below when debugging: how many fall in this category?
  for(i in 1:ncol(adjusted.guide.pos.ll)){
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == -Inf),i] <- min(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] > -Inf), i])
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == Inf),i] <- max(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] < Inf), i])

    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == -Inf),i] <- min(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] > -Inf), i])
    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == Inf),i] <- max(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] < Inf), i])
  }

  total.guide.pos.ll <- rowSums(adjusted.guide.pos.ll)
  total.guide.neg.ll <- rowSums(adjusted.guide.neg.ll)

  ll.guide.ratio <- total.guide.pos.ll - total.guide.neg.ll

  par.summary <- summary_glmm_multiPoisson_pars(par.list)
  write.csv(par.summary, file = paste0(input.specs$dataName,'_glmm_allData_perRepl_multiPois_tmb_Summary.csv') , row.names = F)

  var.covar <- c()
  for(i in 1:length(par.list)){
    var.covar <- rbind(var.covar, par.list[[i]]$model_cov, rep(0, nrow(par.list[[i]]$model_cov)))
  }
  write.csv(var.covar, file = paste0(input.specs$dataName,'_glmm_allData_perRepl_multiPois_Covar_tmb_Summary.csv') , row.names = F)

  out.df <- cbind.data.frame(rawScores = ll.guide.ratio, formatScores = ll.guide.ratio,
    log2_rate_ratio = rep(1, length(ll.guide.ratio)), stringsAsFactors = FALSE)
  #add the info part
  out.df <- cbind(out.df, glmm.info, stringsAsFactors = FALSE)

  return(out.df)

}

summary_glmm_multiPoisson_pars <- function(input.list){

  nr.betas <- 0
  pool.names <- c()
  for(i in 1:length(input.list)){
    nr.betas <- nr.betas + length(input.list[[i]]$beta_neg_pred)
    pool.names <- c(pool.names, input.list[[i]]$guide_names)
  }

  out.df <- data.frame(matrix(ncol = 6, nrow = nr.betas ))
  colnames(out.df) <- c('repl_ID', 'pool_name','orig_beta', 'neg_beta', 'pos_beta', 'beta_diff')

  temp.orig.betas <- c()
  temp.pos.betas <- c()
  temp.neg.betas <- c()
  # temp.tau <- c()
  # temp.theta <- c()
  temp.rownames <- c()
  for(i in 1:length(input.list)){
    temp.name <- paste0('repl_',i)
    temp.rownames <- c(temp.rownames, paste0(temp.name,'_beta_', c(0:(length(input.list[[i]]$beta_neg_pred) - 1) )))
    temp.orig.betas <- c(temp.orig.betas, input.list[[i]]$betas[c(1:length(input.list[[i]]$beta_neg_pred))])
    temp.neg.betas <- c(temp.neg.betas, input.list[[i]]$beta_neg_pred)
    temp.pos.betas <- c(temp.pos.betas,input.list[[i]]$beta_pos_pred)
    # temp.tau <- c(temp.tau, rep(input.list[[i]]$model_tau, length(input.list[[i]]$beta_neg_pred)))
    # temp.theta <- c(temp.theta, rep(input.list[[i]]$model_dispersion, length(input.list[[i]]$beta_neg_pred)))
  }
  beta.diff <- temp.pos.betas - temp.neg.betas

  out.df$repl_ID <- temp.rownames
  out.df$pool_name <- temp.rownames
  out.df$orig_beta <- round(temp.orig.betas,3)
  out.df$neg_beta <- round(temp.neg.betas,3)
  out.df$pos_beta <- round(temp.pos.betas,3)
  out.df$beta_diff <- round(beta.diff,3)
  # out.df$tau <- round(temp.tau,3)
  # out.df$theta <- round(temp.theta,3)
  return(out.df)
}

estimate_multivar_glmm_tmb_poisson_pars <- function(input.counts, input.type, input.specs){

  pool.IDs <- colnames(input.counts)
  guide.nrs <- c(1:nrow(input.counts))
  guide.names <- paste0('guide', guide.nrs)
  nr.pools <- ncol(input.counts)
  pool.nrs <- c(1:nr.pools)
  pool.names <- paste0('pool', pool.nrs)
  positive.types <- c(1:(nr.pools - 1))
  pos.type.label <- c('noPos', paste0('pos', positive.types))
  nr.positives <- length(which(input.type == 'pos'))

  guide.labels <- rep(guide.names, each = length(pool.names))
  pool.labels <- rep(pool.names, length(guide.names))
  type.labels <- c(rep(pos.type.label, nr.positives), rep(rep('noPos', nr.pools), nrow(input.counts) - nr.positives))

  long.df <- cbind.data.frame(observed = c(t(input.counts)), guide_labels = guide.labels,
    pool_labels = pool.labels, type_labels = type.labels, stringsAsFactors = T)

  # tmb.poisson.model <- glmmTMB(observed ~ pool_labels + type_labels + (1|guide_labels)+ (pool_labels-1 | guide_labels),
  #   data = long.df, family='poisson')

  tmb.poisson.model <- glmmTMB(observed ~ pool_labels + type_labels + (pool_labels-1 | guide_labels),
    data = long.df, family='poisson')

  # create the covariance matrix
  # temp.model.sd <- exp(getME(tmb.poisson.model, name = 'theta'))[1:(nr.pools + 1)]
  model.sd <- exp(getME(tmb.poisson.model, name = 'theta'))[1:nr.pools]
  # model.sd <- sqrt(temp.model.sd[2:(nr.pools + 1)]^2 + temp.model.sd[1]^2)  #temp.model.sd[2:(nr.pools + 1)] + temp.model.sd[1]

  model.var <- as.matrix(model.sd) %*% t(as.matrix(model.sd))
  model.corr.info <- VarCorr(tmb.poisson.model)
  model.corr <- attr(model.corr.info$cond$guide_labels, "correlation")  # not sure how to get var-covar directly..
  model.cov <- model.var * model.corr

  nr.beta.neg <- length(unique(long.df$pool_labels))
  nr.beta.pos <- nr.beta.neg - 1
  nr.total.beta <- nr.beta.neg + nr.beta.pos
  tmb.betas.orig <- tmb.poisson.model$fit$par[1:nr.total.beta]
  beta.neg <- tmb.betas.orig[1:nr.beta.neg] + c(0, rep(tmb.betas.orig[1], nr.beta.neg - 1))
  beta.pos <- beta.neg + c(0, tmb.betas.orig[(nr.beta.neg + 1):nr.total.beta])

  if('searchSilencers' %in% names(input.specs)){
    beta.pos <- beta.neg - (beta.pos - beta.neg)  # inverse the effect of the positives
  }

  tmb.neg_logLik <- logLik(tmb.poisson.model)[1]

  return(list( betas = tmb.betas.orig[1:nr.beta.neg], beta_pos_pred = beta.pos,
    beta_neg_pred = beta.neg, model_cov = model.cov, pool_vars = model.sd^2,
    neg_logLik = tmb.neg_logLik, total_model = tmb.poisson.model,
    log_const_one = log(1/sqrt(det(2 * pi * model.cov))), pool_sds = model.sd,
    inverted_sigma = solve(model.cov), pool_IDs = pool.IDs))

}

laplace_multiDimensional <- function(input.counts, input.betas, input.sigma,
  input.sds, input.inverted.sigma, input.const.one, input.const.two){

  f <- function(u, counts, betas, sigma){
    temp.lik <- -0.5 * (t(u) %*% solve(sigma) %*% u)
    for(i in 1:length(counts)){
      temp.lik <- temp.lik + counts[i] * (betas[i] + u[i]) - exp(betas[i] + u[i])
    }
    temp.lik
  }

  opt.out <- optim(par=rep(0, length(input.counts)), counts = input.counts, betas = input.betas, sigma = input.sigma,
    f, method = 'L-BFGS-B', lower = -5*input.sds, upper = 5*input.sds, control=list(fnscale=-1))

  u_0 <- opt.out$par

  inverted.sigma <- input.inverted.sigma

  hessian.matrix <- matrix(0, ncol = length(input.counts), nrow = length(input.counts))
  for(i in 1:length(input.counts)){
    for(j in 1:length(input.counts)){
      if(i != j){
        hessian.matrix[i, j] <- hessian.matrix[j, i] <- -0.5 * inverted.sigma[i, j] -0.5 * inverted.sigma[j, i]
      } else {
        hessian.matrix[i, i] <- -inverted.sigma[i, j] - exp(input.betas[i] + u_0[j])
      }
    }
  }

  log.hessian <- -0.5 * log(det(-hessian.matrix))
  log.dim.const <- log((2 * pi)^(length(input.counts)/2))

  input.const.one + input.const.two + log.dim.const + log.hessian + f(u_0, input.counts, input.betas, input.sigma)
}

poisson_multiDim_laplace_tmb_glmm_analysis <- function(input.counts, input.info, input.specs){

  out.score.list <- list()

  glmm.info <- input.info
  glmm.counts <- input.counts

  pos.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_positiveTraining), ]
  neg.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_negativeTraining), ]

  par.list <- list()

  for(i in 1:length(input.specs$multivar_repl_groups)){
    # estimate sample-specific, type specific paramters
    guide.type <- c(rep('pos', nrow(pos.counts)), rep('neg', nrow(neg.counts)))
    temp.label.counts <- rbind(pos.counts[,input.specs$multivar_repl_groups[[i]]], neg.counts[,input.specs$multivar_repl_groups[[i]]])
    par.list[[i]] <- estimate_multivar_glmm_tmb_poisson_pars(temp.label.counts, guide.type, input.specs)
  }

  guide.pos.ll <- data.frame(matrix(ncol = length(input.specs$multivar_repl_groups), nrow = nrow(glmm.counts)))
  guide.neg.ll <- data.frame(matrix(ncol = length(input.specs$multivar_repl_groups), nrow = nrow(glmm.counts)))

  # browser()

  for(i in 1:nrow(glmm.counts)){
    # print(paste('processed row:', i, 'out of', nrow(glmm.counts), sep = ' '))
    for(j in 1:length(input.specs$multivar_repl_groups)){
      temp.guide.counts <- unlist(glmm.counts[i,input.specs$multivar_repl_groups[[j]]])
      temp.log.const.two <- sum(- lfactorial(temp.guide.counts))
      guide.pos.ll[i,j] <- laplace_multiDimensional(temp.guide.counts, par.list[[j]]$beta_pos_pred,
        par.list[[j]]$model_cov, par.list[[j]]$pool_sds, par.list[[j]]$inverted_sigma,
        par.list[[j]]$log_const_one, temp.log.const.two)
      guide.neg.ll[i,j] <- laplace_multiDimensional(temp.guide.counts, par.list[[j]]$beta_neg_pred,
        par.list[[j]]$model_cov, par.list[[j]]$pool_sds, par.list[[j]]$inverted_sigma,
        par.list[[j]]$log_const_one, temp.log.const.two)
    }
  }

  adjusted.guide.pos.ll <- guide.pos.ll
  adjusted.guide.neg.ll <- guide.neg.ll

  # browser()
  # check below when debugging: how many fall in this category?
  for(i in 1:ncol(adjusted.guide.pos.ll)){
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == -Inf),i] <- min(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] > -Inf), i])
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == Inf),i] <- max(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] < Inf), i])

    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == -Inf),i] <- min(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] > -Inf), i])
    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == Inf),i] <- max(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] < Inf), i])
  }

  total.guide.pos.ll <- rowSums(adjusted.guide.pos.ll)
  total.guide.neg.ll <- rowSums(adjusted.guide.neg.ll)

  ll.guide.ratio <- total.guide.pos.ll - total.guide.neg.ll

  par.summary <- summary_glmm_multiPoisson_pars(par.list)
  write.csv(par.summary, file = paste0(input.specs$dataName,'_glmm_multiPois_tmb_Summary.csv') , row.names = F)

  var.covar <- c()
  for(i in 1:length(par.list)){
    var.covar <- rbind(var.covar, par.list[[i]]$model_cov, rep(0, nrow(par.list[[i]]$model_cov)))
  }
  write.csv(var.covar, file = paste0(input.specs$dataName,'_glmm_multiPois_Covar_tmb_Summary.csv') , row.names = F)

  out.df <- cbind.data.frame(rawScores = ll.guide.ratio, formatScores = ll.guide.ratio,
    log2_rate_ratio = rep(1, length(ll.guide.ratio)), stringsAsFactors = FALSE)
  #add the info part
  out.df <- cbind(out.df, glmm.info, stringsAsFactors = FALSE)

  return(out.df)

}

summary_glmm_poisson_pars <- function(input.list){

  nr.betas <- 0
  pool.names <- c()
  for(i in 1:length(input.list)){
    nr.betas <- nr.betas + length(input.list[[i]]$beta_neg_pred)
    pool.names <- c(pool.names, input.list[[i]]$pool_IDs)
  }

  out.df <- data.frame(matrix(ncol = 7, nrow = nr.betas ))
  colnames(out.df) <- c('repl_ID', 'pool_name','orig_beta', 'neg_beta', 'pos_beta', 'beta_diff', 'tau')

  temp.orig.betas <- c()
  temp.pos.betas <- c()
  temp.neg.betas <- c()
  temp.tau <- c()
  # temp.theta <- c()
  temp.rownames <- c()
  for(i in 1:length(input.list)){
    temp.name <- paste0('repl_',i)
    temp.rownames <- c(temp.rownames, paste0(temp.name,'_beta_', c(0:(length(input.list[[i]]$beta_neg_pred) - 1) )))
    temp.orig.betas <- c(temp.orig.betas, input.list[[i]]$betas[c(1:length(input.list[[i]]$beta_neg_pred))])
    temp.neg.betas <- c(temp.neg.betas, input.list[[i]]$beta_neg_pred)
    temp.pos.betas <- c(temp.pos.betas,input.list[[i]]$beta_pos_pred)
    temp.tau <- c(temp.tau, rep(input.list[[i]]$model_tau, length(input.list[[i]]$beta_neg_pred)))
    # temp.theta <- c(temp.theta, rep(input.list[[i]]$model_dispersion, length(input.list[[i]]$beta_neg_pred)))
  }
  beta.diff <- temp.pos.betas - temp.neg.betas

  out.df$repl_ID <- temp.rownames
  out.df$pool_name <- pool.names
  out.df$orig_beta <- round(temp.orig.betas,3)
  out.df$neg_beta <- round(temp.neg.betas,3)
  out.df$pos_beta <- round(temp.pos.betas,3)
  out.df$beta_diff <- round(beta.diff,3)
  out.df$tau <- round(temp.tau,3)
  # out.df$theta <- round(temp.theta,3)
  return(out.df)
}

estimate_univar_glmm_tmb_poisson_pars <- function(input.counts, input.type, input.specs){

  pool.IDs <- colnames(input.counts)
  guide.nrs <- c(1:nrow(input.counts))
  guide.names <- paste0('guide', guide.nrs)
  nr.pools <- ncol(input.counts)
  pool.nrs <- c(1:nr.pools)
  pool.names <- paste0('pool', pool.nrs)
  positive.types <- c(1:(nr.pools - 1))
  pos.type.label <- c('noPos', paste0('pos', positive.types))
  nr.positives <- length(which(input.type == 'pos'))

  guide.labels <- rep(guide.names, each = length(pool.names))
  pool.labels <- rep(pool.names, length(guide.names))
  type.labels <- c(rep(pos.type.label, nr.positives), rep(rep('noPos', nr.pools), nrow(input.counts) - nr.positives))

  long.df <- cbind.data.frame(observed = c(t(input.counts)), guide_labels = guide.labels,
    pool_labels = pool.labels, type_labels = type.labels, stringsAsFactors = T)

  tmb.poisson.model <- glmmTMB(observed ~ pool_labels + type_labels + (1|guide_labels),
    data = long.df, family='poisson')

  model.tau <- exp(getME(tmb.poisson.model, name = 'theta'))  # SD of the random effect
  #  no idea why this works, no idea how else to get random effect SD

  nr.beta.neg <- length(unique(long.df$pool_labels))
  nr.beta.pos <- nr.beta.neg - 1
  nr.total.beta <- nr.beta.neg + nr.beta.pos
  tmb.betas.orig <- tmb.poisson.model$fit$par[1:nr.total.beta]
  beta.neg <- tmb.betas.orig[1:nr.beta.neg] + c(0, rep(tmb.betas.orig[1], nr.beta.neg - 1))
  beta.pos <- beta.neg + c(0, tmb.betas.orig[(nr.beta.neg + 1):nr.total.beta])

  if('searchSilencers' %in% names(input.specs)){
    beta.pos <- beta.neg - (beta.pos - beta.neg)  # inverse the effect of the positives
  }

  tmb.neg_logLik <- logLik(tmb.poisson.model)[1]

  return(list(betas = tmb.betas.orig[1:nr.beta.neg], beta_pos_pred = beta.pos,
    beta_neg_pred = beta.neg, model_tau = model.tau, neg_logLik = tmb.neg_logLik,
    total_model = tmb.poisson.model, log_const_one = log(1/sqrt(2 * pi * model.tau^2)),
    pool_IDs = pool.IDs))
}

poisson_univariate_integration_function_generator <- function(input.obs, input.beta, input.tau){
  function(x){
    likelihood <- dnorm(x, mean = 0, sd = input.tau)
    for(i in 1:length(input.obs)){
      likelihood <- likelihood * dpois(input.obs[i], lambda = exp(input.beta[i] + x))
    }
    likelihood
  }
}

poisson_univariate_integration_tmb_glmm_analysis <- function(input.counts, input.info, input.specs){

  out.score.list <- list()

  glmm.info <- input.info
  glmm.counts <- input.counts

  pos.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_positiveTraining), ]
  neg.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_negativeTraining), ]

  par.list <- list()

  for(i in 1:length(input.specs$repl_groups)){
    # estimate sample-specific, type specific paramters
    guide.type <- c(rep('pos', nrow(pos.counts)), rep('neg', nrow(neg.counts)))
    temp.label.counts <- rbind(pos.counts[,input.specs$repl_groups[[i]]], neg.counts[,input.specs$repl_groups[[i]]])
    par.list[[i]] <- estimate_univar_glmm_tmb_poisson_pars(temp.label.counts, guide.type, input.specs)
  }

  guide.pos.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(glmm.counts)))
  guide.neg.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(glmm.counts)))
  for(i in 1:nrow(glmm.counts)){
    for(j in 1:length(input.specs$repl_groups)){
      temp.guide.counts <- unlist(glmm.counts[i,input.specs$repl_groups[[j]]])
      temp_pos_integration_function <- poisson_univariate_integration_function_generator(temp.guide.counts,
        par.list[[j]]$beta_pos_pred, par.list[[j]]$model_tau)
      temp_neg_integration_function <- poisson_univariate_integration_function_generator(temp.guide.counts,
        par.list[[j]]$beta_neg_pred, par.list[[j]]$model_tau)
      guide.pos.ll[i,j] <- log(integrate(temp_pos_integration_function, lower = -5*par.list[[j]]$model_tau,
        upper = 5*par.list[[j]]$model_tau)$value[1])
      guide.neg.ll[i,j] <- log(integrate(temp_neg_integration_function, lower = -5*par.list[[j]]$model_tau,
        upper = 5*par.list[[j]]$model_tau)$value[1])
    }
  }

  adjusted.guide.pos.ll <- guide.pos.ll
  adjusted.guide.neg.ll <- guide.neg.ll

  # browser()

  # check below when debugging: how many fall in this category?
  for(i in 1:ncol(adjusted.guide.pos.ll)){
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == -Inf),i] <- min(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] > -Inf), i])
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == Inf),i] <- max(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] < Inf), i])

    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == -Inf),i] <- min(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] > -Inf), i])
    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == Inf),i] <- max(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] < Inf), i])
  }

  total.guide.pos.ll <- rowSums(adjusted.guide.pos.ll)
  total.guide.neg.ll <- rowSums(adjusted.guide.neg.ll)

  ll.guide.ratio <- total.guide.pos.ll - total.guide.neg.ll

  par.summary <- summary_glmm_poisson_pars(par.list)
  write.csv(par.summary, file = paste0(input.specs$dataName,'_glmm_numIntPois_tmb_Summary.csv') , row.names = F)

  out.df <- cbind.data.frame(rawScores = ll.guide.ratio, formatScores = ll.guide.ratio,
    log2_rate_ratio = rep(1, length(ll.guide.ratio)), stringsAsFactors = FALSE)
  #add the info part
  out.df <- cbind(out.df, glmm.info, stringsAsFactors = FALSE)

  return(out.df)

}

laplace_oneDimensional <- function(input.counts, input.betas, input.tau,
  input.const.one, input.const.two){


    f_actual <- function(u, betas, tau, counts){
      temp.lik <- -(u^2)/(2*tau^2)
      for(i in 1:length(betas)){
        temp.lik <- temp.lik + counts[i]*(betas[i] + u) - exp(betas[i] + u)
      }
      temp.lik
    }

    opt.out <- optim(par=0, tau = input.tau, betas = input.betas, counts = input.counts, f_actual,
      method = 'L-BFGS-B', lower = -5 * input.tau, upper = 5 * input.tau, control=list(fnscale=-1))

    u_0 <- opt.out$par

    val = 0.5 * log(2*pi / abs(sum(-exp(input.betas + u_0)) - 1/(input.tau**2))) + f_actual(u_0, input.betas, input.tau, input.counts)

    input.const.one + input.const.two + val

}

poisson_1D_laplace_tmb_glmm_analysis <- function(input.counts, input.info, input.specs){

  out.score.list <- list()

  glmm.info <- input.info
  glmm.counts <- input.counts

  pos.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_positiveTraining), ]
  neg.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_negativeTraining), ]

  par.list <- list()

  for(i in 1:length(input.specs$repl_groups)){
    # estimate sample-specific, type specific paramters
    guide.type <- c(rep('pos', nrow(pos.counts)), rep('neg', nrow(neg.counts)))
    temp.label.counts <- rbind(pos.counts[,input.specs$repl_groups[[i]]], neg.counts[,input.specs$repl_groups[[i]]])
    par.list[[i]] <- estimate_univar_glmm_tmb_poisson_pars(temp.label.counts, guide.type, input.specs)
  }

  guide.pos.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(glmm.counts)))
  guide.neg.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(glmm.counts)))

  # browser()

  for(i in 1:nrow(glmm.counts)){
    # print(paste('processed row:', i, 'out of', nrow(glmm.counts), sep = ' '))
    for(j in 1:length(input.specs$repl_groups)){
      temp.guide.counts <- unlist(glmm.counts[i,input.specs$repl_groups[[j]]])
      temp.log.const.two <- sum(- lfactorial(temp.guide.counts))
      guide.pos.ll[i,j] <- laplace_oneDimensional(temp.guide.counts, par.list[[j]]$beta_pos_pred,
        par.list[[j]]$model_tau, par.list[[j]]$log_const_one, temp.log.const.two)
      guide.neg.ll[i,j] <- laplace_oneDimensional(temp.guide.counts, par.list[[j]]$beta_neg_pred,
        par.list[[j]]$model_tau, par.list[[j]]$log_const_one, temp.log.const.two)
    }
  }

  adjusted.guide.pos.ll <- guide.pos.ll
  adjusted.guide.neg.ll <- guide.neg.ll

  # browser()
  # check below when debugging: how many fall in this category?
  for(i in 1:ncol(adjusted.guide.pos.ll)){
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == -Inf),i] <- min(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] > -Inf), i])
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == Inf),i] <- max(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] < Inf), i])

    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == -Inf),i] <- min(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] > -Inf), i])
    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == Inf),i] <- max(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] < Inf), i])
  }

  total.guide.pos.ll <- rowSums(adjusted.guide.pos.ll)
  total.guide.neg.ll <- rowSums(adjusted.guide.neg.ll)

  ll.guide.ratio <- total.guide.pos.ll - total.guide.neg.ll

  par.summary <- summary_glmm_poisson_pars(par.list)
  write.csv(par.summary, file = paste0(input.specs$dataName,'_glmm_uniPoisson_tmb_Summary.csv') , row.names = F)

  out.df <- cbind.data.frame(rawScores = ll.guide.ratio, formatScores = ll.guide.ratio,
    log2_rate_ratio = rep(1, length(ll.guide.ratio)), stringsAsFactors = FALSE)
  #add the info part
  out.df <- cbind(out.df, glmm.info, stringsAsFactors = FALSE)

  return(out.df)

}

allData_nb_tmb_glmm_analysis <- function(input.counts, input.info, input.specs){

  out.score.list <- list()

  glmm.info <- input.info
  glmm.counts <- input.counts

  pos.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_allData_positiveTraining), ]
  neg.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_allData_negativeTraining), ]

  par.list <- list()

  for(i in 1:length(input.specs$repl_groups)){
    # estimate sample-specific, type specific paramters
    guide.type <- c(rep('pos', nrow(pos.counts)), rep('neg', nrow(neg.counts)))
    temp.label.counts <- rbind(pos.counts[,input.specs$repl_groups[[i]]], neg.counts[,input.specs$repl_groups[[i]]])
    par.list[[i]] <- estimate_glmm_tmb_nb_pars_jointly(temp.label.counts, guide.type, input.specs)
  }

  guide.pos.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(glmm.counts)))
  guide.neg.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(glmm.counts)))
  for(i in 1:nrow(glmm.counts)){
    for(j in 1:length(input.specs$repl_groups)){
      temp.guide.counts <- unlist(glmm.counts[i,input.specs$repl_groups[[j]]])
      temp_pos_integration_function <- nb_integration_function_generator(temp.guide.counts,
        par.list[[j]]$beta_pos_pred, par.list[[j]]$model_tau, par.list[[j]]$model_dispersion)
      temp_neg_integration_function <- nb_integration_function_generator(temp.guide.counts,
        par.list[[j]]$beta_neg_pred, par.list[[j]]$model_tau, par.list[[j]]$model_dispersion)
      guide.pos.ll[i,j] <- log(integrate(temp_pos_integration_function, lower = -5*par.list[[j]]$model_tau,
        upper = 5*par.list[[j]]$model_tau)$value[1])
      guide.neg.ll[i,j] <- log(integrate(temp_neg_integration_function, lower = -5*par.list[[j]]$model_tau,
        upper = 5*par.list[[j]]$model_tau)$value[1])
    }
  }

  adjusted.guide.pos.ll <- guide.pos.ll
  adjusted.guide.neg.ll <- guide.neg.ll

  # check below when debugging: how many fall in this category?
  for(i in 1:ncol(adjusted.guide.pos.ll)){
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == -Inf),i] <- min(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] > -Inf), i])
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == Inf),i] <- max(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] < Inf), i])

    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == -Inf),i] <- min(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] > -Inf), i])
    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == Inf),i] <- max(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] < Inf), i])
  }

  total.guide.pos.ll <- rowSums(adjusted.guide.pos.ll)
  total.guide.neg.ll <- rowSums(adjusted.guide.neg.ll)

  ll.guide.ratio <- total.guide.pos.ll - total.guide.neg.ll

  par.summary <- summary_glmm_nb_pars(par.list)
  write.csv(par.summary, file = paste0(input.specs$dataName,'_glmm_allData_nb_tmb_Summary.csv') , row.names = F)

  out.df <- cbind.data.frame(rawScores = 2*ll.guide.ratio, formatScores = 2*ll.guide.ratio,
    log2_rate_ratio = rep(1, length(ll.guide.ratio)), stringsAsFactors = FALSE)
  #add the info part
  out.df <- cbind(out.df, glmm.info, stringsAsFactors = FALSE)

  return(out.df)

}

nb_tmb_glmm_analysis <- function(input.counts, input.info, input.specs){

  out.score.list <- list()

  glmm.info <- input.info
  glmm.counts <- input.counts

  pos.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_positiveTraining), ]
  neg.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_negativeTraining), ]

  par.list <- list()

  for(i in 1:length(input.specs$repl_groups)){
    # estimate sample-specific, type specific paramters
    guide.type <- c(rep('pos', nrow(pos.counts)), rep('neg', nrow(neg.counts)))
    temp.label.counts <- rbind(pos.counts[,input.specs$repl_groups[[i]]], neg.counts[,input.specs$repl_groups[[i]]])
    par.list[[i]] <- estimate_glmm_tmb_nb_pars_jointly(temp.label.counts, guide.type, input.specs)
  }

  guide.pos.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(glmm.counts)))
  guide.neg.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(glmm.counts)))
  for(i in 1:nrow(glmm.counts)){
    for(j in 1:length(input.specs$repl_groups)){
      temp.guide.counts <- unlist(glmm.counts[i,input.specs$repl_groups[[j]]])
      temp_pos_integration_function <- nb_integration_function_generator(temp.guide.counts,
        par.list[[j]]$beta_pos_pred, par.list[[j]]$model_tau, par.list[[j]]$model_dispersion)
      temp_neg_integration_function <- nb_integration_function_generator(temp.guide.counts,
        par.list[[j]]$beta_neg_pred, par.list[[j]]$model_tau, par.list[[j]]$model_dispersion)
      guide.pos.ll[i,j] <- log(integrate(temp_pos_integration_function, lower = -5*par.list[[j]]$model_tau,
        upper = 5*par.list[[j]]$model_tau)$value[1])
      guide.neg.ll[i,j] <- log(integrate(temp_neg_integration_function, lower = -5*par.list[[j]]$model_tau,
        upper = 5*par.list[[j]]$model_tau)$value[1])
    }
  }

  adjusted.guide.pos.ll <- guide.pos.ll
  adjusted.guide.neg.ll <- guide.neg.ll

  # check below when debugging: how many fall in this category?
  for(i in 1:ncol(adjusted.guide.pos.ll)){
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == -Inf),i] <- min(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] > -Inf), i])
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == Inf),i] <- max(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] < Inf), i])

    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == -Inf),i] <- min(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] > -Inf), i])
    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == Inf),i] <- max(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] < Inf), i])
  }

  total.guide.pos.ll <- rowSums(adjusted.guide.pos.ll)
  total.guide.neg.ll <- rowSums(adjusted.guide.neg.ll)

  ll.guide.ratio <- total.guide.pos.ll - total.guide.neg.ll

  par.summary <- summary_glmm_nb_pars(par.list)
  write.csv(par.summary, file = paste0(input.specs$dataName,'_glmm_nb_tmb_Summary.csv') , row.names = F)

  out.df <- cbind.data.frame(rawScores = ll.guide.ratio, formatScores = ll.guide.ratio,
    log2_rate_ratio = rep(1, length(ll.guide.ratio)), stringsAsFactors = FALSE)
  #add the info part
  out.df <- cbind(out.df, glmm.info, stringsAsFactors = FALSE)

  return(out.df)

}

# this estimates the parameters for the positive and negative controls
# input.type describes if the guide is from a positive pool or not
# input.counts contain first all positive, then all negative guides!
estimate_glmm_tmb_nb_pars_jointly <- function(input.counts, input.type, input.specs){

  pool.IDs <- colnames(input.counts)
  guide.nrs <- c(1:nrow(input.counts))
  guide.names <- paste0('guide', guide.nrs)
  pool.nrs <- c(1:ncol(input.counts))
  pool.names <- paste0('pool', pool.nrs)
  positive.types <- c(1:(ncol(input.counts) - 1))
  pos.type.label <- c('noPos', paste0('pos', positive.types))
  nr.positives <- length(which(input.type == 'pos'))

  guide.labels <- rep(guide.names, each = length(pool.names))
  pool.labels <- rep(pool.names, length(guide.names))
  type.labels <- c(rep(pos.type.label, nr.positives), rep(rep('noPos', ncol(input.counts)), nrow(input.counts) - nr.positives))

  long.df <- cbind.data.frame(observed = c(t(input.counts)), guide_labels = guide.labels,
    pool_labels = pool.labels, type_labels = type.labels, stringsAsFactors = T)

  tmb.nb2.model <- glmmTMB(observed ~ pool_labels + type_labels + (1|guide_labels),
    data = long.df, family=nbinom2)

  model.dispersion <- sigma(tmb.nb2.model) # residual variance, theta corresponding to the NB variance
  model.tau <- exp(getME(tmb.nb2.model, name = 'theta'))  # SD of the random effect
  #  no idea why this works, no idea how else to get random effect SD

  nr.beta.neg <- length(unique(long.df$pool_labels))
  nr.beta.pos <- nr.beta.neg - 1
  nr.total.beta <- nr.beta.neg + nr.beta.pos
  tmb.betas.orig <- tmb.nb2.model$fit$par[1:nr.total.beta]
  beta.neg <- tmb.betas.orig[1:nr.beta.neg] + c(0, rep(tmb.betas.orig[1], nr.beta.neg - 1))
  beta.pos <- beta.neg + c(0, tmb.betas.orig[(nr.beta.neg + 1):nr.total.beta])

  if('searchSilencers' %in% names(input.specs)){
    beta.pos <- beta.neg - (beta.pos - beta.neg)  # inverse the effect of the positives
  }

  tmb.neg_logLik <- logLik(tmb.nb2.model)[1]

  return(list(betas = tmb.betas.orig[1:nr.beta.neg], beta_pos_pred = beta.pos,
    beta_neg_pred = beta.neg, model_tau = model.tau, model_dispersion = model.dispersion,
    neg_logLik = tmb.neg_logLik, pool_IDs = pool.IDs))

}

# summarize the parameter estimated obtained from the controls
summary_glmm_nb_pars <- function(input.list){

  nr.betas <- 0
  pool.names <- c()
  for(i in 1:length(input.list)){
    nr.betas <- nr.betas + length(input.list[[i]]$beta_neg_pred)
    pool.names <- c(pool.names, input.list[[i]]$pool_IDs)
  }

  out.df <- data.frame(matrix(ncol = 8, nrow = nr.betas ))
  colnames(out.df) <- c('repl_ID', 'pool_name','orig_beta', 'neg_beta', 'pos_beta', 'beta_diff', 'tau', 'theta')

  temp.orig.betas <- c()
  temp.pos.betas <- c()
  temp.neg.betas <- c()
  temp.tau <- c()
  temp.theta <- c()
  temp.rownames <- c()
  for(i in 1:length(input.list)){
    temp.name <- paste0('repl1_',i)
    temp.rownames <- c(temp.rownames, paste0(temp.name,'_beta_', c(0:(length(input.list[[i]]$beta_neg_pred) - 1) )))
    temp.orig.betas <- c(temp.orig.betas, input.list[[i]]$betas[c(1:length(input.list[[i]]$beta_neg_pred))])
    temp.neg.betas <- c(temp.neg.betas, input.list[[i]]$beta_neg_pred)
    temp.pos.betas <- c(temp.pos.betas,input.list[[i]]$beta_pos_pred)
    temp.tau <- c(temp.tau, rep(input.list[[i]]$model_tau, length(input.list[[i]]$beta_neg_pred)))
    temp.theta <- c(temp.theta, rep(input.list[[i]]$model_dispersion, length(input.list[[i]]$beta_neg_pred)))
  }
  beta.diff <- temp.pos.betas - temp.neg.betas

  out.df$repl_ID <- temp.rownames
  out.df$pool_name <- pool.names
  out.df$orig_beta <- round(temp.orig.betas,3)
  out.df$neg_beta <- round(temp.neg.betas,3)
  out.df$pos_beta <- round(temp.pos.betas,3)
  out.df$beta_diff <- round(beta.diff,3)
  out.df$tau <- round(temp.tau,3)
  out.df$theta <- round(temp.theta,3)
  return(out.df)
}

# given a set of parameters, create a function capable of computing the likelihood
#  given some observations
nb_integration_function_generator <- function(input.obs, input.beta, input.tau, input.dispersion){
  function(x){
    likelihood <- dnorm(x, mean = 0, sd = input.tau)
    for(i in 1:length(input.obs)){
      likelihood <- likelihood * dnbinom(input.obs[i], mu = exp(input.beta[i] + x), size = input.dispersion)
    }
    likelihood
  }
}

estimate_glmm_pars_byDeviance <- function(input.counts, input.type, input.specs, input.startPool){

  pool.IDs <- colnames(input.counts)
  guide.nrs <- c(1:nrow(input.counts))
  guide.names <- paste0('guide', guide.nrs)
  pool.nrs <- c(1:ncol(input.counts))
  pool.names <- paste0('pool', pool.nrs)
  positive.types <- c(1:(ncol(input.counts) - 1))

  #change to : intercept, pos1, pos
  pos.type.label <- c('interc', paste0('pos', positive.types))
  nr.positives <- length(which(input.type == 'pos'))

  neg.type.label <- c('interc', paste0('neg', positive.types))
  nr.negatives <- length(which(input.type == 'neg'))

  guide.labels <- rep(guide.names, each = length(pool.names))
  pool.labels <- rep(pool.names, length(guide.names))
  type.labels <- c(rep(pos.type.label, nr.positives), rep(neg.type.label, nr.negatives))

  long.df <- cbind.data.frame(observed = c(t(input.counts)), guide_labels = guide.labels,
    pool_labels = pool.labels, type_labels = type.labels, stringsAsFactors = T)

  tmb.nb2.model <- glmmTMB(observed ~ type_labels + (1|guide_labels),
    data = long.df, family=nbinom2)

  model.dispersion <- sigma(tmb.nb2.model) # residual variance, theta corresponding to the NB variance
  model.tau <- exp(getME(tmb.nb2.model, name = 'theta'))  # SD of the random effect
  #  no idea why this works, no idea how else to get random effect SD

  nr.beta.neg <- ncol(input.counts) - 1
  nr.beta.pos <- ncol(input.counts) - 1
  nr.total.beta <- nr.beta.neg + nr.beta.pos
  tmb.betas.orig <- tmb.nb2.model$fit$par[2:(nr.total.beta + 1)]
  beta.neg <- tmb.betas.orig[1:nr.beta.neg] #+ c(0, rep(tmb.betas.orig[1], nr.beta.neg - 1))
  beta.pos <- tmb.betas.orig[(nr.beta.neg + 1):nr.total.beta] #beta.neg + c(0, tmb.betas.orig[(nr.beta.neg + 1):nr.total.beta])

  if('searchSilencers' %in% names(input.specs)){
    beta.pos <- beta.neg - (beta.pos - beta.neg)  # inverse the effect of the positives
  }

  tmb.neg_logLik <- logLik(tmb.nb2.model)[1]

  return(list(betas = tmb.betas.orig[2:(nr.beta.neg + 1)], beta_pos_pred = beta.pos,
    beta_neg_pred = beta.neg, model_tau = model.tau, model_dispersion = model.dispersion,
    neg_logLik = tmb.neg_logLik, pool_IDs = pool.IDs, beta_intercept = tmb.nb2.model$fit$par[1]))

}

nb_deviance_integration_function_generator <- function(input.obs, input.beta, input.tau,
  input.dispersion, input.mean){
  function(x){
    likelihood <- dnorm(x, mean = input.mean, sd = input.tau)
    for(i in 1:length(input.obs)){
      likelihood <- likelihood * dnbinom(input.obs[i], mu = exp(input.beta[i] + x),
        size = input.dispersion)
    }
    likelihood
  }
}

nb_deviance_glmm_analysis <- function(input.counts, input.info, input.specs){

  out.score.list <- list()

  glmm.info <- input.info
  glmm.counts <- input.counts

  pos.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_positiveTraining), ]
  neg.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_negativeTraining), ]

  par.list <- list()

  for(i in 1:length(input.specs$repl_groups)){
    # estimate sample-specific, type specific paramters
    guide.type <- c(rep('pos', nrow(pos.counts)), rep('neg', nrow(neg.counts)))
    temp.label.counts <- rbind(pos.counts[,input.specs$repl_groups[[i]]], neg.counts[,input.specs$repl_groups[[i]]])
    # reduced.temp.label.counts <- temp.label.counts[, c(2:ncol(temp.label.counts))]
    par.list[[i]] <- estimate_glmm_pars_byDeviance(temp.label.counts, guide.type, input.specs) #, temp.label.counts[,1])
  }

  guide.pos.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(glmm.counts)))
  guide.neg.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(glmm.counts)))
  for(i in 1:nrow(glmm.counts)){
    for(j in 1:length(input.specs$repl_groups)){
      temp.guide.counts <- unlist(glmm.counts[i,input.specs$repl_groups[[j]]])
      # pool.guide.counts <- temp.guide.counts[2:length(temp.guide.counts)]
      temp_pos_integration_function <- nb_deviance_integration_function_generator(temp.guide.counts,
        c(0, par.list[[j]]$beta_pos_pred), par.list[[j]]$model_tau, par.list[[j]]$model_dispersion,
        par.list[[j]]$beta_intercept)
      temp_neg_integration_function <- nb_deviance_integration_function_generator(temp.guide.counts,
        c(0, par.list[[j]]$beta_neg_pred), par.list[[j]]$model_tau, par.list[[j]]$model_dispersion,
        par.list[[j]]$beta_intercept)
      guide.pos.ll[i,j] <- log(integrate(temp_pos_integration_function,
        lower = (-5*par.list[[j]]$model_tau + par.list[[j]]$beta_intercept),
        upper = (5*par.list[[j]]$model_tau + par.list[[j]]$beta_intercept) )$value[1])
      guide.neg.ll[i,j] <- log(integrate(temp_neg_integration_function,
        lower = (-5*par.list[[j]]$model_tau + par.list[[j]]$beta_intercept),
        upper = (5*par.list[[j]]$model_tau + par.list[[j]]$beta_intercept) )$value[1])
    }
  }

  adjusted.guide.pos.ll <- guide.pos.ll
  adjusted.guide.neg.ll <- guide.neg.ll

  # check below when debugging: how many fall in this category?
  for(i in 1:ncol(adjusted.guide.pos.ll)){
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == -Inf),i] <- min(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] > -Inf), i])
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == Inf),i] <- max(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] < Inf), i])

    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == -Inf),i] <- min(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] > -Inf), i])
    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == Inf),i] <- max(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] < Inf), i])
  }

  total.guide.pos.ll <- rowSums(adjusted.guide.pos.ll)
  total.guide.neg.ll <- rowSums(adjusted.guide.neg.ll)

  ll.guide.ratio <- total.guide.pos.ll - total.guide.neg.ll

  par.summary <- summary_glmm_nb_deviance_pars(par.list)
  write.csv(par.summary, file = paste0(input.specs$dataName,'_glmm_nb_deviance_Summary.csv') , row.names = F)

  out.df <- cbind.data.frame(rawScores = ll.guide.ratio, formatScores = ll.guide.ratio,
    log2_rate_ratio = rep(1, length(ll.guide.ratio)), stringsAsFactors = FALSE)
  #add the info part
  out.df <- cbind(out.df, glmm.info, stringsAsFactors = FALSE)

  return(out.df)

}

summary_glmm_nb_deviance_pars <- function(input.list){

  nr.betas <- 0
  pool.names <- c()
  for(i in 1:length(input.list)){
    nr.betas <- nr.betas + length(input.list[[i]]$beta_neg_pred) + 1
    pool.names <- c(pool.names, input.list[[i]]$pool_IDs)
  }

  out.df <- data.frame(matrix(ncol = 7, nrow = nr.betas ))
  colnames(out.df) <- c('repl_ID', 'pool_name', 'neg_beta', 'pos_beta', 'beta_diff', 'tau', 'theta')

  temp.pos.betas <- c()
  temp.neg.betas <- c()
  temp.tau <- c()
  temp.theta <- c()
  temp.rownames <- c()
  for(i in 1:length(input.list)){
    temp.name <- paste0('repl1_',i)
    temp.rownames <- c(temp.rownames, 'intercept', paste0(temp.name,'_beta_', c(1:length(input.list[[i]]$beta_neg_pred) )))
    temp.neg.betas <- c(temp.neg.betas, input.list[[i]]$beta_intercept, input.list[[i]]$beta_neg_pred)
    temp.pos.betas <- c(temp.pos.betas, input.list[[i]]$beta_intercept, input.list[[i]]$beta_pos_pred)
    temp.tau <- c(temp.tau, rep(input.list[[i]]$model_tau, length(input.list[[i]]$beta_neg_pred) + 1))
    temp.theta <- c(temp.theta, rep(input.list[[i]]$model_dispersion, length(input.list[[i]]$beta_neg_pred) + 1))
  }
  beta.diff <- temp.pos.betas - temp.neg.betas

  out.df$repl_ID <- temp.rownames
  out.df$pool_name <- pool.names
  out.df$neg_beta <- round(temp.neg.betas,3)
  out.df$pos_beta <- round(temp.pos.betas,3)
  out.df$beta_diff <- round(beta.diff,3)
  out.df$tau <- round(temp.tau,3)
  out.df$theta <- round(temp.theta,3)
  return(out.df)
}

RELICS_summary <- function(input.list){

  nr.betas <- 0
  pool.names <- c()
  for(i in 1:length(input.list)){
    nr.betas <- nr.betas + length(input.list[[i]]$beta_neg_pred) + 1
    pool.names <- c(pool.names, input.list[[i]]$pool_IDs)
  }

  out.df <- data.frame(matrix(ncol = 7, nrow = nr.betas ))
  colnames(out.df) <- c('repl_ID', 'pool_name', 'neg_beta', 'pos_beta', 'beta_diff', 'tau', 'theta')

  temp.pos.betas <- c()
  temp.neg.betas <- c()
  temp.tau <- c()
  temp.theta <- c()
  temp.rownames <- c()
  for(i in 1:length(input.list)){
    temp.name <- paste0('repl1_',i)
    temp.rownames <- c(temp.rownames, 'intercept', paste0(temp.name,'_beta_', c(1:length(input.list[[i]]$beta_neg_pred) )))
    temp.neg.betas <- c(temp.neg.betas, input.list[[i]]$beta_intercept, input.list[[i]]$beta_neg_pred)
    temp.pos.betas <- c(temp.pos.betas, input.list[[i]]$beta_intercept, input.list[[i]]$beta_pos_pred)
    temp.tau <- c(temp.tau, rep(input.list[[i]]$model_tau, length(input.list[[i]]$beta_neg_pred) + 1))
    temp.theta <- c(temp.theta, rep(input.list[[i]]$model_dispersion, length(input.list[[i]]$beta_neg_pred) + 1))
  }
  beta.diff <- temp.pos.betas - temp.neg.betas

  out.df$repl_ID <- temp.rownames
  out.df$pool_name <- pool.names
  out.df$neg_beta <- round(temp.neg.betas,3)
  out.df$pos_beta <- round(temp.pos.betas,3)
  out.df$beta_diff <- round(beta.diff,3)
  out.df$tau <- round(temp.tau,3)
  out.df$theta <- round(temp.theta,3)
  return(out.df)
}

RELICS_function_generator <- function(input.obs, input.beta, input.tau,
  input.dispersion, input.mean){
    # x is ranging around the mean
  function(x){
    likelihood <- dnorm(x, mean = input.mean, sd = input.tau)
    for(i in 1:length(input.obs)){
      likelihood <- likelihood * dnbinom(input.obs[i], mu = exp(input.beta[i] + x),
        size = input.dispersion)
    }
    likelihood
  }
}

# calculates fold change between two groups
# by default a pseudo count of 1 is added everywhere
# if the samples are paired ($foldChangePaired == 'yes'), fold change is calculated
#   per paired sample and then averaged across samples
# if samples are not paired the counts per groups are averaged and one fold change
#   between the average is taken
fold_change_analysis <- function(in_counts,in_info,in_specs){
  spec_contents <- names(in_specs)

  #add 1 to all counts to avoid division by zero.
  #normalize all counts by the sample counts
  pseudo.counts <- in_counts + 1
  norm.counts <- t(apply(pseudo.counts, 1, '/', colSums(pseudo.counts)))
  log10.mean.ratio <- c()
  log2.mean.ratio <- c()

  if(in_specs$foldChangePaired == 'yes'){
    count.ratios <- as.data.frame(norm.counts[,in_specs$new_grp1] / norm.counts[,in_specs$new_grp2], stringsAsFactors = FALSE)
    log10.ratios <- log10(count.ratios)
    log10.mean.ratio <- rowMeans(log10.ratios)
    log2.ratios <- log2(count.ratios)
    log2.mean.ratio <- rowMeans(log2.ratios)
  } else {
    grp1.avgCount <- rowMeans(as.data.frame(norm.counts[,in_specs$new_grp1]))
    grp2.avgCount <- rowMeans(as.data.frame(norm.counts[,in_specs$new_grp2]))
    log10.mean.ratio <- log10(grp1.avgCount / grp2.avgCount)
    log2.mean.ratio <- log2(grp1.avgCount / grp2.avgCount)
  }

  out.df <- cbind.data.frame(rawScores = log10.mean.ratio, formatScores = log10.mean.ratio,
    log2_rate_ratio = log2.mean.ratio, stringsAsFactors = FALSE)
  #add the info part
  out.df <- cbind(out.df, in_info, stringsAsFactors = FALSE)
  # if(ncol(out.df) == 7){
  #   colnames(out.df) <- c('rawScores', 'formatScores', 'log2_rate_ratio', 'chrom', 'label', 'start', 'end')
  # } else {
  #   colnames(out.df) <- c('rawScores', 'formatScores', 'log2_rate_ratio', 'chrom', 'label', 'start')
  # }
  return(out.df)
}

#analyze count data using the zero-inflated negative binaomial model
ZINB_analysis <- function(in_counts,in_info,in_specs){

  #check if any of the columns is completely zero, if so, add pseudo-count of 1 everywhere
  in_colSums <- colSums(in_counts)
  nr_zeroCols <- length(which(in_colSums == 0))
  if(nr_zeroCols > 0){
    print('Pseudo count of 1 added in ZINB analysis due to rows with all zeros')
    in_counts <- in_counts + 1
  }

  spec_contents <- names(in_specs)
  init_disp <- in_specs$InitialDispersion
  init_eta <- in_specs$InitialEta
  count_init_dispersion_est <- rep(init_disp,nrow(in_counts))
  count_init_eta_est <- rep(init_eta,nrow(in_counts))
  count_sizeFactor <- sizeFactor_TC_nrom(in_counts[, c(in_specs$new_grp1,in_specs$new_grp2)])

  ZINB_llr_list <- ZINB_rowWise_LLR(100,in_specs$ZINBminThresh,in_specs$ZINBminThresh,
    in_counts,count_sizeFactor,count_init_dispersion_est,count_init_eta_est,
    in_specs$new_grp1,in_specs$new_grp2)
  #output from above
  #[[1]] vector of log-likelihoods
  #[[2]] data frame of null eta, log-likelihoods
  #[[3]] data frame of null dispersion, rates, log-likelihoods
  #[[4]] data frame of alternative eta, log-likelihoods
  #[[5]] data frame of alternative dispersion, case rates, control rates, log-likelihoods
  ZINB_llr <- ZINB_llr_list[[1]]

  ZINB_llr_filtered <- llr_filtering_switchNeg(ZINB_llr)
  ZINB_pVals <- pchisq(ZINB_llr_filtered[[1]], df = 1, lower.tail = FALSE)

  ZINB.rate.ratio <- log2(ZINB_llr_list[[5]][,2] / ZINB_llr_list[[5]][,3])
  ZINB.formatted.scores <- -log10(ZINB_pVals)
  # ZINB.formatted.scores[ZINB.rate.ratio < 0] <- ZINB.formatted.scores[ZINB.rate.ratio < 0] * -1

  out.df <- cbind.data.frame(rawScores = ZINB_pVals, formatScores = ZINB.formatted.scores,
    log2_rate_ratio = ZINB.rate.ratio, stringsAsFactors = FALSE)
  #add the info part
  out.df <- cbind(out.df, in_info, stringsAsFactors = FALSE)
  # if(ncol(out.df) == 7){
  #   colnames(out.df) <- c('rawScores', 'formatScores', 'log2_rate_ratio', 'chrom', 'label', 'start', 'end')
  # } else {
  #   colnames(out.df) <- c('rawScores', 'formatScores', 'log2_rate_ratio', 'chrom', 'label', 'start')
  # }
  return(out.df)
}

#output: list:
#[[1]] vector of log-likelihoods
#[[2]] data frame of null eta, log-likelihoods
#[[3]] data frame of null dispersion, rates, log-likelihoods
#[[4]] data frame of alternative eta, log-likelihoods
#[[5]] data frame of alternative dispersion, case rates, control rates, log-likelihoods
ZINB_rowWise_LLR <- function(maxIter,eta_LikMin,d_rate_LikMin,in_df,in_sizeFactor,
  in_dispersion_est,in_eta_est,in_case_range,in_control_range){

  out_list <- list()

  #null hypothesis: only one rate
  initial_rates <- row_rate_est(in_df[, c(in_case_range,in_control_range)], in_sizeFactor)
	null_d_rate_1 <- ZINB_obtian_all_d_rates(in_df[, c(in_case_range,in_control_range)],
    initial_rates,in_sizeFactor,in_dispersion_est,in_eta_est)
  null_eta_1 <- ZINB_obtian_all_etas(in_df[, c(in_case_range,in_control_range)],
    initial_rates,in_sizeFactor,in_dispersion_est,in_eta_est)

  null_eta_d_rate <- ZINB_obtain_all_d_rate_eta(maxIter,in_df[, c(in_case_range,in_control_range)],
    null_d_rate_1[,2], in_sizeFactor,null_d_rate_1[,1],null_eta_1[,1],
    null_eta_1[,2],null_d_rate_1[,3],eta_LikMin,d_rate_LikMin)

  print('past null hypothesis')
  #alternative hypothesis: two rates, one dispersion
  case_df <- as.data.frame(in_df[,in_case_range], stringsAsFactors = FALSE)
  control_df <- as.data.frame(in_df[,in_control_range], stringsAsFactors = FALSE)

  case_sizeFactors <- in_sizeFactor[in_case_range]
  control_sizeFactors <- in_sizeFactor[in_control_range]

  initial_case_rate <- row_rate_est(case_df, case_sizeFactors)
  initial_control_rate <- row_rate_est(control_df, control_sizeFactors)

  alternative_rate_d_1 <- ZINB_obtian_all_HA_d_rates(case_df,control_df,initial_case_rate,initial_control_rate,case_sizeFactors,
    control_sizeFactors,in_dispersion_est,in_eta_est)
  alternative_eta_1 <- ZINB_obtian_all_HA_etas(case_df,control_df,initial_case_rate,initial_control_rate,case_sizeFactors,
    control_sizeFactors,in_dispersion_est,in_eta_est)

  alternative_eta_rate_d <- ZINB_obtain_all_HA_d_rate_eta(100,case_df,control_df,alternative_rate_d_1[,2],alternative_rate_d_1[,3],case_sizeFactors,control_sizeFactors,
    alternative_rate_d_1[,1],alternative_eta_1[,1],alternative_eta_1[,2],alternative_rate_d_1[,4],eta_LikMin,d_rate_LikMin)

  lik_diff_h0_hA <- null_eta_d_rate[[2]][,3] - alternative_eta_rate_d[[2]][,4]
  llr_h0_hA <- 2*lik_diff_h0_hA

  out_list[[1]] <- llr_h0_hA
  out_list[[2]] <- null_eta_d_rate[[1]]
  out_list[[3]] <- null_eta_d_rate[[2]]
  out_list[[4]] <- alternative_eta_rate_d[[1]]
  out_list[[5]] <- alternative_eta_rate_d[[2]]

  return(out_list)
}

#estimation of dispersion and rate parameter of one sgRNA given the counts of one sgRNA
ZINB_estimate_d_rate <- function(in_row_counts,in_tc_sj,par,in_eta){
  d <- par[1]
  rate <- par[2]
  S_j <- in_tc_sj

	in_row_counts <- unlist(in_row_counts)
  scaled_rates <- rate*S_j

  zero_scaled_rates <- scaled_rates[which(in_row_counts == 0)]
  zero_eta <- in_eta[which(in_row_counts == 0)]
  non_zero_scaled_rates <- scaled_rates[which(in_row_counts > 0)]
  non_zero_eta <- in_eta[which(in_row_counts > 0)]
  non_zero_obs <- in_row_counts[which(in_row_counts > 0)]
  nr_zeros <- rep(0,(length(in_row_counts)-length(non_zero_obs)))

  zero_lik <- log(zero_eta + (1-zero_eta)*dnbinom(nr_zeros, size = d, mu = zero_scaled_rates))
  non_zero_lik <- log(1-non_zero_eta) + dnbinom(non_zero_obs, size = d, mu = non_zero_scaled_rates, log = TRUE)

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

#obtain rates and dispersion parameters given input counts and etas
ZINB_obtian_all_d_rates <- function(in_df,in_RateVec,in_tcSj,in_d_vec,in_eta_vec){
    out_d <- c()
    out_rate <- c()
    out_values <- c()

    for(i in 1:(length(in_RateVec))){
        optim_row_vec <- in_df[i,]
        temp_optim_par_out <- optim(par=c(in_d_vec[i],in_RateVec[i]),ZINB_estimate_d_rate,
        in_row_counts = optim_row_vec,in_eta = in_eta_vec,
        in_tc_sj = in_tcSj,method = 'L-BFGS-B',lower = c(0.0001,0.0001))
        out_d <- c(out_d,temp_optim_par_out$par[1])
        out_rate <- c(out_rate,temp_optim_par_out$par[2])
        out_values <- c(out_values,temp_optim_par_out$value)
    }
    out_d_rate_values <- cbind(out_d,out_rate,out_values)
    return(out_d_rate_values)
}

#obtain all etas given fixed rates and dispersion parameters
ZINB_obtian_all_etas <- function(in_df,in_RateVec,in_tcSj,in_D_vec,in_etas){
    out_eta <- c()
    out_values <- c()

    for(i in 1:(dim(in_df)[2])){
        optim_col_vec <- unlist(in_df[,i])
        temp_optim_eta_out <- optim(par=in_etas[i],ZINB_estimate_eta,in_col_vec = optim_col_vec,
          in_rate_vec = in_RateVec,in_d_vec = in_D_vec,in_tc_sj = in_tcSj[i],method = 'L-BFGS-B',lower = 0.0001, upper = 0.9999)
        out_eta <- c(out_eta,temp_optim_eta_out$par[1])
        out_values <- c(out_values,temp_optim_eta_out$value)
    }
    out_eta_values <- cbind(out_eta,out_values)
    return(out_eta_values)
}

#estimation of etas given counts and dispersions
ZINB_estimate_eta <- function(in_col_vec,in_rate_vec,in_d_vec,in_tc_sj,par){
  eta <- par[1]

  zero_scaled_rates <- in_rate_vec[which(in_col_vec == 0)]*in_tc_sj
  zero_d <- in_d_vec[which(in_col_vec == 0)]
  non_zero_scaled_rates <- in_rate_vec[which(in_col_vec > 0)]*in_tc_sj
  non_zero_d <- in_d_vec[which(in_col_vec > 0)]
  non_zero_obs <- in_col_vec[which(in_col_vec > 0)]
  nr_zeros <- rep(0,(length(in_col_vec)-length(non_zero_obs)))

  zero_lik <- log(eta + (1-eta)*dnbinom(nr_zeros, size = zero_d, mu = zero_scaled_rates))
  non_zero_lik <- log((1-eta)) + dnbinom(non_zero_obs, size = non_zero_d, mu = non_zero_scaled_rates,log = TRUE)

  lik_vec <- c(zero_lik,non_zero_lik)

  if(length(which(lik_vec == Inf)) > 0){
    print(paste('there are: ',length(which(lik_vec == Inf)),' Inf in lik_vec',sep = ''))
  }
  if(length(which(lik_vec == -Inf)) > 0){
    print(paste('there are: ',length(which(lik_vec == -Inf)),' -Inf in lik_vec',sep = ''))
  }

  result = -sum(lik_vec)
  result
}

#output: list:
#[[1]] results for eta values from latest iteration: eta values | log-likelihood
#[[2]] results for dispersion and rate values from latest iteration: dispersion | rates | log-likelihood values
ZINB_obtain_all_d_rate_eta <- function(maxIter,in_df,in_RateVec,in_tcSj,in_D_vec,in_etas,inEtaLik,inDRateLik,eta_lik_min,rate_d_min){
  out_list <- list()
  maxIterCout <- 1
  cur_eta <- in_etas
  cur_rate <- in_RateVec
  cur_d <- in_D_vec
  cur_eta_lik <- inEtaLik
  cur_d_rate_lik <- inDRateLik
  while(maxIterCout < maxIter){
    print(paste('iteration: ',maxIterCout,sep = ''))
    temp_rate_d_value <- ZINB_obtian_all_d_rates(in_df,cur_rate,in_tcSj,cur_d,cur_eta)
    cur_rate <- temp_rate_d_value[,2]
    cur_d <- temp_rate_d_value[,1]

    temp_eta_value <- ZINB_obtian_all_etas(in_df,cur_rate,in_tcSj,cur_d,cur_eta)
    cur_eta <- temp_eta_value[,1]

    eta_lik_diff <- cur_eta_lik - temp_eta_value[,2]
    rate_d_lik_diff <- cur_d_rate_lik - temp_rate_d_value[,3]

    if(sum(abs(eta_lik_diff)) < eta_lik_min && sum(abs(rate_d_lik_diff)) < rate_d_min){
      print(paste('number of incorrect eta likelihoods: ', length(which(eta_lik_diff < 0)), sep = ''))
      print(paste('number of incorrect rate_d likelihoods: ', length(which(rate_d_lik_diff < 0)), sep = ''))
      out_list[[1]] <- temp_eta_value
      out_list[[2]] <- temp_rate_d_value
      return(out_list)
    }

    cur_eta_lik <- temp_eta_value[,2]
    cur_d_rate_lik <- temp_rate_d_value[,3]
    maxIterCout <- maxIterCout + 1
  }
  print('Did not converge!')
  out_list[[1]] <- temp_eta_value
  out_list[[2]] <- temp_rate_d_value
  return(out_list)
}

#for alternative hypothesis
#output: estimates for dispersion | estimates for case rates | estimates for control rates | log-likelihood values
ZINB_obtian_all_HA_d_rates <- function(in_case_df,in_control_df,in_caseRate,in_controlRate,in_case_tcSj,in_control_tcSj,in_d_vec,in_Etas){
    out_d <- c()
    out_case_rate <- c()
    out_control_rate <- c()
    out_values <- c()

    for(i in 1:(length(in_caseRate))){
        optim_case_row <- in_case_df[i,]
        optim_control_row <- in_control_df[i,]
        temp_optim_par_out <- optim(par=c(in_d_vec[i],in_caseRate[i],in_controlRate[i]),ZINB_estimate_HA_d_rate,in_case_counts = optim_case_row, in_etas = in_Etas,
        in_control_counts = optim_control_row, in_case_tc_sj = in_case_tcSj,in_control_tc_sj = in_control_tcSj,method = 'L-BFGS-B',lower = c(0.0001,0.0001))
        out_d <- c(out_d,temp_optim_par_out$par[1])
        out_case_rate <- c(out_case_rate,temp_optim_par_out$par[2])
        out_control_rate <- c(out_control_rate,temp_optim_par_out$par[3])
        out_values <- c(out_values,temp_optim_par_out$value)
    }
    out_d_rate_values <- cbind(out_d,out_case_rate,out_control_rate,out_values)
    return(out_d_rate_values)
}

#for alternative hypothesis
ZINB_estimate_HA_d_rate <- function(in_case_counts,in_control_counts,in_case_tc_sj,in_control_tc_sj,par,in_etas){
  d <- par[1]
  case_rate <- par[2]
  control_rate <- par[3]

	in_case_row_counts <- unlist(in_case_counts)
  in_control_row_counts <- unlist(in_control_counts)
  scaled_case_rates <- case_rate*in_case_tc_sj
  scaled_control_rates <- control_rate*in_control_tc_sj
  case_etas <- in_etas[1:length(in_case_counts)]
  control_etas <- in_etas[(length(in_case_counts)+1):length(in_etas)]

  zero_scaled_case_rate <- scaled_case_rates[which(in_case_row_counts == 0)]
  zero_case_eta <- case_etas[which(in_case_row_counts == 0)]
  non_zero_case_rate <- scaled_case_rates[which(in_case_row_counts > 0)]
  non_zero_case_eta <- case_etas[which(in_case_row_counts > 0)]
  non_zero_case_obs <- in_case_row_counts[which(in_case_row_counts > 0)]
  nr_case_zeros <- rep(0,(length(in_case_row_counts)-length(non_zero_case_obs)))

  zero_scaled_control_rate <- scaled_control_rates[which(in_control_row_counts == 0)]
  zero_control_eta <- control_etas[which(in_control_row_counts == 0)]
  non_zero_control_rate <- scaled_control_rates[which(in_control_row_counts > 0)]
  non_zero_control_eta <- control_etas[which(in_control_row_counts > 0)]
  non_zero_control_obs <- in_control_row_counts[which(in_control_row_counts > 0)]
  nr_control_zeros <- rep(0,(length(in_control_row_counts)-length(non_zero_control_obs)))

  case_zero_lik <- log(zero_case_eta + (1-zero_case_eta)*dnbinom(nr_case_zeros, size = d, mu = zero_scaled_case_rate))
  case_non_zero_lik <-log(1-non_zero_case_eta) + dnbinom(non_zero_case_obs, size = d, mu = non_zero_case_rate, log = TRUE)

  control_zero_lik <- log(zero_control_eta + (1-zero_control_eta)*dnbinom(nr_control_zeros, size = d, mu = zero_scaled_control_rate))
  control_non_zero_lik <- log(1-non_zero_control_eta) + dnbinom(non_zero_control_obs, size = d, mu = non_zero_control_rate, log = TRUE)

  case_lik_vec <- c(case_zero_lik,case_non_zero_lik)
  control_lik_vec <- c(control_zero_lik,control_non_zero_lik)

  if(length(which(case_lik_vec == Inf)) > 0){
    print(paste('there are: ',length(which(case_lik_vec == Inf)),' Inf in case_lik_vec',par,sep = ''))
    print(paste('par/mu: ',par,sep = ''))
  }
  if(length(which(case_lik_vec == -Inf)) > 0){
    print(paste('there are: ',length(which(case_lik_vec == -Inf)),' -Inf in case_lik_vec',par,sep = ''))
    print(paste('par/mu: ',par,sep = ''))
  }
  if(length(which(control_lik_vec == Inf)) > 0){
    print(paste('there are: ',length(which(control_lik_vec == Inf)),' Inf in control_lik_vec',par,sep = ''))
    print(paste('par/mu: ',par,sep = ''))
  }
  if(length(which(control_lik_vec == -Inf)) > 0){
    print(paste('there are: ',length(which(control_lik_vec == -Inf)),' -Inf in control_lik_vec',par,sep = ''))
    print(paste('par/mu: ',par,sep = ''))
  }

  case_result <- -sum(case_lik_vec)
  control_result <- -sum(control_lik_vec)
  result = case_result + control_result
  result
}

#output: eta values | log-likelihood
ZINB_obtian_all_HA_etas <- function(in_case_df,in_control_df,in_caseRate,in_controlRate,in_case_tcSj,in_control_tcSj,in_D_vec,in_etas){
    out_eta <- c()
    out_values <- c()

    case_etas <- in_etas[1:(dim(in_case_df)[2])]
    control_etas <- in_etas[(dim(in_case_df)[2] + 1):length(in_etas)]

    for(i in 1:(dim(in_case_df)[2])){
        optim_col_vec <- unlist(in_case_df[,i])
        temp_optim_eta_out <- optim(par=case_etas[i],ZINB_estimate_eta,in_col_vec = optim_col_vec,
          in_rate_vec = in_caseRate,in_d_vec = in_D_vec,in_tc_sj = in_case_tcSj[i],method = 'L-BFGS-B',lower = 0.0001, upper = 0.9999)
        out_eta <- c(out_eta,temp_optim_eta_out$par[1])
        out_values <- c(out_values,temp_optim_eta_out$value)
    }
    for(i in 1:(dim(in_control_df)[2])){
        optim_col_vec <- unlist(in_control_df[,i])
        temp_optim_eta_out <- optim(par=control_etas[i],ZINB_estimate_eta,in_col_vec = optim_col_vec,
          in_rate_vec = in_controlRate,in_d_vec = in_D_vec,in_tc_sj = in_control_tcSj[i],method = 'L-BFGS-B',lower = 0.0001, upper = 0.9999)
        out_eta <- c(out_eta,temp_optim_eta_out$par[1])
        out_values <- c(out_values,temp_optim_eta_out$value)
    }
    out_eta_values <- cbind(out_eta,out_values)
    return(out_eta_values)
}

#output: list:
#[[1]] results for eta values from latest iteration: eta values | log-likelihood
#[[2]] results for dispersion and rate values from latest iteration: dispersion | case rates | control rates | log-likelihood values
ZINB_obtain_all_HA_d_rate_eta <- function(maxIter,in_case_df,in_control_df,in_caseRate,in_controlRate,in_case_tcSj,in_control_tcSj,
  in_D_vec,in_etas,inEtaLik,inDRateLik,eta_lik_min,rate_d_min){

  out_list <- list()
  maxIterCout <- 1
  cur_eta <- in_etas
  cur_case_rate <- in_caseRate
  cur_control_rate <- in_controlRate
  cur_d <- in_D_vec
  cur_eta_lik <- inEtaLik
  cur_d_rate_lik <- inDRateLik
  while(maxIterCout < maxIter){
    print(paste('iteration: ',maxIterCout,sep = ''))
    temp_rate_d_value <- ZINB_obtian_all_HA_d_rates(in_case_df,in_control_df,cur_case_rate,cur_control_rate,in_case_tcSj,in_control_tcSj,cur_d,cur_eta)
    cur_case_rate <- temp_rate_d_value[,2]
    cur_control_rate <- temp_rate_d_value[,3]
    cur_d <- temp_rate_d_value[,1]

    temp_eta_value <- ZINB_obtian_all_HA_etas(in_case_df,in_control_df,cur_case_rate,cur_control_rate,in_case_tcSj,in_control_tcSj,cur_d,cur_eta)
    cur_eta <- temp_eta_value[,1]

    eta_lik_diff <- cur_eta_lik - temp_eta_value[,2]
    rate_d_lik_diff <- cur_d_rate_lik - temp_rate_d_value[,4]

    if(sum(abs(eta_lik_diff)) < eta_lik_min && sum(abs(rate_d_lik_diff)) < rate_d_min){
      print(paste('number of incorrect eta likelihoods: ', length(which(eta_lik_diff < 0)), sep = ''))
      print(paste('number of incorrect rate_d likelihoods: ', length(which(rate_d_lik_diff < 0)), sep = ''))
      out_list[[1]] <- temp_eta_value
      out_list[[2]] <- temp_rate_d_value
      return(out_list)
    }

    cur_eta_lik <- temp_eta_value[,2]
    cur_d_rate_lik <- temp_rate_d_value[,4]
    maxIterCout <- maxIterCout + 1
  }
  print('Did not converge!')
  out_list[[1]] <- temp_eta_value
  out_list[[2]] <- temp_rate_d_value
  return(out_list)
}

# used following link to set up data:
#  http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#countmat
DESeq2_analysis <- function(input.counts, input.info, input.specs){

  selected.samples <- input.counts[,c(input.specs$new_grp1, input.specs$new_grp2)]

  select.rowSums <- rowSums(selected.samples)
  nr.zeroRows <- length(which(select.rowSums == 0))
  if(nr.zeroRows > 0){
    print('Warning: DESeq2 method detected all-zero rows. Pseudo count of 1 is added.')
    selected.samples <- selected.samples + 1
  }

  # generate an indexing system to ensure the results can be matched between counts and info
  input.info$index <- paste('row', c(1:nrow(input.info)), sep = '_')
  formatted.counts <- selected.samples
  formattted.colnames <- colnames(formatted.counts)
  row.names(formatted.counts) <- input.info$index

  # create the condition table
  group1.names <- rep('group1', length(input.specs$new_grp1)) #specify the low-pool guides as group2 (even though labelled as group1)
  group2.names <- rep('group2', length(input.specs$new_grp2)) #specify the low-pool guides as group1 (even though labelled as group2)
  # need to add initial pools first, hence group 2 and then group 1
  group.names <- c(group2.names,group1.names)  #edgeR assumes treatment-pool comes first, else ratio will be off
  condition.df <- data.frame(condition = group.names, type = rep(input.specs$DESeq2_type, length(group.names)))
  row.names(condition.df) <- formattted.colnames

  deseq2.data <- DESeqDataSetFromMatrix(countData = formatted.counts,
                              colData = condition.df,
                              design = ~ condition)
  deseq2.analyzed <- DESeq(deseq2.data)
  deseq2.results <- results(deseq2.analyzed)


  output.df <- cbind.data.frame(rawScores =  deseq2.results$pvalue,
    formatScores = -log10(deseq2.results$pvalue) * sign(deseq2.results$log2FoldChange),
    log2_rate_ratio = deseq2.results$log2FoldChange , stringsAsFactors = FALSE)

  # output.df <- cbind.data.frame(rawScores =  deseq2.results$pvalue,
  #   formatScores = -log10(deseq2.results$pvalue),
  #   log2_rate_ratio = deseq2.results$log2FoldChange , stringsAsFactors = FALSE)

  # check if row names need to be adjusted ...

  output.df<- cbind(output.df, input.info, stringsAsFactors = FALSE)

  return(output.df)
}

# first part of CREST analysis
# calculates p-values for the counts of each sample
edgeR_analysis <- function(input.counts, input.info, input.specs){

  edgeR.list <- run_edgeR(input.counts,input.info,input.specs)
  edgeR.ranked.table <- edgeR.list$edgeRtable
  edgeR.filtered.info <- edgeR.list$filteredInfo

  # convert edgeR scores into RRA format
  output.df <- cbind.data.frame(rawScores =  edgeR.ranked.table$PValue, formatScores = edgeR.ranked.table$score,
    log2_rate_ratio = edgeR.ranked.table$logFC, stringsAsFactors = FALSE)

  output.df<- cbind(output.df, edgeR.filtered.info, stringsAsFactors = FALSE)

  # this should happen at another place ...
  # if(ncol(edgeR.filtered.info) == 4){
  #   output.df$end <- edgeR.filtered.info[,4]
  # }

  save_edgeRData(edgeR.list$filteredCounts, edgeR.filtered.info, input.specs)

  return(output.df)
}

#given counts, determine which sgRNAs are differentially observed
# the groups are specified by 'input.specs$new_grp1' and  'input.specs$new_grp2'
#output:
#   logFC: edgeR calculated log fold change
#   logCPM: edgeR calculated log counts per million
#   LR: edgeR calculated log likelihood ratio?
#   PValue: edgeR calculated p-values
#   score: manually added score, -log2 signed p-values based on logFC
#   rank: ranked scores
  #logFC   logCPM           LR     PValue       score rank
  #1 -0.01277404 7.121507 0.0002225714 0.98809694  0.01197447 2937
  #2  1.76918903 7.361801 3.5958577326 0.05792373 -2.84862813  514
run_edgeR <- function(input.counts, input.info, input.specs){
  selected.count.cols <- input.counts[,c(input.specs$new_grp2, input.specs$new_grp1)]
  group1.names <- rep('group1', length(input.specs$new_grp1)) #specify the low-pool guides as group2 (even though labelled as group1)
  group2.names <- rep('group2', length(input.specs$new_grp2)) #specify the low-pool guides as group1 (even though labelled as group2)
  group.names <- c(group1.names, group2.names)  #edgeR assumes treatment-pool comes first, else ratio will be off
  colnames(selected.count.cols) <- group.names
  filtered.matrix <- c()
  filtered.info <- c()
  count.matrix <- as.matrix(selected.count.cols)

  # crestEdgeRfilter is the criterion used by CREST-analysis to remove low-count reads
  if('crestEdgeRfilter' %in% names(input.specs)){
    filtered.rows <- which(rowSums(cpm(selected.count.cols) > 3) >= (ncol(selected.count.cols)/3))  # nor sure what latter part is doing
    filtered.matrix <- count.matrix[filtered.rows,]
    filtered.info <- input.info[filtered.rows,]
  } else {
    filtered.matrix <- count.matrix
    filtered.info <- input.info
  }

  save_edgeR_analysisData(filtered.matrix, filtered.info, input.specs)

  DGE.output <- DGEList(filtered.matrix, group = factor(group.names))  #create DGE list object
  edgeR.design <- model.matrix(~factor(group.names))  # specify there are only 2 groups
  edgeR.dispersion <- estimateDisp(DGE.output, edgeR.design)
  edgeR.fit <- glmFit(edgeR.dispersion, edgeR.design)
  edgeR.lrt <- glmLRT(edgeR.fit, coef=2)
  edgeR.table <- edgeR.lrt$table
  edgeR.table$score <- -log10(edgeR.table$PValue) * sign(edgeR.table$logFC)
  edgeR.score.rank <- rank(edgeR.table$score)  #is this ordering done the right way?
  edgeR.table$rank <- edgeR.score.rank
  return(list(edgeRtable = edgeR.table, filteredInfo = filtered.info, filteredCounts = as.data.frame(filtered.matrix)))
}

#analyze count data using the negative binomial model
# input: count df, info df, specification list
# output: df, each method adds the following to the
#   rawScores, formatScores, log2_rate_ratio, chromosome, labels, sgRNA targets (can be 1 or 2 columns)
#   for more description see: score_calculation
NB_analysis <- function(in_counts,in_info,in_specs){

  #check if any of the columns is completely zero, if so, add pseudo-count of 1 everywhere
  in_rowSums <- rowSums(in_counts)
  nr_zeroRows <- length(which(in_rowSums == 0))
  if(nr_zeroRows > 0){
    print('Warning: NB method detected all-zero rows. Pseudo count of 1 is added.')
    in_counts <- in_counts + 1
  }

  spec_contents <- names(in_specs)
  init_disp <- in_specs$InitialDispersion
  count_init_dispersion_est <- rep(init_disp,nrow(in_counts))

  count_sizeFactor <- sizeFactor_TC_nrom(in_counts[, c(in_specs$subGroup1, in_specs$subGroup2)])

  #NB_llr_list[[1]]: vector of log-likelihoods
  #NB_llr_list[[2]]: data frame of null dispersion, rates, log-likelihoods
  #NB_llr_list[[3]]: data frame of alternative dispersion, case rates, control rates, log-likelihoods
  NB_llr_list <- negativeBinomial_rowWise_LLR(in_counts,count_sizeFactor,count_init_dispersion_est,in_specs$subGroup1,in_specs$subGroup2)
  #extract log-likelihoods
  NB_llr <- NB_llr_list[[1]]
  #remove negative log-likelihoods
  NB_llr_filtered <- llr_filtering_switchNeg(NB_llr)
  #NB p-Values, only using one degree of freedom
  NB_count_pVals <- pchisq(NB_llr_filtered[[1]], df = 1, lower.tail = FALSE)
  NB_rate_ratio <- log2(NB_llr_list[[3]][,2] / NB_llr_list[[3]][,3])
  NB.formatted.scores <- -log10(NB_count_pVals)
  #NB.formatted.scores[NB_rate_ratio < 0] <- NB.formatted.scores[NB_rate_ratio < 0] * -1

  # save negative binomial results
  if('saveNB_LLR' %in% names(in_specs)){
    write.csv(NB_llr_filtered[[1]], 'NB_LLR.csv', row.names = F)
    write.csv(NB_llr_list[[2]], 'NB_nullModel.csv', row.names = F)
    write.csv(NB_llr_list[[3]], 'NB_altModel.csv', row.names = F)
  }

  out.df <- cbind.data.frame(rawScores = NB_count_pVals, formatScores = NB.formatted.scores,
    log2_rate_ratio = NB_rate_ratio, stringsAsFactors = FALSE)
  #add the info part
  out.df <- cbind(out.df, in_info, stringsAsFactors = FALSE)
  # if(ncol(out.df) == 7){
  #   colnames(out.df) <- c('rawScores', 'formatScores', 'log2_rate_ratio', 'chrom', 'label', 'start', 'end')
  # } else {
  #   colnames(out.df) <- c('rawScores', 'formatScores', 'log2_rate_ratio', 'chrom', 'label', 'start')
  # }
  return(out.df)
}

#changes the rows which have a negative log-likelihood to the smallest, positive log-likelihood
#this is based on the assumption that negative log-likelihoods arise from precision errors during numerical optimization
#output: [[1]]: adjusted log-likelihood ratios (for the ones below zero)
#        [[2]]: positions of the negative log-likelihoods for potential follow up
llr_filtering_switchNeg <- function(in_llr){
  out_list <- list()

  neg_llr <- which(in_llr < 0)
  print(paste('Nr. of negative llr: ',length(neg_llr),sep = ''))
  print(paste('Smallest negative llr: ',min(in_llr),sep = ''))
  if(length(neg_llr) > 0){
    temp_neg_llr <- in_llr[-neg_llr]
    min_positiv_llr <- min(temp_neg_llr)
    print(paste('Smallest positive llr: ',min_positiv_llr,sep = ''))
    in_llr[neg_llr] <- min_positiv_llr
  }

  out_list[[1]] <- in_llr
  out_list[[2]] <- neg_llr
  return(out_list)
}

#output: list:
#[[1]] vector of log-likelihoods
#[[2]] data frame of null dispersion, rates, log-likelihoods
#[[3]] data frame of alternative dispersion, case rates, control rates, log-likelihoods
negativeBinomial_rowWise_LLR <- function(in_df,in_sizeFactor,in_dispersion_est,in_case_range,in_control_range){

  out_list = list()

  #null hypothesis: only one rate

  initial_rates <- row_rate_est(in_df[, c(in_case_range, in_control_range)], in_sizeFactor)
	null_d_rate <- NB_obtian_all_d_rates(in_df[, c(in_case_range, in_control_range)],
    initial_rates,in_sizeFactor,in_dispersion_est)
  print('past null hypothesis')
  #alternative hypothesis: two rates, one dispersion
  case_df <- as.data.frame(in_df[,in_case_range], stringsAsFactors = FALSE)
  control_df <- as.data.frame(in_df[,in_control_range], stringsAsFactors = FALSE)

  case_sizeFactors <- in_sizeFactor[in_case_range]
  control_sizeFactors <- in_sizeFactor[in_control_range]

  initial_case_rate <- row_rate_est(case_df, case_sizeFactors)
  initial_control_rate <- row_rate_est(control_df, control_sizeFactors)

  alternative_rates <- NB_obtian_all_HA_d_rates(case_df,control_df,initial_case_rate,initial_control_rate,
    case_sizeFactors,control_sizeFactors,in_dispersion_est)

  lik_diff_h0_hA <- null_d_rate[,3] - alternative_rates[,4]
  llr_h0_hA <- 2*lik_diff_h0_hA

  out_list[[1]] <- llr_h0_hA
  out_list[[2]] <- null_d_rate
  out_list[[3]] <- alternative_rates
  return(out_list)
}

#negative-binomial null-model. Assumes one rate for each guide
#obtain dispersion- and rate-values for each guide
NB_obtian_all_d_rates <- function(in_df,in_RateVec,in_tcSj,in_d_vec){
    out_d <- c()
    out_rate <- c()
    out_values <- c()

    for(i in 1:(length(in_RateVec))){
        optim_row_vec <- in_df[i,]
        temp_optim_par_out <- optim(par=c(in_d_vec[i],in_RateVec[i]),NB_estimate_d_rate,in_row_counts = optim_row_vec,in_tc_sj = in_tcSj,method = 'L-BFGS-B',lower = c(0.0001,0.0001))
        out_d <- c(out_d,temp_optim_par_out$par[1])
        out_rate <- c(out_rate,temp_optim_par_out$par[2])
        out_values <- c(out_values,temp_optim_par_out$value)
    }
    out_d_rate_values <- cbind(out_d,out_rate,out_values)
    return(out_d_rate_values)
}

#maximize likelihood for one rate parameter and one dispersion parameter for a given guide (for a given set of counts)
#scale rate by size factor in order to keep raw counts
NB_estimate_d_rate <- function(in_row_counts,in_tc_sj,par){
  d <- par[1]
  rate <- par[2]
  S_j <- in_tc_sj

	in_row_counts <- unlist(in_row_counts)
  scaled_rates <- rate*S_j

  lik_vec <- dnbinom(in_row_counts, size = d, mu = scaled_rates, log = TRUE)

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

#negative-binomial alternative-model. Assumes one rate for each condition and one dispersion parameter
#obtain one dispersion- and two rate-values for each guide
NB_obtian_all_HA_d_rates <- function(in_case_df,in_control_df,in_caseRate,in_controlRate,in_case_tcSj,in_control_tcSj,in_d_vec){
    out_d <- c()
    out_case_rate <- c()
    out_control_rate <- c()
    out_values <- c()

    for(i in 1:(length(in_caseRate))){
        optim_case_row <- in_case_df[i,]
        optim_control_row <- in_control_df[i,]
        temp_optim_par_out <- optim(par=c(in_d_vec[i],in_caseRate[i],in_controlRate[i]),NB_estimate_HA_d_rate,in_case_counts = optim_case_row,in_control_counts = optim_control_row,
            in_case_tc_sj = in_case_tcSj,in_control_tc_sj = in_control_tcSj,method = 'L-BFGS-B',lower = c(0.0001,0.0001))
        out_d <- c(out_d,temp_optim_par_out$par[1])
        out_case_rate <- c(out_case_rate,temp_optim_par_out$par[2])
        out_control_rate <- c(out_control_rate,temp_optim_par_out$par[3])
        out_values <- c(out_values,temp_optim_par_out$value)
    }
    out_d_rate_values <- cbind(out_d,out_case_rate,out_control_rate,out_values)
    return(out_d_rate_values)
}

#maximize likelihood for two rate parameters and one dispersion parameter for a given guide (for a given set of counts)
#scale rates by size factor in order to keep raw counts
NB_estimate_HA_d_rate <- function(in_case_counts,in_control_counts,in_case_tc_sj,in_control_tc_sj,par){
  d <- par[1]
  case_rate <- par[2]
  control_rate <- par[3]

	in_case_row_counts <- unlist(in_case_counts)
  in_control_row_counts <- unlist(in_control_counts)
  scaled_case_rates <- case_rate*in_case_tc_sj
  scaled_control_rates <- control_rate*in_control_tc_sj

  case_lik_vec <- dnbinom(in_case_row_counts, size = d, mu = scaled_case_rates, log = TRUE)
  control_lik_vec <- dnbinom(in_control_row_counts, size = d, mu = scaled_control_rates, log = TRUE)

  if(length(which(case_lik_vec == Inf)) > 0){
    print(paste('there are: ',length(which(case_lik_vec == Inf)),' Inf in case_lik_vec',par,sep = ''))
    print(paste('par/mu: ',par,sep = ''))
  }
  if(length(which(case_lik_vec == -Inf)) > 0){
    print(paste('there are: ',length(which(case_lik_vec == -Inf)),' -Inf in case_lik_vec',par,sep = ''))
    print(paste('par/mu: ',par,sep = ''))
  }
  if(length(which(control_lik_vec == Inf)) > 0){
    print(paste('there are: ',length(which(control_lik_vec == Inf)),' Inf in control_lik_vec',par,sep = ''))
    print(paste('par/mu: ',par,sep = ''))
  }
  if(length(which(control_lik_vec == -Inf)) > 0){
    print(paste('there are: ',length(which(control_lik_vec == -Inf)),' -Inf in control_lik_vec',par,sep = ''))
    print(paste('par/mu: ',par,sep = ''))
  }

  case_result <- -sum(case_lik_vec)
  control_result <- -sum(control_lik_vec)
  result = case_result + control_result
  result
}

#estimate per-guide rate-values based on initial count data
row_rate_est <- function(in_df, in_sj){
    m_sj <-length(in_sj)
    qij <- rowSums(data.frame(mapply(`/`,in_df,in_sj,SIMPLIFY = FALSE)))/(m_sj)
    return(qij)
}

#obtain size factors via total count normalization
sizeFactor_TC_nrom <- function(in_df){
  col_counts <- colSums(in_df)
  col_sizeFactors <- col_counts/sum(col_counts)
  return(col_sizeFactors)
}

# saving data to be analyzed
# input: count df, info df, specification list
# output: csv tables is specified
#   counts to be analyzed
#   info to be analyzed
#   rows removed in the filtering process
save_analysisData <- function(input.counts, input.info, input.specs){
  if('saveAnalysisCounts' %in% names(input.specs)){
    write.csv(input.counts, file = 'filtered_analysis_counts.csv',row.names = FALSE)
  }
  if('saveAnalysisInfo' %in% names(input.specs)){
    write.csv(input.info, file = 'filtered_analysis_info.csv',row.names = FALSE)
  }
  if('saveAnalysisRowsRemoved' %in% names(input.specs)){
    write.csv(input.specs, file = 'filtered_analysis_rowsRemoved.csv',row.names = FALSE)
  }
}

# saving data to be analyzed
# input: count df, info df, specification list
# output: csv tables is specified
#   counts to be analyzed
#   info to be analyzed
#   rows removed in the filtering process
save_edgeR_analysisData <- function(input.counts, input.info, input.specs){
  if('saveAnalysisCounts' %in% names(input.specs)){
    write.csv(input.counts, file = 'edgeR_filtered_analysis_counts.csv',row.names = FALSE)
  }
  if('saveAnalysisInfo' %in% names(input.specs)){
    write.csv(input.info, file = 'edgeR_filtered_analysis_info.csv',row.names = FALSE)
  }
  if('saveAnalysisRowsRemoved' %in% names(input.specs)){
    write.csv(input.specs, file = 'edgeR_filtered_analysis_rowsRemoved.csv',row.names = FALSE)
  }
}

# saving data to be analyzed
# edgeR has a separate filtering step
# input: count df, info df, specification list
# output: csv tables is specified
#   counts to be analyzed
#   info to be analyzed
#   rows removed in the filtering process
save_edgeRData <- function(input.counts, input.info, input.specs){
  if('saveAnalysisCounts' %in% names(input.specs)){
    write.csv(input.counts, file = 'edgeR_filtered_analysis_counts.csv',row.names = FALSE)
  }
  if('saveAnalysisInfo' %in% names(input.specs)){
    write.csv(input.info, file = 'edgeR_filtered_analysis_info.csv',row.names = FALSE)
  }
  if('saveAnalysisRowsRemoved' %in% names(input.specs)){
    write.csv(input.specs, file = 'edgeR_filtered_analysis_rowsRemoved.csv',row.names = FALSE)
  }
}

# filter counts
# input: count df, info df, specification list
# output: list
#   $counts, contains re-arranged counts
#   $info, contains chromsome and positional information about re-arranged data
#   $specs, contains specs about re-arranged data
#   $rowsRemoved, rows which were removed
# TODO: implement zero rows, beware of unsotred pool
filtering_counts <- function(input.counts, input.info, input.specs){
  out.list <- list()

  #if no filtering is applied
  if(input.specs$dataFilter == 'none'){
    out.list$counts <- input.counts
    out.list$info <- input.info
    out.list$specs <- input.specs
    return(out.list)
  }

  #if zero rows should be removed
  if(input.specs$dataFilter == 'zeroRows'){
    rows.to.remove <- which(rowSums(input.counts) == 0)
    if(length(rows.to.remove) > 0){
      out.list$counts <- input.counts[-rows.to.remove,]
      out.list$info <- input.info[-rows.to.remove,]
      out.list$specs <- input.specs
    } else {
      out.list$counts <- input.counts
      out.list$info <- input.info
      out.list$specs <- input.specs
    }
    return(out.list)
  }

  #if CPMs below indicated shoul dbe removed
  if(input.specs$dataFilter == 'CPMs'){
    cpm.cutoff <- input.specs$cpm
    rows.to.remove <- which(rowSums(cpm(input.counts)) < cpm.cutoff)
    if(length(rows.to.remove) > 0){
      print(paste('Removed ',length(rows.to.remove),' with CPM less than ', cpm.cutoff, sep = ''))
      out.list$counts <- input.counts[-rows.to.remove,]
      out.list$info <- input.info[-rows.to.remove,]
      out.list$specs <- input.specs
    } else {
      out.list$counts <- input.counts
      out.list$info <- input.info
      out.list$specs <- input.specs
    }
    return(out.list)
  }

}

# re-arrange data
# input: count df, info df, specification list
# output: list
#   $counts, contains re-arranged counts
#   $info, contains chromsome and positional information about re-arranged data
#   $specs, contains specs about re-arranged data
# TODO: implement arrangement, use per-chromosome information, beware of unsotred pool
reArrange_data <- function(input.counts, input.info, input.specs){
  output.list <- list()
  if(input.specs$dataArranging == 'none'){
    print('No data rearrangements were made.')
    output.list$counts <- input.counts
    output.list$info <- input.info
    output.list$specs <- input.specs
    return(output.list)
  }
  if(input.specs$dataArranging == 'yes'){
    if('preBinning' %in% names(input.specs)){
      binning.list <- per_sample_binning(input.counts, input.info, input.specs, 'preBinning')
      output.list$counts <- binning.list$counts
      output.list$info <- binning.list$info
      output.list$specs <- input.specs
      return(output.list)
    }
    if('preOverlapBinning' %in% names(input.specs)){
      binning.list <- per_sample_binning(input.counts, input.info, input.specs, 'preOverlapBinning')
      output.list$counts <- binning.list$counts
      output.list$info <- binning.list$info
      output.list$specs <- input.specs
      return(output.list)
    }
  }
}

# binning per chromosome and aggregate information within bin
# input:
#   input.counts: count file
#   input.info: info file
#   input.specs: specification list
#   input.overlap.type: determines whether overlapping or non-overlapping bins are created
# output: list
#   $counts
#   $info
per_sample_binning <- function(input.counts, input.info, input.specs, input.overlap.type){

  all.chrom.rows <- which(input.info$chrom != 'NA')
  all.chrom.counts <- input.counts[all.chrom.rows, ]
  all.chrom.info <- input.info[all.chrom.rows, ]

  if(! 'end' %in% colnames(all.chrom.info)){
    all.chrom.info$end <- all.chrom.info$start
  }

  final.chrom.bin.counts <- c()  # labels associated with each bin
  final.chrom.info <- c()  # name of each bin, used for matching with labels
  final.chrom.bin.means <- c()  # mean number of guides per bin

  all.chroms <- unique(all.chrom.info$chrom)
  for(chrom in all.chroms){
    chrom.counts <- all.chrom.counts[all.chrom.info$chrom == chrom,]
    chrom.info <- all.chrom.info[all.chrom.info$chrom == chrom,]
    range.list <- obtain_data_GRanges(chrom.info, input.specs, rep(0, nrow(chrom.info)))
    chrom.ranges <- range.list$gRanges
    bin.ranges <- c()
    if(input.overlap.type == 'preBinning'){
      bin.range.list <- generate_preBinning_bins(range.list$rangeStartEnd, chrom, input.specs)
      bin.ranges <- bin.range.list$binGRanges
    }
    else if(input.overlap.type == 'preOverlapBinning'){
      bin.range.list <- generate_preBinning_overlBins(range.list$rangeStartEnd, chrom, input.specs)
      bin.ranges <- bin.range.list$binGRanges
    }
    chrom.labels <- obtain_bin_labels(bin.ranges, chrom.info, input.specs, chrom)
    chrom.overlaps <- as.data.frame(findOverlaps(bin.ranges,chrom.ranges,type = 'any'))

    #for each $queryHits, take the rows corresponding to the $subjectHits (x[,2]) and round the mean for each sample in the count df
    chrom.bin.counts <- c()
    if(input.specs$binBy == 'mean'){
      chrom.bin.counts <- t(sapply(split(chrom.overlaps, chrom.overlaps$queryHits), function(x) round(colMeans(chrom.counts[x[,2],]))))
    }
    else if(input.specs$binBy == 'sum'){
      chrom.bin.counts <- t(sapply(split(chrom.overlaps, chrom.overlaps$queryHits), function(x) colSums(chrom.counts[x[,2],])))
    }
    chrom.info <- cbind.data.frame(chrom = rep(chrom, length(unique(chrom.overlaps$queryHits))),
      label = chrom.labels[unique(chrom.overlaps$queryHits)],
      start = bin.range.list$start[unique(chrom.overlaps$queryHits)],
      end = bin.range.list$end[unique(chrom.overlaps$queryHits)], stringsAsFactors = FALSE)
    # not sure about query hits, have feeling might be adding too many labels?

    final.chrom.bin.counts <- rbind(final.chrom.bin.counts, chrom.bin.counts)  # labels associated with each bin
    final.chrom.info <- rbind(final.chrom.info, chrom.info)

    #obtain average number of overlaps per bin (used to obtain average per na-bin)
    chrom.bin.mean <- mean(table(chrom.overlaps$queryHits))
    final.chrom.bin.means <- c(final.chrom.bin.means, chrom.bin.mean)
  }

  if(length(all.chrom.rows) < nrow(input.counts)){
    na.counts <- input.counts[input.info$chrom == 'NA',]
    na.info <- input.info[input.info$chrom == 'NA',]
    na.chrom.bins.overl.list <- na_pre_binning(final.chrom.bin.means, na.counts, na.info, input.specs)
    final.chrom.bin.counts <- rbind(final.chrom.bin.counts, na.chrom.bins.overl.list$naCounts)
    final.chrom.info <- rbind(final.chrom.info, na.chrom.bins.overl.list$naInfo)
  }

  # remove row names as it can give an error in the fold change analysis
  rownames(final.chrom.bin.counts) <- NULL

  return(list(counts = as.data.frame(final.chrom.bin.counts), info = final.chrom.info))
}

# binning of counts of non-targeting guides
# input:
#   input.mean.bin.guides: mean number of guides per bin per chromosome
#   input.na.counts: all counts of non-targeting guides
#   input.na.info: infor about non-targeting guides
#   input.specs: specification list
# output:
#   naCounts: counts after binning
#   naInfo: info after binning
na_pre_binning <- function(input.mean.bin.guides, input.na.counts, input.na.info, input.specs){
  output.na.counts <- c()
  output.na.info <- c()

  bin.mean <- round(mean(input.mean.bin.guides))
  #treat each label within na separately
  all.na.lables <- unique(input.na.info$label)
  for(lab in all.na.lables){
    lab.na.info <- input.na.info[input.na.info$label == lab,]
    lab.na.counts <- input.na.counts[input.na.info$label == lab,]
    # set the range within which NA counts are pooled (ex: all counts from $start[1] to $end[1] are pooled)
    lab.na.bins <- cbind.data.frame(start = seq(1,nrow(lab.na.info),bin.mean),
      end = seq(1,nrow(lab.na.info),bin.mean))
    lab.na.bins$end <- lab.na.bins$end + bin.mean - 1
    lab.na.bins$end[nrow(lab.na.bins)] <- min(lab.na.bins$end[nrow(lab.na.bins)], nrow(lab.na.info))

    # for each group of NA counts (c(x[1]:x[2])), obtain the per-sample mean
    lab.bin.na.counts <- c()
    if(input.specs$binBy == 'mean'){
      lab.bin.na.counts <- t(sapply(apply(lab.na.bins, 1, function(x) c(x[1]:x[2])), function(x) round(colMeans(lab.na.counts[x, ]))))
    }
    else if(input.specs$binBy == 'sum'){
      lab.bin.na.counts <- t(sapply(apply(lab.na.bins, 1, function(x) c(x[1]:x[2])), function(x) colSums(lab.na.counts[x, ])))
    }

    lab.bin.na.info <- cbind.data.frame(chrom = rep('NA',nrow(lab.bin.na.counts)),
      label = rep(lab,nrow(lab.bin.na.counts)), start = rep(NA,nrow(lab.bin.na.counts)),
      end = rep(NA,nrow(lab.bin.na.counts)), stringsAsFactors = FALSE)

    output.na.counts <- rbind(output.na.counts, lab.bin.na.counts)
    output.na.info <- rbind(output.na.info, lab.bin.na.info)
  }
  return(list(naCounts = output.na.counts, naInfo =  output.na.info))
}

# generate a set of non-overlapping bins and also keep track of vectorized start and end positions
# input:
#   input.range, input.chrom, input.specs
# output: list
#   binGRanges: genomicRanges object for the bin ranges
#   start: vectorized start positions of bins
#   end: vectorized end positions of bins
generate_preBinning_bins <- function(input.range, input.chrom, input.specs){
  temp.range.bin.start <- seq(input.range[1], input.range[2], by = input.specs$preBinning)
  temp.range.bin.end <- temp.range.bin.start + input.specs$preBinning

  temp.bin.ranges <- GRanges(seqnames = rep(input.chrom,length(temp.range.bin.start)),
    ranges = IRanges(temp.range.bin.start, temp.range.bin.end))
  temp.bin.ranges$names <- paste(input.chrom, paste(temp.range.bin.start, temp.range.bin.end,
    sep = '_'), sep = '_')
  return(list(binGRanges = temp.bin.ranges, start = temp.range.bin.start, end = temp.range.bin.end))
}

# generate a set of overlapping bins, step size defined by binSteps, and also keep track of vectorized start and end positions
# input:
#   input.range, input.chrom, input.specs
# output: list
#   binGRanges: genomicRanges object for the bin ranges
#   start: vectorized start positions of bins
#   end: vectorized end positions of bins
generate_preBinning_overlBins <- function(input.range, input.chrom, input.specs){
  temp.range.bin.start <- seq(input.range[1], input.range[2], by = input.specs$preOverlapBinning[1])
  temp.range.bin.end <- temp.range.bin.start + input.specs$preOverlapBinning[2]

  temp.bin.ranges <- GRanges(seqnames = rep(input.chrom,length(temp.range.bin.start)),
    ranges = IRanges(temp.range.bin.start, temp.range.bin.end))
  temp.bin.ranges$names <- paste(input.chrom, paste(temp.range.bin.start, temp.range.bin.end,
    sep = '_'), sep = '_')
  return(list(binGRanges = temp.bin.ranges, start = temp.range.bin.start, end = temp.range.bin.end))
}

# extract the samples of interest
# input: count df, info df, specification list
# output: list
#   $counts, contains counts of relevant samples
#   $info, contains chromsome and positional information
#   $specs, contains specs and added $subGroup1, $subGroup2 and subUnsorted (latter if Model_X is used for analysis)
extract_samples <- function(input.counts, input.info, input.specs){

  columns.to.extract <- c()
  # columns.to.extract <- c(input.specs$Group1, input.specs$Group2)
  if('Group1Comb' %in% names(input.specs)){
    input.specs$Group1 <- vector(mode = 'numeric', length = length(input.specs$Group1Comb))
    input.specs$Group2 <- vector(mode = 'numeric', length = length(input.specs$Group2Comb))
    for(i in 1:length(input.specs$Group1Comb)){
      if(grepl('+', input.specs$Group1Comb[i], fixed = T)){
        temp.split <- strsplit(input.specs$Group1Comb[i],'+')[[1]]
        extract.cols <- seq(1,length(temp.split),2)
        temp.cols <- as.numeric(temp.split[extract.cols])
        col.sum <- rowSums(input.counts[, temp.cols])
        added.col.index <- ncol(input.counts) + 1
        input.specs$Group1[i] <- added.col.index
        input.counts[[paste0(temp.cols, collapse = '_')]] <- col.sum
      } else {
        input.specs$Group1[i] <- as.numeric(input.specs$Group1Comb[i])
      }
    }
    for(i in 1:length(input.specs$Group2Comb)){
      if(grepl('+', input.specs$Group2Comb[i], fixed = T)){
        temp.split <- strsplit(input.specs$Group2Comb[i],'+')[[1]]
        extract.cols <- seq(1,length(temp.split),2)
        temp.cols <- as.numeric(temp.split[extract.cols])
        col.sum <- rowSums(input.counts[, temp.cols])
        added.col.index <- ncol(input.counts) + 1
        input.specs$Group2[i] <- added.col.index
        input.counts[[paste0(temp.cols, collapse = '_')]] <- col.sum
      } else {
        input.specs$Group2[i] <- as.numeric(input.specs$Group2Comb[i])
      }
    }
  }

  columns.to.extract <- c(input.specs$Group1, input.specs$Group2)

  # if 'repl_groups' are specified, add them to the columns which have to be extracted
  if('repl_groups' %in% names(input.specs)){
    columns.to.extract <- unique(c(columns.to.extract, as.numeric(unlist(input.specs$repl_groups))))
  }

  if('multivar_repl_groups' %in% names(input.specs)){
    columns.to.extract <- unique(c(columns.to.extract, as.numeric(unlist(input.specs$multivar_repl_groups))))
  }

  # find the mapping to the above...

  # if('Model_X' %in% input.specs$Method) {
  #   if(! 'unsortedGroups' %in% names(input.specs)) stop('Error in "extract_samples": using Model_X without specifying unsorted group(s)')
  #   columns.to.extract <- c(columns.to.extract, input.specs$unsortedGroups)
  # }
  if(max(columns.to.extract) > ncol(input.counts)) stop('Error in "extract_samples": gropus are referencing non-existent columns')
  output.list <- list()
  output.list$counts <- input.counts[,columns.to.extract]
  final.info <- input.info
  if(! all(c('chrom','label','start') %in% colnames(final.info) == TRUE)){
    print("Missing column labels; either 'chrom', 'label' or 'start' ")
  }
  # colnames(final.info)[1:4] <- c('chrom','label','start','end')
  output.list$info <- final.info
  updated.specs <- input.specs

  group1.pos <- which(columns.to.extract %in% input.specs$Group1)
  group2.pos <- which(columns.to.extract %in% input.specs$Group2)

  # adjust the referenced columns from the replicate gropus
  if('repl_groups' %in% names(input.specs)) {
    # need to preserve the ordering to have meaningful intercepts in GLMM
    repl.grp.pos <- lapply(input.specs$repl_groups, function(x){
      out.order <- c()
      for(i in 1:length(x)){
        out.order <- c(out.order, which(columns.to.extract == x[i]))
      }
      out.order
      })
    updated.specs$repl_groups <- repl.grp.pos
  }

  if('multivar_repl_groups' %in% names(input.specs)) {
    # need to preserve the ordering to have meaningful intercepts in GLMM
    mv.repl.grp.pos <- lapply(input.specs$multivar_repl_groups, function(x){
      mv.out.order <- c()
      for(i in 1:length(x)){
        mv.out.order <- c(mv.out.order, which(columns.to.extract == x[i]))
      }
      mv.out.order
      })
    updated.specs$multivar_repl_groups <- mv.repl.grp.pos
  }


  updated.specs$subGroup1 <- group1.pos # c(1:length(input.specs$Group1))
  updated.specs$subGroup2 <- group2.pos # c(1:length(input.specs$Group2)) + length(input.specs$Group1)
  # if('Model_X' %in% input.specs$Method){
  #   updated.specs$subUnsorted <- c(1:length(input.specs$unsortedGroups)) + length(input.specs$Group1) + length(input.specs$Group2)
  # }
  output.list$specs <- updated.specs
  return(output.list)
}

# no longer used
#crest-seq analysis
#workflow was followed as done by Diao, Fang et al., 2017
crest_seq_analysis <- function(input.counts, input.info, input.specs){
  spec_contents <- names(input.specs)

  #rank1 sgRNAs (done using edgeR)
  #what does output look like?
  edgeR.list <- run_edgeR(input.counts,input.info,input.specs)
  edgeR.ranked.table <- edgeR.list$edgeRtable
  edgeR.filtered.info <- edgeR.list$filteredInfo

  # convert edgeR scores into RRA format
  crest.rra.input.df <- cbind.data.frame(rawScores =  edgeR.ranked.table$PValue, formatScores = edgeR.ranked.table$score,
    log2_rate_ratio = edgeR.ranked.table$logFC, chrom = edgeR.filtered.info[,1], label = edgeR.filtered.info[,2],
    start = edgeR.filtered.info[,3], stringsAsFactors = FALSE)
  if(ncol(edgeR.filtered.info) == 4){
    crest.rra.input.df$end <- edgeR.filtered.info[,4]
  }

  rra.out.df <- RRA(crest.rra.input.df, input.specs, 'CRESTanalysis', input.specs$CRESTrra, 'CRESTrra')

  return(rra.out.df)

}

#output: list. Each analysis specification is stored within for donwstream instructions
read_analysis_specs <- function(in_specs_loc, data.dir = NULL){
  raw_specs <- scan(in_specs_loc,what='character')
  out_specs_list <- list()

  data.loc <- ''  #directory where the data is located (used in case data file is only referenced by name and not entrie path)
  if(! is.null(data.dir)){
    data.loc <- data.dir
  }
  #added during analysis:
  # analysis.specs$labelFile

  out_specs_list$InitialDispersion <- 0.5
  out_specs_list$RRAmaxPercent <- 0.1
  out_specs_list$CRISPReffectRange <- 20
  # default initial dispersion and zero-mass proportion estimate:
  out_specs_list$InitialDispersion <- 0.5
  out_specs_list$InitialEta <- 0.5
  out_specs_list$ZINBminThresh <- 0.1
  # default RRA properties:
  out_specs_list$rraPermutNr <- 5
  out_specs_list$rraNaBin <- 10
  out_specs_list$scoreThresh <- 0.3  # threshold used for calclating alpha-RRA, scores below this won't be considered

  # DESeq2 default is for paired end analysis
  out_specs_list$DESeq2_type <- 'paired-end'
  # out_specs_list$dataDir <- getwd()
  # out_specs_list$SaveCountFile <- 'no'
  # out_specs_list$SaveInfoFile <- 'no'
  # out_specs_list$SaveFilteringRemoved <- 'no'
  # out_specs_list$SavePvalues <- 'yes'
  # #default for seed when simulating data is set to on
  # out_specs_list$includeSimSeed <- 'yes'
  # out_specs_list$analysisType <- 'data' #or 'simulation'
  # out_specs_list$PerformDataAnalysis <- 'yes'
  # out_specs_list$PerformSubTiling <- 'no'
  # #default plot colors:
  # out_specs_list$labelNames <- c('chr','pos','neg')
  # out_specs_list$colorLabels <- c(1,2,3)
  # out_specs_list$colorTypes <- c('black','red','green')
  # out_specs_list$pValThresh <- c(0.1,0.05,0.01,0.001)
  # out_specs_list$pValCol <- c('blue','green','orange','red')
  # out_specs_list$positiveLabels <- 2
  # out_specs_list$negativeLabels <- 3
  # #default RRA parameters:
  # out_specs_list$binSize <- 50
  # out_specs_list$RRAgoodsgRNA <- 3
  # out_specs_list$RRAminPeakWidth <- 1000
  #
  # #assuming no genes are to be plotted and chromosomes are to be plotted individually:
  # out_specs_list$plotGenes <- 'no'
  # out_specs_list$PlotSeparateChromsomes <- 'no'
  #
  # #default parameters for model_X:
  # out_specs_list$nrBetCtrl <- 10
  # out_specs_list$alphaMLKStart <- 0.5
  # out_specs_list$betalMLKStart <- 2
  # out_specs_list$probGfunc <- 0.05
  #
  # #default MAGeCK parameters:
  #   # column to use for p-values
  # out_specs_list$mageckTest <- 12
  #   # number of sgRNAs that need to be greater then mean
  # out_specs_list$mageckNrDiff <- 2
  #   # value which max(ctrl_meanVar_cand) (in function:mageckCheckVar) has to be greater than
  # out_specs_list$mageckNumDiff <- 5

  out_specs_list$dataFilter <- 'none'
  out_specs_list$saveGuideScores <- 'yes'
  out_specs_list$createBedGraph <- 'yes'
  out_specs_list$postScoringAnalysis <- 'yes'
  out_specs_list$RELICS_genomeScoring <- 'yes'
  out_specs_list$dataArranging <- 'none'

  for(spec in raw_specs){
    spec_id <- strsplit(spec,':')[[1]][1]
    if(spec_id == 'dataArranging'){
      out_specs_list$pcrDuplicates <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'CountFileLoc'){
      out_specs_list$CountFileLoc <- paste(data.loc, paste(strsplit(spec,':')[[1]][2:length(strsplit(spec,':')[[1]])],collapse = ':'), sep = '')
      print(paste('countFile: ', out_specs_list$CountFileLoc, sep = ''))
    }
    if(spec_id == 'sgRNAInfoFileLoc'){
      out_specs_list$sgRNAInfoFileLoc <- paste(data.loc, paste(strsplit(spec,':')[[1]][2:length(strsplit(spec,':')[[1]])],collapse = ':'), sep = '')
    }
    if(spec_id == 'sgTargetPosFileLoc'){
      out_specs_list$sgTargetPosFileLoc <- paste(data.loc, paste(strsplit(spec,':')[[1]][2:length(strsplit(spec,':')[[1]])],collapse = ':'), sep = '')
    }
    if(spec_id == 'Method'){
      out_specs_list$Method <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'foldChangePaired'){
      out_specs_list$foldChangePaired <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'Group1'){
      out_specs_list$Group1 <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'Group2'){
      out_specs_list$Group2 <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'Group1Comb'){
      out_specs_list$Group1Comb <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'Group2Comb'){
      out_specs_list$Group2Comb <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'unsortedGroups'){
      out_specs_list$unsortedGroups <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'dataArranging'){
      out_specs_list$dataArranging <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'dataFilter'){
      out_specs_list$dataFilter <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'CPMs'){
      out_specs_list$dataFilter <- 'CPMs'
      out_specs_list$cpm <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'dataName'){
      out_specs_list$dataName <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'MAGeCKrra'){
      out_specs_list$MAGeCKrra <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'CRESTrra'){
      out_specs_list$CRESTrra <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'repl_groups'){
      out_specs_list$repl_groups <- lapply(strsplit(strsplit(strsplit(spec,':')[[1]][2],';')[[1]],','), as.numeric)
      names(out_specs_list$repl_groups) <- paste0('repl_',c(1:length(out_specs_list$repl_groups)))
    }
    if(spec_id == 'multivar_repl_groups'){
      out_specs_list$multivar_repl_groups <- lapply(strsplit(strsplit(strsplit(spec,':')[[1]][2],';')[[1]],','), as.numeric)
      names(out_specs_list$multivar_repl_groups) <- paste0('repl_',c(1:length(out_specs_list$multivar_repl_groups)))
    }
    if(spec_id == 'crisprSystem'){
      out_specs_list$crisprSystem <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'crisprEffectRange'){
      out_specs_list$crisprEffectRange <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'deletionSize'){
      out_specs_list$deletionSize <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'RELICS_genomeScoring'){
      out_specs_list$RELICS_genomeScoring <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'noPostScoreAnalysis'){
      out_specs_list$noPostScoreAnalysis <- strsplit(spec,':')[[1]][2]
    }

    #parameters with default settings:
    if(spec_id == 'saveAnalysisCounts'){
      out_specs_list$saveAnalysisCounts <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'saveAnalysisInfo'){
      out_specs_list$saveAnalysisInfo <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'saveAnalysisRowsRemoved'){
      out_specs_list$saveAnalysisRowsRemoved <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'InitialDispersion'){
      out_specs_list$InitialDispersion <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'InitialEta'){
      out_specs_list$InitialEta <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'ZINBminThresh'){
      out_specs_list$ZINBminThresh <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'binSize'){
      out_specs_list$binSize <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'preBinning'){
      out_specs_list$preBinning <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'preOverlapBinning'){
      out_specs_list$preOverlapBinning <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'binBy'){
      out_specs_list$binBy <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'labelFile'){
      out_specs_list$labelFile <- paste(data.loc, strsplit(spec,':')[[1]][2], sep = '')
    }
    if(spec_id == 'labelHierarchy'){
      out_specs_list$labelHierarchy <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'postScoringAnalysis'){
      out_specs_list$postScoringAnalysis <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'postScoringCREST'){
      out_specs_list$postScoringCREST <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'postScoreRRA'){
      out_specs_list$postScoreRRA <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'postScoreAlphaRRA'){
      out_specs_list$postScoreAlphaRRA <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'postScoreSlidingWindow'){
      out_specs_list$postScoreSlidingWindow <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'guidePerSlidingWindow'){
      out_specs_list$guidePerSlidingWindow <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'maxWindowSize'){
      out_specs_list$maxWindowSize <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'postScoringMAGeCK'){
      out_specs_list$postScoringMAGeCK <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'evaluate_perGuide_Performance'){
      out_specs_list$evaluate_perGuide_Performance <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'evaluate_perRegion_Performance'){
      out_specs_list$evaluate_perRegion_Performance <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'evaluate_perElement_Performance'){
      out_specs_list$evaluate_perElement_Performance <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'pos_regions'){
      out_specs_list$pos_regions <- paste(data.loc, paste(strsplit(spec,':')[[1]][2:length(strsplit(spec,':')[[1]])],collapse = ':'), sep = '')
    }
    if(spec_id == 'pos_region_size'){
      out_specs_list$pos_region_size <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'saveGuideScores'){
      out_specs_list$saveGuideScores <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'save_dirichlet_info'){
      out_specs_list$save_dirichlet_info <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'save_glmmUntrained_info'){
      out_specs_list$save_glmmUntrained_info <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'positiveLabels'){
      out_specs_list$positiveLabels <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'negativeLabels'){
      out_specs_list$negativeLabels <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'genePosFileLoc'){
      out_specs_list$plotGenes <- 'yes'
      plotGenesDataLoc <- paste(data.loc, paste(strsplit(spec,':')[[1]][2:length(strsplit(spec,':')[[1]])],collapse = ':'), sep = '')
      out_specs_list$plotGenesData <- read.csv(plotGenesDataLoc,as.is = TRUE,header=TRUE)
    }
    if(spec_id == 'zoomRange'){
      out_specs_list$zoomRange <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
      out_specs_list$zoomChromosomePlot <- 'yes'
    }
    if(spec_id == 'tileZoomRange'){
      out_specs_list$tileZoomRange <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
      out_specs_list$tileZoomChromosomePlot <- 'yes'
    }
    if(spec_id == 'plotAllChrom'){
      out_specs_list$plotAllChrom <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'plotChroms'){
      out_specs_list$plotChroms <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'plotSeparateChroms'){
      out_specs_list$plotSeparateChroms <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'RRAmaxPercent'){
      out_specs_list$RRAmaxPercent <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'CRISPReffectRange'){
      out_specs_list$CRISPReffectRange <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'rraPermutNr'){
      out_specs_list$rraPermutNr <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'rraNaBin'){
      out_specs_list$rraNaBin <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'scoreThresh'){
      out_specs_list$scoreThresh <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'createBedGraph'){
      out_specs_list$createBedGraph <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'saveNB_LLR'){
      out_specs_list$saveNB_LLR <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'postScoreGuideWindowRRA'){
      out_specs_list$postScoreGuideWindowRRA <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'guidePerFixedWindow'){
      out_specs_list$guidePerFixedWindow <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'crestEdgeRfilter'){
      out_specs_list$crestEdgeRfilter <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'simulated_data'){
      out_specs_list$simulated_data <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'glmm_positiveTraining'){
      out_specs_list$glmm_positiveTraining <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'glmm_negativeTraining'){
      out_specs_list$glmm_negativeTraining <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'glmm_untrainedNull'){
      out_specs_list$glmm_untrainedNull <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'dirichlet_null'){
      out_specs_list$dirichlet_null <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'glmm_allData_positiveTraining'){
      out_specs_list$glmm_allData_positiveTraining <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'glmm_allData_negativeTraining'){
      out_specs_list$glmm_allData_negativeTraining <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'glmm_reducedData_positiveTraining'){
      out_specs_list$glmm_reducedData_positiveTraining <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'glmm_reducedData_negativeTraining'){
      out_specs_list$glmm_reducedData_negativeTraining <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'glmm_neg_trainingFraction'){
      out_specs_list$glmm_neg_trainingFraction <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'searchSilencers'){
      out_specs_list$searchSilencers <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'hmmStates'){
      out_specs_list$hmmStates <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'transition_prob'){
      out_specs_list$transition_prob <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'hmmStartStateProbs'){
      out_specs_list$hmmStartStateProbs <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'transition_type'){
      out_specs_list$transition_type <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'dataDir'){
      out_specs_list$dataDir <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'simName'){
      out_specs_list$simName <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'scorePlotting'){
      out_specs_list$scorePlotting <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'nrSims'){
      out_specs_list$nrSims <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'categoryProbs'){
      out_specs_list$categoryProbs <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'nrFreeCores'){
      out_specs_list$nrFreeCores <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'parallelOutfile'){
      out_specs_list$parallelOutfile <- strsplit(spec,':')[[1]][2]
    }

    # if(spec_id == 'dataDir'){
    #   out_specs_list$dataDir <- paste(strsplit(spec,':')[[1]][2:length(strsplit(spec,':')[[1]])],collapse = ':')
    # }
    # if(spec_id == 'enhancerPosFileLoc'){
    #   out_specs_list$enhancerPosFileLoc <- paste(strsplit(spec,':')[[1]][2:length(strsplit(spec,':')[[1]])],collapse = ':')
    # }
    # if(spec_id == 'pcrDuplicates'){
    #   out_specs_list$pcrDuplicates <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'includeSimSeed'){
    #   out_specs_list$includeSimSeed <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'CRISPReffectRange'){
    #   out_specs_list$CRISPReffectRange <- as.numeric(strsplit(spec,':')[[1]][2])
    # }
    # if(spec_id == 'negNonTarget'){
    #   out_specs_list$negNonTarget <- as.numeric(strsplit(spec,':')[[1]][2])
    # }
    # if(spec_id == 'Method'){
    #   out_specs_list$Method <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    # }
    # if(spec_id == 'foldChangePaired'){
    #   out_specs_list$foldChangePaired <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    # }
    # if(spec_id == 'subTileMethod'){
    #   out_specs_list$subTileMethod <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    # }
    # if(spec_id == 'dataName'){
    #   out_specs_list$dataName <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'DataType'){
    #   out_specs_list$DataType <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'PerformSubTiling'){
    #   if(strsplit(spec,':')[[1]][2] == 'yes') out_specs_list$PerformSubTiling <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'PerformDataAnalysis'){
    #   out_specs_list$PerformDataAnalysis <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'AnalyzeSimData'){
    #   out_specs_list$AnalyzeSimData <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'mageckTest'){
    #   mtCol <- strsplit(spec,':')[[1]][2]
    #   if(mtCol == 'low') out_specs_list$mageckTest <- 11 else if(mtCol == 'high') out_specs_list$mageckTest <- 12 else out_specs_list$mageckTest <- 13
    # }
    # if(spec_id == 'RRApath'){
    #   out_specs_list$RRApath <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    # }
    # if(spec_id == 'RRAplotBinLo'){
    #   out_specs_list$RRAplotBinLo <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    # }
    # if(spec_id == 'RRAplotBinPvalue'){
    #   out_specs_list$RRAplotBinPvalue <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    # }
    # if(spec_id == 'RRAplotBinFDR'){
    #   out_specs_list$RRAplotBinFDR <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    # }
    # if(spec_id == 'RRAplotBinGoodSGR'){
    #   out_specs_list$RRAplotBinGoodSGR <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    # }
    # if(spec_id == 'binSize'){
    #   out_specs_list$binSize <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    # }
    # if(spec_id == 'RRAmaxPercent'){
    #   out_specs_list$RRAmaxPercent <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    # }
    # if(spec_id == 'RRAgoodsgRNA'){
    #   out_specs_list$RRAgoodsgRNA <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    # }
    # if(spec_id == 'RRAminPeakWidth'){
    #   out_specs_list$RRAminPeakWidth <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    # }
    # if(spec_id == 'AdditionalAnalysisSteps'){
    #   if(strsplit(spec,':')[[1]][2] == 'yes') out_specs_list$AdditionalAnalysisSteps <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'genePosFileLoc'){
    #   out_specs_list$plotGenes <- 'yes'
    #   plotGenesDataLoc <- paste(strsplit(spec,':')[[1]][2:length(strsplit(spec,':')[[1]])],collapse = ':')
    #   out_specs_list$plotGenesData <- read.csv(plotGenesDataLoc,as.is = TRUE,header=TRUE)
    # }
    # if(spec_id == 'binCellCounts'){
    #   out_specs_list$binCellCounts <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    # }
    # if(spec_id == 'unsortedCells'){
    #   out_specs_list$unsortedCells <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    # }
    # if(spec_id == 'nrBetCtrl'){
    #   out_specs_list$nrBetCtrl <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    # }
    # if(spec_id == 'alphaMLKStart'){
    #   out_specs_list$alphaMLKStart <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    # }
    # if(spec_id == 'betalMLKStart'){
    #   out_specs_list$betalMLKStart <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    # }
    # if(spec_id == 'probGfunc'){
    #   out_specs_list$probGfunc <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    # }
    # if(spec_id == 'SavePvalues'){
    #   out_specs_list$SavePvalues <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'analysisType'){
    #   out_specs_list$analysisType <- strsplit(spec,':')[[1]][2]
    # }
    # #user specification whether the filtered count file is to be saved
    # if(spec_id == 'SaveCountFile'){
    #   out_specs_list$SaveCountFile <- strsplit(spec,':')[[1]][2]
    # }
    # #user specification whether the filtered info file is to be saved
    # if(spec_id == 'SaveInfoFile'){
    #   out_specs_list$SaveInfoFile <- strsplit(spec,':')[[1]][2]
    # }
    # #user specification whether the information about the rows/guides removed is to be saved
    # if(spec_id == 'SaveFilteringRemoved'){
    #   out_specs_list$SaveFilteringRemoved <- strsplit(spec,':')[[1]][2]
    # }
    # #user specification whether the information about the rates, dispersion estimates and log-likelihoods are to be saved
    # if(spec_id == 'saveRateLLR'){
    #   if(strsplit(spec,':')[[1]][2] == 'yes') out_specs_list$saveRateLLR <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'LLRplots'){
    #   if(strsplit(spec,':')[[1]][2] == 'yes') out_specs_list$LLRplots <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'QQPlots'){
    #   if(strsplit(spec,':')[[1]][2] == 'yes') out_specs_list$QQPlots <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'VolcanoPlots'){
    #   if(strsplit(spec,':')[[1]][2] == 'yes') out_specs_list$VolcanoPlots <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'CountRatioPlots'){
    #   if(strsplit(spec,':')[[1]][2] == 'yes') out_specs_list$CountRatioPlots <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'RateRatioPlots'){
    #   if(strsplit(spec,':')[[1]][2] == 'yes') out_specs_list$RateRatioPlots <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'chromosomePlots'){
    #   if(strsplit(spec,':')[[1]][2] == 'yes') out_specs_list$chromosomePlots <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'AUCeval'){
    #   if(strsplit(spec,':')[[1]][2] == 'yes') out_specs_list$AUCeval <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'prAUCeval'){
    #   if(strsplit(spec,':')[[1]][2] == 'yes') out_specs_list$prAUCeval <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'PlotSeparateChromsomes'){
    #   out_specs_list$PlotSeparateChromsomes <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'InitialDispersion'){
    #   out_specs_list$InitialDispersion <- as.numeric(strsplit(spec,':')[[1]][2])
    # }
    # if(spec_id == 'InitialEta'){
    #   out_specs_list$InitialEta <- as.numeric(strsplit(spec,':')[[1]][2])
    # }
    # if(spec_id == 'NoDataFiltering'){
    #   if(strsplit(spec,':')[[1]][2] == 'yes') out_specs_list$NoDataFiltering <- strsplit(spec,':')[[1]][2]
    # }
    # if(spec_id == 'labelNames'){
    #   out_specs_list$labelNames <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    # }
    # if(spec_id == 'colorLabels'){
    #   out_specs_list$colorLabels <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    # }
    # if(spec_id == 'positiveLabels'){
    #   out_specs_list$positiveLabels <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    # }
    # if(spec_id == 'negativeLabels'){
    #   out_specs_list$negativeLabels <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    # }
    # if(spec_id == 'colorTypes'){
    #   out_specs_list$colorTypes <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    # }
    # if(spec_id == 'pValThresh'){
    #   out_specs_list$pValThresh <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    # }
    # if(spec_id == 'ZINBminThresh'){
    #   out_specs_list$ZINBminThresh <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    # }
    # if(spec_id == 'mageckNrDiff'){
    #   out_specs_list$mageckNrDiff <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    # }
    # if(spec_id == 'mageckNumDiff'){
    #   out_specs_list$mageckNumDiff <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    # }
    # if(spec_id == 'pValCol'){
    #   out_specs_list$pValCol <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    # }
    # if(spec_id == 'zoomInPlots'){
    #   out_specs_list$zoomInPlots <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    # }
    # if(spec_id == 'zoomRange'){
    #   out_specs_list$zoomRange <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    # }
  }
  return(out_specs_list)
}

################################################################################
write_specs_file <- function(input.specs.list, input.filename){

  specs.names <- names(input.specs.list)
  for(name in specs.names){
    temp_line <- paste(name,':',paste(input.specs.list[[name]],collapse = ','),sep='')
    write(temp_line,file = paste(input.filename, '.txt', sep = ''),append = TRUE)
  }

}

create_default_specs <- function(input.name){
  default.list <- list()
  default.list$CountFileLoc <- 'counts.csv'
  default.list$sgRNAInfoFileLoc <- 'info.csv'
  default.list$Group1 <- c(1,2,3,4)
  default.list$Group2 <- c(5,6,7,8)
  default.list$positiveLabels <- 'exon'
  default.list$negativeLabels <- 'neg'
  default.list$CRISPReffectRange <- 20
  default.list$genePosFileLoc <- 'genes.csv'
  default.list$Method <- c('NB','edgeR','FoldChange','ZINB', 'NB-GLMM','Viterbi', 'Frwd-Bkwd')
  default.list$foldChangePaired <- 'yes'
  default.list$evaluatePerformance <- 'yes'
  default.list$plotSeparateChroms <- 'yes'
  default.list$labelHierarchy <-  c('chr','neg','exon')
  default.list$dataFilter <- 'none'
  default.list$CRESTrra <- '/iblm/netapp/home/pfiaux/CRISPR_screening_project/crest/bin'
  default.list$rraPermutNr <- 100

  default.list$dataName <- input.name
  default.list$dataArranging <- 'none'

  return(default.list)
}

# analysis.name <- 'u_lo09_med_vs_hi_regLab'  # used as $dataName
# analysis.specs <- create_default_specs(analysis.name)
# analysis.specs$Group1 <- c(12,16)
# analysis.specs$Group2 <- c(10,14)
# analysis.specs$CountFileLoc <- '../formatted_guide_Unique.csv'
# analysis.specs$sgRNAInfoFileLoc <- '../final_guide_regionInfo.csv'
# analysis.specs$genePosFileLoc <- '/iblm/netapp/home/pfiaux/CREST_seq_data/Data/gata3_hg37_formatted_unique_exon_regions.csv'
# analysis.specs$foldChangePaired <- 'no'
# analysis.specs$zoomRange <- '../GATA3_zoom_region1.csv'
# analysis.specs$positiveLabels <- c('exon1','exon2','exon3','exon4','exon5','exon6')
# analysis.specs$labelHierarchy <- c('chr','neg','exon1','exon2','exon3','exon4','exon5','exon6')
# analysis.specs$tileZoomRange <- '../gata3_tileZoom_hg37_formatted_unique_exon_regions.csv'
# analysis.specs$saveGuideScores <- 'yes'
# analysis.specs$analysisDir <- getwd()
# none_no_analysis(analysis.specs)

# analyzing data without any pre- or post-binning
none_no_analysis <- function(input.par.list){
  # create analysis output directory
  analysis.name <- paste0(input.par.list$dataName, '_none_no')
  analysis.dir <- paste0(input.par.list$analysisDir, '/', analysis.name)
  dir.create(file.path(analysis.dir))
  setwd(file.path(analysis.dir))
  input.par.list$dataName <- analysis.name
  write_specs_file(input.par.list, analysis.name)
  analyze_data(paste(analysis.name,'.txt', sep = ''))
  setwd(file.path(input.par.list$analysisDir))
}

# depricated
none_no_multiSimAnalysis <- function(input.par.list){
  # create analysis output directory
  analysis.name <- paste0(input.par.list$dataName, '_none_no')
  analysis.dir <- paste0(input.par.list$analysisDir, '/', analysis.name)
  dir.create(file.path(analysis.dir))
  setwd(file.path(analysis.dir))
  input.par.list$dataName <- analysis.name
  write_specs_file(input.par.list, analysis.name)
  evaluate_multi_simulations(paste(analysis.name,'.txt', sep = ''))
  setwd(file.path(input.par.list$analysisDir))
}

# previously: none_no_multiSimAnalysis
multiSimGuideAnalysis <- function(input.par.list){
  # create analysis output directory
  analysis.name <- paste0(input.par.list$dataName, '_simAnalysis_perGuide')
  analysis.dir <- paste0(input.par.list$analysisDir, '/', analysis.name)
  dir.create(file.path(analysis.dir))
  setwd(file.path(analysis.dir))
  input.par.list$dataName <- analysis.name
  input.par.list$noPostScoreAnalysis <- 'yes'
  write_specs_file(input.par.list, analysis.name)
  evaluate_multi_simulations(paste(analysis.name,'.txt', sep = ''))
  setwd(file.path(input.par.list$analysisDir))
}

# analyzing data without pre-binning but post-binning using 50pb
none_y50_analysis <- function(input.par.list){
  # create analysis output directory
  analysis.name <- paste0(input.par.list$dataName, '_none_y50')
  analysis.dir <- paste0(input.par.list$analysisDir, '/', analysis.name)
  dir.create(file.path(analysis.dir))
  setwd(file.path(analysis.dir))
  input.par.list$dataName <- analysis.name
  input.par.list$postScoringAnalysis <- 'yes'
  # #input.par.list$postScoringCREST <- 'yes'
  input.par.list$postScoreAlphaRRA <- 'yes'
  # input.par.list$binSize <- 50
  input.par.list$postScoreSlidingWindow <- 'yes'
  write_specs_file(input.par.list, analysis.name)
  analyze_data(paste(analysis.name,'.txt', sep = ''))
  setwd(file.path(input.par.list$analysisDir))
}

# analyzing data without pre-binning but post-binning using 50pb
none_y50__multiSimAnalysis <- function(input.par.list){
  # create analysis output directory
  analysis.name <- paste0(input.par.list$dataName, '_none_y50')
  analysis.dir <- paste0(input.par.list$analysisDir, '/', analysis.name)
  dir.create(file.path(analysis.dir))
  setwd(file.path(analysis.dir))
  input.par.list$dataName <- analysis.name
  input.par.list$postScoringAnalysis <- 'yes'
  # input.par.list$postScoringCREST <- 'yes'
  input.par.list$postScoreAlphaRRA <- 'yes'
  input.par.list$binSize <- 50
  input.par.list$postScoreSlidingWindow <- 'yes'
  write_specs_file(input.par.list, analysis.name)
  evaluate_multi_simulations(paste(analysis.name,'.txt', sep = ''))
  setwd(file.path(input.par.list$analysisDir))
}

# previously: none_y50__multiSimAnalysis
multiSimGenomeScoreAnalysis <- function(input.par.list){
  # create analysis output directory
  analysis.name <- paste0(input.par.list$dataName, '_simAnalysis_perGenomeRegion')
  analysis.dir <- paste0(input.par.list$analysisDir, '/', analysis.name)
  dir.create(file.path(analysis.dir))
  setwd(file.path(analysis.dir))
  input.par.list$dataName <- analysis.name
  # input.par.list$postScoringAnalysis <- 'yes'
  # # input.par.list$postScoringCREST <- 'yes'
  # input.par.list$postScoreAlphaRRA <- 'yes'
  # input.par.list$binSize <- 50
  # input.par.list$postScoreSlidingWindow <- 'yes'
  write_specs_file(input.par.list, analysis.name)
  evaluate_multi_simulations(paste(analysis.name,'.txt', sep = ''))
  setwd(file.path(input.par.list$analysisDir))
}

# analyzing data without pre-binning but post-binning using 50pb
none_guide20window_analysis <- function(input.par.list){
  # create analysis output directory
  analysis.name <- paste0(input.par.list$dataName, '_none_guide20window')
  analysis.dir <- paste0(input.par.list$analysisDir, '/', analysis.name)
  dir.create(file.path(analysis.dir))
  setwd(file.path(analysis.dir))
  input.par.list$dataName <- analysis.name
  input.par.list$postScoringAnalysis <- 'yes'
  input.par.list$postScoreGuideWindowRRA <- 'yes'
  input.par.list$guidePerFixedWindow <- 20
  write_specs_file(input.par.list, analysis.name)
  analyze_data(paste(analysis.name,'.txt', sep = ''))
  setwd(file.path(input.par.list$analysisDir))
}

# analyzing data without pre-binning but post-binning using 50pb
none_guide10window_analysis <- function(input.par.list){
  # create analysis output directory
  analysis.name <- paste0(input.par.list$dataName, '_none_guide10window')
  analysis.dir <- paste0(input.par.list$analysisDir, '/', analysis.name)
  dir.create(file.path(analysis.dir))
  setwd(file.path(analysis.dir))
  input.par.list$dataName <- analysis.name
  input.par.list$postScoringAnalysis <- 'yes'
  input.par.list$postScoreGuideWindowRRA <- 'yes'
  input.par.list$guidePerFixedWindow <- 10
  write_specs_file(input.par.list, analysis.name)
  analyze_data(paste(analysis.name,'.txt', sep = ''))
  setwd(file.path(input.par.list$analysisDir))
}

# analyzing data without pre-binning but post-binning using 50pb
none_guide50window_analysis <- function(input.par.list){
  # create analysis output directory
  analysis.name <- paste0(input.par.list$dataName, '_none_guide50window')
  analysis.dir <- paste0(input.par.list$analysisDir, '/', analysis.name)
  dir.create(file.path(analysis.dir))
  setwd(file.path(analysis.dir))
  input.par.list$dataName <- analysis.name
  input.par.list$postScoringAnalysis <- 'yes'
  input.par.list$postScoreGuideWindowRRA <- 'yes'
  input.par.list$guidePerFixedWindow <- 50
  write_specs_file(input.par.list, analysis.name)
  analyze_data(paste(analysis.name,'.txt', sep = ''))
  setwd(file.path(input.par.list$analysisDir))
}

# analyzing data without pre-binning but post-binning using 50pb
none_y200_analysis <- function(input.par.list){
  # create analysis output directory
  analysis.name <- paste0(input.par.list$dataName, '_none_y200')
  analysis.dir <- paste0(input.par.list$analysisDir, '/', analysis.name)
  dir.create(file.path(analysis.dir))
  setwd(file.path(analysis.dir))
  input.par.list$dataName <- analysis.name
  input.par.list$postScoringAnalysis <- 'yes'
  input.par.list$postScoringCREST <- 'yes'
  input.par.list$postScoreAlphaRRA <- 'yes'
  input.par.list$binSize <- 200
  write_specs_file(input.par.list, analysis.name)
  analyze_data(paste(analysis.name,'.txt', sep = ''))
  setwd(file.path(input.par.list$analysisDir))
}

# analyzing data without pre-binning but post-binning using 500pb
none_y500_analysis <- function(input.par.list){
  # create analysis output directory
  analysis.name <- paste0(input.par.list$dataName, '_none_y500')
  analysis.dir <- paste0(input.par.list$analysisDir, '/', analysis.name)
  dir.create(file.path(analysis.dir))
  setwd(file.path(analysis.dir))
  input.par.list$dataName <- analysis.name
  input.par.list$postScoringAnalysis <- 'yes'
  input.par.list$postScoringCREST <- 'yes'
  input.par.list$postScoreAlphaRRA <- 'yes'
  input.par.list$binSize <- 500
  write_specs_file(input.par.list, analysis.name)
  analyze_data(paste(analysis.name,'.txt', sep = ''))
  setwd(file.path(input.par.list$analysisDir))
}

# analyzing data without pre-binning but post-binning using 50pb
nonOverl_700_sum_no_analysis <- function(input.par.list){
  # create analysis output directory
  analysis.name <- paste0(input.par.list$dataName, '_nonOverl_700_sum_no')
  analysis.dir <- paste0(input.par.list$analysisDir, '/', analysis.name)
  dir.create(file.path(analysis.dir))
  setwd(file.path(analysis.dir))
  input.par.list$dataName <- analysis.name
  input.par.list$dataArranging <- 'yes'
  input.par.list$preBinning <- 700
  input.par.list$binBy <- 'sum'
  write_specs_file(input.par.list, analysis.name)
  analyze_data(paste(analysis.name,'.txt', sep = ''))
  setwd(file.path(input.par.list$analysisDir))
}

# analyzing data without pre-binning but post-binning using 50pb
overl_350_700_sum_no_analysis <- function(input.par.list){
  # create analysis output directory
  analysis.name <- paste0(input.par.list$dataName, '_overl_700_350_sum_no')
  analysis.dir <- paste0(input.par.list$analysisDir, '/', analysis.name)
  dir.create(file.path(analysis.dir))
  setwd(file.path(analysis.dir))
  input.par.list$dataName <- analysis.name
  input.par.list$dataArranging <- 'yes'
  input.par.list$preOverlapBinning <- c(350, 700)
  input.par.list$binBy <- 'sum'
  write_specs_file(input.par.list, analysis.name)
  analyze_data(paste(analysis.name,'.txt', sep = ''))
  setwd(file.path(input.par.list$analysisDir))
}

# analyzing data without pre-binning but post-binning using 50pb
none_y1500_analysis <- function(input.par.list){
  # create analysis output directory
  analysis.name <- paste0(input.par.list$dataName, '_none_y1500')
  analysis.dir <- paste0(input.par.list$analysisDir, '/', analysis.name)
  dir.create(file.path(analysis.dir))
  setwd(file.path(analysis.dir))
  input.par.list$dataName <- analysis.name
  input.par.list$postScoringAnalysis <- 'yes'
  input.par.list$postScoringCREST <- 'yes'
  input.par.list$postScoreAlphaRRA <- 'yes'
  input.par.list$binSize <- 1500
  write_specs_file(input.par.list, analysis.name)
  analyze_data(paste(analysis.name,'.txt', sep = ''))
  setwd(file.path(input.par.list$analysisDir))
}

# analyzing data without pre-binning but post-binning using 50pb
none_y3000_analysis <- function(input.par.list){
  # create analysis output directory
  analysis.name <- paste0(input.par.list$dataName, '_none_y3000')
  analysis.dir <- paste0(input.par.list$analysisDir, '/', analysis.name)
  dir.create(file.path(analysis.dir))
  setwd(file.path(analysis.dir))
  input.par.list$dataName <- analysis.name
  input.par.list$postScoringAnalysis <- 'yes'
  input.par.list$postScoringCREST <- 'yes'
  input.par.list$postScoreAlphaRRA <- 'yes'
  input.par.list$binSize <- 3000
  write_specs_file(input.par.list, analysis.name)
  analyze_data(paste(analysis.name,'.txt', sep = ''))
  setwd(file.path(input.par.list$analysisDir))
}

# analyzing data without pre-binning but post-binning using 50pb
none_y5000_analysis <- function(input.par.list){
  # create analysis output directory
  analysis.name <- paste0(input.par.list$dataName, '_none_y5000')
  analysis.dir <- paste0(input.par.list$analysisDir, '/', analysis.name)
  dir.create(file.path(analysis.dir))
  setwd(file.path(analysis.dir))
  input.par.list$dataName <- analysis.name
  input.par.list$postScoringAnalysis <- 'yes'
  input.par.list$postScoringCREST <- 'yes'
  input.par.list$postScoreAlphaRRA <- 'yes'
  input.par.list$binSize <- 5000
  write_specs_file(input.par.list, analysis.name)
  analyze_data(paste(analysis.name,'.txt', sep = ''))
  setwd(file.path(input.par.list$analysisDir))
}
