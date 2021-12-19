#' Importance Sampling (IS) in Bayesian Networks
#'
#' Outputs the normalizing constant given some evidence of a Bayesian network
#'
#' @param BayesNet Bayesian network
#' @param obs list containing the evidence nodes and associated values
#' @return relevant DAG
#' @export
sample.normConst <- function(BayesNet, obs, s_method = "SGS", N_samples = 100, relevantSubNet = TRUE, plot = TRUE)
{ 
  # Importance sampling in Bayesian networks given the evidence
  # Returned is the normalzing constant
  
  # Compute the updated Bayesian network using BP 
  sortUpdatedNet <- FALSE
  if (s_method=="BP-LW"){
    engineTemp <- InferenceEngine(BayesNet)
    engineTemp <- belief.propagation(engineTemp, obs)
    updatedNet <- engineTemp@updated.bn
    sortUpdatedNet <- TRUE
  }
  # Compute the updated Bayesian network using LBP 
  if (s_method=="LBP-LW"){
    engineTemp <- InferenceEngine(BayesNet)
    engineTemp <- loopy_belief.propagation(engineTemp, obs)
    updatedNet <- engineTemp@updated.bn
    sortUpdatedNet <- TRUE
  }
  # Compute the updated Bayesian network using subgroup BP
  if (s_method=="SGS"){
    updatedNet <- subGroup_belief.propagation(BayesNet, obs, group_limit = 15)
    sortUpdatedNet <- TRUE
  }
  
  net <- BayesNet
  lengthBN <- length(net@cpts)
  variablesBN <- variables(net)
  Nparents <- sapply(c(1:lengthBN),function(x) length(dim(net@cpts[[x]]))-1)
  
  observed.vars <- obs$observed.vars
  observed.vals <- obs$observed.vals
  
  # sort CPTs according to dimnames
  nodeNamesCPT <- sapply(c(1:lengthBN),function(x) match(names(dimnames(net@cpts[[x]])),variablesBN))
  for(g0 in 1:lengthBN){
    nodePosition <- match(g0, nodeNamesCPT[[g0]])
    tempPerm <- c(1:length(nodeNamesCPT[[g0]]))
    tempPerm <- c(tempPerm[-nodePosition],tempPerm[nodePosition])
    net@cpts[[g0]] <- aperm(net@cpts[[g0]], perm = tempPerm)
  }
  nodeParents <- sapply(c(1:lengthBN),function(x) match(head(names(dimnames(net@cpts[[x]])),-1),variablesBN))
  
  # repeat the sorting for the updatedNet
  if (sortUpdatedNet){
    # sort CPTs according to dimnames
    nodeNamesCPT2 <- sapply(c(1:lengthBN),function(x) match(names(dimnames(updatedNet@cpts[[x]])),variablesBN))
    for(g1 in 1:lengthBN){
      nodePosition <- match(g1, nodeNamesCPT2[[g1]])
      tempPerm <- c(1:length(nodeNamesCPT2[[g1]]))
      tempPerm <- c(tempPerm[-nodePosition],tempPerm[nodePosition])
      updatedNet@cpts[[g1]] <- aperm(updatedNet@cpts[[g1]], perm = tempPerm)
    }
    nodeParents2 <- sapply(c(1:lengthBN),function(x) match(head(names(dimnames(updatedNet@cpts[[x]])),-1),variablesBN))
  }
  
  # # faster version, but might mix up the ordering
  # nodeParents <- sapply(c(1:lengthBN),function(x) get.allParents(dag(net),x))
  
  # initialize vector of starting variables
  curVec <- rep(1, lengthBN)
  curVec[observed.vars] <- observed.vals
  
  # get topological order for sampling order
  topolOrder <- topo_sort(graph_from_adjacency_matrix(dag(net), "directed"), mode = c("out"))
  
  # get relevant subnetwork
  if (relevantSubNet==FALSE){
    consideredNodes <- c(1:lengthBN)
  }else{
    consideredNodes <- get.allAncestors(dag(net), observed.vars) 
  }
  
  # define range of sampled variables
  vars1 <- intersect(topolOrder, setdiff(consideredNodes,observed.vars))
  vars2 <- intersect(topolOrder, observed.vars)
  
  # intialize partition function
  partitionFuncMean <- 0
  partitionFunc <- rep(0, N_samples)
  partitionFuncTempAll <- rep(0, N_samples)
  
  if(s_method=="GS"){
    nodeChildren <- sapply(c(1:lengthBN),function(x) get.allChildren(dag(net),x))
  }
  
  for (j1 in 1:N_samples){
    
    # reset partition function
    funcP1 <- 0
    funcP2 <- 0
    # funcP1 <- 1
    # funcP2 <- 1
    
    for (j2 in vars1){
      # get probabilities for the sampling
      if(s_method=="FS"){
        # get the prob. of node given the parents
        tempProb <- get.probNode(net, curVec, j2, nodeParents)
        tempProb2 <- tempProb
      }else if(s_method=="GS"){
        # get the prob. of node given the Markov blanket
        probsMB <- get.probNodeMB(net, curVec, j2, nodeParents, nodeChildren, consideredNodes)
        tempProb <- probsMB$probMB
        tempProb2 <- probsMB$tempProbNode
      }else if(s_method %in% c("BP-LW", "LBP-LW", "SGBP-LW","SGS")){
        tempProb <- get.probNode(updatedNet, curVec, j2, nodeParents2)
        tempProb2 <- get.probNode(net, curVec, j2, nodeParents)
      }else{
        stop('s_method is not valid')
      }
      
      # sample
      curVecNode <- sample(1:length(tempProb), size = 1, prob=tempProb)
      curVec[j2] <- curVecNode
      tempProbNode <- tempProb[curVecNode]
      
      # calculate importance function (use log scale to handle small numbers)
      funcP1 <- funcP1 + log(tempProb[curVecNode])
      # funcP1 <- funcP1*mpfr(tempProb[curVecNode], 50)
      # calculate probability function for weighting (use log scale to handle small numbers)
      funcP2 <- funcP2 + log(tempProb2[curVecNode])
      # funcP2 <- funcP2*mpfr(tempProb2[curVecNode], 50)
    }
    
    for (j3 in vars2){
      parentsAndNode <- c(nodeParents[[j3]],j3)
      curVecPar <- rbind(curVec[parentsAndNode])
      tempProb <- net@cpts[[j3]][curVecPar]
      
      # no sampling for evidence nodes
      tempProbNode <- tempProb
      
      # calculate probability function for weighting (use log scale to handle small numbers)
      funcP2 <- funcP2 + log(tempProbNode)
      # funcP2 <- funcP2*mpfr(tempProbNode,50)
    }
    
    # calculate partition function (use log scale to handle small numbers)
    
    partitionFuncTemp <- exp(funcP2 - funcP1)
    # partitionFuncTemp <- as.numeric(funcP2/funcP1)
    
    # caluclate mean of partition function
    partitionFuncMean <- partitionFuncMean*(j1-1)/j1 + partitionFuncTemp/j1
    
    partitionFunc[j1] <- partitionFuncMean
    partitionFuncTempAll[j1] <- partitionFuncTemp
  }
  
  results <- list("NormConst"=partitionFunc[N_samples],"partitionFunc"=partitionFunc,"partitionFuncTempAll"=partitionFuncTempAll)
  
  if (plot){
    plot.normConstCalc(results)
    #plot(1:N_samples, partitionFunc, type="l", col="blue", main="Parition Function", xlab="Sampling Step", ylab="Normalizing Constant")
  }
  
  return(results)
  
}


get.probNode <- function(net, curVec, node, nodeParents)
{
  # get probability of node given vector
  cpts <- cpts(net)
  parents <- nodeParents[[node]]
  if(length(parents)==0){
    tempProb <- cpts[[node]]
  }else{
    curVecPar <- curVec[parents]
    mystring <- paste("cpts[[node]][", paste(as.character(curVecPar), collapse=", "),",]")
    tempProb <- eval(parse(text=mystring))
  }
  
  return(tempProb)
}

get.probNodeMB <- function(net, curVec, node, nodeParents, nodeChildren, consideredNodes)
{
  # get probability of Markov blanket of node (used in Gibbs sampling)
  curVecTemp <- curVec
  
  tempProbNode <- get.probNode(net, curVec, node, nodeParents)
  
  tempChildren <- get.allChildren(dag(net), node)
  
  tempChildren <- intersect(tempChildren, consideredNodes)
  
  probMB <- rep(NA, length(tempProbNode))
  
  for (h2 in 1:length(tempProbNode)){
    curVecTemp[node] <- h2
    ProbChildAll <- 1
    for (h1 in tempChildren){
      tempProbChild <- get.probNode(net, curVecTemp, h1, nodeParents)
      ProbChildAll <- ProbChildAll*tempProbChild[curVec[h1]]
    }
    probMB[h2] <- ProbChildAll*tempProbNode[h2]
  }
  
  
  probMB <- probMB/sum(probMB)
  
  return(list("probMB"=probMB, "tempProbNode"=tempProbNode))
}

plot.normConstCalc <- function(NormCalcList)
{
  # Plot the progress
  
  curVariance <- sapply(c(1:length(NormCalcList$partitionFuncTempAll)), function(x) sd(NormCalcList$partitionFuncTempAll[1:x])/sqrt(x))
  lowBound <- NormCalcList$partitionFunc-curVariance
  upBound <- NormCalcList$partitionFunc+curVariance
  
  df <- data.frame(normC=NormCalcList$partitionFunc,countsS=c(1:length(NormCalcList$partitionFunc)),lowBound=lowBound,upBound=upBound,type="Cool core")
  
  color_list <- '#EE6677'
  
  ggplot(data=df, aes(y = normC, x = countsS))+
    geom_line(aes(y = normC, x = countsS,color=color_list)) + 
    geom_ribbon(aes(ymin=lowBound, ymax=upBound, fill=color_list),alpha=0.5) + 
    theme(legend.position = "none") +
    xlab("Sampling Step") + 
    ylab("Normalizing Constant")
  
}


plot.samplingProgress <- function(sampleResults, exactValue = NULL)
{
  # Plot the progress
  
  lengthMethods <- length(sampleResults)
  
  # set colors
  if(lengthMethods>3){
    N_colors <- lengthMethods
    color_list <- brewer.pal(n = N_colors, name = 'RdYlBu')
  }else{
    N_colors <- 3
    color_list <- brewer.pal(n = N_colors, name = 'RdYlBu')[1:lengthMethods]
  }
  # without Rcolorbrewer:
  # color_list <- c("#EE6677", "#228833", "#4477AA")
  
  dfAll <- NULL
  Smethod <- NULL
  
  for (t0 in 1:length(sampleResults)){
    NormCalcList <- sampleResults[[t0]]
    curVariance <- sapply(c(1:length(NormCalcList$partitionFuncTempAll)), 
                          function(x) sd(NormCalcList$partitionFuncTempAll[1:x])/sqrt(x))
    lowBound <- NormCalcList$partitionFunc-curVariance
    upBound <- NormCalcList$partitionFunc+curVariance
    timeScale <- c(1:length(NormCalcList$partitionFunc))/length(NormCalcList$partitionFunc)*NormCalcList$elapsed_time[[1]]
    
    df <- data.frame(normC=NormCalcList$partitionFunc,
                     countsS=timeScale,
                     lowBound=lowBound,upBound=upBound)
    df$group <- t0
    df$Method <- NormCalcList$method
    
    Smethod <- c(Smethod,NormCalcList$method)
    dfAll <- rbind(dfAll, df)
  }
  
  if(length(exactValue)==0){
    ggplot(data=dfAll, aes(y = normC, x = countsS, group=group))+
      geom_line(aes(y = normC, x = countsS,colour=Method)) + 
      geom_ribbon(aes(ymin=lowBound, ymax=upBound, fill=Method),alpha=0.5) + 
      scale_colour_manual(values=color_list) +
      scale_fill_manual(values=color_list) +
      # guides(fill = FALSE) +
      guides(col = FALSE) +
      xlab("Elapsed Time (s)") + 
      ylab("Normalizing Constant")
  }else{
    ggplot(data=dfAll, aes(y = normC, x = countsS, group=group))+
      geom_line(aes(y = normC, x = countsS,colour=Method)) + 
      geom_ribbon(aes(ymin=lowBound, ymax=upBound, fill=Method),alpha=0.5) + 
      scale_colour_manual(values=color_list) +
      scale_fill_manual(values=color_list) +
      geom_hline(yintercept=exactValue, linetype="dashed", color = "black") +
      # guides(fill = FALSE) +
      guides(col = FALSE) +
      xlab("Elapsed Time (s)") + 
      ylab("Normalizing Constant")
  }
  
  
}

benchmark <- function(BayesNet = NULL, obs = NULL, methods = c("FS","GS", "SGS"), N_samples = 500, N_nodes = 30)
{
  # benchmark different sampling schemes
  
  # sample Bayesian network and observation if missing
  if(length(BayesNet)==0){
    N_nodes <- N_nodes #10
    N_neighbours <- 2 
    nodeDim <- 2
    BayesNet <- randomBN(N_nodes, N_neighbours, nodeDim, visualize = FALSE)
  }
  exactValue <- NULL
  if(length(obs)==0){
    N_obs <- as.integer(N_nodes/2)
    obs <- list("observed.vars" = sample(1:N_nodes, N_obs), "observed.vals" = rep(1,N_obs))
    # exactValue <- (1/nodeDim)^N_obs
  }
  
  #
  sampleResults <- list()
  for (b0 in 1:length(methods)){
    # adjust the number of samples to the method as some sample faster 
    # (for plotting purposes only)
    if(methods[b0]=="GS"){
      N_samplesTemp <- as.integer(N_samples/4)
    }else if(methods[b0]=="FS"){
      N_samplesTemp <- N_samples*2
    }else{
      N_samplesTemp <- N_samples
    }
    
    # perform sampling
    start_time <- Sys.time() # measure computation time
    sampleResults[[b0]] <- sample.normConst(BayesNet, obs, s_method = methods[b0], N_samples = N_samplesTemp)
    end_time <- Sys.time()
    elapsed_time <- end_time - start_time
    sampleResults[[b0]]$elapsed_time <- elapsed_time
    sampleResults[[b0]]$method <- methods[b0]
    
    # if(methods[b0]=="SGS"){
    #   exactValue <- sampleResults[[b0]]$partitionFunc[N_samples]
    # }
  }
  
  # plot the results (i.e. sampling progress)
  if(length(exactValue)==0){
    plot.samplingProgress(sampleResults)
  }else{
    plot.samplingProgress(sampleResults, exactValue)
  }
  
}
