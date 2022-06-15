#' @importFrom stats sd
plot.normConstCalc <- function(NormCalcList)
{
  # Plot the progress
  
  normC=countsS=0
  
  curVariance <- sapply(c(1:length(NormCalcList$partitionFuncTempAll)), function(x) stats::sd(NormCalcList$partitionFuncTempAll[1:x])/sqrt(x))
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

#' @importFrom stats sd
#' @importFrom RColorBrewer brewer.pal
plot.samplingProgress <- function(sampleResults, exactValue = NULL)
{
  # Plot the progress
  
  normC=countsS=group=Method=0
  
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
                          function(x) stats::sd(NormCalcList$partitionFuncTempAll[1:x])/sqrt(x))
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



benchmark <- function(BayesNet = NULL, obs = NULL, methods = c("FS","GS", "SGS"), N_samples = 500, N_nodes = 30, exactValue = FALSE)
{
  
  # benchmark different sampling schemes
  
  # sample Bayesian network and observation if missing
  if(length(BayesNet)==0){
    N_nodes <- N_nodes #10
    N_neighbours <- 2 
    nodeDim <- 2
    BayesNet <- randomBN(N_nodes, N_neighbours, nodeDim, visualize = FALSE)
  }

  if(length(obs)==0){
    N_obs <- as.integer(N_nodes/2)
    obs <- list("observed.vars" = sample(1:N_nodes, N_obs), "observed.vals" = rep(1,N_obs))
  }
  
  # calc exact result
  if(exactValue){
    exactVal <- calcExactInference(BayesNet, obs)
  }
  
  # sample by method
  sampleResults <- list()
  for (b0 in 1:length(methods)){
    # adjust the number of samples to the method as some sample faster 
    # (for plotting purposes only)
    if(methods[b0]=="GS"){
      N_samplesTemp <- as.integer(N_samples/2.5)
    }else if(methods[b0]=="FS"){
      N_samplesTemp <- N_samples*2
    }else if(methods[b0]=="SGS"){
      N_samplesTemp <- N_samples*3
    }else{
      N_samplesTemp <- N_samples
    }
    
    
    # perform sampling
    start_time <- Sys.time() # measure computation time
    sampleResults[[b0]] <- approxInference(BayesNet, obs, s_method = methods[b0], N_samples = N_samplesTemp, returnList=TRUE)
    end_time <- Sys.time()
    elapsed_time <- end_time - start_time
    sampleResults[[b0]]$elapsed_time <- elapsed_time
    sampleResults[[b0]]$method <- methods[b0]
    
    # if(methods[b0]=="SGS"){
    #   exactValue <- sampleResults[[b0]]$partitionFunc[N_samples]
    # }
  }
  
  # plot the results (i.e. sampling progress)
  if(exactValue){
    plot.samplingProgress(sampleResults, exactVal)
  }else{
    plot.samplingProgress(sampleResults)
  }
  
}



# 
# benchmarkExtensive <- function(BayesNet = NULL, obs = NULL, methods = c("FS","GS", "SGS"), N_samples = 500, N_nodes = 30)
# {
#   # benchmark different sampling schemes
#   
#   # sample Bayesian network and observation if missing
#   if(length(BayesNet)==0){
#     N_nodes <- N_nodes #10
#     N_neighbours <- 2 
#     nodeDim <- 2
#     BayesNet <- randomBN(N_nodes, N_neighbours, nodeDim, visualize = FALSE)
#   }
#   exactValue <- NULL
#   if(length(obs)==0){
#     N_obs <- as.integer(N_nodes/2)
#     obs <- list("observed.vars" = sample(1:N_nodes, N_obs), "observed.vals" = rep(1,N_obs))
#     # exactValue <- (1/nodeDim)^N_obs
#   }
#   
#   #
#   sampleResults <- list()
#   for (b0 in 1:length(methods)){
#     # adjust the number of samples to the method as some sample faster 
#     # (for plotting purposes only)
#     if(methods[b0]=="GS"){
#       N_samplesTemp <- as.integer(N_samples/4)
#     }else if(methods[b0]=="FS"){
#       N_samplesTemp <- N_samples*2
#     }else{
#       N_samplesTemp <- N_samples
#     }
#     
#     # perform sampling
#     start_time <- Sys.time() # measure computation time
#     sampleResults[[b0]] <- approxInference(BayesNet, obs, s_method = methods[b0], N_samples = N_samplesTemp, returnList=TRUE)
#     end_time <- Sys.time()
#     elapsed_time <- end_time - start_time
#     sampleResults[[b0]]$elapsed_time <- elapsed_time
#     sampleResults[[b0]]$method <- methods[b0]
#     
#     # if(methods[b0]=="SGS"){
#     #   exactValue <- sampleResults[[b0]]$partitionFunc[N_samples]
#     # }
#   }
#   
#   # plot the results (i.e. sampling progress)
#   if(length(exactValue)==0){
#     plot.samplingProgress(sampleResults)
#   }else{
#     plot.samplingProgress(sampleResults, exactValue)
#   }
#   
# }


benchmarkStudy <- function(BayesNet = NULL, obs = NULL, methods = c("FS","GS", "SGS"), N_samples = 500, N_nodes = 30)
{
  # benchmark different sampling schemes
  
  # sample Bayesian network and observation if missing
  if(length(BayesNet)==0){
    N_nodes <- N_nodes
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
      N_samplesTemp <- as.integer(N_samples/2)
    }else if(methods[b0]=="FS"){
      N_samplesTemp <- N_samples*2
    }else if(methods[b0]=="SGS"){
      N_samplesTemp <- N_samples*3
    }else{
      N_samplesTemp <- N_samples
    }
    
    # perform sampling
    start_time <- Sys.time() # measure computation time
    sampleResults[[b0]] <- approxInference(BayesNet, obs, s_method = methods[b0], N_samples = N_samplesTemp, returnList=TRUE)
    end_time <- Sys.time()
    elapsed_time <- end_time - start_time
    sampleResults[[b0]]$elapsed_time <- elapsed_time
    sampleResults[[b0]]$method <- methods[b0]
    
    # if(methods[b0]=="SGS"){
    #   exactValue <- sampleResults[[b0]]$partitionFunc[N_samples]
    # }
  }
  
  # # plot the results (i.e. sampling progress)
  # if(length(exactValue)==0){
  #   plot.samplingProgress(sampleResults)
  # }else{
  #   plot.samplingProgress(sampleResults, exactValue)
  # }
  
  return(sampleResults)
  
}



benchmarkStudyMain <- function(BayesNet = NULL, obs = NULL, methods = c("FS","GS", "SGS"), N_samples = 500, N_rep = 10, N_nodes = 30, exactValue = TRUE)
{
  # repeat the benchmarking N_rep times

  # calc exact result
  if(exactValue){
    exactVal <- exactInference(BayesNet, obs)
  }

  benchmarkResults <- list()
  for (i in 1:N_rep){
    tempRes <- benchmarkStudy(BayesNet, obs, methods, N_samples, N_nodes)
    benchmarkResults[[i]] <- tempRes
  }



  if(exactValue){
    myResults <- processAndPlot(benchmarkResults,exactVal)

    myResults$simulationDetails <- list("BayesNet"=BayesNet, "obs"=obs, "methods"=methods, "N_samples"=N_samples, "N_rep"=N_rep)

    return(myResults)
  }#else{
    # plot.samplingProgressStudy(benchmarkResults)
  #}

}


#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggpubr ggarrange
processAndPlot <- function(benchmarkResults, exactVal)
{
  normC=countsS=group=Method=0
  
  # Plot the progress
  lengthMethods <- length(benchmarkResults[[1]])
  lengthRes <- length(benchmarkResults)
  
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
  df2All <- NULL
  Smethod <- NULL
  
  # process results for plot
  for (t0 in 1:lengthMethods){
    # NormCalcList <- benchmarkResults[[tt0]][[t0]]
    curMean <- rowMeans(sapply(1:lengthRes, function(x) benchmarkResults[[x]][[t0]]$partitionFunc))
    timeMean <- mean(sapply(1:lengthRes, function(x) benchmarkResults[[x]][[t0]]$elapsed_time[[1]]))
    curSTD <- sqrt(RowVar(sapply(1:lengthRes, function(x) benchmarkResults[[x]][[t0]]$partitionFunc)))
    MSE <- rowMeans(sapply(1:lengthRes, function(x) (benchmarkResults[[x]][[t0]]$partitionFunc-exactVal)^2))
    # RMSE <- sqrt(rowMeans(sapply(1:lengthRes, function(x) (benchmarkResults[[x]][[t0]]$partitionFunc-exactVal)^2)))
    NRMSE <- sqrt(rowMeans(sapply(1:lengthRes, function(x) (benchmarkResults[[x]][[t0]]$partitionFunc-exactVal)^2)))/exactVal
    # curVariance <- sapply(c(1:length(NormCalcList$partitionFuncTempAll)), 
    #                       function(x) sd(NormCalcList$partitionFuncTempAll[1:x])/sqrt(x))

    lowBound <- curMean-curSTD
    upBound <- curMean+curSTD
    timeScale <- c(1:length(curMean))/length(curMean)*timeMean
    
    df <- data.frame(normC=curMean, stdEr=curSTD, MSE=MSE, NRMSE=NRMSE,
                     countsS=timeScale,
                     lowBound=lowBound,upBound=upBound)
    df$group <- t0
    df$Method <- benchmarkResults[[1]][[t0]]$method
    
    timeOfMeasurmentForBoxPlot <- 0.2 # set the fixed time
    
    timePointForBoxPlotRes <- which.min(abs(timeScale-timeOfMeasurmentForBoxPlot))
    
    BoxPlotRes <- sapply(1:lengthRes, function(x) (benchmarkResults[[x]][[t0]]$partitionFunc[timePointForBoxPlotRes]-exactVal)^2)
    BoxPlotResNRMS <- sapply(1:lengthRes, function(x) sqrt((benchmarkResults[[x]][[t0]]$partitionFunc[timePointForBoxPlotRes]-exactVal)^2)/exactVal )
    
    df2 <- data.frame(BoxPlotRes=BoxPlotRes, BoxPlotResNRMS=BoxPlotResNRMS)
    
    df2$group <- t0
    df2$Method <- benchmarkResults[[1]][[t0]]$method
    df2$timePoint <- timePointForBoxPlotRes
    df2$timePointInSec <- timeOfMeasurmentForBoxPlot

    Smethod <- c(Smethod,df$Method)
    dfAll <- rbind(dfAll, df)
    df2All <- rbind(df2All, df2)
  }
  
  # plot
  p1 <- ggplot(data=dfAll, aes(y = normC, x = countsS, group=group))+
    geom_line(aes(y = normC, x = countsS,colour=Method)) + 
    geom_ribbon(aes(ymin=lowBound, ymax=upBound, fill=Method),alpha=0.5) + 
    scale_colour_manual(values=color_list) +
    scale_fill_manual(values=color_list) +
    geom_hline(yintercept=exactVal, linetype="dashed", color = "black") +
    # guides(fill = FALSE) +
    guides(col = FALSE) +
    xlab("Elapsed Time (s)") + 
    ylab("Normalizing Constant")
  
  # plot the MSE
  p2 <- ggplot(data=dfAll, aes(y = MSE, x = countsS, group=group))+
    geom_line(aes(y = MSE, x = countsS,colour=Method)) + 
    scale_colour_manual(values=color_list) +
    # guides(fill = FALSE) +
    # guides(col = FALSE) +
    xlab("Elapsed Time (s)") + 
    ylab("MSE")
  
  # make a boxplot of the results
  p3 <- ggplot(data=df2All, aes(x=Method, y=BoxPlotRes, colour=Method)) + 
    geom_boxplot()+
    geom_jitter(position=position_jitter(0.8))+
    scale_colour_manual(values=color_list)+
    theme(legend.position = "none")+
    ylab(paste0("MSE (", timeOfMeasurmentForBoxPlot, "s)"))
  
  p4 <- ggarrange(
    p2,                # First row with line plot
    # Second row with box and dot plots
    ggarrange(p1, p3, ncol = 2, labels = c("B", "C")),
    nrow = 2,
    labels = "A"       # Label of the line plot
  )
  
  # plot the NRMSE
  p5 <- ggplot(data=dfAll, aes(y = NRMSE, x = countsS, group=group))+
    geom_line(aes(y = NRMSE, x = countsS,colour=Method)) + 
    scale_colour_manual(values=color_list) +
    # guides(fill = FALSE) +
    # guides(col = FALSE) +
    xlab("Elapsed Time (s)") + 
    ylab("NRMSE")
  
  # make a boxplot of the results
  p6 <- ggplot(data=df2All, aes(x=Method, y=BoxPlotResNRMS, colour=Method)) + 
    geom_boxplot()+
    geom_jitter(position=position_jitter(0.8))+
    scale_colour_manual(values=color_list)+
    theme(legend.position = "none")+
    ylab(paste0("NRMSE (", timeOfMeasurmentForBoxPlot, "s)"))
  
  return(list("RAWbenchmarkResults"=benchmarkResults,"results_byTime"=dfAll, "results_final"=df2All, "plot.MSE"=p2, "plot.NC"=p1, "plot.boxplotMSE"=p3, "plot.summary"=p4, "plot.NRMSE"=p5, "plot.boxplotNRMSE"=p6))
  
}


# plot.samplingProgressStudy <- function(benchmarkResults, exactValue = NULL)
# {
#   # Plot the progress
#   
#   lengthMethods <- length(benchmarkResults[[1]])
#   lengthRes <- length(benchmarkResults)
#   
#   # set colors
#   if(lengthMethods>3){
#     N_colors <- lengthMethods
#     color_list <- brewer.pal(n = N_colors, name = 'RdYlBu')
#   }else{
#     N_colors <- 3
#     color_list <- brewer.pal(n = N_colors, name = 'RdYlBu')[1:lengthMethods]
#   }
#   # without Rcolorbrewer:
#   # color_list <- c("#EE6677", "#228833", "#4477AA")
#   
#   dfAll <- NULL
#   df2All <- NULL
#   Smethod <- NULL
#   
#   # process results for plot
#   for (t0 in 1:lengthMethods){
#     # NormCalcList <- benchmarkResults[[tt0]][[t0]]
#     curMean <- rowMeans(sapply(1:lengthRes, function(x) benchmarkResults[[x]][[t0]]$partitionFunc))
#     timeMean <- mean(sapply(1:lengthRes, function(x) benchmarkResults[[x]][[t0]]$elapsed_time[[1]]))
#     curSTD <- sqrt(RowVar(sapply(1:lengthRes, function(x) benchmarkResults[[x]][[t0]]$partitionFunc)))
#     # curVariance <- sapply(c(1:length(NormCalcList$partitionFuncTempAll)), 
#     #                       function(x) sd(NormCalcList$partitionFuncTempAll[1:x])/sqrt(x))
#     lowBound <- curMean-curSTD
#     upBound <- curMean+curSTD
#     timeScale <- c(1:length(curMean))/length(curMean)*timeMean
#     
#     df <- data.frame(normC=curMean, stdEr=curSTD,
#                      countsS=timeScale,
#                      lowBound=lowBound,upBound=upBound)
#     df$group <- t0
#     df$Method <- benchmarkResults[[1]][[t0]]$method
#     
#     timePointForBoxPlotRes <- which.min(abs(timeScale-0.1))
#     
#     BoxPlotRes <- sapply(1:lengthRes, function(x) benchmarkResults[[x]][[t0]]$partitionFunc[timePointForBoxPlotRes])
#     
#     df2 <- data.frame(BoxPlotRes=BoxPlotRes)
#     
#     df2$group <- t0
#     df2$Method <- benchmarkResults[[1]][[t0]]$method
#     df2$timePoint <- timePointForBoxPlotRes
#     
#     Smethod <- c(Smethod,df$Method)
#     dfAll <- rbind(dfAll, df)
#     df2All <- rbind(df2All, df2)
#   }
#   
#   # plot
#   if(length(exactValue)==0){
#     p1 <- ggplot(data=dfAll, aes(y = normC, x = countsS, group=group))+
#       geom_line(aes(y = normC, x = countsS,colour=Method)) + 
#       geom_ribbon(aes(ymin=lowBound, ymax=upBound, fill=Method),alpha=0.5) + 
#       scale_colour_manual(values=color_list) +
#       scale_fill_manual(values=color_list) +
#       # guides(fill = FALSE) +
#       guides(col = FALSE) +
#       xlab("Elapsed Time (s)") + 
#       ylab("Normalizing Constant")
#   }else{
#     ggplot(data=dfAll, aes(y = normC, x = countsS, group=group))+
#       geom_line(aes(y = normC, x = countsS,colour=Method)) + 
#       geom_ribbon(aes(ymin=lowBound, ymax=upBound, fill=Method),alpha=0.5) + 
#       scale_colour_manual(values=color_list) +
#       scale_fill_manual(values=color_list) +
#       geom_hline(yintercept=exactValue, linetype="dashed", color = "black") +
#       # guides(fill = FALSE) +
#       guides(col = FALSE) +
#       xlab("Elapsed Time (s)") + 
#       ylab("Normalizing Constant")
#   }
#   
#   # plot
#   if(length(exactValue)==0){
#     p2 <- ggplot(data=dfAll, aes(y = stdEr, x = countsS, group=group))+
#       geom_line(aes(y = stdEr, x = countsS,colour=Method)) + 
#       scale_colour_manual(values=color_list) +
#       # guides(fill = FALSE) +
#       # guides(col = FALSE) +
#       xlab("Elapsed Time (s)") + 
#       ylab("Standard Error")
#   }else{
#     ggplot(data=dfAll, aes(y = normC, x = countsS, group=group))+
#       geom_line(aes(y = normC, x = countsS,colour=Method)) + 
#       geom_ribbon(aes(ymin=lowBound, ymax=upBound, fill=Method),alpha=0.5) + 
#       scale_colour_manual(values=color_list) +
#       scale_fill_manual(values=color_list) +
#       geom_hline(yintercept=exactValue, linetype="dashed", color = "black") +
#       # guides(fill = FALSE) +
#       guides(col = FALSE) +
#       xlab("Elapsed Time (s)") + 
#       ylab("Standard Error")
#   }
#   
#   # make a boxplot of the results
#   
#   p3 <- ggplot(df2All, aes(x=Method, y=BoxPlotRes, colour=Method)) + 
#     geom_boxplot()+
#     scale_colour_manual(values=color_list)+
#     theme(legend.position = "none")
#   
#   ggarrange(
#     p2,                # First row with line plot
#     # Second row with box and dot plots
#     ggarrange(p1, p3, ncol = 2, labels = c("B", "C")),
#     nrow = 2,
#     labels = "A"       # Label of the line plot
#   )
#   
#   # ggarrange(
#   #   p2,                # First row with line plot
#   #   # Second row with box and dot plots
#   #   ggarrange(p1, p3, ncol = 2), 
#   #   nrow = 2 
#   # ) 
#   
# }

RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

plot.BNwithObs <- function(myBN, myObs){
  # plot the Bayes net with observations in grey
  nodeCol <- rep('white',num.nodes(myBN))
  nodeCol[myObs$observed.vars] <- rep('grey',length(myObs$observed.vars))

  plot.BN(myBN, node.col = nodeCol)
}

#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics par
plot.BNwithSubGroups <- function(myBN, myObs, visualizeAll=FALSE){
  # plot BN with subgroups in color and corresponding observations in grey
  
  mySubGroups <- get.allSubGroups(myBN@dag, myObs$observed.vars) 
  
  
  N_colors <- length(mySubGroups$allChiEvidence)
  color_list <- brewer.pal(n = (N_colors+2), name = 'Greys')[2:(N_colors+1)]
  
  color_list2 <- brewer.pal(n = N_colors, name = 'RdYlBu')
  
  nodeCol <- rep('white',num.nodes(myBN))
  
  for (yy in 1:length(mySubGroups$allSubGroups)){
    
    nodeCol[mySubGroups$allSubGroups[[yy]]] <- rep(color_list2[yy],length(mySubGroups$allSubGroups[[yy]]))
  }
  
  for (yy in 1:length(mySubGroups$allChiEvidence)){
    
    nodeCol[mySubGroups$allChiEvidence[[yy]]] <- rep(color_list[yy],length(mySubGroups$allChiEvidence[[yy]]))
  }
  
  # nodeCol[myObs$observed.vars] <- rep('grey',length(myObs$observed.vars))
  
  plot.BN(myBN, node.col = nodeCol)
  
  if(visualizeAll){
    dim2 <- ceiling(sqrt(length(mySubGroups$allSubGroups)))
    dim1 <- ceiling(length(mySubGroups$allSubGroups)/dim2)
    storeMar <- par("mar")
    par(mfrow=c(dim1,dim2),mar=c(1,1,1,1))
    for (jj in 1:length(mySubGroups$allSubGroups)){
      # graphviz.plot(allSubBNs[[jj]])
      # graphviz.plot(allSubGroups[[jj]], highlight = list(nodes = allEvidenceSubGroups[[jj]], col = "black", fill = "grey"))
      # subDAG <- as.matrix(as_adjacency_matrix(induced_subgraph(
        # graph_from_adjacency_matrix(DAG,mode="directed"),c(mySubGroups$allEvidenceSubGroups[[jj]],mySubGroups$allSubGroups[[jj]]))))
      subBN <- get.subBN(myBN, c(mySubGroups$allEvidenceSubGroups[[jj]],mySubGroups$allSubGroups[[jj]]))
      plot.BN(subBN)
    }
    par(mfrow = c(1,1))
    par(mar=storeMar)
  }
}

#' Benchmark Inference Methods
#'
#' Outputs benchmark results of different inference methods
#'
#' @param N_var number of variables
#' @param N_Obs number of observations
#' @param N_nets number of Bayes nets
#' @param N_rep number of repetitions
#' @param N_samples number of samples
#' @param samplingMethods benchmarked sampling methods
#' @param DAGmethod DAG method
#' @param N_neighbours number of neighbours
#' @param nodeDim category size
#' @return benchmark results
#' @export
benchmarkMultipleNets <- function(N_var, N_Obs, N_nets, N_rep, N_samples=100, samplingMethods=c("GS","SGS", "LBP"),DAGmethod="er", N_neighbours=2, nodeDim=2){
  
  # benchmark multiple Bayes nets of same size
  
  resultTable <<- list()
  
  for (o1 in 1:N_nets){
    
    print(paste0("Benchmarking Bayes net #", o1, " of ", N_nets))
    
    start_time <- Sys.time() # measure computation time

    # Simulate a random BN and observations
    myBN <- randomBN(N_var, method=DAGmethod, N_neighbours=N_neighbours, nodeDim=nodeDim)#, N_neighbours=2) #, uniformCPTs = FALSE)
    myObs <- list("observed.vars" = sample(1:N_var, N_Obs), "observed.vals" = rep(1,N_Obs))
    
    # perform benchmark
    resultTable[[o1]] <- benchmarkStudyMain(myBN,myObs,methods=samplingMethods,N_rep=N_rep,N_samples=N_samples)
    
    end_time <- Sys.time()
    elapsed_time <- end_time - start_time
    print(elapsed_time)
    
  }
  
  return(resultTable)
  
}


# ProcessResultTable <- function(myResultTable){
#   
#   averageResults_byTime <- myResultTable[[1]]$results_byTime
#   averageResults_final <- myResultTable[[1]]$results_final
#   
#   lengthMethods <- length(myResultTable[[1]]$simulationDetails$methods)
#   
#   averageResults_byTime$MSE <- rowMeans(sapply(1:length(myResultTable), function(x) myResultTable[[x]]$results_byTime$MSE))
#   averageResults_byTime$countsS <- rowMeans(sapply(1:length(myResultTable), function(x) myResultTable[[x]]$results_byTime$countsS))
#   averageResults_final$MSE <- rowMeans(sapply(1:length(myResultTable), function(x) myResultTable[[x]]$results_final$BoxPlotRes))
#   
#   # set colors
#   if(lengthMethods>3){
#     N_colors <- lengthMethods
#     color_list <- brewer.pal(n = N_colors, name = 'RdYlBu')
#   }else{
#     N_colors <- 3
#     color_list <- brewer.pal(n = N_colors, name = 'RdYlBu')[1:lengthMethods]
#   }
#   
#   # plot the MSE
#   p2 <- ggplot(data=averageResults_byTime, aes(y = MSE, x = countsS, group=group))+
#     geom_line(aes(y = MSE, x = countsS,colour=Method)) + 
#     scale_colour_manual(values=color_list) +
#     # guides(fill = FALSE) +
#     # guides(col = FALSE) +
#     xlab("Elapsed Time (s)") + 
#     ylab("MSE") # + xlim(-0.01, 0.2)  + ylim(-0.1e-11, 0.75e-11)
#   
#   # make a boxplot of the results
#   p3 <- ggplot(data=averageResults_final, aes(x=Method, y=BoxPlotRes, colour=Method)) + 
#     geom_boxplot()+
#     scale_colour_manual(values=color_list)+
#     theme(legend.position = "none")+
#     ylab(paste0("MSE (", averageResults_final$timePointInSec, "s)"))
#   
#   p4 <- ggplot(data=averageResults_final, aes(x=Method, y=BoxPlotRes, colour=Method)) + 
#     geom_boxplot()+
#     scale_colour_manual(values=tail(color_list,2))+
#     theme(legend.position = "none")+
#     ylab(paste0("MSE (", averageResults_final$timePointInSec, "s)"))+
#     scale_x_discrete(limits=c("LBP","SGS"))
#   
#   print(p2)
#   
#   print(p3)
#   
#   print(p4)
# }


#' @importFrom utils tail
ProcessResultTableNRMSE <- function(myResultTable, returnResults=FALSE){
  
  NRMSE=countsS=group=Method=BoxPlotResNRMS=0
  averageResults_byTime <- myResultTable[[1]]$results_byTime
  averageResults_final <- myResultTable[[1]]$results_final
  
  lengthMethods <- length(myResultTable[[1]]$simulationDetails$methods)
  
  averageResults_byTime$NRMSE <- rowMeans(sapply(1:length(myResultTable), function(x) myResultTable[[x]]$results_byTime$NRMSE))
  averageResults_byTime$countsS <- rowMeans(sapply(1:length(myResultTable), function(x) myResultTable[[x]]$results_byTime$countsS))
  
  averageResults_final <- myResultTable[[1]]$results_final
  for (h1 in 2:length(myResultTable)){
    averageResults_final <- rbind(averageResults_final,myResultTable[[h1]]$results_final)
  }
  # averageResults_final$BoxPlotResNRMS <- sapply(1:length(myResultTable), function(x) myResultTable[[x]]$results_final$BoxPlotResNRMS)
  
  # set colors
  if(lengthMethods>3){
    N_colors <- lengthMethods
    color_list <- brewer.pal(n = N_colors, name = 'RdYlBu')
  }else{
    N_colors <- 3
    color_list <- brewer.pal(n = N_colors, name = 'RdYlBu')[1:lengthMethods]
  }
  
  if(returnResults==TRUE){
    return(list("averageResults_byTime"=averageResults_byTime, "averageResults_final"=averageResults_final))
  }
  
  # plot the MSE
  p2 <- ggplot(data=averageResults_byTime, aes(y = NRMSE, x = countsS, group=group))+
    geom_line(aes(y = NRMSE, x = countsS,colour=Method)) + 
    scale_colour_manual(values=color_list) +
    # guides(fill = FALSE) +
    # guides(col = FALSE) +
    xlab("Elapsed Time (s)") + 
    ylab("NRMSE") # + xlim(-0.01, 0.2)  + ylim(-0.1e-11, 0.75e-11)
  
  # make a boxplot of the results
  p3 <- ggplot(data=averageResults_final, aes(x=Method, y=BoxPlotResNRMS, colour=Method)) + 
    geom_boxplot()+
    geom_jitter(position=position_jitter(0.8))+
    scale_colour_manual(values=color_list)+
    theme(legend.position = "none")+
    ylab(paste0("NRMSE (", averageResults_final$timePointInSec, "s)"))
  
  p4 <- ggplot(data=averageResults_final, aes(x=Method, y=BoxPlotResNRMS, colour=Method)) + 
    geom_boxplot()+
    geom_jitter(position=position_jitter(0.8))+
    scale_colour_manual(values=utils::tail(color_list,2))+
    theme(legend.position = "none")+
    ylab(paste0("NRMSE (", averageResults_final$timePointInSec, "s)"))+
    scale_x_discrete(limits=c("LBP","SGS"))
  
  print(p2)
  
  print(p3)
  
  print(p4)
  
}


# pathfinderNet <- readRDS("/Users/frbayer/Downloads/pathfinder.rds")
# 
# hailfinder <- readRDS("/Users/frbayer/Downloads/hailfinder.rds")
# 
# 
# pathfinderNet
# 
# tempBN <- convert_bnlearn(hailfinder)

# convert_bnlearn <- function(BNconvert){
#   # convert Bayes net of type "bnlearn" to type "SubGroupSeparation"
#   
#   lengthBNC <- length(BNconvert)
#   
#   # create new BN
#   tempBN <- BN()
#   name(tempBN) <- "converted Bayes net"
#   num.nodes(tempBN) <- lengthBNC
#   variables(tempBN) <- names(BNconvert)
#   discreteness(tempBN) <- rep(TRUE, lengthBNC)
#   node.sizes(tempBN) <- sapply(1:lengthBNC, function(x) unname(dim(BNconvert[[x]]$prob)[1]))
#   cpts(tempBN) <- sapply(1:lengthBNC, function(x) BNconvert[[x]]$prob)
#   dag(tempBN) <- amat(BNconvert)
#   wpdag(tempBN) <- matrix(0, lengthBNC, lengthBNC)
#   
#   return(tempBN)
# }


#' Visualize benchmark results
#'
#' Outputs plots of benchmark results for different inference methods
#'
#' @param resultTable1 result table 1 of function benchmarkMultipleNets()
#' @param resultTable2 result table 2 of function benchmarkMultipleNets()
#' @param resultTable3 result table 3 of function benchmarkMultipleNets()
#' @param resultTable4 result table 4 of function benchmarkMultipleNets()
#' @param labelBP label
#' @param labelHead label head
#' @param fileName file name
#' @param width figure width
#' @param height figure height
#' @param niceBoxPlot if TRUE, only a fraction of points (150) are plotted in boxplot
#' @param retP if TRUE, a plot is returned for summary plot, otherwise no ruturn
#' @return plot of benchmark results
#' 
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices cairo_pdf
#' @importFrom utils tail
#' 
#' @export
makeAllPlots <- function(resultTable1, resultTable2, resultTable3, resultTable4, labelBP = c("1","2","3","4"), labelHead = c("1","2","3","4"), fileName = "Test", width = 7, height = 4.1, niceBoxPlot=FALSE, retP=FALSE){
  # function to create the final plots from the benchmark studies
  
  dimBN=BoxPlotResNRMS=Method=NRMSE=countsS=group=NULL
  
  ## FinalPlots:
  restuls1 <- ProcessResultTableNRMSE(resultTable1, returnResults = TRUE)
  restuls2 <- ProcessResultTableNRMSE(resultTable2, returnResults = TRUE)
  restuls3 <- ProcessResultTableNRMSE(resultTable3, returnResults = TRUE)
  restuls4 <- ProcessResultTableNRMSE(resultTable4, returnResults = TRUE)
  
  restuls1$averageResults_final$dimBN <- labelBP[1]
  restuls2$averageResults_final$dimBN <- labelBP[2]
  restuls3$averageResults_final$dimBN <- labelBP[3]
  restuls4$averageResults_final$dimBN <- labelBP[4]
  
  resAllDims <- rbind(restuls1$averageResults_final,restuls2$averageResults_final,restuls3$averageResults_final,restuls4$averageResults_final)
  
  # if FS is removed
  resAllDims2 <- resAllDims[!resAllDims$Method=='FS',]
  # resAllDims2$Method[resAllDims2$Method=="LBP (IS)"] <- "LBP"
  
  # set colors
  color_list <- brewer.pal(n = 4, name = 'RdYlBu')
  
  # change LBP name for plotting
  resAllDims2$Method[resAllDims2$Method=="LBP"] <- "LBP-IS"
  
  # # create subset of data for plotting with "geom_jitter"
  numb_points <- 150
  
  if (niceBoxPlot==TRUE){
    dftest1 = resAllDims2[resAllDims2$Method=="LBP-IS"&resAllDims2$dimBN==labelBP[1],][1:numb_points,]
    dftest2 = resAllDims2[resAllDims2$Method=="SGS"&resAllDims2$dimBN==labelBP[1],][1:numb_points,]
    dftest3 = resAllDims2[resAllDims2$Method=="GS"&resAllDims2$dimBN==labelBP[1],][1:numb_points,]
    resAllDims2_subset1 <- merge(merge(dftest1,dftest2,  all = TRUE),dftest3,  all = TRUE)
    
    dftest1 = resAllDims2[resAllDims2$Method=="LBP-IS"&resAllDims2$dimBN==labelBP[2],][1:numb_points,]
    dftest2 = resAllDims2[resAllDims2$Method=="SGS"&resAllDims2$dimBN==labelBP[2],][1:numb_points,]
    dftest3 = resAllDims2[resAllDims2$Method=="GS"&resAllDims2$dimBN==labelBP[2],][1:numb_points,]
    resAllDims2_subset2 <- merge(merge(dftest1,dftest2,  all = TRUE),dftest3,  all = TRUE)
    
    dftest1 = resAllDims2[resAllDims2$Method=="LBP-IS"&resAllDims2$dimBN==labelBP[3],][1:numb_points,]
    dftest2 = resAllDims2[resAllDims2$Method=="SGS"&resAllDims2$dimBN==labelBP[3],][1:numb_points,]
    dftest3 = resAllDims2[resAllDims2$Method=="GS"&resAllDims2$dimBN==labelBP[3],][1:numb_points,]
    resAllDims2_subset3 <- merge(merge(dftest1,dftest2,  all = TRUE),dftest3,  all = TRUE)
    
    dftest1 = resAllDims2[resAllDims2$Method=="LBP-IS"&resAllDims2$dimBN==labelBP[4],][1:numb_points,]
    dftest2 = resAllDims2[resAllDims2$Method=="SGS"&resAllDims2$dimBN==labelBP[4],][1:numb_points,]
    dftest3 = resAllDims2[resAllDims2$Method=="GS"&resAllDims2$dimBN==labelBP[4],][1:numb_points,]
    resAllDims2_subset4 <- merge(merge(dftest1,dftest2,  all = TRUE),dftest3,  all = TRUE)
    
    resAllDims2_subset <- merge(merge(merge(resAllDims2_subset1,resAllDims2_subset2, all = TRUE),resAllDims2_subset3,all = TRUE), resAllDims2_subset4,  all = TRUE)
  }
  
  dodge <- position_dodge(width = 0.4)
  
  if(niceBoxPlot==TRUE){
    # make a boxplot of the results
    p1 <- ggplot(data=resAllDims2, aes(x=factor(dimBN, levels = c(labelBP[1],labelBP[2],labelBP[3],labelBP[4])), y=BoxPlotResNRMS, fill=Method, colour=Method), alpha = 0.3) + 
      geom_boxplot(outlier.shape = NA, alpha = 0.3) + ## THIS IS FOR JITTER, otherwise geom_boxplot()+
      # geom_jitter(position=position_jitter(0.42), cex=0.4, data=resAllDims2_subset, position = dodge)+
      geom_point(pch = 21,data=resAllDims2_subset, position = position_jitterdodge(0.15),cex=0.4, alpha = 0.3)+ ## THIS IS FOR JITTER
      scale_colour_manual(values=color_list[c(1,3,4)])+
      scale_fill_manual(values=color_list[c(1,3,4)])+ # values=c("white","white","white"))+
      labs(y = "NRMSE", x="Dimensions")+
      ylim(0,1.5)+
      theme_minimal()+
      theme(legend.position = "bottom")
  }else{
    # make a boxplot of the results
    p1 <- ggplot(data=resAllDims2, aes(x=factor(dimBN, level = c(labelBP[1],labelBP[2],labelBP[3],labelBP[4])), y=BoxPlotResNRMS, fill=Method, colour=Method), alpha = 0.3) + 
      geom_boxplot(outlier.shape = NA, alpha = 0.3) + ## THIS IS FOR JITTER, otherwise geom_boxplot()+
      # geom_jitter(position=position_jitter(0.42), cex=0.4, data=resAllDims2_subset, position = dodge)+
      geom_point(pch = 21,data=resAllDims2, position = position_jitterdodge(0.15),cex=0.4, alpha = 0.3)+ ## THIS IS FOR JITTER
      scale_colour_manual(values=color_list[c(1,3,4)])+
      scale_fill_manual(values=color_list[c(1,3,4)])+ # values=c("white","white","white"))+
      labs(y = "NRMSE", x="Dimensions")+
      ylim(0,1.5)+
      theme_minimal()+
      theme(legend.position = "bottom")
  }
  
  print(p1)
  
  # if FS is removed
  plotData <- restuls1$averageResults_byTime
  plotData <- plotData[!plotData$Method=='FS',]
  plotData$Method[plotData$Method=="LBP"] <- "LBP-IS"
  
  minLims <- min(c(tail(plotData$countsS[plotData$Method=="LBP-IS"],1),tail(plotData$countsS[plotData$Method=="GS"],1), tail(plotData$countsS[plotData$Method=="SGS"],1)))
  
  p2 <- ggplot(data=plotData, aes(y = NRMSE, x = countsS, group=group))+
    geom_line(aes(y = NRMSE, x = countsS,colour=Method)) + 
    scale_colour_manual(values=color_list[c(1,3,4)]) +
    lims(x=c(0,minLims))+
    # guides(fill = FALSE) +
    # guides(col = FALSE) +
    xlab("Elapsed Time") + 
    theme_minimal()+
    theme(legend.direction = "horizontal")+
    ylab("NRMSE") # + xlim(-0.01, 0.2)  + ylim(-0.1e-11, 0.75e-11)
  
  print(p2)
  
  p2s <- ggplot(data=plotData, aes(y = NRMSE, x = countsS, group=group))+
    geom_line(aes(y = NRMSE*sqrt(countsS), x = countsS,colour=Method)) +
    scale_colour_manual(values=color_list[c(1,3,4)]) +
    lims(x=c(0,minLims))+
    # guides(fill = FALSE) +
    # guides(col = FALSE) +
    xlab("Elapsed Time") +
    theme_minimal()+
    ylab(expression(NRMSE*sqrt(t))) # + xlim(-0.01, 0.2)  + ylim(-0.1e-11, 0.75e-11)
  
  print(p2s)
  
  p2ss <- ggplot(data=plotData, aes(y = NRMSE, x = countsS, group=group))+
    geom_line(aes(y = log(NRMSE), x = log(countsS),colour=Method)) +
    scale_colour_manual(values=color_list[c(1,3,4)]) +
    # lims(x=c(0,minLims))+
    # guides(fill = FALSE) +
    # guides(col = FALSE) +
    # scale_x_continuous(trans='log2') +
    # scale_y_continuous(trans='log2') +
    xlab("log(Elapsed Time)") +
    theme_minimal()+
    ylab("log(NRMSE)") # + xlim(-0.01, 0.2)  + ylim(-0.1e-11, 0.75e-11)
  
  print(p2ss)
  
  # if FS is removed
  plotData <- restuls2$averageResults_byTime
  plotData <- plotData[!plotData$Method=='FS',]
  plotData$Method[plotData$Method=="LBP"] <- "LBP-IS"
  
  minLims <- min(c(tail(plotData$countsS[plotData$Method=="LBP-IS"],1),tail(plotData$countsS[plotData$Method=="GS"],1), tail(plotData$countsS[plotData$Method=="SGS"],1)))
  
  p3 <- ggplot(data=plotData, aes(y = NRMSE, x = countsS, group=group))+
    geom_line(aes(y = NRMSE, x = countsS,colour=Method)) + 
    scale_colour_manual(values=color_list[c(1,3,4)]) +
    lims(x=c(0,minLims))+
    # guides(fill = FALSE) +
    # guides(col = FALSE) +
    xlab("Elapsed Time") + 
    theme_minimal()+
    ylab("NRMSE") # + xlim(-0.01, 0.2)  + ylim(-0.1e-11, 0.75e-11)
  
  print(p3)
  
  p3s <- ggplot(data=plotData, aes(y = NRMSE*sqrt(countsS), x = countsS, group=group))+
    geom_line(aes(y = NRMSE*sqrt(countsS), x = countsS,colour=Method)) +
    scale_colour_manual(values=color_list[c(1,3,4)]) +
    lims(x=c(0,minLims))+
    # guides(fill = FALSE) +
    # guides(col = FALSE) +
    xlab("Elapsed Time") +
    theme_minimal()+
    ylab(expression(NRMSE*sqrt(t))) # + xlim(-0.01, 0.2)  + ylim(-0.1e-11, 0.75e-11)
  
  print(p3s)
  
  p3ss <- ggplot(data=plotData, aes(y = NRMSE, x = countsS, group=group))+
    geom_line(aes(y = log(NRMSE), x = log(countsS),colour=Method)) +
    scale_colour_manual(values=color_list[c(1,3,4)]) +
    # lims(x=c(0,minLims))+
    # guides(fill = FALSE) +
    # guides(col = FALSE) +
    # scale_x_continuous(trans='log2') +
    # scale_y_continuous(trans='log2') +
    xlab("log(Elapsed Time)") +
    theme_minimal()+
    ylab("log(NRMSE)") # + xlim(-0.01, 0.2)  + ylim(-0.1e-11, 0.75e-11)
  
  print(p3ss)
  
  # if FS is removed
  plotData <- restuls3$averageResults_byTime
  plotData <- plotData[!plotData$Method=='FS',]
  plotData$Method[plotData$Method=="LBP"] <- "LBP-IS"
  
  minLims <- min(c(tail(plotData$countsS[plotData$Method=="LBP-IS"],1),tail(plotData$countsS[plotData$Method=="GS"],1), tail(plotData$countsS[plotData$Method=="SGS"],1)))
  
  p4 <- ggplot(data=plotData, aes(y = NRMSE, x = countsS, group=group))+
    geom_line(aes(y = NRMSE, x = countsS,colour=Method)) + 
    scale_colour_manual(values=color_list[c(1,3,4)]) +
    lims(x=c(0,minLims))+
    # guides(fill = FALSE) +
    # guides(col = FALSE) +
    xlab("Elapsed Time") + 
    theme_minimal()+
    ylab("NRMSE") # + xlim(-0.01, 0.2)  + ylim(-0.1e-11, 0.75e-11)
  
  print(p4)
  
  p4s <- ggplot(data=plotData, aes(y = NRMSE*sqrt(countsS), x = countsS, group=group))+
    geom_line(aes(y = NRMSE*sqrt(countsS), x = countsS,colour=Method)) +
    scale_colour_manual(values=color_list[c(1,3,4)]) +
    lims(x=c(0,minLims))+
    # guides(fill = FALSE) +
    # guides(col = FALSE) +
    xlab("Elapsed Time") +
    theme_minimal()+
    ylab(expression(NRMSE*sqrt(t))) # + xlim(-0.01, 0.2)  + ylim(-0.1e-11, 0.75e-11)
  
  print(p4s)
  
  p4ss <- ggplot(data=plotData, aes(y = NRMSE, x = countsS, group=group))+
    geom_line(aes(y = log(NRMSE), x = log(countsS),colour=Method)) +
    scale_colour_manual(values=color_list[c(1,3,4)]) +
    # lims(x=c(0,minLims))+
    # guides(fill = FALSE) +
    # guides(col = FALSE) +
    # scale_x_continuous(trans='log2') +
    # scale_y_continuous(trans='log2') +
    xlab("log(Elapsed Time)") +
    theme_minimal()+
    ylab("log(NRMSE)") # + xlim(-0.01, 0.2)  + ylim(-0.1e-11, 0.75e-11)
  
  print(p4ss)
  
  # if FS is removed
  plotData <- restuls4$averageResults_byTime
  plotData <- plotData[!plotData$Method=='FS',]
  plotData$Method[plotData$Method=="LBP"] <- "LBP-IS"
  
  minLims <- min(c(tail(plotData$countsS[plotData$Method=="LBP-IS"],1),tail(plotData$countsS[plotData$Method=="GS"],1), tail(plotData$countsS[plotData$Method=="SGS"],1)))
  
  p5 <- ggplot(data=plotData, aes(y = NRMSE, x = countsS, group=group))+
    geom_line(aes(y = NRMSE, x = countsS,colour=Method)) + 
    scale_colour_manual(values=color_list[c(1,3,4)]) +
    lims(x=c(0,minLims))+
    # guides(fill = FALSE) +
    # guides(col = FALSE) +
    xlab("Elapsed Time") + 
    theme_minimal()+
    ylab("NRMSE") # + xlim(-0.01, 0.2)  + ylim(-0.1e-11, 0.75e-11)
  
  print(p5)
  
  p5s <- ggplot(data=plotData, aes(y = NRMSE*sqrt(countsS), x = countsS, group=group))+
    geom_line(aes(y = NRMSE*sqrt(countsS), x = countsS,colour=Method)) +
    scale_colour_manual(values=color_list[c(1,3,4)]) +
    lims(x=c(0,minLims))+
    # guides(fill = FALSE) +
    # guides(col = FALSE) +
    xlab("Elapsed Time") +
    theme_minimal()+
    ylab(expression(NRMSE*sqrt(t))) # + xlim(-0.01, 0.2)  + ylim(-0.1e-11, 0.75e-11)
  
  print(p5s)
  
  p5ss <- ggplot(data=plotData, aes(y = NRMSE, x = countsS, group=group))+
    geom_line(aes(y = log(NRMSE), x = log(countsS),colour=Method)) +
    scale_colour_manual(values=color_list[c(1,3,4)]) +
    # lims(x=c(0,minLims))+
    # guides(fill = FALSE) +
    # guides(col = FALSE) +
    # scale_x_continuous(trans='log2') +
    # scale_y_continuous(trans='log2') +
    xlab("log(Elapsed Time)") +
    theme_minimal()+
    ylab("log(NRMSE)") # + xlim(-0.01, 0.2)  + ylim(-0.1e-11, 0.75e-11)
  
  print(p5ss)
  
  # store figures
  
  # plot 1
  # png("~/Desktop/Figure1.png", width = 10, height = 10, units = 'cm', res = 300)
  # grid.arrange(p1)
  # dev.off()
  cairo_pdf(paste0("~/Desktop/Figure", fileName,"BoxPlot.pdf"), width = width, height = height)
  grid.arrange(p1)
  dev.off()
  
  # tikz(paste0("~/Desktop/Figure", fileName,"BoxPlot.tex"), width = width, height = height)
  # grid.arrange(p1)
  # dev.off()
  
  p_grid <- grid.arrange(p2+ggtitle(labelHead[1])+theme(legend.position = "none"),p3+ggtitle(labelHead[2])+theme(legend.position = "none"),
                         p4+ggtitle(labelHead[3])+theme(legend.position = "none"),p5+ggtitle(labelHead[4])+theme(legend.position = "none"))
  legend <- cowplot::get_legend(p2)
  
  # p6 <- cowplot::plot_grid(p_grid, legend,ncol = 2, rel_heights = c(1, 1), rel_widths = c(1,0.15))
  p6 <- cowplot::plot_grid(p_grid, legend, nrow = 2, rel_heights = c(1, 0.1), rel_widths = c(1,1))
  
  # plot 6
  # png("~/Desktop/Figure6.png", width = 10, height = 10, units = 'cm', res = 300)
  # grid.arrange(p6)
  # dev.off()
  cairo_pdf(paste0("~/Desktop/Figure", fileName,".pdf"), width = width, height = height)
  grid.arrange(p6)
  dev.off()
  
  # tikz(paste0("~/Desktop/Figure", fileName,".tex"), width = width, height = height)
  # grid.arrange(p6)
  # dev.off()
  
  p_grids <- grid.arrange(p2s+ggtitle(labelHead[1])+theme(legend.position = "none"),p3s+ggtitle(labelHead[2])+theme(legend.position = "none"),
                          p4s+ggtitle(labelHead[3])+theme(legend.position = "none"),p5s+ggtitle(labelHead[4])+theme(legend.position = "none"))
  legend <- cowplot::get_legend(p2)
  
  # p6s <- cowplot::plot_grid(p_grids, legend,ncol = 2, rel_heights = c(1, 1), rel_widths = c(1,0.15))
  p6s <- cowplot::plot_grid(p_grids, legend,nrow = 2, rel_heights = c(1, 0.1), rel_widths = c(1,1))
  
  # plot 6
  # png("~/Desktop/Figure6.png", width = 10, height = 10, units = 'cm', res = 300)
  # grid.arrange(p6s)
  # dev.off()
  cairo_pdf(paste0("~/Desktop/Figure", fileName,"_time.pdf"), width = width, height = height)
  grid.arrange(p6s)
  dev.off()
  
  p_grids <- grid.arrange(p2ss+ggtitle(labelHead[1])+theme(legend.position = "none"),p3ss+ggtitle(labelHead[2])+theme(legend.position = "none"),
                          p4ss+ggtitle(labelHead[3])+theme(legend.position = "none"),p5ss+ggtitle(labelHead[4])+theme(legend.position = "none"))
  # legend2 <- cowplot::get_legend(p2s)
  
  # p6ss <- cowplot::plot_grid(p_grids, legend2,ncol = 2, rel_heights = c(1, 1), rel_widths = c(1,0.15))
  p6ss <- cowplot::plot_grid(p_grids, legend,nrow = 2, rel_heights = c(1, 0.1), rel_widths = c(1,1))
  
  # plot 6
  # png("~/Desktop/Figure6.png", width = 10, height = 10, units = 'cm', res = 300)
  # grid.arrange(p6s)
  # dev.off()
  cairo_pdf(paste0("~/Desktop/Figure", fileName,"_log.pdf"), width = width, height = height)
  grid.arrange(p6ss)
  dev.off()
  
  if(retP==TRUE) return(p1)
}

