#' @rdname loopy_belief.propagation
#' @aliases loopy_belief.propagation,InferenceEngine
setMethod("loopy_belief.propagation",
          c("InferenceEngine"),
          function(ie, observations = NULL, return.potentials = FALSE){
            {
              ###############################
              # moved inside in order to eliminate a NOTE in R CMD check
              proc.order <- function(node, from, adj)
              {
                # Recursive method to compute order of the message passing in the upward step.
                #
                # node : current node
                # from : (local) root
                # adj  : adjacency matrix
                neighbours <- setdiff(which(adj[node,] > 0), from)
                
                if (length(neighbours) > 0)
                {
                  for (n in neighbours) {
                    proc.order(n, node, adj)
                    parents.list <<- c(parents.list, node)
                  }
                }
                
                process.order <<- c(process.order, node)
              }
              
              ##############################
              
              # if (missing(net))
              # {
              net <- bn(ie)
              # }
              
              if (missing(observations))
              {
                obs <- observations(ie)
                observed.vars <- obs$observed.vars
                observed.vals <- obs$observed.vals
              }
              else
              {
                observed.vars <- observations[[1]]
                observed.vals <- observations[[2]]
                obs <- unique.observations(observed.vars, observed.vals)
                observed.vars <- obs$observed.vars
                observed.vals <- obs$observed.vals
                observations(ie) <- list(observed.vars, observed.vals)
              }
              
              num.nodes  <- num.nodes(net)
              num.cliqs  <- num.nodes(net) ### adjusted for LBP
              
              # cliques contains the variables that compose each clique
              cliques <- getLoopyGroupCliques(net@dag)
              
              # cliques <- lapply(1:num.cliqs, function(x) names(dimnames(engine@bn@cpts[[x]]))) ### adjusted for LBP
              # cliques <- lapply(cliques, function(xxx) sapply(xxx, function(x) unname(which(sapply(engine@bn@variables, function(y) x %in% y))) ) ) ### adjusted for LBP
              # deleteNodes <- c()
              # for(ii in 1:num.cliqs){
              #   if (length(cliques[[ii]])==1){
              #     deleteNodes <- c(deleteNodes,ii)
              #   }
              # }
              # cliques[deleteNodes] <- NULL
              
              ctree      <- dag(net) ### adjusted for LBP
              
              cpts       <- cpts(net)
              
              variables  <- variables(net)
              
              node.sizes <- node.sizes(net)
              
              dim.vars   <- lapply(1:num.nodes,
                                   function(x)
                                     #as.list(
                                     match(
                                       c(unlist(
                                         names(dimnames(cpts[[x]])), F, F
                                       )),
                                       variables
                                     )
                                   #)
              )
              
              quantiles <- quantiles(net)
              
              is.discrete <- discreteness(net)
              
              # discretize continuous observed variables
              if (inherits(observed.vars, "character")) # hope that the user gives coherent input...
                observed.vars <- sapply(observed.vars, function(x) which(x == variables))
              
              observed.vals.disc <- c()
              if (length(observed.vars) > 0) {
                for (i in 1:length(observed.vars)) {
                  ovr <- observed.vars[i]
                  ovl <- observed.vals[i]
                  if (is.discrete[ovr]) {
                    observed.vals.disc <- c(observed.vals.disc, ovl)
                  } else {
                    q <- quantiles[[ovr]]
                    cv <- cut(ovl, q, labels=FALSE, include.lowest=TRUE)
                    if (is.na(cv)) {
                      if (ovl <= q[1]) cv <- 1
                      else cv <- node.sizes[ovr]
                    }
                    observed.vals.disc <- c(observed.vals.disc, cv)
                  }
                }
              }
              
              # the following lines are necessary due to a sub-group specific issue:
              # if a subgroup is taken, sometimes a variables parent might be cut off
              # this can only happen if nodes that are observed 
              # here, we adjust the respective CPTs accordingly
              repeat_dimVars <- FALSE
              for (u0 in 1:num.nodes){
                if(!is.na(match(NA,dim.vars[[u0]]))){
                  repeat_dimVars <- TRUE
                  tempIntOut <- which(dim.vars[[u0]] %in% NA)
                  
                  pot <- cpts[[u0]]
                  dms <- names(dimnames(cpts[[u0]]))
                  tempDimNames <- names(dimnames(cpts[[u0]]))
                  tempDimms <- tempDimNames[tempIntOut]
                  for (u1 in 1:length(tempIntOut)){
                    
                    out <- marginalize(pot, dms, tempDimms[u1])
                    
                    pot <- out$potential
                    dms <- out$vars
                    
                  }
                  
                  dms <- match(dms, variables)
                  
                  pot <- as.array(pot)
                  dmns <- list(NULL)
                  
                  for (i in 1:length(dms))
                  {
                    dmns[[i]] <- c(1:node.sizes[dms[i]])
                  }
                  dimnames(pot)        <- dmns
                  names(dimnames(pot)) <- as.list(variables[dms])
                  cpts[[u0]] <- pot
                  
                }
              }
              if(repeat_dimVars == TRUE){
                dim.vars   <- lapply(1:num.nodes,
                                     function(x)
                                       #as.list(
                                       match(
                                         c(unlist(
                                           names(dimnames(cpts[[x]])), F, F
                                         )),
                                         variables
                                       )
                                     #)
                )
              }
              
              
              # potentials is a list containing the probability tables of each clique
              potentials <- as.list(rep(as.list(c(1)), num.cliqs))
              
              # dimensions.contained contains the variables that effectively compose the cpt
              # currently contained in each node of the clique tree.
              # After last round it will match corresponding clique.
              dimensions.contained <- lapply(1:num.cliqs, function(x) as.list(c(NULL)))
              
              
              # choose as root (one among) the clique(s) whose connected edges have the highest overall sum
              
              # root <- which.max(rowSums(ctree)) ### not needed for LBP
              
              # Assign factors to a cluster graph
              # Construct initial potentials:
              # initial potentials are conditional or joint probability tables, depending on the initial BN
              # and how we assign each CPT to the cliques.
              #
              # The probabilities are stored as multidimensional arrays, with the convention that
              # the variables in the tables are in alphanumerical order.
              # E.g.: the network A -> C <- D -> B, whose junction tree has two nodes (and tables) ACD and BD.
              #
              # We construct a table for a clique this way:
              # - start with an empty table (NULL)
              # - whenever a {C,J}PT is assigned to that clique, its variables are ordered
              # - then, we control if the variables of the table we're inserting are already present in the clique table:
              #   - if no, we can multiply the two tables, ensuring the variables are ordered after that
              #   - if yes, we multiply the two tables
              # If the clique is currently empty, just add the cpt.
              # We have, however, to maintain the order of the variables in the probability table.
              
              
              # for (cpt in 1:num.nodes)
              # {
              #   # find a suitable clique for the CPT of node `cpt` (a clique that contains node `cpt` and all of its parents)
              #   #                 target.clique <- which.min(lapply(1:num.cliqs,
              #   #                                                   function(x){
              #   #                                                     length(
              #   #                                                       which(unlist(
              #   #                                                         is.element(
              #   #                                                           c(unlist(dim.vars[[cpt]])),
              #   #                                                           c(cliques[[x]])
              #   #                                                         )
              #   #                                                       ) == FALSE) == 0)
              #   #                                                   }
              #   #                 ))
              #   
              #   # TODO: PROFILING: is there anything better?
              #   target.clique <- which.min(lapply(1:num.cliqs,
              #                                     function(x){
              #                                       length(which(!is.na(match(
              #                                         c(unlist(dim.vars[[cpt]])),
              #                                         c(cliques[[x]])
              #                                       )) == FALSE))
              #                                     }))
              #   
              #   # get the variables currently contained in the selected clique
              #   ds <- unlist(dimensions.contained[[target.clique]], F, F)
              #   
              #   if (length(ds) == 0)
              #   {
              #     # if current clique is empty, just insert the cpt
              #     out <- sort.dimensions(cpts[[cpt]], dim.vars[[cpt]])
              #     potentials[[target.clique]]           <- out$potential
              #     dimensions.contained[[target.clique]] <- out$vars
              #   }
              #   else
              #   {
              #     # multiply current prob. table for the already inserted prob. table
              #     out <- mult(potentials[[target.clique]],
              #                 dimensions.contained[[target.clique]],
              #                 cpts[[cpt]],
              #                 dim.vars[[cpt]],
              #                 node.sizes)
              #     potentials[[target.clique]]           <- out$potential
              #     dimensions.contained[[target.clique]] <- out$vars
              #   }
              #   
              # }
              
              for (cpt in 1:num.cliqs) ### replaces code above for LBP 
              {
                # print("===============")
                # print(dim.vars[[cpt]])
                # print(cpts[[cpt]])
                out <- sort.dimensions(cpts[[cpt]], dim.vars[[cpt]])
                potentials[[cpt]] <- out$potential
                dimensions.contained[[cpt]] <- out$vars
                # print("out")
                # print(out)
              }
              
              # INCORPORATE EVIDENCE
              # If there are any observed variables, insert the knowledge.
              # Each observation is inserted by setting to zero all of the combinations that
              # do not match the observation. This is done this way:
              # - find a clique that contains the observed variable
              # - create a probability table for the observed variable, with only one non-zero entry,
              #   the one corresponding to the observed value
              # - multiply the table of the clique for the newly created table
              # - normalize after belief propagation
              
              if (length(observed.vars) > 0)
              {
                # observed.vars <- c(unlist(observed.vars, F, F))
                
                if (inherits(observed.vars, "character")) # hope that the user gives coherent input...
                  observed.vars <- sapply(observed.vars, function(x) which(x == variables))
                
                for (var in 1:length(observed.vars))
                {
                  obs.var <- observed.vars[var]
                  obs.val <- observed.vals.disc[var]
                  
                  if (obs.val <= 0 || obs.val > node.sizes[obs.var])
                  {
                    message(cat("Variable", obs.var, "cannot take value", obs.val, ", skipping..."))
                  }
                  else
                  {
                    # look for MULTIPLE clique containing the variable
                    target.cliques <- which(sapply(lapply(1:num.cliqs,
                                                          function(x) {
                                                            which(is.element(
                                                              unlist(dimensions.contained[[x]]),
                                                              obs.var
                                                            ) == TRUE)
                                                          }
                    ), function(x) !is.integer0(x)
                    )
                    )
                    
                    ### adjusted for LBP
                    for(kk in 1:length(target.cliques))
                    {
                      target.clique <- target.cliques[kk]
                      tmp <- rep(0, node.sizes[obs.var])
                      tmp[obs.val] <- 1
                      out <- mult(potentials[[target.clique]],
                                  dimensions.contained[[target.clique]],
                                  tmp,
                                  obs.var,
                                  node.sizes)
                      potentials[[target.clique]]           <- out$potential
                      dimensions.contained[[target.clique]] <- out$vars
                    }
                  }
                }
              }
              
              
              # COMPUTE PROCESSING ORDER
              # (not important for LBP, though it should not contain nodes
              # witout parents)
              process.order <- c()
              parents.list  <- list() 
              hh <- 1
              topolOrder <- as.vector(topo_sort(graph_from_adjacency_matrix(dag(net), "directed"), mode = c("in")))
              for (ll in 1:num.nodes){
                curnode <- topolOrder[ll]
                curparent <- get.allParents(dag(net), curnode)
                if(!is.integer0(curparent)){
                  parents.list[[hh]] <- curparent
                  process.order <- c(process.order,curnode)
                  hh <- hh+1
                }
              }
              # proc.order(root, c(), ctree)
              
              #print("PROCESS ORDER")
              #print(process.order)
              #print("PARENTS LIST")
              #print(parents.list)
              #print("network")
              #print(net)
              
              # N_lbp <- 3
              
              # for (aa in 1:N_lbp){
                
                # check whether message passing can be performed
                if (length(process.order) > 1) {
                  
                  # MESSAGE PASSING FROM LEAVES TO ROOT
                  
                  # msg.pots contains the prob.tables for the messages,
                  # while msg.vars contains the corresponding variables
                  msg.pots <- NULL
                  msg.vars <- NULL
                  for (ii0 in 1:length(process.order))
                  {
                    msg.pots[[ii0]] <- as.list(c(NULL))
                    msg.vars[[ii0]] <- as.list(c(NULL))
                    
                    # msg.pots[[ii0]][[1]] <- as.list(c(NULL)) 
                    # msg.vars[[ii0]][[1]] <- as.list(c(NULL))
                    
                    for (l1 in 1:length(parents.list[[ii0]])){
                      msg.pots[[ii0]][[l1]] <- as.list(c(NULL))
                      msg.vars[[ii0]][[l1]] <- as.list(c(NULL))
                    }
                    
                  }
                  
                  # For each clique (excluding the root) compute the message by marginalizing
                  # the variables not in the separator, then store the message and multiply it
                  # for the cpt contained in the neighbour clique, overwriting the corresponding
                  # potential and associated variables.
                  for (clique in 1:(length(process.order)))
                  {
                    # print("clique")
                    # print(clique)
                    # print("parents.list[clique]")
                    # print(parents.list[clique])
                    # print("cliques[[parents.list[clique]]]")
                    
                    for (jj in 1:length(parents.list[[clique]]))
                    {
                      curparent <- parents.list[[clique]][jj]
                      
                      # print(curparent)
                      
                      out <- compute.message(potentials[[process.order[clique]]],
                                             dimensions.contained[[process.order[clique]]],
                                             cliques[[process.order[clique]]],
                                             cliques[[curparent]],
                                             node.sizes)
                      msg.pots[[clique]][[jj]] <- out$potential
                      msg.vars[[clique]][[jj]] <- out$vars
                      
                      bk <- potentials[[curparent]]
                      bkd <- dimensions.contained[[curparent]]
                      out <- mult(potentials[[curparent]],
                                  dimensions.contained[[curparent]],
                                  msg.pots[[clique]][[jj]],
                                  msg.vars[[clique]][[jj]],
                                  node.sizes)
                      potentials[[curparent]]           <- out$potential
                      dimensions.contained[[curparent]] <- out$vars
                    }
                  }
                  
                  # Upward step is thus completed. Now go backward from root to leaves.
                  # This step is done by taking the CPT of the root node and dividing it (for each child)
                  # by the message received from the corresponding child, then marginalize the variables
                  # not in the separator and pass the new nessage to the child. As the messages computed
                  # in the upward step are not needed anymore after the division, they can be overwritten.
                  # Then multiply the cpt of the child for the message computed, and iterate by treating
                  # each (internal) node as root.
                  
                  # In contrast to junction tree alogrithm, the messages computed in the upward step are 
                  # still needed after division and therefore "msg.potsOld" is introduced to prevent 
                  # overwriting
                  
                  msg.potsOld <- msg.pots
                  msg.varsOld <- msg.vars
                  
                  for (clique in (length(process.order)):1)
                  {                
                    # print("clique")
                    # print(clique)
                    # print("parents.list[clique]")
                    # print(parents.list[clique])
                    # print("cliques[[parents.list[clique]]]")
                    
                    for (jj in 1:length(parents.list[[clique]]))
                    {
                      curparent <- parents.list[[clique]][jj]
                      # 
                      # print(curparent)
                      # print(process.order[clique])
                      # print("-------------------")
                      
                      out <- divide(potentials[[curparent]],
                                    dimensions.contained[[curparent]],
                                    msg.potsOld[[clique]][[jj]],
                                    msg.varsOld[[clique]][[jj]],
                                    node.sizes)
                      msg.pots[[clique]][[jj]] <- out$potential
                      msg.vars[[clique]][[jj]] <- out$vars
                      
                      out  <- compute.message(msg.pots[[clique]][[jj]],
                                              msg.vars[[clique]][[jj]],
                                              cliques[[curparent]],
                                              cliques[[process.order[clique]]],
                                              node.sizes
                      )
                      msg.pots[[clique]][[jj]] <- out$potential
                      msg.vars[[clique]][[jj]] <- out$vars
                      
                      out <- mult(potentials[[process.order[clique]]],
                                  dimensions.contained[[process.order[clique]]],
                                  msg.pots[[clique]][[jj]],
                                  msg.vars[[clique]][[jj]],
                                  node.sizes)
                      potentials[[process.order[clique]]] <- out$potential
                      dimensions.contained[[process.order[clique]]] <- out$vars
                    }
                  }
                  
                  # Finally, normalize and add dimension names and return the potentials computed (will be all JPTs).
                  for (x in 1:num.cliqs) {
                    s <- sum(potentials[[x]])
                    potentials[[x]] <- potentials[[x]] / s
                    dmns <- list(NULL)
                    for (i in 1:length(dimensions.contained[[x]]))
                    {
                      dmns[[i]] <- c(1:node.sizes[dimensions.contained[[x]][[i]]])
                    }
                    dimnames(potentials[[x]])        <- dmns
                    names(dimnames(potentials[[x]])) <- as.list(variables[c(unlist(dimensions.contained[[x]]))])
                  }
                  
                } # end if (length(process.order) > 1)
                
              # } # end for (aa in 1:N_lbp)
              
              if (return.potentials)
                return(potentials)
              
              ###################
              # Now create new BN with updated beliefs
              ###################
              #
              # To buil new conditional probability tables for the original network,
              # starting from the updated beliefs, we proceed this way:
              # for each node of the BN:
              # - if it is an observed variable, construct a prob.table containing only
              #   one variable (the one of the node) with one only non-zero (in fact, 1)
              #   entry, the one corresponding to the observed value
              # - if it is a non-observed variable, find a clique that contains the variables
              #   of the original cpt (it must exist, because of the moralization - we have
              #   already used this property), sum out the possible other variables introduced
              #   by the triangulation, and divide the JPT by the JPT of the parent nodes
              #   (e.g.: if we start from P(ABCD) and want P(C|A,B), we sum out D, then we
              #   obtain P(AB) by summing out C, and then we compute P(ABC)/P(AB) = P(C|A,B)).
              #   Easy peasy.
              
              
              
              nbn <- BN()
              name(nbn)         <- name(net)
              num.nodes(nbn)    <- num.nodes(net)
              variables(nbn)    <- variables(net)
              node.sizes(nbn)   <- node.sizes(net)
              discreteness(nbn) <- discreteness(net)
              dag(nbn)          <- dag(net)
              wpdag(nbn)        <- wpdag(net)
              scoring.func(nbn) <- scoring.func(net)
              struct.algo(nbn)  <- struct.algo(net)
              quantiles(nbn)    <- quantiles(net)
              
              ncpts <- NULL # lapply(1:num.nodes, function(x) as.list(c(NULL)))
              
              for (node in 1:num.nodes)
              {
                # faster, but result does not change. While debugging, better keep this out...
                mpos <- match(node, observed.vars)
                if (!is.na(mpos))
                  # works also when observed.vars == c()
                  # in that case, the `else` branch will be the chosen one for every variable
                {
                  ncpts[[node]] <- array(rep(0, node.sizes[observed.vars[mpos]]),
                                         c(node.sizes[observed.vars[mpos]]))
                  ncpts[[node]][observed.vals.disc[mpos]] <- 1
                  dimnames(ncpts[[node]]) <- list(c(1:node.sizes[observed.vars[mpos]]))
                  names(dimnames(ncpts[[node]])) <- as.list(variables[observed.vars[mpos]])
                }
                else
                {
                  # dnode <- unlist(dim.vars[[node]], F, F)
                  
                  #                   target.clique <- which.min(lapply(1:num.cliqs,
                  #                                                     function(x){
                  #                                                       length(
                  #                                                         which(c(
                  #                                                           is.element(
                  #                                                             dnode,
                  #                                                             unlist(dimensions.contained[[x]])
                  #                                                           )
                  #                                                         ) == FALSE) == 0)
                  #                                                     }
                  #                   ))
                  
                  target.clique <- node
                  # target.clique <- which.min(lapply(1:num.cliqs,
                  #                                   function(x){
                  #                                     length(which(!is.na(match(
                  #                                       dnode,
                  #                                       unlist(dimensions.contained[[x]])
                  #                                     )) == FALSE))
                  #                                   }))
                  
                  pot <- potentials[[target.clique]]
                  dms <- c(unlist(dimensions.contained[[target.clique]]))
                  vs  <- c(unlist(dim.vars[[node]]))
                  
                  #                   others <- setdiff(dms,vs)
                  #                   for (var in others)
                  #                   {
                  #                     out <- marginalize(pot, dms, var)
                  #                     pot <- out$potential
                  #                     dms <- out$vars
                  #                   }
                  remaining <- match(vs, dms)
                  dms <- dms[remaining]
                  pot <- apply(pot, remaining, sum)
                  pot <- pot / sum(pot)
                  
                  #                   cat(node, " ", dms,"\n")
                  #                   print(pot)
                  
                  if (length(dms) > 1)
                  {
                    pot.bak <- pot
                    dms.bak <- dms
                    # readLines(file("stdin"),1)
                    #                     out <- marginalize(pot, dms, node)
                    #                     pot <- out$potential
                    #                     dms <- out$vars
                    remaining <- (1:length(dms))[-which(dms == node)]
                    dms <- dms[remaining]
                    pot <- apply(pot, remaining, sum)
                    pot <- pot / sum(pot)
                    out <- divide(pot.bak,
                                  dms.bak,
                                  as.array(pot),
                                  dms, # out$vars,
                                  node.sizes)
                    
                    pot <- out$potential
                    #pot <- pot / sum(pot)
                    dms <- out$vars
                  }
                  pot <- as.array(pot)
                  dmns <- list(NULL)
                  for (i in 1:length(dms))
                  {
                    dmns[[i]] <- c(1:node.sizes[dms[i]])
                  }
                  dimnames(pot)        <- dmns
                  names(dimnames(pot)) <- as.list(variables[dms])
                  ncpts[[node]] <- pot
                  #                   print(pot)
                  #                   readLines(file("stdin"),1)
                }
                
              }
              
              cpts(nbn) <- ncpts
              
              updated.bn(ie) <- nbn
              jpts(ie) <- potentials
              
              # print(cliques)
              
              return(ie)
            }
          })

is.integer0 <- function(x)
{
  # check if variable is equal to integer(0)
  # added for LBP
  is.integer(x) && length(x) == 0L
}

get.allParents <- function(DAG, node)
{
  # get all parents of a node for DAG
  which(DAG[,node] > 0)
}

get.allChildren <- function(DAG, node)
{
  # get all children of a node for DAG
  which(DAG[node,] > 0)
}



