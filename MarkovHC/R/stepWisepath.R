#' @export
# extract stepwise path points and transition points
stepWisepath = function(MarkovObject=None,
                        MarkovLevel=None,
                        stepBasin=None){
  #path point
  MarkovHCpathPoint <- c()
  transitionPoint <- list()
  for(i in 1:(length(stepBasin)-1)){
    MarkovHCpath <- findTransitionPath(MarkovObject = MarkovObject,
                                       level = MarkovLevel,
                                       basinA = stepBasin[i],
                                       basinB = stepBasin[i+1])

    for(j in MarkovHCpath[[2]]){
      MarkovHCpathPoint <- c(MarkovHCpathPoint, MarkovObject$hierarchicalStructure[[1]]$basinPoints[[j]])
    }

    #transition point
    for(k in 1:length(MarkovHCpath[[2]])){
      if((MarkovHCpath[[2]][k] %in% MarkovObject$hierarchicalStructure[[MarkovLevel]]$graphvertex_basins[[stepBasin[i]]])&(!(MarkovHCpath[[2]][k+1] %in% MarkovObject$hierarchicalStructure[[MarkovLevel]]$graphvertex_basins[[stepBasin[i]]]))){
        transitionPointtemp <- c(MarkovObject$hierarchicalStructure[[1]]$basinPoints[[MarkovHCpath[[2]][k]]],
                                 MarkovObject$hierarchicalStructure[[1]]$basinPoints[[MarkovHCpath[[2]][k+1]]])
      }
    }

    transitionPoint[[length(transitionPoint)+1]] <- transitionPointtemp
    names(transitionPoint)[length(transitionPoint)] <- paste('cp',stepBasin[i],stepBasin[i+1],sep = '')
  }
  return(list(MarkovHCpathPoint, transitionPoint))
}
