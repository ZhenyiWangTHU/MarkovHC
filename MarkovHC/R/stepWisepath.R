#' Extract stepwise path points and transition points
#'
#' Function \code{stepWisepath} The function finds the samples on the
#' stepwise transition path.
#' @param MarkovObject The output of the function, \code{MarkovHC}.
#' @param MarkovLevel An integer value indicates the customized level.
#' @param stepBasin An integer vector indicates the customized basins on the transition path.
#' @details This function finds the samples on the transition path along the customized basins
#' in stepBasin parameter.
#' @return  This function returns a list consists of two components:
#' MarkovHCpathPoint is a vector contains the indices of samples along the transition path.
#'
#' transitionPoint is a list contains the indices of samples identified as critical points
#' along the transition path.
#'
#' @author Zhenyi Wang wangzy17@mails.tsinghua.edu.cn
#' @export

stepWisepath = function(MarkovObject=NULL,
                        MarkovLevel=NULL,
                        stepBasin=NULL){
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
    for(k in 1:(length(MarkovHCpath[[2]])-1)){
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
