#' Recommend levels with possible biological meaning
#'
#' Function \code{energyGap_selection} The function recommends levels
#' with possible biological meaning.
#' @param MarkovObject The output of the function, \code{MarkovHC}.
#' @param m An integer value. A (local) peak is defined as
#' a point such that m points either side of it has a lower or equal value to it.
#' @return  This function plots the energy gap along the increasing levels,
#' red points indicate recommend levels with possible biological meaning,
#' and a red circle highlight the level that may with an optimal cluster number.
#'
#' A list consist of two components are returned. p is a vector contains the recommend
#' levels with possible biological meaning. recommend_level is the level that may
#' with an optimal cluster number.
#'
#' @author Zhenyi Wang wangzy17@mails.tsinghua.edu.cn
#' @export

energyGap_selection = function(MarkovObject=NULL,
                               m=3){
  C_cut_gap <- diff(MarkovObject$midResults$C_cut_seq)
  model <- piecewise.linear(x=1:length(C_cut_gap), y=C_cut_gap)
  recommend_level <- model$change.point
  # find peaks
  p <- find_peaks(C_cut_gap, m = m)
  #temp <- which.min(abs(p-recommend_level))
  #recommend_level <- p[temp]
  recommend_level <- ceiling(recommend_level)
  #
  C_cut_gap <- c(0, C_cut_gap)
  recommend_level <- recommend_level+1
  p <- p+1
  #
  plot(C_cut_gap, xlab='levels', ylab='energy gap')
  lines(C_cut_gap)
  points(p, C_cut_gap[p], col = 'red', pch = 19)
  points(recommend_level, C_cut_gap[recommend_level], col = 'red', pch = 1,cex=3,lwd=3)
  points(p[length(p)], C_cut_gap[p[length(p)]], col = 'red', pch = 1,cex=3,lwd=3)
  #return(list(recommend_level,p))
  print('levels with possible biological meaning:')
  print(p)
  print('the level may with an optimal cluster number is among:')
  print(paste("levels:from",recommend_level,"to",p[length(p)]))
}
