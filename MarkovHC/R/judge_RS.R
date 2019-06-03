#' Get a vector to judge recurrent state
#'
#' Function '\code{judge_RS}' gets a vector helping to judge recurrent state based on a transition matrix.
#' @param P A transition matrix.
#' @details Assume that we have an input transition matrix P. We can first find out all
#' eigen values and eigen vectors of P^T. As we know that if a state i is a recurrent state,
#' then there must be a eigen value 1 with a corresponding eigen vector whose i th entry
#' is non-zero. If a state i is a transient state, then for all eigen values 1, the i th entry
#' of the corresponding eigen vector must be 0. According to the above conclusions, we can easily
#' judge whether a state is recurrent or not.
#' @return
#'     \tabular{ll}{
#'     \code{RS} \tab A vector used to judge recurrency. If the i th entry is non-zero, then state i is a recurrent state. \cr
#'     }
#' @author Zhenyi Wang
#' @examples
#' P<-matrix(c(1,0,0,0,0.5,0,0,0.5,1),nrow=3)
#' RS<-judge_RS(P)

judge_RS<-function(P=NULL){
  RSmatrix<-eigen(t(P),symmetric = FALSE)
  #the eigen value is close to 1 enough
  eigen_value_one<-which(abs(RSmatrix$values-1)<1e-10)
  if(length(eigen_value_one)==1){
    return(abs(Re(RSmatrix$vectors[,eigen_value_one])))
  }else{
    RS<-apply(abs(Re(RSmatrix$vectors[,eigen_value_one])),1,sum)
  }
  return(RS)
}
