#' List all the observable partial capture histories
#'
#' This function returns a list of all the observable partial capture histories which can be recorded in a discrete-time capture-recapture setting with \eqn{t} consecutive trapping occasions. The observable partial capture histories are \eqn{2^t-1}
#'
#' @usage list.historylabels(t,t.max=15)
#'
#' @param t a positive integer representing the total number of trapping occasions
#' @param t.max a positive integer representing upper bound on the total number of trapping occasions allowed.

#' @details For obvious computing/memory reasons t is not allowed to be arbitrarily large. With \code{t.max=15} there are 32767 possible partial capture histories. If \code{t>t.max} the function stops with an error message.
#'
#' @return
#'
#' A list of all the observable partial capture histories which can be recorded in a discrete-time capture-recapture setting with t consecutive trapping occasions. If \code{t>t.max} the function stops with an error message.
#'
#' @author Danilo Alunni Fegatelli and Luca Tardella
#'
#' @seealso \code{\link{partition.ch}}
#'
#' @examples
#'
#' list.historylabels(t=4)
#'
#' @export

list.historylabels = function(t, t.max=15){

  if(t>t.max) stop(paste("argument t cannot be greater than t.max=",t.max,sep=""))

  ph.list=vector(mode="list",length=0)

  for(j in 1:(t-1)){
    fac=vector(mode="list",length=j)
    lapply(fac[1:j], function(x) factor(c(0,1))) -> fac[1:j]
    matindex=expand.grid(fac)
    matindex=as.data.frame(t(matindex))
    lapply(matindex[1:ncol(matindex)],function(x) factor(x,levels=c("0","1"))) -> matindex[1:ncol(matindex)]
    ph.list=c(ph.list,as.list(matindex))
  }

  ph.list=c("",ph.list)

  out=lapply(ph.list,paste,collapse="")

  return(out)
}
