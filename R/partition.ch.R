#' Partition of partial capture histories according to equivalence classes of numerical quantification corresponding to supplied intervals
#'
#' All the possible partial capture histories observable during a capture-recapture experiment with \eqn{t} sampling occasions can be partitioned according to numerical values corresponding to some meaningful covariate (quantification of binary sequences corresponding to partial capture histories). Each subset of the partition corresponds to all partial capture histories which returns numerical values of the quantification within one of the intervals represented by two consecutive values in the optional argument vector \code{breaks}.
#'
#' @usage partition.ch(quantify.ch.fun, t, breaks, include.lowest = T,
#' type = c("list", "index"), ...)
#'
#' @param quantify.ch.fun a function which returns a numerical value for each possible partial capture history
#' @param t an integer. \code{t} is number of trapping occasions
#' @param breaks a vector of numerical values which are used as bounds for the interval of numerical values corresponding to partial capture histories that belongs to the same partition
#' @param include.lowest a logical, indicating if an x[i] equal to the lowest (or highest, when  right = FALSE) breaks value should be included
#' @param type a character string. It can be either  \code{"list"} or  \code{"index"}. See examples.
#' @param \\dots additional arguments to be passed to \code{quantify.ch.fun}
#'
#' @details It is useful in conjunction with \code{LBRecap.custom.part}. See examples.
#'
#' @return
#'
#' If the argument \code{type="list"} a list is returned. If \code{type="index"} a numerical index corresponding to the numeric integer equivalent of the consecutive interval according to the convention used in objects of class \code{factor}
#'
#' @references REFERENCES
#'
#' @author Danilo Alunni Fegatelli and Luca Tardella
#'
#' @seealso \code{\link{LBRecap.custom.part}}, \code{\link{BBRecap.custom.part}}
#'
#' @examples
#'
#' data(mouse)
#' head(mouse)
#' t=ncol(mouse)
#'
#' Mc1.partition=partition.ch(quantify.ch.fun=quant.binary,t=t,breaks=c(0,0.5,1))
#' Mc1.partition
#'
#' mod.Mc1.cust=BBRecap.custom.part(mouse,partition=Mc1.partition)
#' mod.Mc1.cust
#'
#' mod.Mc1.easy=BBRecap(mouse,mod="Mc",markov.ord=1,output="complete")
#'
#' mod.Mc1.easy$N.hat.RMSE
#' mod.Mc1.easy$HPD.N
#' mod.Mc1.easy$log.marginal.likelihood
#'
#' # the two functions give the same results!
#'
#'
#' @export

partition.ch=function(quantify.ch.fun,t,breaks,include.lowest=T,type=c("list","index"), ... ){

  ch.list=unlist(list.historylabels(t))
  names(ch.list)=ch.list
  names(ch.list)[1]="empty"
  ch.quantified=sapply(ch.list,FUN=quantify.ch.fun, ... )
partition.index=cut(ch.quantified,breaks=breaks,
                      include.lowest=include.lowest)


  if(type[1]=="list"){
    out=split(unlist(ch.list),partition.index)
  }else{
    out=as.numeric(partition.index)
    names(out)=names(ch.list)
  }

  return(out)

}

