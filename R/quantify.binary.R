#' Quantification of binary capture histories
#'
#' The \code{quant.binary} family of functions allow to quantify binary capture histories (partial or complete) in terms of a meaningful quantity which can be interpreted as a possibly meaningful behavioral covariate (like memory persistence of previous capture history)
#'
#' @param x either a character string or a numeric vector exclusively made by binary entries 0 or 1.
#' @param markov.ord a positive integer representing the order of the Markovian structure which one is willing to reproduce with suitable partition of the unit interval and the quantification of capture history standardized in the unit interval
#' @details For a more detailed description of instances of meaningful behavioral covariates see Alunni Fegatelli and Tardella (2012) and Alunni Fegatelli (2013)[PhD Thesis]
#' @return For \code{quant.binary} it returns a numeric value within the unit interval [0,1] and \cr \code{quant.binary.markov}. For \code{quant.binary.integer} it returns an integer value.
#' @examples
#'
#' ## Example of quantification with character input
#'
#' capt.hist="0110"
#'
#' quant.binary(capt.hist)
#' quant.binary.markov(capt.hist,markov.ord=2)
#' quant.binary.integer(capt.hist)
#' quant.binary.counts(capt.hist)
#' quant.binary.counts.integer(capt.hist)
#'
#' ## Example of quantification with numeric input
#'
#' ch=c(0,1,1,0)
#' quant.binary(ch)
#' quant.binary.markov(ch,markov.ord=2)
#' quant.binary.integer(ch)
#' quant.binary.counts(ch)
#' quant.binary.counts.integer(ch)
#'
#' @keywords Partial_capture_history_quantification
#'
#' @name quant
NULL

#' @rdname quant
#' @export
quant.binary.integer=function(x){

  if(is.character(x)){
    binarystring=x
    bh=as.numeric(unlist(strsplit(binarystring,split="")))
  }else{
    bh=x
  }
  if(length(setdiff(unique(bh),c(0,1)) )>0){
    stop('Binary history contains non binary 0/1 content')
  }

  n=length(bh)

  if(n==0){
    out=0
  }else{
  out=sum(bh*(2^(0:(n-1))))
}

  return(out)

}

#' @rdname quant
#' @export
quant.binary=function(x){

  if(is.character(x)){
    binarystring=x
    bh=as.numeric(unlist(strsplit(binarystring,split="")))
  }else{
    bh=x
  }
  if(length(setdiff(unique(bh),c(0,1)) )>0){
    stop('Binary history contains non binary 0/1 content')
  }

  n=length(bh)
  if(n==0){
    out=0
  }else{

  num=sum(bh*(2^(0:(n-1))))
  den=sum(rep(1,n)*(2^(0:(n-1))))

  out=num/den
}

  return(out)

}

#' @rdname quant
#' @export
quant.binary.markov=function(x,markov.ord){

  if(is.character(x)){
    binarystring=x
    bh=as.numeric(unlist(strsplit(binarystring,split="")))
  }else{
    bh=x
  }
  if(length(setdiff(unique(bh),c(0,1)) )>0){
    stop('Binary history contains non binary 0/1 content')
  }

  bh=c(rep(0,markov.ord-1),bh)
  n=length(bh)

  if(n==0){
    out=0
  }else{

  num=sum(bh*(2^(0:(n-1))))
  den=sum(rep(1,n)*(2^(0:(n-1))))

  out=num/den
}

  return(out)

}

#' @rdname quant
#' @export
quant.binary.counts=function(x){

  if(is.character(x)){
    binarystring=x
    bh=as.numeric(unlist(strsplit(binarystring,split="")))
  }else{
    bh=x
  }
  if(length(setdiff(unique(bh),c(0,1)) )>0){
    stop('Binary history contains non binary 0/1 content')
  }

  n=length(bh)
  if(n==0){
    out=0
  }else{

  num=sum(bh)
  den=length(bh)

  out=num/den
}
  return(out)

}

#' @rdname quant
#' @export
quant.binary.counts.integer=function(x){

  if(is.character(x)){
    binarystring=x
    bh=as.numeric(unlist(strsplit(binarystring,split="")))
  }else{
    bh=x
  }
  if(length(setdiff(unique(bh),c(0,1)) )>0){
    stop('Binary history contains non binary 0/1 content')
  }

  n=length(bh)
  if(n==0){
    out=0
  }else{

  out=sum(bh)

}
  return(out)

}
