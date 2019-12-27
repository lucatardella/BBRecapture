#' Rissanen's universal prior for integers
#'
#' It returns (up to normalizing constant) the mass assigned to each positive integer or a vector of integers by Rissanen's universal prior for positive integers
#'
#' @usage rissanen(n)
#'
#' @param n a vector of positive integers
#'
#' @details Rissanen's universal prior on positive integers is one of the default options for eliciting a noninformative prior distribution on the unknown population size \eqn{N}. It is a proper prior with tails of the order between \eqn{1/N} and \eqn{1/N^2}
#'
#' @return
#'
#' The mass assigned to each positive integer in the input vector of integers \code{n}  by Rissanen's universal prior for positive integers
#' \deqn{Q(n)=2^{-\log^*(n)} \qquad n>0}
#' where
#' \eqn{\log^*(x)= \log(x)+\log( \log (x))  + \log( \log (\log (x)))....}
#' where the sum involves only the non-negative terms.
#' Notice that masses are not normalized hence they do not add to one,
#' but to a finite positive real constant
#' \deqn{c = \sum_{n=1}^\infty Q(n)}
#'
#' @references Rissanen, J. (1983) A universal prior for integers and estimation by minimum description length. Ann. Statist. 11, no. 2, 416-431
#'
#' @author Danilo Alunni Fegatelli and Luca Tardella
#'
#' @examples
#' # Notice that masses are not normalized hence they do not add to one,
#' # but to a finite positive real constant c
#'
#' rissanen(1:5)
#'
#' @export
#'
#' @importFrom stats as.formula binomial coef density glm logLik median na.omit optimize qnorm quantile
#' @importFrom graphics abline axis plot points
#' @importFrom HI arms
#' @importFrom lme4 findbars fixef glmer ranef

rissanen=function(n) {

	n=ceiling(n)
	res=0
	a=log(n,2)
        aa=max(a)
	while((2^aa)>1) {
		res=res+a
                a=log(a,2)*(log(a,2)>0)
		aa=log(aa,2)
}

        ind=which(is.na(res))

        if(any(ind)){
        temp=c()

res.temp=res
n.single=n[ind]

for(j in 1:length(ind)){

       n.new=n.single[j]
	n.new=ceiling(n.new)
	res=0
	a=log(n.new,2)
        aa=max(a)

	while((2^aa)>1) {
		res=res+a
                a=log(a,2)*(log(a,2)>0)
		aa=log(aa,2)
              }

temp=c(temp,res)

     }

res=res.temp
res[ind]=temp

     }

		return(2^(-res))
}


