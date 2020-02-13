#' Bayesian inference for behavioural effect models based on a partition of the set of all partial capture histories
#'
#' Bayesian inference for a general model framework based on the capture probabilities conditioned on each possible partial capture history.  As suggested in Alunni Fegatelli and Tardella (2012) the conditional approach originally proposed in Farcomeni (2011)  [saturated reparameterization] is reviewed in terms of partitions into equivalence classes of conditional probabilities. In this function the user can directly provide the model as a partition.
#'
#' @usage BBRecap.custom.part (data,last.column.count=FALSE, partition, neval = 1000,
#'    by.incr = 1, prior.N = c("Rissanen", "Uniform", "one.over.N", "one.over.N2"),
#'    output = c("base", "complete"))
#'
#' @param data can be one of the following: \enumerate{
#'   \item an \eqn{M} by \eqn{t} binary matrix/data.frame. In this case the input is interpreted as a matrix whose rows contain individual capture histories for all \eqn{M} observed units
#'   \item a matrix/data.frame with \eqn{(t+1)} columns. The first \eqn{t} columns contain binary entries corresponding to capture occurrences, while the last column contains non negative integers corresponding to frequencies. This format is allowed only when \code{last.column.count} is set to \code{TRUE}
#'   \item a \eqn{t}-dimensional array or table representing the counts of the \eqn{2^t} contingency table of binary outcomes
#' }
#' \eqn{M} is the number of units captured at least once and \eqn{t} is the number of capture occasions.
#' @param last.column.count a logical. In the default case \code{last.column.count=FALSE} each row of the input argument \code{data} represents the complete capture history for each observed unit. When \code{last.column.count} is set to \code{TRUE} in each row the first \eqn{t} entries represent one of the observed complete capture histories and the last entry in the last column is the number of observed units with that capture history
#' @param partition list. \code{partition} represents a partition of the set of all partial capture histories.
#' @param neval a positive integer. \code{neval} is the number of values of the population size \eqn{N} where the posterior is evaluated starting from \eqn{M}. The default value is \code{neval}=1000.
#' @param by.incr a positive integer. \code{by.incr} represents the increment on the sequence of possible population sizes \eqn{N} where the posterior is evaluated. The default value is \code{by.incr}=1. The use of \code{by.incr}>1 is discouraged unless the range of \eqn{N} values of interest is very large
#' @param prior.N a character. \code{prior.N} is the label for the prior distribution for \eqn{N}. When \code{prior.N} is set to \code{"Rissanen"} (default) the Rissanen prior is used as a prior on \eqn{N}. This distribution was first proposed in Rissanen 1983 as a universal prior on integers. \code{prior.N="Uniform"} stands for a prior on \eqn{N} proportional to a constant value. \code{prior.N="one.over.N"} stands for a prior on \eqn{N} proportional to \eqn{1/N}. \code{prior.N="one.over.N2"} stands for a prior on \eqn{N} proportional to \eqn{1/N^2}.
#' @param output a character. \code{output} selects the kind of output from a very basic summary info on the posterior output (point and interval estimates for the unknown \eqn{N}) to more complete details
#'
#' @details
#' Uniform prior distribution is considered for the nuisance parameters.
#'
#' @return
#'
#' \item{Prior}{prior distribution for \code{N}.}
#' \item{N.hat.mean}{posterior mean for \eqn{N}}
#' \item{N.hat.median}{posterior median for \eqn{N}}
#' \item{N.hat.mode}{posterior mode for \eqn{N}}
#' \item{N.hat.RMSE}{minimizer of a specific loss function connected with the Relative Mean Square Error.}
#' \item{HPD.N}{\eqn{95 \%} highest posterior density interval estimate for \eqn{N}.}
#' \item{log.marginal.likelihood}{log marginal likelihood.}
#' \item{N.range}{values of N considered.}
#' \item{posterior.N}{values of the posterior distribution for each N considered}
#' \item{partition}{partition of the set H}
#'
#' @references
#'
#' Farcomeni A. (2011) Recapture models under equality constraints for the conditional capture probabilities. Biometrika 98(1):237--242
#'
#' Alunni Fegatelli, D. and Tardella, L. (2012) Improved inference on capture recapture models with behavioural effects. Statistical Methods & Applications Applications Volume 22, Issue 1, pp 45-66 10.1007/s10260-012-0221-4
#'
#' Alunni Fegatelli D. (2013) New methods for capture-recapture modelling with behavioural response and individual heterogeneity. http://hdl.handle.net/11573/918752
#'
#' @author Danilo Alunni Fegatelli and Luca Tardella
#'
#' @seealso \code{\link{partition.ch}}, \code{\link{LBRecap.custom.part}}, \code{\link{BBRecap}}
#'
#' @examples
#'
#' data(greatcopper)
#' partition.Mc1=partition.ch(quant.binary,t=ncol(greatcopper),breaks=c(0,0.5,1))
#' mod.Mc1=BBRecap.custom.part(greatcopper,partition=partition.Mc1)
#' str(mod.Mc1)
#'
#' @keywords Behavioural_models Bayesian_inference

#' @export
BBRecap.custom.part=function(
data,
last.column.count=FALSE,
partition,
neval=1000,
by.incr=1,
prior.N=c("Rissanen","Uniform","one.over.N","one.over.N2"),
output=c("base","complete")){

 prior.N=match.arg(prior.N)
 output=match.arg(output)

 prior=switch(prior.N,
 Rissanen="Rissanen.prior",
 Uniform="Uniform.prior",
 one.over.N="1overN.prior",
 one.over.N2="1overN^2.prior"
 )

  if(!(any(c("data.frame","matrix","array","table") %in% class(data)))){
       stop("input data must be a data.frame or a matrix object or an array")
  }

 data.matrix=data

   if( !("matrix" %in% class(data)) & any(c("array","table") %in% class(data))){

   	n.occ=length(dim(data))
mm=matrix(ncol=n.occ,nrow=2^n.occ)

for(i in 1:2^n.occ){
mm[i,]=as.numeric(intToBits(i-1)>0)[1:n.occ]
}

data.vec=c(data)
data[1]=0
dd=c()
for(i in 1:length(data)){

	dd=c(dd,rep(mm[i,],data.vec[i]))

		}

data.matrix=matrix(dd,ncol=length(dim(data)),byrow=T)

   }


   if(!("matrix" %in% class(data))){
       data.matrix=as.matrix(data)
       }

  if(last.column.count){

if (any(data[,ncol(data)]<0)){
    stop("Last column must contain non negative frequencies/counts")
}

  	data=as.matrix(data)
  	data.matrix=matrix(ncol=(ncol(data)-1))

for(i in 1:nrow(data)){

d=rep(data[i,1:(ncol(data)-1)],(data[i,ncol(data)]))
dd=matrix(d,ncol=(ncol(data)-1),byrow=T)
data.matrix=rbind(data.matrix,dd)

}


data.matrix=data.matrix[-1,]

  }

  if(any(data.matrix!=0 & data.matrix!=1)) stop("data must be a binary matrix")

  if(sum(apply(data.matrix,1,sum)==0)){
 	warning("input data argument contains rows with all zeros: these rows will be removed and ignored")

        data.matrix=data.matrix[apply(data.matrix,1,sum)!=0,]
    }

 t=ncol(data.matrix)
 M=nrow(data.matrix)

	p.p=c()

	for(r in 1:length(partition)){
	pp=unique(partition[[r]])
	p.p=c(p.p,pp)
	}

 	p.p=sort(unique(p.p))

 	cond1=!all(sort(unlist(list.historylabels(ncol(data.matrix))))==p.p)

 	cond2=max(table(unlist(partition)))>1

 	if(cond1 | cond2)
 		stop("The input argument 'partition' does not represent a partition of the set of all partial capture histories")


 	prior.distr.N=function(x){

 	if(prior=="Rissanen.prior") {out=(rissanen(x))}
 	if(prior=="Uniform.prior") {out=(1/(x^0))}
 	if(prior=="1overN.prior") {out=(1/(x^1))}
 	if(prior=="1overN^2.prior") {out=(1/(x^2))}

 	return(out)
 	}

	partial=pch(data.matrix)

	mm1=matrix(NA,ncol=ncol(data.matrix),nrow=nrow(data.matrix))
	mm0=matrix(NA,ncol=ncol(data.matrix),nrow=nrow(data.matrix))
	cc=c()
	n.obs.1=c()
	n.obs.0=c()
	n.unobs=c()

	log.post.N=c()
 	post.N=c()
 	vv=c()
  	prior.inv.const=0

	for(k in 1:length(partition)){

	for(i in 1:nrow(data.matrix)){
	for(j in 1:ncol(data.matrix)){
	mm1[i,j]=any(partition[[k]]==partial[i,j]) & data.matrix[i,j]==1
	mm0[i,j]=any(partition[[k]]==partial[i,j]) & data.matrix[i,j]==0
	}
	}
	n.obs.0[k]=sum(mm0)
	n.obs.1[k]=sum(mm1)
	}

	for(k in 1:length(partition)){

	for(j in 1:ncol(data.matrix)){
	cc[j]=any(partition[[k]]==pch(matrix(rep(0,ncol(data.matrix)),nrow=1))[j])
	}
	n.unobs[k]=sum(cc)
	}

  for(l in 1:neval){

	val=seq(M,((M+(neval-1)*by.incr)),by=by.incr)
  	nn.0=n.obs.0+(n.unobs*(l-1)*by.incr)
  	nn.1=n.obs.1

 prior.inv.const=prior.inv.const+ prior.distr.N((M+l-1))

 log.post.N[l]=lchoose((sum(nn.1)+sum(nn.0))/t,M)+sum(lbeta((nn.1+1),(nn.0+1)))+log(prior.distr.N((M+l-1)))

	}

	l.max=max(log.post.N)

	for(k in 1:neval){
	vv[k]<-exp(log.post.N[k]-l.max)
	}

	ss=sum(vv)
	post.N=vv/ss

#	return(post.N)

ord=order(post.N,decreasing=T)
p.max=ord[1]
mode.N=val[p.max]

mean.N<-round(sum(post.N*val))

median.N=M

g<-c()
for(k in 1:neval)	{
g=c(g,post.N[k])
if (sum(g)<=0.5) median.N=median.N+1
			}

funzioneRMSE=function(x){
	sum((((x/val)-1)^2)*post.N)
	}

estimate.N=round(optimize(funzioneRMSE,c(min(val),max(val)))$minimum)

### Credible Set ###

alpha<-0.05

g=0
d=0

aa=c()

ordine<-order(post.N,decreasing=T)

w<-val
w<-w[ordine]
p<-post.N
p<-p[ordine]

for(k in 1:neval)		{
if (g<(1-alpha)) 	{g=g+p[k]
  		 	 d=d+1}
				}
aa<-w[1:d]
inf.lim<-min(aa)
sup.lim<-max(aa)

log.marg.likelihood=log(sum(exp(log.post.N-max(log.post.N)-log(prior.inv.const))))+max(log.post.N)

out=switch(output,
base=
list(Prior.N=prior.N,N.hat.RMSE=estimate.N,HPD.N.=c(inf.lim,sup.lim),log.marginal.likelihood=log.marg.likelihood),
complete=
list(Prior.N=prior.N,N.hat.mean=mean.N,N.hat.median=median.N,N.hat.mode=mode.N,N.hat.RMSE=estimate.N,HPD.N=c(inf.lim,sup.lim),log.marginal.likelihood=log.marg.likelihood,
N.range=val,posterior.N=post.N,Partition=partition))

  return(out)

}
