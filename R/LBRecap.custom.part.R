#' Unconditional likelihood inference for behavioural effect models based on an ad-hoc partition of the set of all partial capture histories
#'
#' Unconditional likelihood inference for a general model framework based on the capture probabilities conditioned on each possible partial capture history. As suggested in Alunni Fegatelli and Tardella (2012) the conditional approach originally proposed in Farcomeni (2011)  [saturated reparameterization] is reviewed in terms of partitions into equivalence classes of conditional probabilities. In this function the user can directly provide the model as a partition.
#'
#' @usage LBRecap.custom.part (data,last.column.count=FALSE, partition, neval = 1000,
#'    by.incr = 1, output = c("base", "complete"))
#'
#' @param data can be one of the following:
#' \enumerate{
#'   \item an \eqn{M} by \eqn{t} binary matrix/data.frame
#'   \item a matrix/data.frame with \eqn{(t+1)} columns according to the value of \cr \code{last.column.count}
#'   \item a \eqn{t}-dimensional array or table representing the counts of the \eqn{2^t} contingency table of binary outcomes
#'   \eqn{M} is the number of units captured at least once and \eqn{t} is the number of capture occasions.
#'   }
#' @param last.column.count a logical. In the default case \code{last.column.count=FALSE} each row of \code{data} represents the complete capture history for each observed unit. When code{last.column.count=TRUE} in each row the first \eqn{t} entries represent one of the possible observed complete capture histories and the last entry (last column) is the number of observed units with that capture history


#' @param partition list. \code{partition} represents a partition of the set of all partial capture histories.
#' @param neval a positive integer. \code{neval} is the number of values evaluated for the population size N. The default value is \code{neval}=1000.
#' @param by.incr a positive integer. \code{by.incr} represents the increment on the sequence of evaluated values for \eqn{N}. The default value is \code{by.incr}=1.
#' @param output character. \code{output} select the kind of output from a very basic summary info on the posterior output (point and interval estimates for the unknown \code{N}) to more complete details.
#'
#' @details The unconditional likelihood is evaluated by means of \code{glm/glmer} for each value of the \code{N} parameter and it is then maximized.
#'
#' @return
#'  (if \code{output="complete"}) the function \code{LBRecap} returns a list of:
#'  \enumerate{
#'   \item{N.hat}{unconditional maximum likelihood estimate for \eqn{N}}
#'   \item{CI}{interval estimate for \eqn{N}}
#'   \item{pH.hat}{point estimate of nuisance parameters (conditional probabilities)}
#'   \item{AIC}{Akaike information criterion.}
#'   \item{L.Failure}{Likelihood Failure condition}
#'   \item{N.range}{sequence of \eqn{N} values considered}
#'   \item{log.lik}{values of the log-likelihood distribution for each \eqn{N} considered}
#'   \item{partitions}{list of subsets of partial capture histories corresponding to equivalence classes}
#'   }
#'
#' @references
#' Alunni Fegatelli, D. and Tardella, L. (2016), Flexible behavioral capture–recapture modeling. Biometrics, 72(1):125-135. doi:10.1111/biom.12417
#'
#' Alunni Fegatelli, D. and Tardella, L. (2012) Improved inference on capture recapture models with behavioural effects. Statistical Methods & Applications Applications Volume 22, Issue 1, pp 45-66 10.1007/s10260-012-0221-4
#'
#' Farcomeni A. (2011) Recapture models under equality constraints for the conditional capture probabilities. Biometrika 98(1):237--242
#'
#' @author Danilo Alunni Fegatelli and Luca Tardella
#'
#' @seealso \code{\link{partition.ch}}, \code{\link{BBRecap.custom.part}}, \code{\link{LBRecap}}
#'
#' @examples
#' data(greatcopper)
#' partition.Mc1=partition.ch(quant.binary,t=ncol(greatcopper),breaks=c(0,0.5,1))
#' mod.Mc1=LBRecap.custom.part(greatcopper,partition=partition.Mc1)
#' str(mod.Mc1)
#'
#' @keywords Behavioural_models Unconditional_MLE
#'
#' @export
LBRecap.custom.part=function(
data,
last.column.count=FALSE,
partition,
neval=1000,
by.incr=1,
output=c("base","complete")){

 output=match.arg(output)

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

	partial=pch(data.matrix)

	mm1=matrix(NA,ncol=ncol(data.matrix),nrow=nrow(data.matrix))
	mm0=matrix(NA,ncol=ncol(data.matrix),nrow=nrow(data.matrix))
	cc=c()
	n.obs.1=c()
	n.obs.0=c()
	n.unobs=c()
 	p.1=matrix(ncol = length(partition), nrow = neval)

 	log.lik=c()

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

	failure.conditions=F

	if(sum(as.logical(n.unobs))==1){failure.conditions=T}

  for(l in 1:neval){

	val=seq(M,((M+(neval-1)*by.incr)),by=by.incr)
  	nn.0=n.obs.0+(n.unobs*(l-1)*by.incr)
  	nn.1=n.obs.1

p.1[l,]=nn.1/(nn.0+nn.1)

n.00=nn.0[nn.0>0]
n.10=nn.1[nn.0>0]

n.01=nn.0[nn.1>0]
n.11=nn.1[nn.1>0]

log.lik[l]=lchoose((sum(nn.1)+sum(nn.0))/t,M)+(sum(n.11*log(n.11/(n.11+n.01)))+(sum(n.00*log(n.00/(n.10+n.00)))))

	}

 AIC=-2*max(log.lik)+2*(length(partition)+1)

 u=which(log.lik==max(log.lik))-1

 N.hat=M+(u*by.incr)

 pH.hat=p.1[which(log.lik==max(log.lik)),]
 names(pH.hat)=paste("pH",1:length(pH.hat),sep="")

  ## CONFIDENCE INTERVAL (normal approximation)

	z2=(qnorm(0.975))^2

 Nplus=N.hat

 	if(N.hat<(M+neval*(by.incr))){
 i=1

 while((2*(max(log.lik)-log.lik[(which(log.lik==(max(log.lik)))+i)])<z2 & (N.hat+i*(by.incr))< M+neval*(by.incr)))
 {
 	i=i+1
 	}

 	 Nplus=N.hat+i*by.incr
 	}

 Nminus=N.hat

 	if(N.hat>M){
 i=1

 while(2*(max(log.lik)-log.lik[(which(log.lik==(max(log.lik)))-i)])<z2 &(N.hat-i*(by.incr))>M){
 	i=i+1
 	}
 	 Nminus=N.hat-i*by.incr
	}


 out=switch(output,
 base=list(N.hat=N.hat,CI=c(Nminus,Nplus),AIC=AIC,L.Failure=FALSE),
 complete=list(N.hat=N.hat,CI=c(Nminus,Nplus),pH.hat=pH.hat,AIC=AIC,L.Failure=FALSE,log.lik=log.lik,N=seq(M,neval*by.incr+M,by=by.incr)[-neval+1],Partition=partition))

  if(failure.conditions==T){
  	n1=n.obs.1[which(n.unobs!=0)]
  	n0=n.obs.0[which(n.unobs!=0)]
  if((M*(t-1)-n0<=((M-1)*(t-1))/2-1) & n1==M) {
  out=switch(output,
  base=list(N.hat=NA,CI=c(NA,NA),AIC=-Inf,L.Failure=TRUE),
  complete=list(N.hat=NA,CI=c(NA,NA),AIC=-Inf,L.Failure=TRUE,log.lik=log.lik,N=seq(M,neval*by.incr+M,by=by.incr)[-neval+1],Partition=partition))
		 }
	 }

  return(out)
}
