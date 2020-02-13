#' Standard behavioural and time effect models via unconditional (complete) likelihood approach
#'
#' Comparative point and interval estimates for the population size \eqn{N} obtained fitting many alternative behavioural and time effect capture-recapture models.
#' AIC index is reported for each alternative model.
#'
#' @usage LBRecap.all(data, last.column.count=FALSE, neval=1000, by.incr=1,
#'    which.mod=c("all","standard"), sort=c("default","AIC"))
#'
#' @param data can be one of the following:
#' \enumerate{
#'   \item an \eqn{M} by \eqn{t} binary matrix/data.frame
#'   \item a matrix/data.frame with \eqn{(t+1)} columns according to the value of \cr \code{last.column.count}
#'   \item a \eqn{t}-dimensional array or table representing the counts of the \eqn{2^t} contingency table of binary outcomes
#'   \eqn{M} is the number of units captured at least once and \eqn{t} is the number of capture occasions.
#'   }
#' @param last.column.count a logical. In the default case \code{last.column.count=FALSE} each row of \code{data} represents the complete capture history for each observed unit. When code{last.column.count=TRUE} in each row the first \eqn{t} entries represent one of the possible observed complete capture histories and the last entry (last column) is the number of observed units with that capture history
#' @param neval a positive integer. \code{neval} is the number of alternative values of the population size N where the likelihood is evaluated and then maximized. They run from the minimum value M and they are increased by \code{by.incr} (see below the description of the \code{by.incr} argument). The default value is \code{neval}=1000.
#' @param by.incr a positive integer. \code{by.incr} represents the increment on the sequence of evaluated values for \eqn{N}. The default value is \code{by.incr}=1.
#' @param which.mod a character. \code{which.mod} selects which models are fitted and compared. In the default setting \code{which.mod}=\code{"all"} all alternative models are fitted including new behavioural models based on alternative meaningful covariates (see Details). When \code{which.mod}=\code{"standard"} the function only fits classical behavioural models with either enduring effects as in \eqn{M_b}, \eqn{M_{c_1b}}, \eqn{M_{c_2b}} or ephemeral effects as in purely Markovian \eqn{M_{c_1}} and \eqn{M_{c_2}}
#' @param sort character. \code{sort} selects the order of models.
#'
#' @details The available models are: \eqn{M_0}, \eqn{M_b}, \eqn{M_t}, \eqn{M_{c_1}}, \eqn{M_{c_1b}}, \eqn{M_{c_2}}, \eqn{M_{c_2b}}, \eqn{M_{mc}},\eqn{M_{mc_{int}}}, \eqn{M_{mc_{count}}} and \eqn{M_{mc_{count.int}}}.
#' This function \code{LBRecap.all} can be computing intensive for high values of \code{neval}.
#'
#' @return
#'
#' A \code{data.frame} with one row corresponding to each model and the following columns:
#'
#' model: model considered
#'
#' npar: number of parameters
#'
#' AIC: Akaike's information criterion
#'
#' Nhat: estimate of population size
#'
#' Ninf: lower \eqn{95 \%} confidence limit
#'
#' Nsup: upper \eqn{95 \%} confidence limit
#'
#' @references
#' Alunni Fegatelli, D. and Tardella, L. (2016), Flexible behavioral capture–recapture modeling. Biometrics, 72(1):125-135. doi:10.1111/biom.12417
#'
#' Alunni Fegatelli D. (2013) New methods for capture-recapture modelling with behavioural response and individual heterogeneity. http://hdl.handle.net/11573/918752
#'
#' Alunni Fegatelli D., Tardella L. (2012) Improved inference on capture recapture models with behavioural effects. Statistical Methods & Applications Applications Volume 22, Issue 1, pp 45-66 10.1007/s10260-012-0221-4
#'
#' Farcomeni A. (2011) Recapture models under equality constraints for the conditional capture probabilities. Biometrika 98(1):237--242
#'
#' Otis D. L., Burnham K. P., White G. C, Anderson D. R. (1978) Statistical Inference From Capture Data on Closed Animal Populations, Wildlife Monographs.
#'
#' Yang H.C., Chao A. (2005) Modeling animals behavioral response by Markov chain models for capture-recapture experiments, Biometrics 61(4), 1010-1017
#'
#' @author Danilo Alunni Fegatelli and Luca Tardella
#'
#' @seealso \code{\link{LBRecap}}, \code{\link{BBRecap.all}}
#'
#' @examples
#' \dontrun{
#'  data(greatcopper)
#'  LBRecap.all(greatcopper)
#' }
#' @keywords Behavioural_models Unconditional_MLE
#' @export
LBRecap.all=function(data, last.column.count=FALSE, neval=1000, by.incr=1,which.mod=c("all","standard"),sort=c("default","AIC")){

    ss=sort[1]
    mm=which.mod[1]

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

  if(last.column.count==TRUE){
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
        #### BETTER USING data.frame and subset(...)
        data.matrix=data.matrix[apply(data.matrix,1,sum)!=0,]
    }

 t=ncol(data.matrix)
 M=nrow(data.matrix)

	mod.M0=LBRecap(data.matrix,mod="M0",neval=neval,by.incr=by.incr)
	mod.Mb=LBRecap(data.matrix,mod="Mb",neval=neval,by.incr=by.incr)
	mod.Mlogistic=LBRecap(data.matrix,neval=neval,by.incr=by.incr)
	mod.Mc1=LBRecap(data.matrix,mod="Mc",markov.ord=1,neval=neval,by.incr=by.incr)
	mod.Mc2=LBRecap(data.matrix,mod="Mc",markov.ord=2,neval=neval,by.incr=by.incr)
	mod.Mc1b=LBRecap(data.matrix,mod="Mcb",markov.ord=1,neval=neval,by.incr=by.incr)
	mod.Mc2b=LBRecap(data.matrix,mod="Mcb",markov.ord=2,neval=neval,by.incr=by.incr)
        mod.Mt=LBRecap(data.matrix,mod="Mt",neval=neval)

output=data.frame()
model=c("M0","Mb","Mlogistic","Mc1","Mc2","Mc1b","Mc2b","Mt")
npar=c(2,3,3,3,5,4,6,t+1)
AIC=c(mod.M0$AIC,mod.Mb$AIC,mod.Mlogistic$AIC,mod.Mc1$AIC,mod.Mc2$AIC,mod.Mc1b$AIC,mod.Mc2b$AIC,mod.Mt$AIC)
Nhat=c(mod.M0$N.hat,mod.Mb$N.hat,mod.Mlogistic$N.hat,mod.Mc1$N.hat,mod.Mc2$N.hat,mod.Mc1b$N.hat,mod.Mc2b$N.hat,mod.Mt$N.hat)
Ninf=c(mod.M0$CI[1],mod.Mb$CI[1],mod.Mlogistic$CI[1],mod.Mc1$CI[1],mod.Mc2$CI[1],mod.Mc1b$CI[1],mod.Mc2b$CI[1],mod.Mt$CI[1])
Nsup=c(mod.M0$CI[2],mod.Mb$CI[2],mod.Mlogistic$CI[2],mod.Mc1$CI[2],mod.Mc2$CI[2],mod.Mc1b$CI[2],mod.Mc2b$CI[2],mod.Mt$CI[2])

         output=data.frame(Model=model,npar=npar,AIC=AIC,Nhat=Nhat,Ninf=Ninf,Nsup=Nsup)

if(which.mod=="all"){

	mod.Mlogistic.integer=LBRecap(data.matrix,neval=neval,by.incr=by.incr,mbc.function="integer")
	mod.Mlogistic.counts=LBRecap(data.matrix,neval=neval,by.incr=by.incr,mbc.function="counts")
	mod.Mlogistic.counts.integer=LBRecap(data.matrix,neval=neval,by.incr=by.incr,mbc.function="counts.integer")



model=c("M0","Mb","Mlogistic","Mc1","Mc2","Mc1b","Mc2b","Mt","Mlogistic.integer","Mlogistic.counts","Mlogistic.counts.integer")
npar=c(2,3,3,3,5,4,6,t+1,rep(3,3))
AIC=c(mod.M0$AIC,mod.Mb$AIC,mod.Mlogistic$AIC,mod.Mc1$AIC,mod.Mc2$AIC,mod.Mc1b$AIC,mod.Mc2b$AIC,mod.Mt$AIC,mod.Mlogistic.integer$AIC,mod.Mlogistic.counts$AIC,mod.Mlogistic.counts.integer$AIC)
Nhat=c(mod.M0$N.hat,mod.Mb$N.hat,mod.Mlogistic$N.hat,mod.Mc1$N.hat,mod.Mc2$N.hat,mod.Mc1b$N.hat,mod.Mc2b$N.hat,mod.Mt$N.hat,mod.Mlogistic.integer$N.hat,mod.Mlogistic.counts$N.hat,mod.Mlogistic.counts.integer$N.hat)
Ninf=c(mod.M0$CI[1],mod.Mb$CI[1],mod.Mlogistic$CI[1],mod.Mc1$CI[1],mod.Mc2$CI[1],mod.Mc1b$CI[1],mod.Mc2b$CI[1],mod.Mt$CI[1],mod.Mlogistic.integer$CI[1],mod.Mlogistic.counts$CI[1],mod.Mlogistic.counts.integer$CI[1])
Nsup=c(mod.M0$CI[2],mod.Mb$CI[2],mod.Mlogistic$CI[2],mod.Mc1$CI[2],mod.Mc2$CI[2],mod.Mc1b$CI[2],mod.Mc2b$CI[2],mod.Mt$CI[2],mod.Mlogistic.integer$CI[2],mod.Mlogistic.counts$CI[2],mod.Mlogistic.counts.integer$CI[2])

output=data.frame(Model=model,npar=npar,AIC=AIC,Nhat=Nhat,Ninf=Ninf,Nsup=Nsup)


}


#if(ss=="npar"){
#    output=output[order(npar),]
#}


if(ss=="AIC"){
    output=output[order(AIC),]
}

print(paste("neval =",neval,sep=" "))
print(paste("M =",M,sep=" "))
return(output)

}
