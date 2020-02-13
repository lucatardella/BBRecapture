#' Comparative Bayesian analysis of alternative flexible behavioural and time effect models
#'
#' Comparative point and interval estimates for the population size \eqn{N} obtained fitting many alternative behavioural and time effect capture-recapture models.
#' Log marginal likelihood is reported for each alternative model.
#'
#' @usage BBRecap.all(data, last.column.count=FALSE, neval=1000, by.incr=1,nsim=10000,
#' burnin=round(nsim/10),nsim.ML=500,burnin.ML=round(nsim.ML/10), num.t = 50,
#' prior.N = c("Rissanen","Uniform","one.over.N","one.over.N2"),
#' which.mod=c("all","standard"), sort=c("default","log.ML"))
#'
#' @param data can be one of the following: \enumerate{
#'   \item an \eqn{M} by \eqn{t} binary matrix/data.frame. In this case the input is interpreted as a matrix whose rows contain individual capture histories for all \eqn{M} observed units
#'   \item a matrix/data.frame with \eqn{(t+1)} columns. The first \eqn{t} columns contain binary entries corresponding to capture occurrences, while the last column contains non negative integers corresponding to frequencies. This format is allowed only when \code{last.column.count} is set to \code{TRUE}
#'   \item a \eqn{t}-dimensional array or table representing the counts of the \eqn{2^t} contingency table of binary outcomes
#' }
#' \eqn{M} is the number of units captured at least once and \eqn{t} is the number of capture occasions.
#' @param last.column.count a logical. In the default case \code{last.column.count=FALSE} each row of the input argument \code{data} represents the complete capture history for each observed unit. When \code{last.column.count} is set to \code{TRUE} in each row the first \eqn{t} entries represent one of the observed complete capture histories and the last entry in the last column is the number of observed units with that capture history
#' @param neval a positive integer. \code{neval} is the number of values of the population size \eqn{N} where the posterior is evaluated starting from \eqn{M}. The default value is \code{neval}=1000.
#' @param by.incr a positive integer. \code{nsim} is the number of iterations for the Metropolis-within-Gibbs algorithm which allows the approximation of the posterior. It is considered only if \code{mod} is \code{"linear.logistic"} or \code{"Msubjective"}. In the other cases closed form evaluation of the posterior is available up to a proportionality constant. The default value is \code{nsim=10000.}
#' @param nsim a positive integer. \code{nsim} is the number of iterations for the Metropolis-within-Gibbs algorithm which allows the approximation of the posterior. It is considered only if \code{mod} is \code{"linear.logistic"} or \code{"Msubjective"}. In the other cases closed form evaluation of the posterior is available up to a proportionality constant. The default value is \code{nsim=10000.}
#' @param burnin a positive integer. \code{burnin} is the initial number of MCMC samples discarded. It is considered only if \code{mod} is \code{"linear.logistic"}  or \code{"Msubjective"}. The default value for \code{burnin} is \code{round(nsim/10).}
#' @param nsim.ML a positive integer. Whenever MCMC is needed \code{nsim.ML} is the number of iterations used in the marginal likelihood estimation procedure via power posterior method of Friel and Pettit (2008)
#' @param burnin.ML a positive integer. Whenever MCMC is needed \code{burnin.ML} is the initial number of samples discarded for marginal likelihood estimation via power-posterior approach. The default value is \code{burnin.ML} is \code{round(nsim/10)}.
#' @param num.t a positive integer. Whenever MCMC is needed \code{num.t} is the number of powers used in the power posterior approximation method for the marginal likelihood evaluation.
#' @param prior.N a character. \code{prior.N} is the label for the prior distribution for \eqn{N}. When \code{prior.N} is set to \code{"Rissanen"} (default) the Rissanen prior is used as a prior on \eqn{N}. This distribution was first proposed in Rissanen 1983 as a universal prior on integers. \code{prior.N="Uniform"} stands for a prior on \eqn{N} proportional to a constant value. \code{prior.N="one.over.N"} stands for a prior on \eqn{N} proportional to \eqn{1/N}. \code{prior.N="one.over.N2"} stands for a prior on \eqn{N} proportional to \eqn{1/N^2}.
#' @param which.mod a character. \code{which.mod} selects which models are fitted and compared. In the default setting \code{which.mod}=\code{"all"} all alternative models are fitted including new behavioural models based on alternative meaningful covariates (see Details). When \code{which.mod}=\code{"standard"} the function only fits classical behavioural models with either enduring effects as in \eqn{M_b}, \eqn{M_{c_1b}}, \eqn{M_{c_2b}} or ephemeral effects as in purely Markovian \eqn{M_{c_1}} and \eqn{M_{c_2}}
#' @param sort character. \code{sort} selects the order of models.
#'
#' @details
#' The available models are: \eqn{M_0}, \eqn{M_b}, \eqn{M_t}, \eqn{M_{c_1}}, \eqn{M_{c_1b}}, \eqn{M_{c_2}}, \eqn{M_{c_2b}}, \eqn{M_{mc}}, \eqn{M_{mc_{int}}}, \eqn{M_{mc_{count}}} and \eqn{M_{mc_{count.int}}}.
#' This function \code{BBRecap.all} can be computing intensive for high values of \code{neval} and \code{nsim}.
#'
#' @return A \code{data.frame} with one row corresponding to each model and the following columns:
#'
#' A dataframe with one row corresponding to each model and the following columns:
#'
#' model: model considered
#'
#' npar: number of parameters
#'
#' log.marginal.likelihood: log marginal likelihood
#'
#' Nhat: estimate of population size
#'
#' Ninf: lower \eqn{95 \%} highest posterior density interval
#'
#' Nsup: upper \eqn{95 \%} highest posterior density interval
#'
#' @references
#' Otis D. L., Burnham K. P., White G. C, Anderson D. R. (1978) Statistical Inference From Capture Data on Closed Animal Populations, Wildlife Monographs.
#'
#' Yang H.C., Chao A. (2005) Modeling animals behavioral response by Markov chain models for capture-recapture experiments, Biometrics 61(4), 1010-1017
#'
#' N. Friel and A. N. Pettitt. Marginal likelihood estimation via power posteriors. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 70(3):589, 607--2008
#'
#' Farcomeni A. (2011) Recapture models under equality constraints for the conditional capture probabilities. Biometrika 98(1):237--242
#'
#' Alunni Fegatelli, D. and Tardella, L. (2012) Improved inference on capture recapture models with behavioural effects. Statistical Methods & Applications Applications Volume 22, Issue 1, pp 45-66 10.1007/s10260-012-0221-4
#'
#' Alunni Fegatelli D. (2013) New methods for capture-recapture modelling with behavioural response and individual heterogeneity. http://hdl.handle.net/11573/918752
#'
#' @author Danilo Alunni Fegatelli and Luca Tardella
#'
#' @seealso \code{\link{BBRecap}}
#'
#' @examples
#'
#' \dontrun{
#' data(hornedlizard)
#' BBRecap.all(hornedlizard,neval=200)
#' }
#'
#' @keywords Behavioural_models Bayesian_inference
#' @export
BBRecap.all=function(data, last.column.count=FALSE, neval=1000, by.incr=1,nsim=10000,
  burnin=round(nsim/10),nsim.ML=500,burnin.ML=round(nsim.ML/10), num.t = 50,
  prior.N = c("Rissanen","Uniform","one.over.N","one.over.N2"),
  which.mod=c("all","standard"), sort=c("default","log.ML")){

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

	mod.M0=BBRecap(data.matrix,mod="M0",neval=neval,by.incr=by.incr,output="complete")
	mod.Mb=BBRecap(data.matrix,mod="Mb",neval=neval,by.incr=by.incr,output="complete")
	mod.Mlogistic=BBRecap(data.matrix,neval=neval,by.incr=by.incr,output="complete.ML",nsim=nsim,nsim.ML=nsim.ML,burnin=burnin,burnin.ML=burnin.ML,num.t=num.t,prior.N=prior.N)
	mod.Mc1=BBRecap(data.matrix,mod="Mc",markov.ord=1,neval=neval,by.incr=by.incr,output="complete")
	mod.Mc2=BBRecap(data.matrix,mod="Mc",markov.ord=2,neval=neval,by.incr=by.incr,output="complete")
	mod.Mc1b=BBRecap(data.matrix,mod="Mcb",markov.ord=1,neval=neval,by.incr=by.incr,output="complete")
	mod.Mc2b=BBRecap(data.matrix,mod="Mcb",markov.ord=2,neval=neval,by.incr=by.incr,output="complete")
    mod.Mt=BBRecap(data.matrix,mod="Mt",neval=neval,output="complete")

output=data.frame()
model=c("M0","Mb","Mlogistic","Mc1","Mc2","Mc1b","Mc2b","Mt")
npar=c(2,3,3,3,5,4,6,t+1)
log.marginal.likelihood=c(mod.M0$log.marginal.likelihood,mod.Mb$log.marginal.likelihood,mod.Mlogistic$log.marginal.likelihood,mod.Mc1$log.marginal.likelihood,mod.Mc2$log.marginal.likelihood,mod.Mc1b$log.marginal.likelihood,mod.Mc2b$log.marginal.likelihood,mod.Mt$log.marginal.likelihood)
Nhat.RMSE=c(mod.M0$N.hat.RMSE,mod.Mb$N.hat.RMSE,mod.Mlogistic$N.hat.RMSE,mod.Mc1$N.hat.RMSE,mod.Mc2$N.hat.RMSE,mod.Mc1b$N.hat.RMSE,mod.Mc2b$N.hat.RMSE,mod.Mt$N.hat.RMSE)
Nhat.mean=c(mod.M0$N.hat.mean,mod.Mb$N.hat.mean,mod.Mlogistic$N.hat.mean,mod.Mc1$N.hat.mean,mod.Mc2$N.hat.mean,mod.Mc1b$N.hat.mean,mod.Mc2b$N.hat.mean,mod.Mt$N.hat.mean)
Nhat.median=c(mod.M0$N.hat.median,mod.Mb$N.hat.median,mod.Mlogistic$N.hat.median,mod.Mc1$N.hat.median,mod.Mc2$N.hat.median,mod.Mc1b$N.hat.median,mod.Mc2b$N.hat.median,mod.Mt$N.hat.median)
Nhat.mode=c(mod.M0$N.hat.mode,mod.Mb$N.hat.mode,mod.Mlogistic$N.hat.mode,mod.Mc1$N.hat.mode,mod.Mc2$N.hat.mode,mod.Mc1b$N.hat.mode,mod.Mc2b$N.hat.mode,mod.Mt$N.hat.mode)
Ninf=c(mod.M0$HPD.N[1],mod.Mb$HPD.N[1],mod.Mlogistic$HPD.N[1],mod.Mc1$HPD.N[1],mod.Mc2$HPD.N[1],mod.Mc1b$HPD.N[1],mod.Mc2b$HPD.N[1],mod.Mt$HPD.N[1])
Nsup=c(mod.M0$HPD.N[2],mod.Mb$HPD.N[2],mod.Mlogistic$HPD.N[2],mod.Mc1$HPD.N[2],mod.Mc2$HPD.N[2],mod.Mc1b$HPD.N[2],mod.Mc2b$HPD.N[2],mod.Mt$HPD.N[2])

         output=data.frame(Model=model,npar=npar,log.marginal.likelihood=log.marginal.likelihood,Nhat.RMSE=Nhat.RMSE,Nhat.mean=Nhat.mean,Nhat.median=Nhat.median,Nhat.mode=Nhat.mode,Ninf=Ninf,Nsup=Nsup)

if(which.mod=="all"){

	mod.Mlogistic.integer=BBRecap(data.matrix,neval=neval,by.incr=by.incr,mbc.function="integer",output="complete.ML",nsim=nsim,nsim.ML=nsim.ML,burnin=burnin,burnin.ML=burnin.ML,num.t=num.t,prior.N=prior.N)
	mod.Mlogistic.counts=BBRecap(data.matrix,neval=neval,by.incr=by.incr,mbc.function="counts",output="complete.ML",nsim=nsim,nsim.ML=nsim.ML,burnin=burnin,burnin.ML=burnin.ML,num.t=num.t,prior.N=prior.N)
	mod.Mlogistic.counts.integer=BBRecap(data.matrix,neval=neval,by.incr=by.incr,mbc.function="counts.integer",output="complete.ML",nsim=nsim,nsim.ML=nsim.ML,burnin=burnin,burnin.ML=burnin.ML)



model=c("M0","Mb","Mlogistic","Mc1","Mc2","Mc1b","Mc2b","Mt","Mlogistic.integer","Mlogistic.counts","Mlogistic.counts.integer")
npar=c(2,3,3,3,5,4,6,t+1,rep(3,3))
log.marginal.likelihood=c(mod.M0$log.marginal.likelihood,mod.Mb$log.marginal.likelihood,mod.Mlogistic$log.marginal.likelihood,mod.Mc1$log.marginal.likelihood,mod.Mc2$log.marginal.likelihood,mod.Mc1b$log.marginal.likelihood,mod.Mc2b$log.marginal.likelihood,mod.Mt$log.marginal.likelihood,mod.Mlogistic.integer$log.marginal.likelihood,mod.Mlogistic.counts$log.marginal.likelihood,mod.Mlogistic.counts.integer$log.marginal.likelihood)
Nhat.RMSE=c(mod.M0$N.hat.RMSE,mod.Mb$N.hat.RMSE,mod.Mlogistic$N.hat.RMSE,mod.Mc1$N.hat.RMSE,mod.Mc2$N.hat.RMSE,mod.Mc1b$N.hat.RMSE,mod.Mc2b$N.hat.RMSE,mod.Mt$N.hat.RMSE,mod.Mlogistic.integer$N.hat.RMSE,mod.Mlogistic.counts$N.hat.RMSE,mod.Mlogistic.counts.integer$N.hat.RMSE)
Nhat.mean=c(mod.M0$N.hat.mean,mod.Mb$N.hat.mean,mod.Mlogistic$N.hat.mean,mod.Mc1$N.hat.mean,mod.Mc2$N.hat.mean,mod.Mc1b$N.hat.mean,mod.Mc2b$N.hat.mean,mod.Mt$N.hat.mean,mod.Mlogistic.integer$N.hat.mean,mod.Mlogistic.counts$N.hat.mean,mod.Mlogistic.counts.integer$N.hat.mean)
Nhat.median=c(mod.M0$N.hat.median,mod.Mb$N.hat.median,mod.Mlogistic$N.hat.median,mod.Mc1$N.hat.median,mod.Mc2$N.hat.median,mod.Mc1b$N.hat.median,mod.Mc2b$N.hat.median,mod.Mt$N.hat.median,mod.Mlogistic.integer$N.hat.median,mod.Mlogistic.counts$N.hat.median,mod.Mlogistic.counts.integer$N.hat.median)
Nhat.mode=c(mod.M0$N.hat.mode,mod.Mb$N.hat.mode,mod.Mlogistic$N.hat.mode,mod.Mc1$N.hat.mode,mod.Mc2$N.hat.mode,mod.Mc1b$N.hat.mode,mod.Mc2b$N.hat.mode,mod.Mt$N.hat.mode,mod.Mlogistic.integer$N.hat.mode,mod.Mlogistic.counts$N.hat.mode,mod.Mlogistic.counts.integer$N.hat.mode)
Ninf=c(mod.M0$HPD.N[1],mod.Mb$HPD.N[1],mod.Mlogistic$HPD.N[1],mod.Mc1$HPD.N[1],mod.Mc2$HPD.N[1],mod.Mc1b$HPD.N[1],mod.Mc2b$HPD.N[1],mod.Mt$HPD.N[1],mod.Mlogistic.integer$HPD.N[1],mod.Mlogistic.counts$HPD.N[1],mod.Mlogistic.counts.integer$HPD.N[1])
Nsup=c(mod.M0$HPD.N[2],mod.Mb$HPD.N[2],mod.Mlogistic$HPD.N[2],mod.Mc1$HPD.N[2],mod.Mc2$HPD.N[2],mod.Mc1b$HPD.N[2],mod.Mc2b$HPD.N[2],mod.Mt$HPD.N[2],mod.Mlogistic.integer$HPD.N[2],mod.Mlogistic.counts$HPD.N[2],mod.Mlogistic.counts.integer$HPD.N[2])

output=data.frame(Model=model,npar=npar,log.marginal.likelihood=log.marginal.likelihood,Model=model,npar=npar,log.marginal.likelihood=log.marginal.likelihood,Nhat.RMSE=Nhat.RMSE,Nhat.mean=Nhat.mean,Nhat.median=Nhat.median,Nhat.mode=Nhat.mode,Ninf=Ninf,Nsup=Nsup)


}


#if(ss=="npar"){
#    output=output[order(npar),]
#}


if(ss=="log.ML"){
    output=output[order(log.marginal.likelihood),]
}

print(paste("neval =",neval,sep=" "))
print(paste("M =",M,sep=" "))
return(output)

}
