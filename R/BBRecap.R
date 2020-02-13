#' Bayesian inference for capture-recapture analysis with emphasis on behavioural effect modelling
#'
#' Bayesian inference for a large class of discrete-time capture-recapture models under closed population with special emphasis on behavioural effect modelling including also the \emph{meaningful behavioral covariate} approach proposed in Alunni Fegatelli (2013) [PhD thesis]. Many of the standard classical models such as \eqn{M_0}, \eqn{M_b}, \eqn{M_{c_1}}, \eqn{M_t} or \eqn{M_{bt}} can be regarded as particular instances of the aforementioned approach. Other flexible alternatives can be fitted through a careful choice of a meaningful behavioural covariate and a possible partition of its admissible range.
#'
#' @usage
#' BBRecap (data,last.column.count=FALSE, neval = 1000, by.incr = 1,
#' mbc.function = c("standard","markov","counts","integer","counts.integer"),
#' mod = c("linear.logistic", "M0", "Mb", "Mc", "Mcb", "Mt", "Msubjective.cut",
#'         "Msubjective"), nsim = 5000, burnin = round(nsim/10),
#'         nsim.ML = 1000, burnin.ML = round(nsim.ML/10), num.t = 50,
#'         markov.ord=NULL, prior.N = c("Rissanen","Uniform","one.over.N","one.over.N2"),
#'         meaningful.mat.subjective = NULL, meaningful.mat.new.value.subjective = NULL,
#'         z.cut=NULL, output = c("base", "complete", "complete.ML"))
#'
#' @param data can be one of the following:
#' \enumerate{
#'   \item an \eqn{M} by \eqn{t} binary matrix/data.frame. In this case the input is interpreted as a matrix whose rows contain individual capture histories for all \eqn{M} observed units
#'   \item a matrix/data.frame with \eqn{(t+1)} columns. The first \eqn{t} columns contain binary entries corresponding to capture occurrences, while the last column contains non negative integers corresponding to frequencies. This format is allowed only when \code{last.column.count} is set to \code{TRUE}
#'   \item a \eqn{t}-dimensional array or table representing the counts of the \eqn{2^t} contingency table of binary outcomes
#' }
#' \eqn{M} is the number of units captured at least once and \eqn{t} is the number of capture occasions.
#' @param last.column.count a logical. In the default case \code{last.column.count=FALSE} each row of the input argument \code{data} represents the complete capture history for each observed unit. When \code{last.column.count} is set to \code{TRUE} in each row the first \eqn{t} entries represent one of the observed complete capture histories and the last entry in the last column is the number of observed units with that capture history
#' @param neval a positive integer. \code{neval} is the number of values of the population size \eqn{N} where the posterior is evaluated starting from \eqn{M}. The default value is \code{neval}=1000.
#' @param by.incr a positive integer. \code{nsim} is the number of iterations for the Metropolis-within-Gibbs algorithm which allows the approximation of the posterior. It is considered only if \code{mod} is \code{"linear.logistic"} or \code{"Msubjective"}. In the other cases closed form evaluation of the posterior is available up to a proportionality constant. The default value is \code{nsim=10000.}


#' @param mbc.function a character string with possible entries (see Alunni Fegatelli (2013) for further details)
#' \enumerate{
#'   \item \code{"standard"} meaningful behavioural covariate in [0,1] obtained through the normalized binary representation of integers relying upon partial capture history
#'   \item \code{"markov"} slight modification of \code{"standard"} providing consistency with arbitrary Markov order models when used in conjunction with the options \code{"Msubjective"} and \code{z.cut}.
#'   \item \code{"counts"} covariate in [0,1] obtained by normalizing the integer corresponding to the sum of binary entries i.e. the number of previous captures
#'   \item \code{"integer"} un-normalized integer corresponding to the binary entries of the partial capture history
#'   \item \code{"counts.integer"} un-normalized covariate obtained as the sum of binary entries i.e. the number of previous captures
#'   }
#' @param  mod a character. \code{mod} represents the behavioural model considered for the analysis. \code{mod="linear.logistic"} is the model proposed in Alunni Fegatelli (2013) based on the meaningful behavioural covariate. \code{mod="M0"} is the most basic model where no effect is considered and all capture probabilities are the same. \code{mod="Mb"} is the classical behavioural model where the capture probability varies only once when first capture occurs. Hence it represents an \emph{enduring} effect to capture. \code{mod="Mc"} is the \emph{ephemeral} behavioural Markovian model originally introduced in Yang and Chao (2005) and subsequently extended in Farcomeni (2011) and reviewed in Alunni Fegatelli and Tardella (2012) where capture probability depends only on the capture status (captured or uncaptured) in the previous \code{k=markov.ord} occasions. \code{mod="Mcb"} is an extension of Yang and Chao's model (2005); it considers both \emph{ephemeral} and \emph{enduring} effect to capture. \code{mod="Mt"} is the standard temporal effect with no behavioural effect. \code{mod="Msubjective.cut"} is an alternative behavioural model obtained through a specific cut on the meaningful behavioural covariate interpreted as memory effect. \code{mod="Msubjective"} is a customizable (subjective) behavioural model within the linear logistic model framework requiring the specification of the two additional arguments: the first one is \code{meaningful.mat.subjective} and contains an \eqn{M} by \eqn{t} matrix of ad-hoc meaningful covariates depending on previous capture history; the second one is \code{meaningful.mat.new.value.subjective} and contains a vector of length \eqn{t} corresponding to meaningful covariates for a generic uncaptured unit. The default value for \code{mod} is \code{"linear.logistic".}
#' @param nsim a positive integer. \code{nsim} is the number of iterations for the Metropolis-within-Gibbs algorithm which allows the approximation of the posterior. It is considered only if \code{mod} is \code{"linear.logistic"} or \code{"Msubjective"}. In the other cases closed form evaluation of the posterior is available up to a proportionality constant. The default value is \code{nsim=10000.}
#' @param burnin a positive integer. \code{burnin} is the initial number of MCMC samples discarded. It is considered only if \code{mod} is \code{"linear.logistic"}  or \code{"Msubjective"}. The default value for \code{burnin} is \code{round(nsim/10).}
#' @param nsim.ML a positive integer. Whenever MCMC is needed \code{nsim.ML} is the number of iterations used in the marginal likelihood estimation procedure via power posterior method of Friel and Pettit (2008)
#' @param burnin.ML a positive integer. Whenever MCMC is needed \code{burnin.ML} is the initial number of samples discarded for marginal likelihood estimation via power-posterior approach. The default value is \code{burnin.ML} is \code{round(nsim/10)}.
#' @param num.t a positive integer. Whenever MCMC is needed \code{num.t} is the number of powers used in the power posterior approximation method for the marginal likelihood evaluation.
#' @param markov.ord a positive integer. \code{markov.ord} is the order of Markovian model \eqn{M_c} or \eqn{M_{cb}}. It is considered only if \code{mod="Mc"} or \code{mod="Mcb"}.
#' @param prior.N a character. \code{prior.N} is the label for the prior distribution for \eqn{N}. When \code{prior.N} is set to \code{"Rissanen"} (default) the Rissanen prior is used as a prior on \eqn{N}. This distribution was first proposed in Rissanen 1983 as a universal prior on integers. \code{prior.N="Uniform"} stands for a prior on \eqn{N} proportional to a constant value. \code{prior.N="one.over.N"} stands for a prior on \eqn{N} proportional to \eqn{1/N}. \code{prior.N="one.over.N2"} stands for a prior on \eqn{N} proportional to \eqn{1/N^2}.
#' @param meaningful.mat.subjective \code{M x t} matrix containing  numerical covariates to be used for a customized logistic model approach
#' @param meaningful.mat.new.value.subjective \code{1 x t} numerical vector corresponding to auxiliary covariate to be considered for unobserved unit
#' @param z.cut numeric vector. \code{z.cut} is a vector containing the cut point for the memory effect covariate. It is considered only if \code{mod="Msubjective.cut"}
#' @param output a character. \code{output} selects the kind of output from a very basic summary info on the posterior output (point and interval estimates for the unknown \eqn{N}) to more complete details including MCMC simulations for all parameters in the model when appropriate.
#'
#' @details
#' Independent uniform distributions are considered as default prior for the nuisance parameters. If \code{model="linear.logistic"} or \code{model="Msubjective"} and \code{output="complete.ML"}  the marginal likelihood estimation is performed through the \emph{power posteriors method} suggested in Friel and Pettit (2008). In that case the \code{BBRecap} procedure is computing intensive for high values of \code{neval} and \code{nsim}.
#'
#' @return
#'   \item{Model}{model considered}
#'   \item{Prior}{prior distribution for \eqn{N}}
#'   \item{N.hat.mean}{posterior mean for \eqn{N}}
#'   \item{N.hat.median}{posterior median for \eqn{N}}
#'   \item{N.hat.mode}{posterior mode for \eqn{N}}
#'   \item{N.hat.RMSE}{minimizer of a specific loss function connected with the Relative Mean Square Error}
#'   \item{HPD.N}{\eqn{95 \%} highest posterior density interval estimate for \eqn{N}}
#'   \item{log.marginal.likelihood}{log marginal likelihood}
#'   \item{N.range}{values of N considered}
#'   \item{posterior.N}{values of the posterior distribution for each N considered}
#'   \item{z.matrix}{meaningful behavioural covariate matrix for the observed data}
#'   \item{vec.cut}{cut point used to set up meaningful partitions the set of the partial capture histories according to the value of the value of the meaningful behavioural covariate}
#'   \item{N.vec}{simulated values from the posterior marginal distribution of N}
#'   \item{mean.a0}{posterior mean of the parameter a0}
#'   \item{hpd.a0}{highest posterior density interval estimate of the parameter a0}
#'   \item{a0.vec}{simulated values from the posterior marginal distribution of a0}
#'   \item{mean.a1}{posterior mean of the parameter a1}
#'   \item{hpd.a1}{highest posterior density interval estimate of the parameter a1}
#'   \item{a1.vec}{simulated values from the posterior marginal distribution of a1}
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
#' @seealso \code{\link{BBRecap.custom.part}}, \code{\link{LBRecap}}
#'
#' @examples
#'
#' \dontrun{
#' data(greatcopper)
#'
#' mod.Mb=BBRecap(greatcopper,mod="Mb")
#' str(mod.Mb)
#' }
#'
#' @keywords Behavioural_models Bayesian_inference
#' @export
BBRecap=function(
  data,
  last.column.count=FALSE,
  neval=1000,
  by.incr=1,
  mbc.function=c("standard","markov","counts","integer","counts.integer"),
  mod=c("linear.logistic","M0","Mb","Mc","Mcb","Mt","Msubjective.cut","Msubjective"),
  nsim=5000,
  burnin=round(nsim/10),
  nsim.ML=1000,
  burnin.ML=round(nsim.ML/10),
  num.t=50,
  markov.ord=NULL,
  prior.N=c("Rissanen","Uniform","one.over.N","one.over.N2"),
  meaningful.mat.subjective=NULL,
  meaningful.mat.new.value.subjective=NULL,
  z.cut=NULL,
  output=c("base","complete","complete.ML")){

          log.marg.likelihood=NULL

  mod=match.arg(mod)
  prior.N=match.arg(prior.N)
  output=match.arg(output)
  mbc.function=match.arg(mbc.function)

  model=switch(mod,
               linear.logistic="mod.linear.logistic",
               M0="mod.M0",
               Mb="mod.Mb",
               Mc="mod.Mc",
               Mcb="mod.Mcb",
               Mt="mod.Mt",
               Msubjective.cut="mod.Msubjective.cut",
               Msubjective="mod.Msubjective")

  prior=switch(prior.N,
               Rissanen="Rissanen.prior",
               Uniform="Uniform.prior",
               one.over.N="1overN.prior",
               one.over.N2="1overN^2.prior")

  qb.function=switch(mbc.function,
  					standard="qb.standard",
  					markov="qb.markov",
  					integer="qb.integer",
  					counts="qb.count",
  					counts.integer="qb.count.integer")

  if((mod!="Msubjective.cut") & length(z.cut)>0){
 	warning("z.cut without mod='Msubjective.cut' will be IGNORED!")
}
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
data.matrix=rbind(dd,data.matrix)

}

data.matrix=data.matrix[1:(nrow(data.matrix)-1),]

  }

  if(any(data.matrix!=0 & data.matrix!=1)) stop("data must be a binary matrix")

   if(sum(apply(data.matrix,1,sum)==0)){
 	warning("input data argument contains rows with all zeros: these rows will be removed and ignored")

        data.matrix=data.matrix[apply(data.matrix,1,sum)!=0,]
    }

  t=ncol(data.matrix)
  M=nrow(data.matrix)
  logistic=0

  if(model=="mod.Mc" | model=="mod.Mcb"){

    if (!is.numeric(markov.ord) || markov.ord != round(markov.ord) || markov.ord <= 0 || markov.ord > t)
      stop("The Markov order must be a non-negative integer smaller than t!")

    mbc.function="markov"
    z.cut=seq(0,1,by=1/(2^markov.ord))
    z.cut[1]=-0.1
    if(model=="mod.Mcb"){
      z.cut=c(z.cut[1],0,z.cut[-1])
    }
  }

if(qb.function=="qb.standard") recodebinary=quant.binary
if(qb.function=="qb.integer") recodebinary=quant.binary.integer
if(qb.function=="qb.markov") recodebinary=quant.binary.markov
if(qb.function=="qb.count") recodebinary=quant.binary.counts
if(qb.function=="qb.count.integer") recodebinary=quant.binary.counts.integer

#####
  v=c()

 meaningful.mat=matrix(NA,ncol=ncol(data.matrix),nrow=M)

  for(i in 1:ncol(data.matrix)){
    for(j in 1:M){
      v=data.matrix[j,1:i]
      v=paste(v,collapse="")
      if(qb.function=="qb.markov"){
      		meaningful.mat[j,i]=recodebinary(v,markov.ord)}
      else{meaningful.mat[j,i]=recodebinary(v)}
    }
  }

  meaningful.mat=cbind(rep(0,M),meaningful.mat[,1:t-1])
  meaningful.mat.new.value=rep(0,t*by.incr)

  if(model=="mod.M0"){
          z.cut=c(-0.1,1)
        }

  if(model=="mod.Mb"){
          z.cut=c(-0.1,0,1)
        }

  if(model=="mod.Mt"){
          z.cut=seq(0,t)
          meaningful.mat=matrix(rep(1:t),nrow=M,ncol=t,byrow=T)
          meaningful.mat.new.value=1:t
          }

  if(model=="mod.Msubjective"){
    if(is.matrix(meaningful.mat.subjective)){
      if(ncol(meaningful.mat.subjective)!=ncol(data.matrix)||
           nrow(meaningful.mat.subjective)!=nrow(data.matrix))
        stop("meaningful.mat.subjective has different dimensions from data matrix")
      meaningful.mat=meaningful.mat.subjective
    }
    else{ stop("meaningful.mat.subjective must be a numerical matrix") }
    if(is.numeric(meaningful.mat.new.value.subjective)){
      meaningful.mat.vec.new=meaningful.mat.new.value.subjective
    }
  }

  meaningful.mat.vec=c(meaningful.mat)
  y=c(data.matrix)

  nn.1=c()
  nn.0=c()

  coeff.vec1=c()
  coeff.vec2=c()

  log.post.N=c()
  post.N=c()
  vv=c()

  prior.inv.const=0

  prior.distr.N=function(x){

    if(prior=="Rissanen.prior") {out=(rissanen(x))}
    if(prior=="Uniform.prior") {out=(1/(x^0))}
    if(prior=="1overN.prior") {out=(1/(x^1))}
    if(prior=="1overN^2.prior") {out=(1/(x^2))}

    return(out)
  }

  if(model=="mod.linear.logistic" | model=="mod.Msubjective"){logistic=1}

  if(logistic==0){

            if(output=="complete.ML"){
        output="complete"
       }

    for(k in 1:neval){
 meaningful.mat.new=rbind(meaningful.mat,matrix(rep(meaningful.mat.new.value,(k>1)),nrow=((k-1) *by.incr),ncol=t,byrow=T))
 dati.new=rbind(data.matrix,matrix(0,nrow=((k-1)*by.incr),ncol=t))
 meaningful.mat.vec.new=c(meaningful.mat.new)
 y.new=c(dati.new)

        for(j in 1:(length(z.cut)-1)){

          if(j==1){nn.0[j]=sum(meaningful.mat.vec.new>=z.cut[j] &
          	meaningful.mat.vec.new<=z.cut[j+1] & y.new==0)}
          else{nn.0[j]=sum(meaningful.mat.vec.new>z.cut[j] &
          	meaningful.mat.vec.new<=z.cut[j+1] & y.new==0)}
          if(j==1){nn.1[j]=sum(meaningful.mat.vec.new>=z.cut[j] &
          	meaningful.mat.vec.new<=z.cut[j+1] & y.new==1)}
          else{nn.1[j]=sum(meaningful.mat.vec.new>z.cut[j] &
          	meaningful.mat.vec.new<=z.cut[j+1] & y.new==1)}
        }

        prior.inv.const=prior.inv.const+ prior.distr.N((M+k-1))

        log.post.N[k]=lchoose((sum(nn.1)+sum(nn.0))/t,M)+sum(lbeta((nn.1+1),(nn.0+1)))+log(prior.distr.N((M+k-1)))

      }

     # closed k cycle

      l.max=max(log.post.N)

      for(k in 1:neval){
        vv[k]<-exp(log.post.N[k]-l.max)
      }

      ss=sum(vv)
      post.N=vv/ss

	  if(output=="complete"){
      if(model=="mod.Mt"){
      	Exp.beta=matrix(ncol=(length(z.cut)-1),nrow=neval)
      	for(i in 1:neval){
      	Exp.beta[i,]=(nn.1+1)/(nn.1+nn.0-(neval-i)+2)
 	 	}
      	pH.post.mean=apply((post.N*Exp.beta),2,sum)
      }
      else{
      	Exp.beta=matrix(ncol=(length(z.cut)-1),nrow=neval)
      	for(i in 1:neval){
      	Exp.beta[i,]=(nn.1+1)/(nn.1+nn.0+2)
 	 	}
  	 	rr=sort(0:(neval-1)*t,decreasing=T)
 	 	Exp.beta[,1]=Exp.beta[,1]*(nn.1[1]+nn.0[1]+2)/(nn.1[1]+nn.0[1]-rr+2)
      	pH.post.mean=apply((post.N*Exp.beta),2,sum)
      	}

    names(pH.post.mean)=paste("pH",1:(length(z.cut)-1),sep="")
      }

    val=seq(M,((M+(neval-1)*by.incr)),by=by.incr)

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

  }

  ## ARMS

  if(model=="mod.linear.logistic" | model=="mod.Msubjective"){logistic=1}

  if(logistic==1){

    print("implementation via ARMS")

    ### LOG-PRIOR N

    log.prior.N=function(x){

      if(prior=="Rissanen.prior") {out=log(rissanen(x))}
      if(prior=="Uniform.prior") {out=0}
      if(prior=="1overN.prior") {out=log(1/(x^1))}
      if(prior=="1overN^2.prior") {out=log(1/(x^2))}

      return(out)
    }

    ### LOG-LIKELIHOOD

    if(logistic==1){

      log.likelihood=function(NN,aa0,aa1){

        z=matrix(0,ncol=t,nrow=NN)
        x=matrix(0,ncol=t,nrow=NN)

        z[(1:M),]=meaningful.mat
        x[(1:M),]=data.matrix
        l = x*log(exp(aa0+aa1*z)/(1+exp(aa0+aa1*z)))-(1-x)*log((1+exp(aa0+aa1*z)))
        out=sum(l)	+lchoose(NN,M)
        return(out)
      }
      ###

      N.vec=c(M,rep(0,(nsim-1)))
      a0.vec=c(0,rep(0,(nsim-1)))
      a1.vec=c(0,rep(0,(nsim-1)))

      ### MODIFIED-LOG-LIKELIHOOD

      log.likelihood.t=function(NN,aa0,aa1,tt){

        z=matrix(0,ncol=t,nrow=NN)
        x=matrix(0,ncol=t,nrow=NN)

        z[(1:M),]=meaningful.mat
        x[(1:M),]=data.matrix
        l = x*log(exp(aa0+aa1*z)/(1+exp(aa0+aa1*z)))-(1-x)*log((1+exp(aa0+aa1*z)))
        out=tt*(sum(l)	+lchoose(NN,M))
        return(out)
      }
      ###

      for(g in 2:(nsim+1)){

        N.vec[g]=arms(
          N.vec[g-1],
          function(N) log.prior.N(N)+log.likelihood(N,a0.vec[g-1],a1.vec[g-1]),
          function(N) (N>=M)*(N<=M+neval),
          1
        )
        N.vec[g]=round(N.vec[g])

        a0.vec[g]=arms(
          a0.vec[g-1],
          function(a0) log.likelihood(N.vec[g],a0,a1.vec[g-1]),
          function(a0) (a0>-5)*(a0<=5),
          1
        )

        a1.vec[g]=arms(
          a1.vec[g-1],
          function(a1) log.likelihood(N.vec[g],a0.vec[g],a1),
          function(a1) (a1>-5)*(a1<=5),
          1
        )

      }

      mean.N=round(mean(N.vec[-(1:burnin)]))
      mean.a0=mean(a0.vec[-(1:burnin)])
      mean.a1=mean(a1.vec[-(1:burnin)])
      median.N=median(N.vec[-(1:burnin)])

      mode.N=as.numeric(names(sort(table(N.vec[-(1:burnin)]))[length(table(N.vec[-(1:burnin)]))]))

      val=seq(M,(M+neval))

      post.N=(table(c(N.vec,M:(M+neval)))-1)/sum(N.vec)

      funzioneRMSE=function(x){
        sum((((x/val)-1)^2)*post.N)
      }

      estimate.N=round(optimize(funzioneRMSE,c(min(val),max(val)))$minimum)
      HPD.interval=function(x.vec){
        obj=density(x=x.vec,from=min(x.vec),to=max(x.vec),n=diff(range(x.vec))+1,kernel="rectangular")
        density.estimation.N=obj$y
        names(density.estimation.N)=obj$x
        ordered.x.vec=x.vec[order(density.estimation.N[as.character(x.vec)],decreasing=T)]
        HPD.interval=range(ordered.x.vec[1:round(quantile(1:length(x.vec),prob=c(0.95)))])
        return(HPD.interval)
      }

      HPD.N=HPD.interval(N.vec[-(1:burnin)])
      HPD.a0=HPD.interval(a0.vec[-(1:burnin)])
      HPD.a1=HPD.interval(a1.vec[-(1:burnin)])

      if(output=="complete.ML"){

        output="complete"

        a=seq(0.0016,1,length=num.t)

        a=a^4

        l=matrix(ncol=length(a),nrow=nsim.ML)

        for(j in 1:length(a)){

  print(paste("ML evaluation via power-posterior method of Friel & Pettit (2008) ; t =",j,"of",length(a)))

          N.vec.ML=c(M,rep(0,(nsim-1)))
          a0.vec.ML=c(0,rep(0,(nsim-1)))
          a1.vec.ML=c(0,rep(0,(nsim-1)))

          for(g in 2:(nsim.ML+1)){

            N.vec.ML[g]=arms(
              N.vec.ML[g-1],
              function(N) log.prior.N(N)+log.likelihood.t(N,a0.vec.ML[g-1],a1.vec.ML[g-1],a[j]),
              function(N) (N>=M)*(N<=M+neval),
              1
            )

            N.vec.ML[g]=round(N.vec.ML[g])

            a0.vec.ML[g]=arms(
              a0.vec.ML[g-1],
              function(a0) log.likelihood.t(N.vec.ML[g],a0,a1.vec.ML[g-1],a[j]),
              function(a0) (a0>-5)*(a0<=5),
              1
            )

            a1.vec.ML[g]=arms(
              a1.vec.ML[g-1],
              function(a1) log.likelihood.t(N.vec.ML[g],a0.vec.ML[g],a1,a[j]),
              function(a1) (a1>-5)*(a1<=5),
              1
            )

            l[g-1,j]=log.likelihood(N.vec.ML[g],a0.vec.ML[g],a1.vec.ML[g])

          }
        }

        ll=apply(l,2,mean)

        v=c()
        for(k in 1:(length(a)-1)){
          v[k]=(a[k+1]-a[k])*(ll[k+1]+ll[k])/2
        }

        log.marg.likelihood=sum(v)

      }
    }
  }

  if(logistic==0){
    if(mod=="M0" | mod=="Mb" | mod=="Mt" | mod=="Mc" | mod=="Mcb"){
    out=switch(output,
               base=
                 list(Model=model,prior.N=prior.N,N.hat.RMSE=estimate.N,HPD.N=c(inf.lim,sup.lim)),
                 complete=
                 list(Model=model,prior.N=prior.N,N.hat.mean=mean.N,N.hat.median=median.N,N.hat.mode=mode.N,N.hat.RMSE=estimate.N,HPD.N=c(inf.lim,sup.lim),pH.post.mean=pH.post.mean,                    N.range=val,posterior.N=post.N,z.matrix=meaningful.mat,vec.cut=z.cut,log.marginal.likelihood=log.marg.likelihood)
                 )
  }
 else{
      out=switch(output,
               base=
                 list(Model=model,mbc.function=mbc.function,prior.N=prior.N,N.hat.RMSE=estimate.N,HPD.N=c(inf.lim,sup.lim)),
                 complete=
                 list(Model=model,mbc.function=mbc.function,prior.N=prior.N,N.hat.mean=mean.N,N.hat.median=median.N,N.hat.mode=mode.N,N.hat.RMSE=estimate.N,HPD.N=c(inf.lim,sup.lim),pH.post.mean=pH.post.mean, N.range=val,posterior.N=post.N,z.matrix=meaningful.mat,vec.cut=z.cut,log.marginal.likelihood=log.marg.likelihood
                 ) )}
}

  if(logistic==1) {
   if(mod=="Msubjective"){
    out=switch(output,
               base=
                 list(Model=model,prior.N=prior.N,N.hat.RMSE=estimate.N,HPD.N=HPD.N),
               complete=
                 list(Model=model,prior.N=prior.N,N.hat.mean=mean.N,N.hat.median=median.N,N.hat.mode=mode.N,N.hat.RMSE=estimate.N,HPD.N=HPD.N,N.vec=N.vec,mean.a0=mean.a0,HPD.a0=HPD.a0,a0.vec=a0.vec,mean.a1=mean.a1,HPD.a1=HPD.a1,a1.vec=a1.vec,z.matrix=meaningful.mat,log.marginal.likelihood=log.marg.likelihood))
    }

    else{
    	out=switch(output,
               base=
                 list(Model=model,mbc.function=mbc.function,prior.N=prior.N,N.hat.RMSE=estimate.N,HPD.N=HPD.N),
               model=
                 list(Model=model,mbc.function=mbc.function,prior.N=prior.N,N.hat.mean=mean.N,N.hat.median=median.N,N.hat.mode=mode.N,N.hat.RMSE=estimate.N,HPD.N=HPD.N,N.vec=N.vec,mean.a0=mean.a0,HPD.a0=HPD.a0,a0.vec=a0.vec,mean.a1=mean.a1,HPD.a1=HPD.a1,a1.vec=a1.vec,z.matrix=meaningful.mat),
               complete=
                 list(Model=model,mbc.function=mbc.function,prior.N=prior.N,N.hat.mean=mean.N,N.hat.median=median.N,N.hat.mode=mode.N,N.hat.RMSE=estimate.N,HPD.N=HPD.N,N.vec=N.vec,mean.a0=mean.a0,HPD.a0=HPD.a0,a0.vec=a0.vec,mean.a1=mean.a1,HPD.a1=HPD.a1,a1.vec=a1.vec,z.matrix=meaningful.mat,log.marginal.likelihood=log.marg.likelihood))}

  }

  return(out)

}
