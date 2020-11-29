## ----echo=FALSE---------------------------------------------------------------
library("mvtnorm")


## ----child = "latexmacros.Rmd"------------------------------------------------




## ----child = "Sections/prelim.Rmd"--------------------------------------------

## -----------------------------------------------------------------------------
    set.seed(1234)
    xs=rnorm(1000)
    mean(xs) ## should be 0
    var(xs) ## should be 1
    mean(xs<=0) ## P(X \le 0) should be 1/2
    plot(density(xs))
    x<-seq(-3,3,len=10000)
    lines(x,dnorm(x),col="red")
    quantile(xs,prob=c(0.5,0.975)) ## should be 0 and 1.96


## -----------------------------------------------------------------------------
nt=20
probs=rep(0,nt+1) ## including the probability of dry at time 0
thisprob=0;
probs[1]=thisprob
for (i in 1:nt) {
    thisprob=thisprob*0.75+(1-thisprob)*0.4
    probs[i+1]=thisprob
}
plot(0:nt,probs,pch="x",xlab="n",ylab="P(Xn=dry)",ylim=c(0,1))
abline(h=8/13,col="red")


## ----echo=FALSE---------------------------------------------------------------
P1=matrix(nrow=4,ncol=4,byrow=TRUE,data=c(0.5, 0.5, 0, 0,
                                         0.3, 0.7, 0, 0,
					 0, 0, 0.8, 0.2,
                                         0, 0.0, 0, 1))
P2=matrix(nrow=4,ncol=4,byrow=TRUE,data=c(0.0, 0.5, 0.0, 0.5,
                                         0.3, 0.0, 0.7, 0,
					 0, 0.25, 0.0, 0.75,
                                         0.8, 0.0, 0.2, 0))
nt=10
ps1=rep(0,nt+1); ps2=ps1
x=c(1,0,0,0)
ps1[1]=x[1]
for (i in 1:nt) {
  x=t(t(x) %*% P1)
  ps1[i+1]=x[1]
}
x=c(1,0,0,0)
ps2[1]=x[1]
for (i in 1:nt) {
  x=t(t(x) %*% P2)
  ps2[i+1]=x[1]
}

par(mfrow=c(1,3))
plot(0:nt,ps1,ylim=c(0,1),xlab="n",ylab="P(X=1)",type="b",main="Mx 1 started at X0=1")
plot(0:nt,rep(0,nt+1),ylim=c(0,1),xlab="n",ylab="P(X=1)",type="b",main="Mx 1 started at X0=3")
plot(0:nt,ps2,ylim=c(0,1),xlab="n",ylab="P(X=1)",type="b",main="Mx 2 started at X0=1")


## -----------------------------------------------------------------------------
set.seed(5431)
x=10
nt=200
rho=0.8
xs=rep(0,nt+1)
xs[1]=x
for (i in 1:nt) {
    x=rho*x+sqrt(1-rho^2)*rnorm(1)
    xs[i+1]=x
}
par(mfrow=c(1,2))
plot(0:nt,xs,xlab="t",ylab="x")
plot(0:nt,xs,type="l",xlab="t",ylab="x")



## ----child = "Sections/MH1.Rmd"-----------------------------------------------

## ----echo=TRUE----------------------------------------------------------------
library("mvtnorm")


## ----child = "../../Generic_R/generic_indepMVN.Rmd", echo=TRUE----------------

## -----------------------------------------------------------------------------
## MVN independence sampler
## qV has to be a matrix,
## in dimensions d=1 it is a 1x1 matrix.
## ... are any extra arguments, these are
## passed straight through to the log.pi function.

indep.samp.MVN<-function(nits,theta0,log.pi,qmode,qV,...) {
    d=length(theta0)
    A=t(chol(qV))
    store=matrix(nrow=nits+1,ncol=d+1)
    psis=matrix(nrow=nits,ncol=d)
	
    nacc=0
    theta.curr=theta0
    log.pi.curr=log.pi(theta.curr,...)
    log.q.curr=dmvnorm(theta.curr,qmode,qV,log=TRUE) 
    store[1,]=c(theta.curr,log.pi.curr)

    for (i in 1:nits) {
        psi=qmode+A%*%rnorm(d); psis[i,]=psi
        log.pi.prop=log.pi(psi,...)
        log.q.prop=dmvnorm(psi,qmode,qV,log=TRUE) 
        log.alpha=log.pi.prop+log.q.curr-log.pi.curr-log.q.prop
        if (log(runif(1))<log.alpha) {
            theta.curr=psi
            log.pi.curr=log.pi.prop
            log.q.curr=log.q.prop
            nacc=nacc+1
        }
        store[i+1,]=c(theta.curr,log.pi.curr)
    }
    return(list(acc=nacc/nits,store=store,psis=psis))
}



## ----echo=TRUE----------------------------------------------------------------
## N(0,I) log density
log.dens.N0I<-function(theta) {
    return(sum(dnorm(theta,log=TRUE)))			    
}


## ----echo=TRUE----------------------------------------------------------------
    set.seed(12341)
    mcoutIS1=indep.samp.MVN(500,c(0),log.dens.N0I,c(0),matrix(nrow=1,ncol=1,data=4))
    par(mfrow=c(2,2))
    xs=seq(-6,6,len=500)
    plot(xs,dnorm(xs),lwd=2,type="l",ylab="Density",xlab=" ")
    lines(xs,dnorm(xs,sd=2),lwd=2,col="red",lty=2)
    plot(0:100,mcoutIS1$store[1:101,1],type="l",xlab="iteration",ylab="theta")
    plot(density(mcoutIS1$store[1:101,1]),main="Density from chain",xlim=c(-6,6),xlab="")
    allpts=c(mcoutIS1$store[1:101,1],mcoutIS1$psis[1:100])
    lo=min(allpts); hi=max(allpts)
    plot(0:100,mcoutIS1$store[1:101,1],type="l",xlab="iteration",ylab="theta",ylim=c(lo,hi))
    points(1:100,mcoutIS1$psis[1:100],pch="x",col="red")


## ----child = "../../Generic_R/generic_rwm.Rmd", echo=TRUE---------------------

## -----------------------------------------------------------------------------
rwm=function(nits,theta0,log.pi,lambda,qV,...) {
    d=length(theta0)
    if (is.null(qV)) {
        A=diag(d)
    }
    else {
        A=t(chol(qV))
    }
    
    store=matrix(nrow=nits+1,ncol=d+1)
    psis=matrix(nrow=nits,ncol=d)
	
    nacc=0
    theta.curr=theta0
    log.pi.curr=log.pi(theta.curr,...)
    store[1,]=c(theta.curr,log.pi.curr)
    for (i in 1:nits) {
        psi=theta.curr+lambda*A%*%rnorm(d); psis[i,]=psi
        log.pi.prop=log.pi(psi,...)
        log.alpha=log.pi.prop-log.pi.curr
        if (log(runif(1))<log.alpha) {
            theta.curr=psi
            log.pi.curr=log.pi.prop
            nacc=nacc+1
        }
        store[i+1,]=c(theta.curr,log.pi.curr)
    }
    
    return(list(acc=nacc/nits,store=store,psis=psis))
}




## ----echo=TRUE----------------------------------------------------------------
    set.seed(12341)
    mcoutRWM1=rwm(1000,c(0),log.dens.N0I,2,NULL)
    par(mfrow=c(2,2))
    xs=seq(-6,6,len=500)
    plot(xs,dnorm(xs),lwd=2,type="l",ylab="Density",xlab=" ")
    lines(xs,dnorm(xs,1,sd=2),lwd=2,col="red",lty=2)
    points(c(1),c(0),pch="+")
    plot(0:100,mcoutRWM1$store[1:101,1],type="l",xlab="iteration",ylab="theta")
    plot(density(mcoutRWM1$store[1:101,1]),main="Density from chain",xlim=c(-6,6),xlab="")
    allpts=c(mcoutRWM1$store[1:101,1],mcoutRWM1$psis)
    lo=min(allpts); hi=max(allpts)
    plot(0:100,mcoutRWM1$store[1:101,1],type="l",xlab="iteration",ylab="theta",ylim=c(lo,hi))
    points(1:100,mcoutRWM1$psis[1:100],pch="x",col="red")
    library(coda)
    essRWM1=effectiveSize(mcoutRWM1$store[,1])


## ----echo=TRUE----------------------------------------------------------------
    set.seed(123412)
    mcoutRWM2=rwm(500,c(20),log.dens.N0I,2,NULL)
    par(mfrow=c(1,2))
    xs=seq(-6,6,len=500)
    plot(0:100,mcoutRWM2$store[1:101,1],type="l",xlab="iteration",ylab="theta",main="Trace plot")
    brks=c(-5,-2,-1,0,1,2,5,10,15,20)
    hist(mcoutRWM2$store[1:101,1],prob=TRUE,main="Histogram of chain",xlab="",ylim=c(0,max(dnorm(xs))),breaks=brks)
    lines(xs,dnorm(xs,0,sd=1),lwd=2,col="blue",lty=2)


## ----echo=TRUE----------------------------------------------------------------
    par(mfrow=c(1,2))
    hist(mcoutRWM2$store[32:101,1],prob=TRUE,main="Histogram after burn in",xlab="")
    lines(xs,dnorm(xs,0,sd=1),lwd=2,col="blue",lty=2)
    hist(mcoutRWM2$store[32:501,1],prob=TRUE,main="Histogram (500) after burn in",xlab="")
    lines(xs,dnorm(xs,0,sd=1),lwd=2,col="blue",lty=2)


## ----echo=TRUE----------------------------------------------------------------
    set.seed(123412)
    mcoutRWM3=rwm(1000,c(0),log.dens.N0I,0.3,NULL)
    par(mfrow=c(2,2))
    xs=seq(-6,6,len=500)
    plot(0:500,mcoutRWM3$store[1:501,1],type="l",xlab="iteration",ylab="theta",main="RWM Trace plot")
    hist(mcoutRWM3$store[1:501,1],prob=TRUE,main="Histogram of RWM chain",xlab="",xlim=c(-2,2))
    lines(xs,dnorm(xs,0,sd=1),lwd=2,col="blue",lty=2)
    mcoutIS1b=indep.samp.MVN(500,c(0),log.dens.N0I,c(0),matrix(nrow=1,ncol=1,data=100))
    plot(0:500,mcoutIS1b$store[,1],type="l",xlab="iteration",ylab="theta",main="IS Trace plot")
    hist(mcoutIS1b$store[,1],prob=TRUE,main="Histogram of IS chain",xlab="",xlim=c(-2,2))
    lines(xs,dnorm(xs,0,sd=1),lwd=2,col="blue",lty=2)    


## ----echo=TRUE----------------------------------------------------------------
    par(mfrow=c(2,2))
    acf(mcoutRWM2$store[31:500,1],main="lambda=2 (RWM)",lag.max=25)
    acf(mcoutRWM3$store[1:501,1],main="lambda=0.3 (RWM)",lag.max=25)
    acf(mcoutIS1$store[31:500,1],main="sigmaq=2 (IS)",lag.max=25)
    acf(mcoutIS1b$store[,1],main="sigmaq=10 (IS)",lag.max=25)
    essRWM2=effectiveSize(mcoutRWM2$store[31:500,1])
    essRWM3=effectiveSize(mcoutRWM3$store[1:501,1])
    essIS1=effectiveSize(mcoutIS1$store[,1])
    essIS1b=effectiveSize(mcoutIS1b$store[,1])


## ----echo=TRUE----------------------------------------------------------------
    set.seed(123413)
    log.dens.logistic01=function(x) {
        return(sum(dlogis(x,log=TRUE)))
    }   
    mcoutIS2=indep.samp.MVN(1000,c(7),log.dens.logistic01,c(0),matrix(nrow=1,ncol=1,data=2))
    par(mfrow=c(2,2))
    xs=seq(-7,7,len=200)
    plot(xs,dlogis(xs),lwd=2,type="l",ylim=c(0,max(c(dlogis(xs),dnorm(xs,sd=sqrt(2))))),xlab="theta",ylab="Density")
    lines(xs,dnorm(xs,sd=sqrt(2)),col="red",lwd=2)
    plot(0:1000,mcoutIS2$store[,1],type="l",ylab="theta",xlab="iteration")
    a=density(mcoutIS2$store[,1])
    plot(a,lty=2,ylim=c(0,max(c(a$y,dlogis(xs)))),col="blue",main="",xlab="",lwd=2)
    lines(xs,dlogis(xs),lwd=2)
    xs=seq(-7,7,len=200)
    concat=c(dlogis(xs),dnorm(xs,sd=sqrt(2)))
    plot(xs,dlogis(xs,log=TRUE),lwd=2,type="l",ylim=log(c(min(concat),max(concat))),xlab="theta",ylab="Log density")
    lines(xs,dnorm(xs,sd=sqrt(2),log=TRUE),col="red",lwd=2)


## ----child = "../../Generic_R/generic_indept10.Rmd", echo=TRUE----------------

## -----------------------------------------------------------------------------
## MVt10 independence sampler
indep.samp.t10=function(nits,theta0,log.pi,qmode,qV,...) {
    d=length(theta0)
    store=matrix(nrow=nits+1,ncol=d+1)
    psis=matrix(nrow=nits,ncol=d)
	
    nacc=0
    theta.curr=theta0
    log.pi.curr=log.pi(theta.curr,...)
    log.q.curr=dmvt(theta.curr,qmode,qV,df=10) ## defaults to log
    store[1,]=c(theta.curr,log.pi.curr)

    for (i in 1:nits) {
        psi=qmode+as.vector(rmvt(1,qV,df=10)); psis[i,]=psi
        log.pi.prop=log.pi(psi,...)
        log.q.prop=dmvt(psi,qmode,qV,df=10) ## defaults to log
        log.alpha=log.pi.prop+log.q.curr-log.pi.curr-log.q.prop
        if (log(runif(1))<log.alpha) {
            theta.curr=psi
            log.pi.curr=log.pi.prop
            log.q.curr=log.q.prop
            nacc=nacc+1
        }
        store[i+1,]=c(theta.curr,log.pi.curr)
    }
    
    return(list(acc=nacc/nits,store=store,psis=psis))
}



## ----echo=TRUE, fig.height=7--------------------------------------------------
    mcoutIS3=indep.samp.t10(1000,c(7),log.dens.logistic01, c(0), matrix(nrow=1,ncol=1,data=c(2)))

    par(mfrow=c(2,2))
    xs=seq(-7,7,len=200)
    dq=1/sqrt(2)*dt(xs/sqrt(2),df=10)
    dpi=dlogis(xs)
    plot(xs,dpi,lwd=2,type="l",ylim=c(0,max(c(dpi,dq))),xlab="theta",ylab="Density")
    lines(xs,dq,col="red",lwd=2)
    plot(0:1000,mcoutIS3$store[,1],type="l",ylab="theta",xlab="iteration")
    a=density(mcoutIS3$store[,1])
    plot(a,lty=2,ylim=c(0,max(c(a$y,dpi))),col="blue",main="",xlab="",lwd=2)
    lines(xs,dpi,lwd=2)
    xs=seq(-7,7,len=200)
    concat=c(dpi,dq)
    plot(xs,log(dpi),lwd=2,type="l",ylim=log(c(min(concat),max(concat))),xlab="theta",ylab="Log density")
    lines(xs,log(dq),col="red",lwd=2)


## ----echo=TRUE, fig.height=3.5------------------------------------------------
    par(mfrow=c(1,2))
    plot(0:50,mcoutIS3$store[1:51,1],type="l",ylab="theta",xlab="iteration")
    acf(mcoutIS3$store[6:1001,1],main="")
    ess.IS3=effectiveSize(mcoutIS3$store[6:1001])


## ----echo=TRUE, fig.height=7--------------------------------------------------
    set.seed(76543)
    mcoutRWM4=rwm(1000,c(0),log.dens.N0I,10,NULL)
    par(mfrow=c(3,2))
    plot(0:1000,mcoutRWM1$store[,1],type="l",ylab="theta",xlab="iteration",main="lambda=2")
    acf(mcoutRWM1$store[,1],main="lambda=2")
    plot(0:1000,mcoutRWM3$store[,1],type="l",ylab="theta",xlab="iteration",main="lambda=0.3")
    acf(mcoutRWM3$store[,1],main="lambda=0.3")
    plot(0:1000,mcoutRWM4$store[,1],type="l",ylab="theta",xlab="iteration",main="lambda=10")
    acf(mcoutRWM4$store[,1],main="lambda=10")


## ----echo=TRUE, fig.height=7--------------------------------------------------
    log.dens.N110100=function(x) {
      d=length(x)
      pows=0:(d-1)
      sds=10^(pows/2)
      return(sum(dnorm(x,0,sds,log=TRUE)))
    }
    set.seed(423)
    x0=rep(0,3)
    D=diag(c(1,10,100))
    mcoutRWM5=rwm(1000,x0,log.dens.N110100,3,NULL)
    mcoutRWM6=rwm(1000,x0,log.dens.N110100,1.4,D)
    par(mfrow=c(3,2))
    for (i in 1:3) {
    plot(0:1000,mcoutRWM5$store[,i],type="l",ylab=paste("theta",i),xlab="iteration",main="lambda^2 I")
    plot(0:1000,mcoutRWM6$store[,i],type="l",ylab=paste("theta",i),xlab="iteration",main="lambda^2 D")
    }



## ----child = "Sections/MHB.Rmd"-----------------------------------------------

## ----echo=TRUE----------------------------------------------------------------
library(mvtnorm)
library(coda)
set.seed(981237)


## ----child = "../../Generic_R/generic_rwm.Rmd", echo=TRUE---------------------

## -----------------------------------------------------------------------------
rwm=function(nits,theta0,log.pi,lambda,qV,...) {
    d=length(theta0)
    if (is.null(qV)) {
        A=diag(d)
    }
    else {
        A=t(chol(qV))
    }
    
    store=matrix(nrow=nits+1,ncol=d+1)
    psis=matrix(nrow=nits,ncol=d)
	
    nacc=0
    theta.curr=theta0
    log.pi.curr=log.pi(theta.curr,...)
    store[1,]=c(theta.curr,log.pi.curr)
    for (i in 1:nits) {
        psi=theta.curr+lambda*A%*%rnorm(d); psis[i,]=psi
        log.pi.prop=log.pi(psi,...)
        log.alpha=log.pi.prop-log.pi.curr
        if (log(runif(1))<log.alpha) {
            theta.curr=psi
            log.pi.curr=log.pi.prop
            nacc=nacc+1
        }
        store[i+1,]=c(theta.curr,log.pi.curr)
    }
    
    return(list(acc=nacc/nits,store=store,psis=psis))
}




## ----child = "../../Generic_R/generic_indept10.Rmd", echo=TRUE----------------

## -----------------------------------------------------------------------------
## MVt10 independence sampler
indep.samp.t10=function(nits,theta0,log.pi,qmode,qV,...) {
    d=length(theta0)
    store=matrix(nrow=nits+1,ncol=d+1)
    psis=matrix(nrow=nits,ncol=d)
	
    nacc=0
    theta.curr=theta0
    log.pi.curr=log.pi(theta.curr,...)
    log.q.curr=dmvt(theta.curr,qmode,qV,df=10) ## defaults to log
    store[1,]=c(theta.curr,log.pi.curr)

    for (i in 1:nits) {
        psi=qmode+as.vector(rmvt(1,qV,df=10)); psis[i,]=psi
        log.pi.prop=log.pi(psi,...)
        log.q.prop=dmvt(psi,qmode,qV,df=10) ## defaults to log
        log.alpha=log.pi.prop+log.q.curr-log.pi.curr-log.q.prop
        if (log(runif(1))<log.alpha) {
            theta.curr=psi
            log.pi.curr=log.pi.prop
            log.q.curr=log.q.prop
            nacc=nacc+1
        }
        store[i+1,]=c(theta.curr,log.pi.curr)
    }
    
    return(list(acc=nacc/nits,store=store,psis=psis))
}



## ----echo=TRUE----------------------------------------------------------------
    load(file="../Data/MCMC.rData")


## ----echo=TRUE----------------------------------------------------------------
log.pi.cens.gamma.fn=function(theta,xs) {
  alpha=2
  sstar=5
  A=1;B=1
  beta=theta[1]

  if (beta<=0) {
    log.prior=-Inf; log.like1=0; log.like2=0
  }
  else {
    uncensored=(xs<sstar)
    xsu=xs[uncensored]
    ncensored=length(xs)-length(xsu)

    log.prior=dgamma(beta,A,rate=B,log=TRUE)
    log.like1=sum(dgamma(xsu,alpha,rate=beta,log=TRUE))
    log.like2=ncensored*pgamma(sstar,alpha,rate=beta,lower.tail=FALSE,log.p=TRUE)
  }
  
  return(log.prior+log.like1+log.like2)
}



## ----echo=TRUE----------------------------------------------------------------
nits=2000
theta0=1
lambda=.15

mcoutRWMcens=rwm(nits,theta0,log.pi.cens.gamma.fn,lambda,NULL,xsgamcens)

par(mfrow=c(2,2))
hist(xsgamcens,xlab="",main="Data and histogram")
rug(xsgamcens)
plot(0:nits,mcoutRWMcens$store[,1],type="l",xlab="iteration",ylab="beta")
disc=1:50
acf(mcoutRWMcens$store[-disc,1],main="beta")
plot(density(mcoutRWMcens$store[-disc,1]),main="",xlab="beta")
essRWMcensbeta=effectiveSize(mcoutRWMcens$store[-disc,1])
posmeds=signif(median(mcoutRWMcens$store[,1]),digits=3)


## ----echo=TRUE----------------------------------------------------------------
log.pi.Gauss.fn=function(theta,xs) {
  lambda0=0; omega0=1/10000
  alpha0=1; gamma0=1
  mu=theta[1]; tau=theta[2]

  if (tau<=0) {
    log.prior=-Inf; log.like=0
  }	
  else  {
    log.prior=dnorm(mu,lambda0,sd=1/sqrt(omega0),log=TRUE)+dgamma(tau,alpha0,rate=gamma0,log=TRUE)
    log.like=sum(dnorm(xs,mu,1/sqrt(tau),log=TRUE))
  }
  
  return(log.prior+log.like)
}


## ----echo=TRUE----------------------------------------------------------------
ctl=list(fnscale=-1)
fit1=optim(c(1,1),log.pi.Gauss.fn,control=ctl,hessian=TRUE,xs=xsnorm)
qmode=fit1$par
qV= -solve(fit1$hessian)


## -----------------------------------------------------------------------------
qmode
qV


## ----fig.height=10,echo=TRUE--------------------------------------------------
nits=1000
mcoutIndepGauss = indep.samp.t10(nits,qmode,log.pi.Gauss.fn,qmode,qV,xsnorm)
par(mfrow=c(4,2))
plot(0:nits,mcoutIndepGauss$store[,1],type="l",xlab="iteration",ylab="mu")
plot(0:nits,mcoutIndepGauss$store[,2],type="l",xlab="iteration",ylab="tau")
acf(mcoutIndepGauss$store[,1],main="mu")
acf(mcoutIndepGauss$store[,2],main="tau")
plot(density(mcoutIndepGauss$store[,1]),main="",xlab="mu")
plot(density(mcoutIndepGauss$store[,2]),main="",xlab="tau")
plot(mcoutIndepGauss$store[,1],mcoutIndepGauss$store[,2],xlab="mu",ylab="tau")
plot(density(xsnorm),xlab="",main="Data and density estimate")
rug(xsnorm)
essmu=effectiveSize(mcoutIndepGauss$store[,1])
esstau=effectiveSize(mcoutIndepGauss$store[,2])
posmeds=signif(apply(mcoutIndepGauss$store[,1:2],2,median),digits=3)


## ----child = "../../Generic_R/log_pi_logreg.Rmd", echo=TRUE-------------------

## -----------------------------------------------------------------------------
log.pi.logistic.reg.fn=function(betas,ys,X) {
  ## Independent N(0,100) priors for the betas
  sds=rep(10,length(betas))
  log.prior=sum(dnorm(betas,0,sds,log=TRUE))

  etas= X %*% betas
  ps=exp(etas)/(1+exp(etas))
  log.like=sum(ys*log(ps)+(1-ys)*log(1-ps))

  return(log.prior+log.like)
}




## ----echo=TRUE----------------------------------------------------------------
ctl=list(fnscale=-1)
fit2=optim(c(1,1,1),log.pi.logistic.reg.fn,control=ctl,hessian=TRUE,ys=LogisticReg$ys,X=LogisticReg$X)
qmodeLR=fit2$par
qVLR= -solve(fit2$hessian)


## -----------------------------------------------------------------------------
qmodeLR
qVLR


## ----echo=TRUE, fig.height=8--------------------------------------------------
nits=2000
theta0=qmodeLR
lambda=1.3
mcoutRWMlogreg=rwm(nits,theta0,log.pi.logistic.reg.fn,lambda,qVLR,ys=LogisticReg$ys,X=LogisticReg$X)

par(mfrow=c(3,2))
plot(0:nits,mcoutRWMlogreg$store[,1],type="l",xlab="iteration",ylab="beta1")
plot(density(mcoutRWMlogreg$store[,1]),main="",xlab="beta1")
plot(0:nits,mcoutRWMlogreg$store[,2],type="l",xlab="iteration",ylab="beta2")
plot(density(mcoutRWMlogreg$store[,2]),main="",xlab="beta2")
plot(0:nits,mcoutRWMlogreg$store[,3],type="l",xlab="iteration",ylab="beta3")
plot(density(mcoutRWMlogreg$store[,3]),main="",xlab="beta3")

essb1=effectiveSize(mcoutRWMlogreg$store[,1])
essb2=effectiveSize(mcoutRWMlogreg$store[,2])
essb3=effectiveSize(mcoutRWMlogreg$store[,3])
posmeds=signif(apply(mcoutRWMlogreg$store[,1:3],2,median),digits=3)


## ----echo=TRUE----------------------------------------------------------------
round(var(mcoutRWMlogreg$store[,1:3]),dig=4)


## ----echo=TRUE, fig.height=4--------------------------------------------------
    plot(density(mixture),xlab="",main="Data and density estimate")
    rug(mixture)


## ----echo=TRUE----------------------------------------------------------------
log.pi.Gauss.mixture.fn=function(thetas,ys) {
  omega0=1/100; alpha0=1; gamma0=1; a0=1; b0=1
  mu1=thetas[1]; mu2=thetas[2]
  tau1=thetas[3]; tau2=thetas[4]
  p=thetas[5]
  n=length(ys)

  if ((p<=0) || (p>=1) || (tau1<=0) || (tau2<=0)) {
    return(-Inf)
  }

  log.pri1=sum(dnorm(c(mu1,mu2),0,sqrt(1/omega0),log=TRUE))
  log.pri2=sum(dgamma(c(tau1,tau2),alpha0,rate=gamma0,log=TRUE))
  log.pri3=dbeta(p,a0,b0,log=TRUE)

  log.like1= n*log(p)+n*log(tau1)/2-tau1/2*sum((ys-mu1)^2)
  diffs=tau1/2*(ys-mu1)^2-tau2/2*(ys-mu2)^2
  log.like2= sum(log(1+(1-p)/p*sqrt(tau2/tau1)*exp(diffs)))

  return(log.pri1+log.pri2+log.pri3+log.like1+log.like2)
}


## ----echo=TRUE----------------------------------------------------------------
ctl=list(fnscale=-1)
fit3=optim(c(0,3,1,1,.5),log.pi.Gauss.mixture.fn,control=ctl,hessian=TRUE,ys=mixture)
qmode=fit3$par
qV= -solve(fit3$hessian)


## -----------------------------------------------------------------------------
round(qmode,dig=4)
round(qV,dig=4)


## ----echo=TRUE,fig.height=10--------------------------------------------------
nits=2000

mcoutIndepMixture = indep.samp.t10(nits,qmode,log.pi.Gauss.mixture.fn,qmode,qV,ys=mixture)

par(mfrow=c(5,2))
plot(0:nits,mcoutIndepMixture$store[,1],type="l",xlab="iteration",ylab="mu1")
plot(density(mcoutIndepMixture$store[,1]),main="",xlab="mu1")
plot(0:nits,mcoutIndepMixture$store[,2],type="l",xlab="iteration",ylab="mu2")
plot(density(mcoutIndepMixture$store[,2]),main="",xlab="mu2")
plot(0:nits,mcoutIndepMixture$store[,3],type="l",xlab="iteration",ylab="tau1")
plot(density(mcoutIndepMixture$store[,3]),main="",xlab="tau1")
plot(0:nits,mcoutIndepMixture$store[,4],type="l",xlab="iteration",ylab="tau2")
plot(density(mcoutIndepMixture$store[,4]),main="",xlab="tau2")
plot(0:nits,mcoutIndepMixture$store[,5],type="l",xlab="iteration",ylab="p")
plot(density(mcoutIndepMixture$store[,5]),main="",xlab="p")

essIM=rep(0,5)
for (i in 1:5) {
  essIM[i]=effectiveSize(mcoutIndepMixture$store[,i])
}
posmeds=signif(apply(mcoutIndepMixture$store[,1:5],2,median),digits=3)


## ----echo=TRUE----------------------------------------------------------------
fit4=optim(c(3,0,1,1,.5),log.pi.Gauss.mixture.fn,control=ctl,hessian=TRUE,ys=mixture)
qmode=fit4$par
qV= -solve(fit4$hessian)


## -----------------------------------------------------------------------------
round(qmode,dig=4)


## ----echo=TRUE,fig.height=4---------------------------------------------------
mcoutIndepMixtureB = indep.samp.t10(nits,qmode,log.pi.Gauss.mixture.fn,qmode,qV,ys=mixture)
par(mfrow=c(1,2))
plot(density(mcoutIndepMixtureB$store[,1]),main="",xlab="mu1")
plot(density(mcoutIndepMixtureB$store[,2]),main="",xlab="mu2")


## ----echo=TRUE, fig.height=5, fig.width=5-------------------------------------
npt=100
mus=seq(-1,4,len=npt)
logposts=matrix(nrow=npt,ncol=npt)
for (i in 1:npt) {
  for (j in 1:npt) {
    logposts[i,j]=log.pi.Gauss.mixture.fn(c(mus[i],mus[j],1.5,1.5,.5),mixture)
  }
}
contour(mus,mus,logposts,xlab="mu1",ylab="mu2")


## ----echo=TRUE,fig.height=7---------------------------------------------------
rwm.cpt=function(nits,theta0,log.pi,lambdas,...) {
    d=length(theta0)
    if (is.null(qV)) {
       A=diag(d)
    }
    else {
       A=t(chol(qV))
    }
    store=matrix(nrow=nits+1,ncol=d+1)
	
    naccs=rep(0,d) ## one acceptance rate per component
    theta.curr=theta0
    log.pi.curr=log.pi(theta.curr,...)
    store[1,]=c(theta.curr,log.pi.curr)
	
    for (i in 1:nits) {
        for (j in 1:d) {
              psi=theta.curr
              psi[j]=theta.curr[j]+lambdas[j]*rnorm(1)
              log.pi.prop=log.pi(psi,...)
              log.alpha=log.pi.prop-log.pi.curr
              if (log(runif(1))<log.alpha) {
                theta.curr=psi
                log.pi.curr=log.pi.prop
                naccs[j]=naccs[j]+1
              }
          }
          store[i+1,]=c(theta.curr,log.pi.curr)
      }
      return(list(accs=naccs/nits,store=store))
}


nits=2000
theta0=qmodeLR
lambdas=c(0.6,1.0,0.7)
mcoutRWMlogregcpt=rwm.cpt(nits,theta0,log.pi.logistic.reg.fn,lambdas,ys=LogisticReg$ys,X=LogisticReg$X)

par(mfrow=c(3,2))
plot(0:nits,mcoutRWMlogregcpt$store[,1],type="l",xlab="iteration",ylab="beta1")
plot(density(mcoutRWMlogregcpt$store[,1]),main="",xlab="beta1")
plot(0:nits,mcoutRWMlogregcpt$store[,2],type="l",xlab="iteration",ylab="beta2")
plot(density(mcoutRWMlogregcpt$store[,2]),main="",xlab="beta2")
plot(0:nits,mcoutRWMlogregcpt$store[,3],type="l",xlab="iteration",ylab="beta3")
plot(density(mcoutRWMlogregcpt$store[,3]),main="",xlab="beta3")

essLRcpt=rep(0,3)
for (i in 1:3) {
  essLRcpt[i]=effectiveSize(mcoutRWMlogregcpt$store[,i])
}


## ----echo=TRUE----------------------------------------------------------------
par(mfrow=c(1,2))
boxplot(LogisticReg$X[,3]~LogisticReg$X[,2], xlab="x2",ylab="x3")
plot(mcoutRWMlogregcpt$store[,1],mcoutRWMlogregcpt$store[,2],xlab="beta2",ylab="beta3",cex=.2)
abline(v=0.25,lty=2,col="red")


## ----echo=TRUE----------------------------------------------------------------
TT=length(StochVol$ys)
plot(1:TT,StochVol$ys,type="l",xlab="t",ylab="y")


## ----echo=TRUE----------------------------------------------------------------
plot(1:TT,StochVol$xs,type="b",xlab="t",ylab="x")


## ----echo=TRUE----------------------------------------------------------------
qVSV=matrix(nrow=3,ncol=3,data=c(0.068,0.0007,0.014,
                                 0.00007,0.024,0.016,
                                 0.014,0.016,0.22),byrow=TRUE)
qVSV


## ----child = "../../Generic_R/log_pi_SV.Rmd", echo=TRUE-----------------------

## -----------------------------------------------------------------------------
## Log posterior for stochastic volatility model
log.pi.SV.fn<-function(theta,zs,ys) {
    mu=theta[1]; logsigma=theta[2]; logitrho=theta[3]
    lpri1=dnorm(mu,0,10,log=TRUE)
    lpri2=dnorm(logsigma,0,10,log=TRUE)
    lpri3=dnorm(logitrho,0,1,log=TRUE)

    sigma=exp(logsigma)
    rho=exp(logitrho)/(1+exp(logitrho))
    TT=length(zs)
    
    xs=zs ## get vector of right length, with right first element
    for (t in 2:TT) {
        xs[t]=rho*xs[t-1]+sqrt(1-rho*rho)*zs[t]
    }
    
    ll1=sum(dnorm(ys,0,exp(mu+sigma*xs),log=TRUE))
    ll2=sum(dnorm(zs,0,1,log=TRUE))

    return(lpri1+lpri2+lpri3+ll1+ll2)
}



## ----echo=TRUE----------------------------------------------------------------
## See Thinning in the final chapter of these notes
## for an explanation of the "thin" argument
SV.mcmc<-function(nits,theta0,zs0,lambdas,qV,ys,thin=100) {
    d=length(theta0)
    TT=length(zs0)
    
    naccs=rep(0,length(lambdas))
    store=matrix(nrow=nits/thin,ncol=d+2)
    jstore=1

    A=t(chol(qV))
    
    theta.curr=theta0
    zs.curr=zs0
    
    log.pi.curr=log.pi.SV.fn(theta.curr,zs.curr,ys)
    
    for (i in 1:nits) {

        psi=theta.curr+lambdas[1]*A%*%rnorm(d)
        log.pi.prop=log.pi.SV.fn(psi,zs.curr,ys)
        log.alpha=log.pi.prop-log.pi.curr
        if (log(runif(1))<log.alpha) {
            theta.curr=psi
            log.pi.curr=log.pi.prop
            naccs[1]=naccs[1]+1
        }

        logitrho=theta.curr[3]
        rho=exp(logitrho)/(1+exp(logitrho))

        epsilons=rnorm(TT)
        zs.prop=zs.curr+lambdas[2]*epsilons ## T-dimensional RWM
        log.pi.prop=log.pi.SV.fn(theta.curr,zs.prop,ys)

        log.alpha=log.pi.prop-log.pi.curr
        if (log(runif(1))<log.alpha) {
            zs.curr=zs.prop
            log.pi.curr=log.pi.prop
            naccs[2]=naccs[2]+1
        }

        if (i/thin==floor(i/thin)) {
            store[jstore,]=c(theta.curr,zs.curr[90],log.pi.curr)
            jstore=jstore+1
        }
    }
    
    return(list(accs=naccs/nits,store=store))
}


## ----echo=TRUE,fig.height=9---------------------------------------------------
set.seed(17231)

TT=length(StochVol$ys) ## Never name a variable T  !!!
theta0=c(0,0,log(.95/.05))
zs0=rep(0,TT)

lambdas=c(0.4,0.12)
nits=100000

mcoutSV=SV.mcmc(nits,theta0,zs0,lambdas,qVSV,StochVol$ys)

par(mfrow=c(4,2))
its=(1:(nits/100))*100

mus=mcoutSV$store[,1]
sigmas=exp(mcoutSV$store[,2])
rhos=exp(mcoutSV$store[,3])/(1+exp(mcoutSV$store[,3]))

burn=1:200

plot(its,mcoutSV$store[,1],type="l",xlab="iteration",ylab="mu")
plot(density(mus[-burn]),main="",xlab="mu")
plot(its,mcoutSV$store[,2],type="l",xlab="iteration",ylab="log sigma")
plot(density(sigmas[-burn]),main="",xlab="sigma")
plot(its,mcoutSV$store[,3],type="l",xlab="iteration",ylab="logit rho")
plot(density(rhos[-burn]),main="",xlab="rho")
plot(its,mcoutSV$store[,4],type="l",xlab="iteration",ylab="z90")
plot(density(mcoutSV$store[-burn,4]),main="",xlab="z90")


ESS.SV=effectiveSize(as.mcmc(mcoutSV$store[-burn,1:4]))



## ----child = "Sections/Gibbs.Rmd"---------------------------------------------

## ----echo=TRUE----------------------------------------------------------------
load("../Data/MCMC.rData")
library("coda")
set.seed(1235789)  ## for reproducibility of notes


## ----echo=TRUE----------------------------------------------------------------
GibbsGauss<-function(nits,mu0,xs) {
  lambda0=0; omega0=1/10000; alpha0=1; gamma0=1

  ## The updates depend on a few, fixed summary statistics.
  ## Calculate these once at the start rather than.
  ## recalculate them every iteration.

  n=length(xs); S=sum(xs)

  store=matrix(nrow=nits+1,ncol=2)

  mu.curr=mu0
  store[1,1]=mu.curr
  store[1,2]=1/var(xs) ## ballpark, to make tau trace plot look nice 

  for (i in 1:nits) {
     SS=sum((xs-mu.curr)^2)
     tau.curr=rgamma(1,alpha0+n/2,rate=gamma0+0.5*SS)

     tmp1=omega0+n*tau.curr
     tmp2=omega0*lambda0+tau.curr*S
     mu.curr=rnorm(1, tmp2/tmp1, sd=1/sqrt(tmp1))

     store[i+1,]=c(mu.curr,tau.curr)
  }

  return(list(store=store))
}


nits=1000
mu0=mean(xsnorm)

mcoutGibbsGauss=GibbsGauss(nits,mu0,xsnorm)

par(mfrow=c(2,2))
plot(0:nits,mcoutGibbsGauss$store[,1],type="l",xlab="iteration",ylab="mu")
plot(density(mcoutGibbsGauss$store[,1]),xlab="mu",main="")
plot(0:nits,mcoutGibbsGauss$store[,2],type="l",xlab="iteration",ylab="tau")
plot(density(mcoutGibbsGauss$store[,2]),xlab="tau",main="")

essGibbsGauss=effectiveSize(as.mcmc(mcoutGibbsGauss$store[-1,1:2]))
posmeds=signif(apply(mcoutGibbsGauss$store[-1,1:2],2,median),digits=3)


## ----echo=TRUE----------------------------------------------------------------
rCensoredGamma<-function(n,alpha,beta,lo) {
    uu=runif(n)*(1-pgamma(lo,alpha,beta))
    return(qgamma(1-uu,alpha,beta))
}


## ----child = "../../Generic_R/Gibbs4Gamma.Rmd", echo=TRUE---------------------

## -----------------------------------------------------------------------------
## Call the vector of precisely observed x values xsexact
gibbs.censored<-function(nits,beta0,xs) {
    sstar=5
    A0=1;B0=1
    alpha=2

    n=length(xs)
    censs=sort(xs) ## ascending order
    ncens=sum(xs==sstar)
    nexact=n-ncens
    xsexact=censs[1:nexact] ## vector of all the known (fixed) values

    store=matrix(nrow=nits+1,ncol=2) ## store the first latent variable too

    beta.curr=beta0
    store[1,]=c(beta.curr,sstar) 

    for (i in 1:nits) {

        ## Simulate from the conditional distribution of each latent variable
        latents=rCensoredGamma(ncens,alpha,beta.curr,sstar)
        ## Could have used a for loop, but since they are conditionally 
        ## independent given beta, we can simulate them as a vector.
        
        ## Simulate from the conditional distribution of beta
        A1=A0+2*n; B1=B0+sum(c(xsexact,latents))
        beta.curr=rgamma(1,A1,rate=B1)
        
        store[i+1,]=c(beta.curr,latents[1])
    }

    return(list(store=store))
}



## ----echo=TRUE, fig.height=7--------------------------------------------------
nits=5000
mcoutGibbsCens=gibbs.censored(nits,1,xsgamcens)

par(mfrow=c(2,2))
plot(0:nits,mcoutGibbsCens$store[,1],type="l",xlab="iteration",ylab="beta")
plot(density(mcoutGibbsCens$store[,1]),xlab="beta",main="")
plot(0:nits,mcoutGibbsCens$store[,2],type="l",xlab="iteration",ylab="z1")
plot(density(mcoutGibbsCens$store[,2]),xlab="z1",main="")

essGibbsCens=effectiveSize(as.mcmc(mcoutGibbsCens$store[-1,1:2]))
posmeds=signif(apply(mcoutGibbsCens$store[-1,1:2],2,median),digits=4)


## ----echo=TRUE----------------------------------------------------------------
GibbsGaussMixture<-function(nits,mu10,mu20,tau10,tau20,p0,xs) {
  lambda0=0; omega0=1/100; alpha0=1; gamma0=1; a0=1; b0=1

  xs=sort(xs) ## we will store the first, some middle and last z
              ## this will be more meaningful if the xs are sorted

  n=length(xs)

  store=matrix(nrow=nits+1,ncol=9)
  zstostore=c(1, floor(n/3), floor(2*n/3), n)

  mu1.curr=mu10; mu2.curr=mu20; tau1.curr=tau10; tau2.curr=tau20;
  p.curr=p0
  
  zs=rep(1,n) ## create storage
  store[1,1:5]=c(mu1.curr,mu2.curr,tau1.curr,tau2.curr,p.curr)
  store[1,6:9]=zs[zstostore]

  for (i in 1:nits) {

     ## Update the latent variables
     ## could have used a for loop, but can vectorise as they are
     ## indpendent of each other conditional on the parameters
     
     g1s=p.curr*sqrt(tau1.curr)*exp(-0.5*tau1.curr*(xs-mu1.curr)^2)
     g2s=(1-p.curr)*sqrt(tau2.curr)*exp(-0.5*tau2.curr*(xs-mu2.curr)^2)
     prob1s=g1s/(g1s+g2s)
     us=runif(n)
     zs=(us<=prob1s)*1 + (us>prob1s)*2

     xs1=xs[zs==1]; xs2=xs[zs==2]
     n1=length(xs1); n2=length(xs2)

     ## Update p
     
     p.curr=rbeta(1,a0+n1,b0+n2)

     ## Update mu1 and tau1
     S=sum(xs1)
     SS=sum((xs1-mu1.curr)^2)
     tau1.curr=rgamma(1,alpha0+n1/2,rate=gamma0+0.5*SS)

     tmp1=omega0+n1*tau1.curr
     tmp2=omega0*lambda0+tau1.curr*S
     mu1.curr=rnorm(1, tmp2/tmp1, sd=1/sqrt(tmp1))

     ## Update mu2 and tau2
     S=sum(xs2)
     SS=sum((xs2-mu2.curr)^2)
     tau2.curr=rgamma(1,alpha0+n2/2,rate=gamma0+0.5*SS)

     tmp1=omega0+n2*tau2.curr
     tmp2=omega0*lambda0+tau2.curr*S
     mu2.curr=rnorm(1, tmp2/tmp1, sd=1/sqrt(tmp1))

     store[i+1,1:5]=c(mu1.curr,mu2.curr,tau1.curr,tau2.curr,p.curr)
     store[i+1,6:9]=zs[zstostore]
   }

  return(list(store=store))
}


## ----fig.height=10, echo=TRUE-------------------------------------------------
nits=2000
mu10=mean(mixture); mu20=mu10
tau10=1/var(mixture); tau20=tau10
p0=1/2

mcoutGibbsMix=GibbsGaussMixture(nits,mu10,mu20,tau10,tau20,p0,mixture) 

par(mfrow=c(5,2))
plot(0:nits,mcoutGibbsMix$store[,1],type="l",xlab="iteration",ylab="mu1")
plot(density(mcoutGibbsMix$store[,1]),main="",xlab="mu1")
plot(0:nits,mcoutGibbsMix$store[,2],type="l",xlab="iteration",ylab="mu2")
plot(density(mcoutGibbsMix$store[,2]),main="",xlab="mu2")
plot(0:nits,mcoutGibbsMix$store[,3],type="l",xlab="iteration",ylab="tau1")
plot(density(mcoutGibbsMix$store[,3]),main="",xlab="tau1")
plot(0:nits,mcoutGibbsMix$store[,4],type="l",xlab="iteration",ylab="tau2")
plot(density(mcoutGibbsMix$store[,4]),main="",xlab="tau2")
plot(0:nits,mcoutGibbsMix$store[,5],type="l",xlab="iteration",ylab="p")
plot(density(mcoutGibbsMix$store[,5]),main="",xlab="p")

essGM=effectiveSize(as.mcmc(mcoutGibbsMix$store[-(1:51),1:5]))
posmeds=signif(apply(mcoutGibbsMix$store[-(1:51),1:5],2,median),digits=3)


## ----echo=TRUE----------------------------------------------------------------
par(mfrow=c(2,2))
plot(0:nits,mcoutGibbsMix$store[,6],type="p",cex=0.1,xlab="iteration",ylab="z1")
plot(0:nits,mcoutGibbsMix$store[,7],type="p",cex=0.1,xlab="iteration",ylab="zmid1")
plot(0:nits,mcoutGibbsMix$store[,8],type="p",cex=0.1,xlab="iteration",ylab="zmid2")
plot(0:nits,mcoutGibbsMix$store[,9],type="p",cex=0.1,xlab="iteration",ylab="zend")


## ----echo=TRUE----------------------------------------------------------------
## Create data
s1<-c(60,62,58,60,59,58)
s2<-c(55,52,54,57,52,54,58,54)
s3<-c(54,56,56,56,56)
s4<-c(61,59,58,61,59,63)
schooldata=list(s1=s1,s2=s2,s3=s3,s4=s4)
alls=c(s1,s2,s3,s4)
sbars=c(mean(s1),mean(s2),mean(s3),mean(s4))
theta0=c(mean(alls),1/var(s1),1/var(sbars)) ## sensible initial values

HLM2.Gibbs<-function(nits,theta0,ys) {
  lambda0=50; omega0=1/10000; a0=1; b0=1; c0=1; d0=1
  m=length(ys) ## length of list = number of schools
  d=length(theta0) ## three parameters
  store=matrix(nrow=nits+1,ncol=d+m) ## store parameters and latent variables

  theta.curr=theta0
  ns=rep(0,m); ydots=rep(0,m); xs.curr=rep(0,m)
  
  for (j in 1:m) {
    sj=ys[[j]]
    ns[j]=length(sj)
    ydots[j]=sum(sj)
    xs.curr[j]=mean(sj)  ## just to stop trace plots looking silly
  }
  ndot=sum(ns)
  
  store[1,]=c(theta.curr,xs.curr)

  for (i in 1:nits) {
    mu=theta.curr[1]; tau=theta.curr[2]; nu=theta.curr[3]

    ## Latent variance updates
    ## Can use vector updates as conditional on theta the xs are independent
    vrs=1/(tau+nu*ns)  ## vector of variances
    ees=(tau*mu+nu*ydots)*vrs  ## vector of expectations
    xs.curr=ees+sqrt(vrs)*rnorm(m)

    ## Update mu
    S=sum(xs.curr)
    vmu=1/(omega0+m*tau)
    emu=(omega0*lambda0+tau*S)*vmu
    mu=emu+sqrt(vmu)*rnorm(1)

    ## Update tau
    SS=sum((xs.curr-mu)^2)
    tau=rgamma(1,a0+m/2,rate=b0+SS/2)

    ## Update nu
    SSS=0
    for (j in 1:m) {
      sj=ys[[j]]
      SSS=SSS+sum((sj-xs.curr[j])^2)
    }
    nu=rgamma(1,c0+ndot/2,rate=d0+SSS/2)

    theta.curr=c(mu,tau,nu)
    store[i+1,]=c(theta.curr,xs.curr)
  }
  
  return(list(store=store))
}

nits=2000
mcoutHLM=HLM2.Gibbs(nits,theta0,schooldata)
essHLM=effectiveSize(as.mcmc(mcoutHLM$store))
posmeds=signif(apply(mcoutHLM$store[,1:3],2,median),digits=3)


## ----fig.height=6, echo=TRUE--------------------------------------------------

par(mfrow=c(3,2))
plot(0:nits,mcoutHLM$store[,1],type="l",xlab="iteration",ylab="mu")
plot(density(mcoutHLM$store[,1]),xlab="mu",main="")
plot(0:nits,mcoutHLM$store[,2],type="l",xlab="iteration",ylab="tau")
plot(density(mcoutHLM$store[,2]),xlab="tau",main="")
plot(0:nits,mcoutHLM$store[,3],type="l",xlab="iteration",ylab="nu")
plot(density(mcoutHLM$store[,3]),xlab="nu",main="")


## ----fig.height=8, echo=TRUE--------------------------------------------------
par(mfrow=c(4,2))
for (i in 1:4) {
  plot(0:nits,mcoutHLM$store[,i+3],type="l",xlab="iteration",ylab=paste("x",i,sep=""))
  plot(density(mcoutHLM$store[,i+3]),xlab=paste("x",i,sep=""),main="")
}


## ----echo=TRUE, fig.height=4--------------------------------------------------
set.seed(1237)
mu=10; tau=1
xs=10+rnorm(4)/sqrt(tau)
ysA=matrix(nrow=5,ncol=4,data=rnorm(20))
nuA=10; nuB=.05
ysB=ysA/sqrt(nuB)
ysA=ysA/sqrt(nuA)
for (i in 1:4) {
    ysA[,i]=ysA[,i]+xs[i]
    ysB[,i]=ysB[,i]+xs[i]	
}
par(mfrow=c(1,2))
ylo=min(c(as.vector(ysA),as.vector(ysB)))
yhi=max(c(as.vector(ysA),as.vector(ysB)))
plot(rep(1:4,5),as.vector(t(ysA)),pch="x",ylim=c(ylo,yhi),xlab="Group",ylab="y",main="Scenario A")
points(1:4,xs,col="red",pch=5,cex=1)
plot(rep(1:4,5),as.vector(t(ysB)),pch="x",ylim=c(ylo,yhi),xlab="Group",ylab="y",main="Scenario B")
points(1:4,xs,col="red",pch=5,cex=1)



## ----child = "Sections/MH2.Rmd"-----------------------------------------------

## ----echo=TRUE----------------------------------------------------------------
library(mvtnorm)
library(coda)
set.seed(981237)


## ----child = "../../Generic_R/generic_rwm.Rmd", echo=TRUE---------------------

## -----------------------------------------------------------------------------
rwm=function(nits,theta0,log.pi,lambda,qV,...) {
    d=length(theta0)
    if (is.null(qV)) {
        A=diag(d)
    }
    else {
        A=t(chol(qV))
    }
    
    store=matrix(nrow=nits+1,ncol=d+1)
    psis=matrix(nrow=nits,ncol=d)
	
    nacc=0
    theta.curr=theta0
    log.pi.curr=log.pi(theta.curr,...)
    store[1,]=c(theta.curr,log.pi.curr)
    for (i in 1:nits) {
        psi=theta.curr+lambda*A%*%rnorm(d); psis[i,]=psi
        log.pi.prop=log.pi(psi,...)
        log.alpha=log.pi.prop-log.pi.curr
        if (log(runif(1))<log.alpha) {
            theta.curr=psi
            log.pi.curr=log.pi.prop
            nacc=nacc+1
        }
        store[i+1,]=c(theta.curr,log.pi.curr)
    }
    
    return(list(acc=nacc/nits,store=store,psis=psis))
}




## ----child = "../../Generic_R/generic_indept10.Rmd", echo=TRUE----------------

## -----------------------------------------------------------------------------
## MVt10 independence sampler
indep.samp.t10=function(nits,theta0,log.pi,qmode,qV,...) {
    d=length(theta0)
    store=matrix(nrow=nits+1,ncol=d+1)
    psis=matrix(nrow=nits,ncol=d)
	
    nacc=0
    theta.curr=theta0
    log.pi.curr=log.pi(theta.curr,...)
    log.q.curr=dmvt(theta.curr,qmode,qV,df=10) ## defaults to log
    store[1,]=c(theta.curr,log.pi.curr)

    for (i in 1:nits) {
        psi=qmode+as.vector(rmvt(1,qV,df=10)); psis[i,]=psi
        log.pi.prop=log.pi(psi,...)
        log.q.prop=dmvt(psi,qmode,qV,df=10) ## defaults to log
        log.alpha=log.pi.prop+log.q.curr-log.pi.curr-log.q.prop
        if (log(runif(1))<log.alpha) {
            theta.curr=psi
            log.pi.curr=log.pi.prop
            log.q.curr=log.q.prop
            nacc=nacc+1
        }
        store[i+1,]=c(theta.curr,log.pi.curr)
    }
    
    return(list(acc=nacc/nits,store=store,psis=psis))
}



## ----echo=TRUE----------------------------------------------------------------
set.seed(1234756)


## ----child = "../../Generic_R/log_pi_logreg.Rmd", echo=TRUE-------------------

## -----------------------------------------------------------------------------
log.pi.logistic.reg.fn=function(betas,ys,X) {
  ## Independent N(0,100) priors for the betas
  sds=rep(10,length(betas))
  log.prior=sum(dnorm(betas,0,sds,log=TRUE))

  etas= X %*% betas
  ps=exp(etas)/(1+exp(etas))
  log.like=sum(ys*log(ps)+(1-ys)*log(1-ps))

  return(log.prior+log.like)
}




## ----echo=TRUE----------------------------------------------------------------
load(file="../Data/MCMC.rData")

ctl=list(fnscale=-1)
fit=optim(c(1,1,1),log.pi.logistic.reg.fn,control=ctl,hessian=TRUE,ys=LogisticReg$ys,X=LogisticReg$X)
qmodeLR=fit$par
qVLR= -solve(fit$hessian)


## ----echo=TRUE----------------------------------------------------------------
nits=2000
d=length(qmodeLR)
lambda=1.3
theta0=qmodeLR+rnorm(d) 
mcLR1=rwm(nits,theta0,log.pi.logistic.reg.fn,lambda,qVLR,ys=LogisticReg$ys,X=LogisticReg$X)
theta0=qmodeLR+rnorm(d) 
mcLR2=rwm(nits,theta0,log.pi.logistic.reg.fn,lambda,qVLR,ys=LogisticReg$ys,X=LogisticReg$X)
theta0=qmodeLR+rnorm(d) 
mcLR3=rwm(nits,theta0,log.pi.logistic.reg.fn,lambda,qVLR,ys=LogisticReg$ys,X=LogisticReg$X)


## ----echo=TRUE, fig.height=7--------------------------------------------------
combinedchains=mcmc.list(as.mcmc(mcLR1$store),as.mcmc(mcLR2$store),as.mcmc(mcLR3$store))
plot(combinedchains) ## trace plots on top of each other


## ----echo=TRUE----------------------------------------------------------------
gelman.diag(combinedchains)
gelman.plot(combinedchains)


## ----child = "../../Generic_R/log_pi_SV.Rmd", echo=TRUE-----------------------

## -----------------------------------------------------------------------------
## Log posterior for stochastic volatility model
log.pi.SV.fn<-function(theta,zs,ys) {
    mu=theta[1]; logsigma=theta[2]; logitrho=theta[3]
    lpri1=dnorm(mu,0,10,log=TRUE)
    lpri2=dnorm(logsigma,0,10,log=TRUE)
    lpri3=dnorm(logitrho,0,1,log=TRUE)

    sigma=exp(logsigma)
    rho=exp(logitrho)/(1+exp(logitrho))
    TT=length(zs)
    
    xs=zs ## get vector of right length, with right first element
    for (t in 2:TT) {
        xs[t]=rho*xs[t-1]+sqrt(1-rho*rho)*zs[t]
    }
    
    ll1=sum(dnorm(ys,0,exp(mu+sigma*xs),log=TRUE))
    ll2=sum(dnorm(zs,0,1,log=TRUE))

    return(lpri1+lpri2+lpri3+ll1+ll2)
}



## ----echo=TRUE----------------------------------------------------------------
## Gradient wrt the z values
grad.log.pi.SV.fn<-function(theta,zs,ys) {
    mu=theta[1]; logsigma=theta[2]; logitrho=theta[3]

    sigma=exp(logsigma)
    rho=exp(logitrho)/(1+exp(logitrho))
    TT=length(zs)
    
    xs=rep(0,TT); xs[1]=zs[1]
    for (t in 2:TT) {
        xs[t]=rho*xs[t-1]+sqrt(1-rho*rho)*zs[t]
    }

    g=rep(0,TT) ## gradient vector
    rhopowers=0:(TT-1)
    rhopowers=rho^(rhopowers) ## (1, rho, rho^2, ..., rho^{TT-1} )
    ## gfx is gradient wrt x of log f(y|x,\theta)
    gfx = -sigma+sigma*ys^2*exp(-2*mu-2*sigma*xs) ##

    g[1]=sum(rhopowers*gfx)
    for (t in 2:TT) {
      tkeep=t:TT
      g[t]=sqrt(1-rho^2)*sum(rhopowers[tkeep-t+1]*gfx[tkeep])
    }

## Add in grad of log f(z|theta)

    g=g-zs
    
    return(g)
}


## ----echo=TRUE----------------------------------------------------------------
    TT=length(StochVol$ys)


## ----echo=TRUE----------------------------------------------------------------
## Check gradient function gives the same as numerical differentiation
## *Always* worth doing when the calculations are messy.
theta0=c(0,0,log(.95/.05))
zs0=rnorm(TT)/10 ## so not zero in case special

g1=grad.log.pi.SV.fn(theta0,zs0,StochVol$ys)
g2=rep(0,TT)
eps=.001
for (t in 1:TT) {
  zsp=zs0; zsm=zs0
  zsp[t]=zsp[t]+eps; zsm[t]=zsm[t]-eps
  g2[t]=(log.pi.SV.fn(theta0,zsp,StochVol$ys)-log.pi.SV.fn(theta0,zsm,StochVol$ys))/(2*eps)
}
print(g1[1:5])
print(g2[1:5]) ## perfect match


## ----echo=TRUE----------------------------------------------------------------
qVSV=matrix(nrow=3,ncol=3,data=c(0.068,0.0007,0.014,
                                 0.00007,0.024,0.016,
                                 0.014,0.016,0.22),byrow=TRUE)
qVSV


## ----echo=TRUE----------------------------------------------------------------
SV.mcmc2<-function(nits,theta0,zs0,lambdas,qV,ys) {
    d=length(theta0)
    TT=length(zs0)
    thin=10
    
    naccs=rep(0,length(lambdas))
    store=matrix(nrow=nits/thin,ncol=d+2)
    jstore=1

    A=t(chol(qV))
    
    theta.curr=theta0
    zs.curr=zs0
    Vz=lambdas[2]^2*diag(TT)
    
    log.pi.curr=log.pi.SV.fn(theta.curr,zs.curr,ys)
        
    for (i in 1:nits) {

        psi=theta.curr+lambdas[1]*A%*%rnorm(d)
        log.pi.prop=log.pi.SV.fn(psi,zs.curr,ys)
        log.alpha=log.pi.prop-log.pi.curr
        if (log(runif(1))<log.alpha) {
            theta.curr=psi
            log.pi.curr=log.pi.prop
            naccs[1]=naccs[1]+1
        }

        logitrho=theta.curr[3]
        rho=exp(logitrho)/(1+exp(logitrho))

        glp.curr=grad.log.pi.SV.fn(theta.curr,zs.curr,ys)
	mean.from.curr=zs.curr+lambdas[2]^2/2*glp.curr
        epsilons=rnorm(TT)
        zs.prop=mean.from.curr+lambdas[2]*epsilons 
        log.pi.prop=log.pi.SV.fn(theta.curr,zs.prop,ys)

        glp.prop=grad.log.pi.SV.fn(theta.curr,zs.prop,ys)
	mean.from.prop=zs.prop+lambdas[2]^2/2*glp.prop
	lqfromcurr=dmvnorm(zs.prop,mean.from.curr,Vz,log=TRUE)
	lqfromprop=dmvnorm(zs.curr,mean.from.prop,Vz,log=TRUE)

	log.alpha=log.pi.prop+lqfromprop-log.pi.curr-lqfromcurr
	if (log(runif(1))<log.alpha) {
            zs.curr=zs.prop
	    log.pi.curr=log.pi.prop
	    naccs[2]=naccs[2]+1
        }

        if (i/thin==floor(i/thin)) {
            store[jstore,]=c(theta.curr,zs.curr[90],log.pi.curr)
            jstore=jstore+1
        }
    }
    
    return(list(accs=naccs/nits,store=store))
}


## ----echo=TRUE, fig.height=8--------------------------------------------------
lambdas=c(0.4,0.25) ## lambda[2] was 0.12 for RWM
nitsMALASV=30000 # CHANGE BACK TO 30000
thin=10

mcoutSV2=SV.mcmc2(nitsMALASV,theta0,zs0,lambdas,qVSV,StochVol$ys)

par(mfrow=c(4,2))
its=(1:(nitsMALASV/thin))*thin

mus=mcoutSV2$store[,1]
sigmas=exp(mcoutSV2$store[,2])
rhos=exp(mcoutSV2$store[,3])/(1+exp(mcoutSV2$store[,3]))

burn=1:50

plot(its,mcoutSV2$store[,1],type="l",xlab="iteration",ylab="mu")
plot(density(mus[-burn]),main="",xlab="mu")
plot(its,mcoutSV2$store[,2],type="l",xlab="iteration",ylab="log sigma")
plot(density(sigmas[-burn]),main="",xlab="sigma")
plot(its,mcoutSV2$store[,3],type="l",xlab="iteration",ylab="logit rho")
plot(density(rhos[-burn]),main="",xlab="rho")
plot(its,mcoutSV2$store[,4],type="l",xlab="iteration",ylab="z90")
plot(density(mcoutSV2$store[-burn,4]),main="",xlab="z90")


ESS.SV2=effectiveSize(as.mcmc(mcoutSV2$store[-burn,1:4]))


## ----echo=TRUE----------------------------------------------------------------
load("dimension.Rdata")
par(mfrow=c(1,2))
plot(log10(dimds),dimacc[,1],xlab="log10(d)",ylab="acc",pch=1,ylim=c(0,1))
points(log10(dimds),dimacc[,2],pch=2,col=2)
points(log10(dimds),dimacc[,3],pch=3,col=3)
legend(1.0,1.0,legend=c("IS","RWM","MALA"),col=1:3,pch=1:3)
abline(h=0.574,lty=2)
abline(h=0.234,lty=2)
elo=min(dimESSmean); ehi=max(dimESSmean)
plot(log10(dimds),log10(dimESSmean[,1]),ylim=log10(c(elo,ehi)),xlab="log10(d)",ylab="log10(meanESS)")
points(log10(dimds),log10(dimESSmean[,2]),pch=2,col=2)
points(log10(dimds),log10(dimESSmean[,3]),pch=3,col=3)
legend(0.2,2.5,legend=c("IS","RWM","MALA"),col=1:3,pch=1:3)


## ----child = "../../Generic_R/generic_rwm.Rmd", echo=TRUE---------------------

## -----------------------------------------------------------------------------
rwm=function(nits,theta0,log.pi,lambda,qV,...) {
    d=length(theta0)
    if (is.null(qV)) {
        A=diag(d)
    }
    else {
        A=t(chol(qV))
    }
    
    store=matrix(nrow=nits+1,ncol=d+1)
    psis=matrix(nrow=nits,ncol=d)
	
    nacc=0
    theta.curr=theta0
    log.pi.curr=log.pi(theta.curr,...)
    store[1,]=c(theta.curr,log.pi.curr)
    for (i in 1:nits) {
        psi=theta.curr+lambda*A%*%rnorm(d); psis[i,]=psi
        log.pi.prop=log.pi(psi,...)
        log.alpha=log.pi.prop-log.pi.curr
        if (log(runif(1))<log.alpha) {
            theta.curr=psi
            log.pi.curr=log.pi.prop
            nacc=nacc+1
        }
        store[i+1,]=c(theta.curr,log.pi.curr)
    }
    
    return(list(acc=nacc/nits,store=store,psis=psis))
}




## ----child = "../../Generic_R/generic_mala.Rmd", echo=TRUE--------------------

## -----------------------------------------------------------------------------
mala=function(nits,theta0,log.pi,grad.log.pi,lambda,qV,...) {
    d=length(theta0)
    if (is.null(qV)) {
        qV=diag(d)
    }
    
    store=matrix(nrow=nits+1,ncol=d+1)
	
    nacc=0
    theta.curr=theta0
    log.pi.curr=log.pi(theta.curr,...)
    grad.log.pi.curr=grad.log.pi(theta.curr,...)
    mu.from.theta=theta.curr+lambda^2/2*qV%*%grad.log.pi.curr

    store[1,]=c(theta.curr,log.pi.curr)
    
    for (i in 1:nits) {
        psi=as.vector(rmvnorm(1,mu.from.theta,lambda^2*qV))

        log.pi.prop=log.pi(psi,...)
        grad.log.pi.prop=grad.log.pi(psi,...)
        mu.from.psi=psi+lambda^2/2*qV%*%grad.log.pi.prop
        lq.theta.to.psi=dmvnorm(psi,mu.from.theta,lambda^2*qV,log=TRUE)
        lq.psi.to.theta=dmvnorm(theta.curr,mu.from.psi,lambda^2*qV,log=TRUE)

        log.alpha=log.pi.prop+lq.psi.to.theta-log.pi.curr-lq.theta.to.psi
        if (log(runif(1))<log.alpha) {
            theta.curr=psi
            log.pi.curr=log.pi.prop
            grad.log.pi.curr=grad.log.pi.prop
            mu.from.theta=mu.from.psi
            nacc=nacc+1
        }
        store[i+1,]=c(theta.curr,log.pi.curr)
    }
    
    return(list(acc=nacc/nits,store=store))
}



## ----echo=TRUE----------------------------------------------------------------
log.pi.light<-function(theta){
  return(-theta^4/4)
}
log.pi.heavy<-function(theta){
  return(-log(1+theta^2))
}
grad.log.pi.light<-function(theta){
  return(-theta^3)
}
grad.log.pi.heavy<-function(theta){
  return(-2*theta/(1+theta^2))
}

theta0=100
nits=5000
lamrwm=2.38; lammala=1.65

mcRWMlight=rwm(nits,theta0,log.pi.light,lamrwm,NULL)
mcRWMheavy=rwm(nits,theta0,log.pi.heavy,lamrwm,NULL)
mcMALAlight=mala(nits,theta0,log.pi.light,grad.log.pi.light,lammala,NULL)
mcMALAheavy=mala(nits,theta0,log.pi.heavy,grad.log.pi.heavy,lammala,NULL)

par(mfrow=c(2,2))
plot(0:nits,mcRWMlight$store[,1],xlab="Iteration",ylab="theta",main="RWM for pi light",type="l")
plot(0:nits,mcRWMheavy$store[,1],xlab="Iteration",ylab="theta",main="RWM for pi heavy",type="l")
plot(0:nits,mcMALAlight$store[,1],xlab="Iteration",ylab="theta",main="MALA for pi light",type="l")
plot(0:nits,mcMALAheavy$store[,1],xlab="Iteration",ylab="theta",main="MALA for pi heavy",type="l")


## ----echo=TRUE----------------------------------------------------------------
xs=seq(-2.5,2.5,len=1000)
plot(xs,-dnorm(xs,log=TRUE),type="l",xlab="theta",ylab="U")
points(1,-dnorm(1,log=TRUE),pch=19,col=1,cex=1.5)



## ----child = "Sections/App.Rmd"-----------------------------------------------

## ----echo=FALSE---------------------------------------------------------------
library(mvtnorm)
library(coda)
set.seed(981237)


## ----child = "../../Generic_R/generic_rwm.Rmd", echo=FALSE--------------------

## -----------------------------------------------------------------------------
rwm=function(nits,theta0,log.pi,lambda,qV,...) {
    d=length(theta0)
    if (is.null(qV)) {
        A=diag(d)
    }
    else {
        A=t(chol(qV))
    }
    
    store=matrix(nrow=nits+1,ncol=d+1)
    psis=matrix(nrow=nits,ncol=d)
	
    nacc=0
    theta.curr=theta0
    log.pi.curr=log.pi(theta.curr,...)
    store[1,]=c(theta.curr,log.pi.curr)
    for (i in 1:nits) {
        psi=theta.curr+lambda*A%*%rnorm(d); psis[i,]=psi
        log.pi.prop=log.pi(psi,...)
        log.alpha=log.pi.prop-log.pi.curr
        if (log(runif(1))<log.alpha) {
            theta.curr=psi
            log.pi.curr=log.pi.prop
            nacc=nacc+1
        }
        store[i+1,]=c(theta.curr,log.pi.curr)
    }
    
    return(list(acc=nacc/nits,store=store,psis=psis))
}




## ----child = "../../Generic_R/generic_indept10.Rmd", echo=FALSE---------------

## -----------------------------------------------------------------------------
## MVt10 independence sampler
indep.samp.t10=function(nits,theta0,log.pi,qmode,qV,...) {
    d=length(theta0)
    store=matrix(nrow=nits+1,ncol=d+1)
    psis=matrix(nrow=nits,ncol=d)
	
    nacc=0
    theta.curr=theta0
    log.pi.curr=log.pi(theta.curr,...)
    log.q.curr=dmvt(theta.curr,qmode,qV,df=10) ## defaults to log
    store[1,]=c(theta.curr,log.pi.curr)

    for (i in 1:nits) {
        psi=qmode+as.vector(rmvt(1,qV,df=10)); psis[i,]=psi
        log.pi.prop=log.pi(psi,...)
        log.q.prop=dmvt(psi,qmode,qV,df=10) ## defaults to log
        log.alpha=log.pi.prop+log.q.curr-log.pi.curr-log.q.prop
        if (log(runif(1))<log.alpha) {
            theta.curr=psi
            log.pi.curr=log.pi.prop
            log.q.curr=log.q.prop
            nacc=nacc+1
        }
        store[i+1,]=c(theta.curr,log.pi.curr)
    }
    
    return(list(acc=nacc/nits,store=store,psis=psis))
}



## ----echo=FALSE---------------------------------------------------------------
set.seed(12347)


