#' Some Bivariate Copulas
#' @name copula2d
#' @description Parametric bivariate copulas,  densities, and random number generators
# constructive-asymmetric-2d-copula of Wu (2014)
#' @param u,v vectors of same length at which the copula and its density is evaluated
#' @param n number of random vectors to be generated
#' @param lambda a vector of three mixing proportions which sum to one
#' @param copula the name of a copula to be called or a base copula for 
#'    construncting asymmetric copula(see Details)
#' @param ... the parameter(s) of \code{copula}, \code{theta} for most of the models, and
#'    \code{df}, the degrees of freedom if \code{copula='t'}, or \code{m} if \code{copula='nakagami'}
#' @details The names of available copulas are \code{'amh'} (Ali-Mikhai-Haq), \code{'bern'} (Bernstein polynomial model),
#'   \code{'clayton'}(Clayton), \code{'exponential'} (Exponential), \code{'fgm'}(Farlie–Gumbel–Morgenstern),
#'    \code{'frank'} (Frank), \code{'gauss'} (Gaussian), \code{'gumbel'} (Gumbel),
#'   \code{'indep'} (Independence), \code{'joe'} (Joe), \code{'nakagami'} (Nakagami-m), \code{'plackett'} (Plackett), 
#'   \code{'t'} (Student's t).
#'  \code{d2dcop.asym}, etc, calculate the constructive assymmetric copula of Wu (2014)
#' using base \code{copula} \eqn{C_{\theta}} with mixing proportions \eqn{p=(\lambda_1,\lambda_2,\lambda_3)} and 
#'   parameter values \eqn{\theta=(\theta_1,\theta_2,\theta_3)}: 
#'  \eqn{\lambda_0C_{\theta_0}(u,v)+\lambda_1[v-C_{\theta_1}(1-u,v)]+\lambda_2[u-C_{\theta_2}(u,1-v)]}. 
#'   If \code{copula='t'} or \code{'nakagami'}, \code{df} or \code{m} must be also given.
#' @return a vector of copula ot its density values evaluated at \code{(u,v)} 
#'     or an \code{n x 2} matrix of the generated observations
#' @references 
#' Nelsen, R. B. (1999). An Introduction to Copulas. Springer Series in Statistics. New York: Springer.
#'  Wu, S. (2014). Construction of asymmetric copulas and its application 
#'     in two-dimensional reliability modelling. 
#'     European Journal of Operational Research 238 (2), 476–485.
#' @importFrom stats qnorm dt qt
#' @export
d2dcop.asym<-function(u,v, lambda, copula='clayton', ...){
  para<-list(...)
  theta<-para$theta
  if(copula!='t' && copula!='nakagami'){
    cop<-lambda[1]*dcopula(u,v, copula, theta=theta[1])+lambda[2]*dcopula(1-u,v, copula, theta=theta[2])
    cop<-cop+lambda[3]*dcopula(u,1-v, copula, theta=theta[3])
  }
  else if(copula=='t'){
    # copula='t'
    df<-para$df
    if(length(df)<3) df<-rep(df,3)
    cop<-lambda[1]*dcop.t(u,v,theta=theta[1], df=df[1])+lambda[2]*dcop.t(1-u,v,theta=theta[2], df=df[2])
    cop<-cop+lambda[3]*dcop.t(u,1-v,theta=theta[3], df=df[3])
  }
  else{
    #copula="nakagami"
    m<-para$m
    if(length(m)<3) m<-rep(m,3)
    cop<-lambda[1]*dcop.nakagami(u,v,theta=theta[1], m=m[1])+lambda[2]*dcop.nakagami(1-u,v,theta=theta[2], m=m[2])
    cop<-cop+lambda[3]*dcop.nakagami(u,1-v,theta=theta[3], m=m[3])
  }
    
  cop
}
#' @rdname copula2d
#' @export
p2dcop.asym<-function(u,v, lambda, copula='clayton', ...){
  para<-list(...)
  theta<-para$theta
  if(copula!='t' && copula!='nakagami'){
    cop<-lambda[1]*pcopula(u,v, copula, theta=theta[1])+lambda[2]*(v-pcopula(1-u,v, copula, theta=theta[2]))
    cop<-cop+lambda[3]*(u-pcopula(u,1-v, copula, theta=theta[3]))
  }
  else if(copula=='t'){
    # copula='t'
    df<-para$df
    if(length(df)<3) df<-rep(df,3)
    cop<-lambda[1]*pcop.t(u,v,theta=theta[1], df=df[1])+lambda[2]*(v-pcop.t(1-u,v,theta=theta[2], df=df[2]))
    cop<-cop+lambda[3]*(u-pcop.t(u,1-v,theta=theta[3], df=df[3]))
  }
  else{
    #copula="nakagami"
    m<-para$m
    if(length(m)<3) m<-rep(m,3)
    cop<-lambda[1]*pcop.nakagami(u,v,theta=theta[1], m=m[1])+lambda[2]*pcop.nakagami(1-u,v,theta=theta[2], m=m[2])
    cop<-cop+lambda[3]*pcop.nakagami(u,1-v,theta=theta[3], m=m[3])
  }
  cop
}
#' @rdname copula2d
#' @export
r2dcop.asym<-function(n, lambda, copula='clayton', ...){
  para<-list(...)
  theta<-para$theta
  N<-sample(1:3, n, replace=TRUE, prob=lambda)
  out<-NULL
  if(copula!='t' && copula!='nakagami'){
    out<-rbind(out,rcopula(sum(N==1), copula, theta=theta[1]))
    out<-rbind(out, t(c(1,0)+c(-1,1)*t(rcopula(sum(N==2), copula, theta=theta[2]))))
    out<-rbind(out, t(c(0,1)+c(1,-1)*t(rcopula(sum(N==3), copula, theta=theta[3]))))
  }
  else if(copula=='t'){
    # copula='t'
    df<-para$df
    if(length(df)<3) df<-rep(df,3)
    out<-rbind(out,rcop.t(sum(N==1), theta=theta[1], df=df[1]))
    out<-rbind(out, t(c(1,0)+c(-1,1)*t(rcop.t(sum(N==2), theta=theta[2], df=df[2]))))
    out<-rbind(out, t(c(0,1)+c(1,-1)*t(rcop.t(sum(N==3), theta=theta[3], df=df[3]))))
  }
  else {
    #copula="nakagami"
    m<-para$m
    if(length(df)<3) df<-rep(df,3)
    out<-rbind(out,rcop.nakagami(sum(N==1), theta=theta[1], m=m[1]))
    out<-rbind(out, t(c(1,0)+c(-1,1)*t(rcop.nakagami(sum(N==2), theta=theta[2], m=m[2]))))
    out<-rbind(out, t(c(0,1)+c(1,-1)*t(rcop.nakagami(sum(N==3), theta=theta[3], m=m[3]))))
  }
  out
}

# unified calls
#' @rdname copula2d
#' @export
dcopula<-function(u, v, copula, ...){
    copula.arg.check(copula,...)
    switch(copula,
        bern=dcop.bern(u,v,...),
        amh=dcop.amh(u,v,...), 
        clayton=dcop.clayton(u,v,...), 
        exponential=dcop.exp(u,v,...), 
        fgm=dcop.fgm(u,v,...), 
        frank=dcop.frank(u,v,...), 
        gauss=dcop.gauss(u,v,...), 
        gumbel=dcop.gumbel(u,v,...),
        indep=dcop.indep(u,v,...), 
        joe=dcop.joe(u,v,...), 
        nakagami=dcop.nakagami(u,v,...), 
        plackett=dcop.plackett(u,v,...), 
        t=dcop.t(u,v,...)
    )
}   

#' @rdname copula2d
#' @export
pcopula<-function(u, v, copula, ...){
    copula.arg.check(copula,...)
    switch(copula,
        bern=pcop.bern(u,v,...),
        amh=pcop.amh(u,v,...), 
        clayton=pcop.clayton(u,v,...), 
        exponential=pcop.exp(u,v,...), 
        fgm=pcop.fgm(u,v,...), 
        frank=pcop.frank(u,v,...), 
        gauss=pcop.gauss(u,v,...), 
        gumbel=pcop.gumbel(u,v,...),
        indep=pcop.indep(u,v,...), 
        joe=pcop.joe(u,v,...), 
        nakagami=pcop.nakagami(u,v,...), 
        plackett=pcop.plackett(u,v,...), 
        t=pcop.t(u,v,...)
    )
}   

#' @rdname copula2d
#' @export
rcopula<-function(n, copula, ...){
    copula.arg.check(copula,...)
    switch(copula,
        bern=rcop.bern(n,...),
        amh=rcop.amh(n,...), 
        clayton=rcop.clayton(n,...), 
        exponential=rcop.exp(n,...), 
        fgm=rcop.fgm(n,...), 
        frank=rcop.frank(n,...), 
        gauss=rcop.gauss(n,...), 
        gumbel=rcop.gumbel(n,...),
        indep=rcop.indep(n,...), 
        joe=rcop.joe(n,...), 
        nakagami=rcop.nakagami(n,...), 
        plackett=rcop.plackett(n,...), 
        t=rcop.t(n,...)
    )
}   

####################################
# Checking parameter(s) 
####################################
#' @keywords internal
#' @noRd
copula.arg.check<-function(copula,...){
    para<-list(...)
    test<-
    switch(copula,
        bern=(any(para$p < 0) || sum(para$p)>1 || length(para$p)!=prod(para$m+1)),
        amh=(para$theta < -1 || para$theta>1), 
        clayton=(para$theta< -1 || para$theta==0), 
        exponential=(para$theta< 0 || para$theta>=1), 
        fgm=(abs(para$theta)>1), 
        frank=(para$theta< -1 || para$theta==0), 
        gauss=(para$theta<= -1 || para$theta>=1), 
        gumbel=(para$theta< 1),
        indep=FALSE, 
        joe=(para$theta<=1), 
        nakagami=(para$theta< 0 || para$theta>=1 || para$m<=0), 
        plackett=(para$theta<=0 || para$theta==1), 
        t=(abs(para$theta)>=1 || para$df<=0)
    )
    if(test) stop("Invalid paramter(s).")
}   

#########################################################
# Copula densities, distribution functions, and 
#  random number generators to be called
#########################################################
# Bernstein-polynomial-copula 
#' @keywords internal
#' @noRd
dcop.bern<-function(u,v, p, m){
  dmixmvbeta(cbind(u,v), p, m, interval=cbind(0:1,0:1))
}
#' @keywords internal
#' @noRd
pcop.bern<-function(u, v, p, m){
  pmixmvbeta(cbind(u,v), p, m, interval=cbind(0:1,0:1))
}
#' @keywords internal
#' @noRd

rcop.bern<-function(n, p, m){
  rmixmvbeta(n, p, m, interval=cbind(0:1,0:1))
}
#########################################################
# AMH-copula 
#' @keywords internal
#' @noRd
dcop.amh<-function(u,v, theta){
  t<-theta
  a<-1-t*(1-u)*(1-v); 
  ((1-t)*a+2*t*u*v)/a^3
}
#' @keywords internal
#' @noRd
pcop.amh<-function(u, v, theta=.8){
  u*v/(1-theta*(1-u)*(1-v))
}
#' @keywords internal
#' @noRd

rcop.amh<-function(n, theta=.8){
  v<-runif(n)
  u<-runif(n)
  for(i in 1:n)
    u[i]<-qcop.amh.cond(u[i], v[i], theta)
  cbind(u,v)
}
##############################################
# clayton-copula
#' @keywords internal
#' @noRd
dcop.clayton<-function(u,v, theta){
  t<-theta
  # (u^t+v^t-(u*v)^t>0)*(1+t)*(u*v)^t/(u^t+v^t-(u*v)^t)^(2+1/t)
  (u^(-t)+v^(-t)>1)*(1+t)*(u*v)^(-t-1)/(u^(-t)+v^(-t)-1)^(2+1/t)
}
#' @keywords internal
#' @noRd
pcop.clayton<-function(u, v, theta=.8){
  if(theta>=-1 && theta!=0) (pmax(u^(-theta)+v^(-theta)-1, 0))^(1/theta)
  else NULL
}

#' @keywords internal
#' @noRd
rcop.clayton<-function(n, theta=.8){
  v<-runif(n)
  if(theta==-1) u<-1-v
  else if(theta>-1 && theta!=0){    
    u<-runif(n)
    for(i in 1:n)
      u[i]<-qcop.clayton.cond(u[i], v[i], theta)
  }
  else u<-NULL
  cbind(u,v)
}
##############################################
# exponential-copula
#' @keywords internal
#' @noRd
dcop.exp<-function(u,v, theta=.8){
  t<-theta
  if(t==0) rep(1, length(u))
  else ((1-u)*(1-v))^(t/(1-t))*besselI(2*sqrt(t*log(1-u)*log(1-v))/(1-t), nu=0)/(1-t)
}
#' @importFrom stats pchisq
#' @keywords internal
#' @noRd
#pcop.exp<-function(u, v, theta=.8){
#  if(theta==0) u*v
#  else {
#    x<--log(1-u)/(1-theta)
#    y<--log(1-v)/(1-theta)
#    tmpf<-function(s, x, t) exp(-s)*besselI(2*sqrt(t*x*s), nu=0)
#    n<-length(u)
#    out<-rep(0,n)
#    for(i in 1:n) {
#        out[i]<-1+(1-v[i])*(exp(-x[i])*integrate(tmpf, 0, theta*y[i], x=x[i], t=1)$value-1)
#        out[i]<-out[i]-exp(-x[i])*integrate(tmpf, 0, y[i], x=x[i], t=theta)$value
#    }
#    out
#  }
#}
## the R function for noncentral chisq is inaccurate when ncp is large
pcop.exp<-function(u, v, theta=.8){
  eps<-.Machine$double.eps^.5
  if(theta==0) u*v
  else {
    x<--log(1-u)/(1-theta)
    y<--log(1-v)/(1-theta)
    n<-length(u)
    out<-rep(0,n)
    for(i in 1:n) {
        if(1-v[i]<eps && 1-u[i]<eps) out[i]<-1
        else if (1-v[i]<eps) out[i]<-1-(1-u[i])*pchisq(2*y[i], 2, 2*theta*x[i])
        else if (1-u[i]<eps) out[i]<-1+(1-v[i])*(pchisq(2*theta*y[i], 2, 2*x[i])-1)
        else{
            out[i]<-1+(1-v[i])*(pchisq(2*theta*y[i], 2, 2*x[i])-1)
            out[i]<-out[i]-(1-u[i])*pchisq(2*y[i], 2, 2*theta*x[i])
        }
    }
    out
  }
}
#u<-seq(0,1,length=20)
#pcop.exp(u, v=u, theta=.8)
#pcop.exp1(u, v=u, theta=.8)
#pcop.exp(u, v=u, theta=.8)
#pcop.exp2(u, v=u, theta=.8)
#' @keywords internal
#' @noRd
rcop.exp<-function(n, theta=.8){
  v<-runif(n)
  u<-runif(n)
  if(theta>0) {    
    for(i in 1:n)
      u[i]<-qcop.exp.cond(u[i], v[i], theta)
  }
  cbind(u,v)
}

##########################################
# fgm-copula
#' @keywords internal
#' @noRd
dcop.fgm<-function(u,v, theta) 1+theta*(1-2*u)*(1-2*v)
#' @keywords internal
#' @noRd
pcop.fgm<-function(u, v, theta=.8){
  u*v*(1+theta*(1-u)*(1-v))
}

#' @keywords internal
#' @noRd
rcop.fgm<-function(n, theta=.8){
  v<-runif(n)
  u<-runif(n)
  for(i in 1:n)
    u[i]<-qcop.fgm.cond(u[i], v[i], theta)
  cbind(u,v)
}
########################################
# frank-copula
#' @keywords internal
#' @noRd
dcop.frank<-function(u,v, theta) {
  t<-theta
  t*(1-exp(-t))*exp(-t*(u+v))/((exp(-t*u)-1)*(exp(-t*v)-1)+exp(-t)-1)^2
}
#' @keywords internal
#' @noRd
pcop.frank<-function(u, v, theta=.8){
  -log(1+(exp(-theta*u)-1)*(exp(-theta*v)-1)/(exp(-theta)-1))/theta
}

#' @keywords internal
#' @noRd
rcop.frank<-function(n, theta=.8){
  v<-runif(n)
  u<-runif(n)
  for(i in 1:n)
    u[i]<-qcop.frank.cond(u[i], v[i], theta)
  cbind(u,v)
}
##########################################
# gaussian-copula
#' @keywords internal
#' @noRd
dcop.gauss<-function(u, v, theta){
  t<-theta
  exp(-(t^2*(qnorm(u)^2+qnorm(v)^2)-2*t*qnorm(u)*qnorm(v))/(2*(1-t^2)))/sqrt(1-t^2)
}
#' @importFrom mnormt pmnorm 
#' @importFrom stats qnorm rnorm
#' @keywords internal
#' @noRd
pcop.gauss<-function(u, v, theta=.8){
   pmnorm(cbind(qnorm(u), qnorm(v)), mean=c(0,0), varcov=matrix(c(1,theta, theta,1),2,2))
}

#' @keywords internal
#' @noRd
rcop.gauss<-function(n, theta=.8){
  v<-runif(n)
  u<-runif(n)
  for(i in 1:n)
    u[i]<-qcop.gauss.cond(u[i], v[i], theta)
  cbind(u,v)
}
############################################
# gumbel-copula
#' @keywords internal
#' @noRd
dcop.gumbel<-function(u,v, theta){
  t<-theta
  a<-(-log(u))^t+(-log(v))^t; 
  (log(u)*log(v))^(t-1)*(a^(1/t)+t-1)/exp(a^(1/t))/(u*v*a^(2-1/t))
}
#' @keywords internal
#' @noRd
pcop.gumbel<-function(u, v, theta=.8){
   exp(-((-log(u))^theta+(-log(v))^theta)^(1/theta))
}

#' @keywords internal
#' @noRd
rcop.gumbel<-function(n, theta=.8){
  v<-runif(n)
  u<-runif(n)
  for(i in 1:n)
    u[i]<-qcop.gumbel.cond(u[i], v[i], theta)
  cbind(u,v)
}
#############################################
# indepence-copula
#' @keywords internal
#' @noRd
dcop.indep<-function(u,v){
  dunif(u)*dunif(v)
}
#' @keywords internal
#' @noRd
pcop.indep<-function(u, v, theta=.8){
   punif(u)*punif(v)
}

#' @keywords internal
#' @noRd
rcop.indep<-function(n, theta=.8){
  v<-runif(n)
  u<-runif(n)
  cbind(u,v)
}
###########################################
#  joe-copula
#' @keywords internal
#' @noRd
dcop.joe<-function(u,v, theta){
  t<-theta
  a<-(1-u)^t+(1-v)^t-((1-u)*(1-v))^t; 
  ((1-u)*(1-v))^(t-1)*(a+t-1)/(a^(2-1/t))
}
#' @keywords internal
#' @noRd
pcop.joe<-function(u, v, theta=.8){
   1-((1-u)^theta+(1-v)^theta-((1-u)*(1-v))^theta)^(1/theta)
}

#' @keywords internal
#' @noRd
rcop.joe<-function(n, theta=.8){
  v<-runif(n)
  u<-runif(n)
  for(i in 1:n)
    u[i]<-qcop.joe.cond(u[i], v[i], theta)
  cbind(u,v)
}

#########################################
# Nakagami-m-copula
#' @keywords internal
#' @noRd
dcop.nakagami<-function(u,v, theta=.8, m=2){
  t<-theta
  x<-qgamma(u,m)
  y<-qgamma(v,m) 
  xyt<-sqrt(t*x*y)
  out<-gamma(m)*exp(-t*(x+y)/(1-t))*besselI(2*xyt/(1-t), nu=m-1)/xyt^(m-1)/(1-t)
  out
}
#dcop.nakagami(u<-seq(0,1,length=20),u)

#' @keywords internal
#' @noRd
pcop.nakagami<-function(u, v, theta=.8, m=2){
  t<-theta
  cn<-(1-t)^m/gamma(m)
  qm<-cbind(qgamma(u,m)/(1-t),qgamma(v,m)/(1-t))
  tmpf<-function(k)
    pgamma(qm[,1],m+k, log.p=TRUE)+pgamma(qm[,2],m+k, log.p=TRUE) 
  out<-cn*exp(lgamma(m)+tmpf(0))
  eps<-1e-10
  del<-cn*t*exp(lgamma(m+1)+tmpf(1))
  k<-1
  while(max(del, na.rm =TRUE)>eps){
    out<-out+del
    k<-k+1
    del<-cn*t^k*exp(lgamma(m+k)+tmpf(k)-lgamma(k+1))
    #cat(max(del, na.rm =TRUE),"\n")
  }
  out
}
#pcop.nakagami(u<-seq(0,1,length=20),u)

#' @keywords internal
#' @noRd
rcop.nakagami<-function(n, theta=.8, m=2){
  v<-runif(n)
  u<-runif(n)
  for(i in 1:n)
    u[i]<-qcop.nakagami.cond(u[i], v[i], theta, m)
  cbind(u,v)
}

#########################################
# plackett-copula
#' @keywords internal
#' @noRd
dcop.plackett<-function(u,v, theta){
  t<-theta
  t*(1+(u-2*u*v+v)*(t-1))/((1+(t-1)*(u+v))^2-4*t*(t-1)*u*v)^1.5
}
#' @keywords internal
#' @noRd
pcop.plackett<-function(u, v, theta=.8){
    a<-1+(theta-1)*(u+v)
   .5*(a-sqrt(a^2-4*u*v*theta*(theta-1)))/(theta-1)
}

#' @keywords internal
#' @noRd
rcop.plackett<-function(n, theta=.8){
  v<-runif(n)
  u<-runif(n)
  for(i in 1:n)
    u[i]<-qcop.plackett.cond(u[i], v[i], theta)
  cbind(u,v)
}

#############################
# student-t-colula
#' @keywords internal
#' @noRd
dcop.t<-function(u,v, theta, df){
  r<-theta
  d<-df
  tu<-qt(u,d);tv<-qt(v,d); 
  (1+(tu^2-2*r*tu*tv+tv^2)/(d*(1-r^2)))^(-(d+2)/2)/(2*pi*sqrt(1-r^2)*dt(tu,d)*dt(tv,d))
}

#' @importFrom mnormt pmt 
#' @importFrom stats qt rt
#' @keywords internal
#' @noRd
pcop.t<-function(u, v, theta=.8, df=4){
   pmt(cbind(qt(u, df), qt(v, df)), mean=c(0,0), S=matrix(c(1,theta, theta,1),2,2), df)
}


#' @keywords internal
#' @noRd
rcop.t<-function(n, theta=.8, df=4){
  v<-runif(n)
  u<-runif(n)
  for(i in 1:n)
    u[i]<-qcop.t.cond(u[i], v[i], theta, df=4)
  cbind(u,v)
}

###############################################
# Bivariate Normal
###############################################
rbivarnorm<-function(n, mean=c(0,0), sd=c(1,1), rho=.4){
  y<-rnorm(n, mean=mean[2], sd=sd[2])
  x<-rnorm(n, mean=mean[1]+sd[1]*rho*(y-mean[2])/sd[2], sd=sqrt(1-rho^2)*sd[1])
  cbind(x,y)
}
dbivarnorm<-function(x, mean=c(0,0), sd=c(1,1), rho=.4){
  x<-t((t(x)-mean)/sd)
  out<-0.5*exp(-0.5*(x[,1]^2-2*rho*x[,1]*x[,2]+x[,2]^2)/(1-rho^2))
  out/(pi*sd[1]*sd[2]*sqrt(1-rho^2))
}
##################################################### 
#####################################################
# Conditional copula cdf, pdf, quatile, 
#  and random number generator
#####################################################
#' Some Parametric Conditional Bivariate Copulas   
#' @name copula2d.cond
#' @description Density, distribution function, quantile function and 
#'   random generation for conditional copula \eqn{C(u|V=v)} of \eqn{U} given \eqn{V=v}
#'    related to parametric bivariate copula \eqn{C(u,v)=P(U\le u, V\le v)}.
#' @param u vector of \eqn{U} values at which the copula density is evaluated
#' @param v a given value of \eqn{V} under which the conditional copula 
#'    and its density is evaluated
#' @param n number of observations to be generated from conditional copula \eqn{C(u|V=v)}.
#' @param p a vector of probabilities
#' @param copula the name of a copula density to be called (see Details)
#' @param ... the parameter(s) of \code{copula}
#' @details the names of available copulas are \code{'amh'} (Ali-Mikhai-Haq), \code{'bern'} (Bernstein polynomial model),
#'   \code{'clayton'}(Clayton), \code{'exponential'} (Exponential),  
#'   \code{'fgm'}(Farlie–Gumbel–Morgenstern),  \code{'frank'} (Frank), 
#'   \code{'gauss'} (Gaussian), \code{'gumbel'} (Gumbel), \code{'indep'} (Independence), 
#'   \code{'joe'} (Joe), \code{'nakagami'} (Nakagami-m), \code{'plackett'} (Plackett), 
#'   \code{'t'} (Student's t).
#' @return a vector of copula density values evaluated at \code{u} gvien \code{V=v}
#'     or a vector of \code{n} generated \code{u} values from conditional copula \eqn{C(u|V=v)}.
#' @export

dcopula.cond<-function(u, v, copula, ...){
    copula.arg.check(copula,...)
    switch(copula,
        bern=dcop.bern.cond(u,v,...),
        amh=dcop.amh.cond(u,v,...), 
        clayton=dcop.clayton.cond(u,v,...), 
        exponential=dcop.exp.cond(u,v,...), 
        fgm=dcop.fgm.cond(u,v,...), 
        frank=dcop.frank.cond(u,v,...), 
        gauss=dcop.gauss.cond(u,v,...), 
        gumbel=dcop.gumbel.cond(u,v,...),
        indep=dcop.indep.cond(u,v,...), 
        joe=dcop.joe.cond(u,v,...), 
        nakagami=dcop.nakagami.cond(u,v,...), 
        plackett=dcop.plackett.cond(u,v,...), 
        t=dcop.t.cond(u,v,...)
    )
}   

#' @rdname copula2d.cond
#' @export
pcopula.cond<-function(u, v, copula, ...){
    copula.arg.check(copula,...)
    switch(copula,
        bern=pcop.bern.cond(u,v,...),
        amh=pcop.amh.cond(u,v,...), 
        clayton=pcop.clayton.cond(u,v,...), 
        exponential=pcop.exp.cond(u,v,...), 
        fgm=pcop.fgm.cond(u,v,...), 
        frank=pcop.frank.cond(u,v,...), 
        gauss=pcop.gauss.cond(u,v,...), 
        gumbel=pcop.gumbel.cond(u,v,...),
        indep=pcop.indep.cond(u,v,...), 
        joe=pcop.joe.cond(u,v,...), 
        nakagami=pcop.nakagami.cond(u,v,...), 
        plackett=pcop.plackett.cond(u,v,...), 
        t=pcop.t.cond(u,v,...)
    )
}   
#' @rdname copula2d.cond
#' @export
qcopula.cond<-function(p, v, copula, ...){
    copula.arg.check(copula,...)
    switch(copula,
        bern=qcop.bern.cond(p,v,...),
        amh=qcop.amh.cond(p,v,...), 
        clayton=qcop.clayton.cond(p,v,...), 
        exponential=qcop.exp.cond(p,v,...), 
        fgm=qcop.fgm.cond(p,v,...), 
        frank=qcop.frank.cond(p,v,...), 
        gauss=qcop.gauss.cond(p,v,...), 
        gumbel=qcop.gumbel.cond(p,v,...),
        indep=qcop.indep.cond(p,v,...), 
        joe=qcop.joe.cond(p,v,...), 
        nakagami=qcop.nakagami.cond(p,v,...), 
        plackett=qcop.plackett.cond(p,v,...), 
        t=qcop.t.cond(p,v,...)
    )
}   

#' @rdname copula2d.cond
#' @export
rcopula.cond<-function(n, v, copula, ...){
    copula.arg.check(copula,...)
    switch(copula,
        bern=rcop.bern.cond(n,v,...),
        amh=rcop.amh.cond(n,v,...), 
        clayton=rcop.clayton.cond(n,v,...), 
        exponential=rcop.exp.cond(n,v,...), 
        fgm=rcop.fgm.cond(n,v,...), 
        frank=rcop.frank.cond(n,v,...), 
        gauss=rcop.gauss.cond(n,v,...), 
        gumbel=rcop.gumbel.cond(n,v,...),
        indep=rcop.indep.cond(n,v,...), 
        joe=rcop.joe.cond(n,v,...), 
        nakagami=rcop.nakagami.cond(n,v,...), 
        plackett=rcop.plackett.cond(n,v,...), 
        t=rcop.t.cond(n,v,...)
    )
}   
##################################################
# Internal functions
#############################
# Bernstein
#############################
#' @keywords internal
#' @noRd
dcop.bern.cond<-function(u, v=.5, p, m){
  dcop.bern(u, v, p, m)
}
#' @keywords internal
#' @noRd
pcop.bern.cond<-function(u, v=.5, p, m){
  P<-matrix(p, nrow=m[1]+1) # rowSums=1/([m[1]+1)
  pv<-dmixbeta(v, rowSums(P))
  pmixbeta(u, pv)
}
#' @keywords internal
#' @noRd
qcop.bern.cond<-function(u, v=.5, p, m){
  P<-matrix(p, nrow=m[1]+1) # rowSums=1/([m[1]+1)
  pv<-dmixbeta(v, rowSums(P))
  qmixbeta(u, pv)
}
#' @keywords internal
#' @noRd
rcop.bern.cond<-function(n, v=.5, p, m){
  P<-matrix(p, nrow=m[1]+1) # rowSums=1/([m[1]+1)
  pv<-dmixbeta(v, rowSums(P))
  rmixbeta(n, pv)
}

#############################
# AMH
#############################
#' @keywords internal
#' @noRd
dcop.amh.cond<-function(u, v=.5, theta=.8){
  dcop.amh(u, v, theta)
}
#' @keywords internal
#' @noRd
pcop.amh.cond<-function(u, v=.5, theta=.8){
  u*(1-theta*(1-u))/(1-theta*(1-u)*(1-v))^2
}
#' @keywords internal
#' @noRd
qcop.amh.cond<-function(p, v=.5, theta=.8){
  if(theta==0) p
  else{
      A<-(p*theta*(1-v)^2-1)*theta
      B<-1+theta-2*p*theta*(1-v)
      C<-p-1
      q<-apply(round((-B+sqrt(B^2-4*A*C)%o%c(-1,1))/(2*A),10), 1, function(x) x[x>=0 & x<=1])
      1-q
  }
}
#qcop.amh.cond1<-function(p, v=.5, theta=.8){
#   if(theta==0) p
#   else{
#      A<-(p*theta*(1-v)^2-1)*theta
#      B<-1+theta-2*p*theta*(1-v)
#      C<-p-1
#      q<-round((-B+sign(theta)*sqrt(B^2-4*A*C))/(2*A),10)
#      1-q
#   }
#}
#u<-seq(0,1,length=20)
#qcop.amh.cond(u)
#qcop.amh.cond1(u)
#' @keywords internal
#' @noRd
rcop.amh.cond<-function(n, v=.5, theta=.8){
  qcop.amh.cond(runif(n), v, theta)
}


#############################
# Clayton
#############################

#' @keywords internal
#' @noRd
dcop.clayton.cond<-function(u, v=.5, theta=.8){
  if(theta==-1) 1*(u==1-v)
  else if(theta>-1 && theta!=0) dcop.clayton(u, v, theta)
  else NULL
}
#' @keywords internal
#' @noRd
pcop.clayton.cond<-function(u, v=.5, theta=.8){
  if(theta==-1) 1*(u>=1-v)
  else if(theta>-1 && theta!=0) 1/(v*(u^(-theta)+v^(-theta)-1)^(1/theta))^(1+theta)
  else NULL
}
#' @keywords internal
#' @noRd
qcop.clayton.cond<-function(p, v=.5, theta=.8){
  if(theta==-1) (1-v)*(p>0)
  else if(theta>-1 && theta!=0)  (1-v^(-theta)*(1-p^(-theta/(1+theta))))^(-1/theta)
  else NULL
}
#' @keywords internal
#' @noRd
rcop.clayton.cond<-function(n, v=.5, theta=.8){
  if(theta==-1) (1-v)*rep(1,n)
  else if(theta>-1 && theta!=0){
    p<-runif(n)
    (1-v^(-theta)*(1-p^(-theta/(1+theta))))^(-1/theta)
  }
  else NULL
}

#############################
# Exponential
#############################

#' @keywords internal
#' @noRd
dcop.exp.cond<-function(u, v=.5, theta=.8){
  dcop.exp(u, v, theta)
}

#' @importFrom stats dchisq qchisq
#' @keywords internal
#' @noRd
#pcop.exp.cond<-function(u, v=.5, theta=.8){
#  if(theta==0) u
#  else {
#    x<--log(1-u)/(1-theta)
#    y<--log(1-v)/(1-theta)
#    tmpf<-function(s) exp(-s)*besselI(2*sqrt(theta*y*s), nu=0)
#    n<-length(u)
#    out<-rep(0,n)
#    for(i in 1:n) {
#        out[i]<-(1-v)^(theta/(1-theta))*integrate(tmpf, 0, x[i])$value
#    }
#    out
#  }
#}
# the R function for noncentral chisq is inaccurate when ncp is large
#pcop.exp.cond1<-function(u, v=.5, theta=.8){
#  if(theta==0) u
#  else {
#    x<--log(1-u)/(1-theta)
#    y<--log(1-v)/(1-theta)
#    out<-pchisq(2*theta*y, 2, 2*x, FALSE)  
#    out<-out+2*(theta*(1-v)*dchisq(2*theta*y,2,2*x)-(1-u)*dchisq(2*y,2,2*theta*x))/(1-theta)/(1-v) 
#    out
#  }
#}
pcop.exp.cond<-function(u, v=.5, theta=.8){
  if(theta==0) u
  else {
    x<--log(1-u)/(1-theta)
    y<--log(1-v)/(1-theta)
    out<-pchisq(2*x,2,2*theta*y) 
    out
  }
}
#u<-seq(0,1,length=20)
#pcop.exp.cond(u, theta=.1)
#pcop.exp.cond2(u, theta=.1)
#' @keywords internal
#' @noRd
#qcop.exp.cond<-function(p, v=.5, theta=.8){
#  np<-length(p)
#  a<-min(p[p!=0:1])/2
#  b<-(1+max(p[p!=0:1]))/2
#  q<-rep(0,np)
#  for(i in 1:np){
#    if(p[i]==0 || p[i]==1) q[i]<-p[i]
#    else q[i]<-uniroot(function(u) pcop.exp.cond(u, v=v, theta=theta)-p[i], c(a,b))$root
#  }
#  q
#}
qcop.exp.cond<-function(p, v=.5, theta=.8){
  y<--log(1-v)/(1-theta)
  q<-1-exp(-.5*(1-theta)*qchisq(p, 2, 2*theta*y))
  q
}
#u<-seq(0,1,length=20)
#qcop.exp.cond(u, theta=.1)
#qcop.exp.cond1(u, theta=.1)

#' @keywords internal
#' @noRd
#rcop.exp.cond<-function(n, v=.5, theta=.8){
#  p<-runif(n)
#  if(theta==0) p
#  else  {
#    qcop.exp.cond(p, v, theta)
#  } 
#}
rcop.exp.cond<-function(n, v=.5, theta=.8){
  p<-runif(n)
  y<--log(1-v)
  if(theta==0) p
  else  {
    x<-(1-theta)*qchisq(p, 2, 2*theta*y/(1-theta))/2
    1-exp(-x)
  } 
}
#par(mfrow=c(1,2))
#hist(rcop.exp.cond(1000), freq=FALSE)
#hist(rcop.exp.cond1(1000), freq=FALSE)


#############################
# FGM
#############################

#' @keywords internal
#' @noRd
dcop.fgm.cond<-function(u, v=.5, theta=.8){
  dcop.fgm(u, v, theta)
}
#' @keywords internal
#' @noRd
pcop.fgm.cond<-function(u, v=.5, theta=.8){
  u*(1+theta*(1-u)*(1-2*v)) 
}
#' @keywords internal
#' @noRd
qcop.fgm.cond<-function(p, v=.5, theta=.8){
  if(theta==0) p
  else{
      A<- theta*(1-2*v) 
      B<--(1+A)
      C<-p 
      apply(round((-B+sqrt(B^2-4*A*C)%o%c(-1,1))/(2*A),10), 1, function(x) x[x>=0 & x<=1])
  }
}
#' @keywords internal
#' @noRd
rcop.fgm.cond<-function(n, v=.5, theta=.8){
  qcop.fgm.cond(runif(n), v, theta)
}


#############################
# Frank
#############################

#' @keywords internal
#' @noRd
dcop.frank.cond<-function(u, v=.5, theta=.8){
  dcop.frank(u, v, theta)
}
#' @keywords internal
#' @noRd
pcop.frank.cond<-function(u, v=.5, theta=.8){
  exp(-theta*v)/(exp(-theta*v)-1+(exp(-theta)-1)/(exp(-theta*u)-1))
}
#' @keywords internal
#' @noRd
qcop.frank.cond<-function(p, v=.5, theta=.8){
  -log(1+p*(exp(-theta)-1)/(p+(1-p)*exp(-theta*v)))/theta
}
#' @keywords internal
#' @noRd
rcop.frank.cond<-function(n, v=.5, theta=.8){
  qcop.frank.cond(runif(n), v, theta)
}



#############################
# Gaussian
#############################

#' @keywords internal
#' @noRd
dcop.gauss.cond<-function(u, v=.5, theta=.8){
  dcop.gauss(u, v, theta)
}
#' @keywords internal
#' @noRd
pcop.gauss.cond<-function(u, v=.5, theta=.8){
  pnorm((qnorm(u)-theta*qnorm(v))/sqrt(1-theta^2))
}
#' @keywords internal
#' @noRd
qcop.gauss.cond<-function(p, v=.5, theta=.8){
  pnorm(sqrt(1-theta^2)*qnorm(p)+theta*qnorm(v))
}
#' @keywords internal
#' @noRd
rcop.gauss.cond<-function(n, v=.5, theta=.8){
  pnorm(sqrt(1-theta^2)*rnorm(n)+theta*qnorm(v))
}


#############################
# Gumbel
#############################

#' @keywords internal
#' @noRd
dcop.gumbel.cond<-function(u, v=.5, theta=1.8){
  dcop.gumbel(u, v, theta)
}
#' @keywords internal
#' @noRd
pcop.gumbel.cond<-function(u, v=.5, theta=1.8){
  e<-((-log(u))^theta+(-log(v))^theta)
  (-log(v))^(theta-1)*exp(-e^(1/theta))/e^(1-1/theta)/v
}
#' @keywords internal
#' @noRd
qcop.gumbel.cond<-function(p, v=.5, theta=1.8){
  np<-length(p)
  a<-min(p[p!=0:1])/2
  b<-(1+max(p[p!=0:1]))/2
  q<-rep(0,np)
  for(i in 1:np){
    if(p[i]==0 || p[i]==1) q[i]<-p[i]
    else q[i]<-uniroot(function(u) pcop.gumbel.cond(u, v=v, theta=theta)-p[i], c(a,b))$root
  }
  q
}
#' @keywords internal
#' @noRd
rcop.gumbel.cond<-function(n, v=.5, theta=1.8){
  qcop.gumbel.cond(runif(n), v, theta)
}



#############################
# Independence
#############################

#' @keywords internal
#' @noRd
dcop.indep.cond<-function(u, v=.5){
  dunif(u)
}
#' @keywords internal
#' @noRd
pcop.indep.cond<-function(u, v=.5){
  punif(u)
}
#' @keywords internal
#' @noRd
qcop.indep.cond<-function(p, v=.5){
  qunif(p)
}
#' @keywords internal
#' @noRd
rcop.indep.cond<-function(n, v=.5){
  runif(n)
}


#############################
# Joe
#############################

#' @keywords internal
#' @noRd
dcop.joe.cond<-function(u, v=.5, theta=1.8){
  dcop.joe(u, v, theta)
}
#' @keywords internal
#' @noRd
pcop.joe.cond<-function(u, v=.5, theta=1.8){
  (1-v)^(theta-1)*(1-(1-u)^theta)/((1-u)^theta+(1-v)^theta-((1-u)*(1-v))^theta)^(1-1/theta)
}
#' @keywords internal
#' @noRd
qcop.joe.cond<-function(p, v=.5, theta=1.8){
  np<-length(p)
  a<-min(p[p!=0:1])/2
  b<-(1+max(p[p!=0:1]))/2
  q<-rep(0,np)
  for(i in 1:np){
    if(p[i]==0 || p[i]==1) q[i]<-p[i]
    else q[i]<-uniroot(function(u) pcop.joe.cond(u, v=v, theta=theta)-p[i], c(a,b))$root
  }
  q
}
#' @keywords internal
#' @noRd
rcop.joe.cond<-function(n, v=.5, theta=1.8){
  qcop.joe.cond(runif(n), v, theta)
}

#############################
# Nakagami-m
#############################

#' @keywords internal
#' @noRd
dcop.nakagami.cond<-function(u, v=.5, theta=.8, m=2){
  dcop.nakagami(u, v, theta, m)
}
#' @importFrom stats pgamma qgamma
#' @keywords internal
#' @noRd
pcop.nakagami.cond<-function(u, v=.5, theta=.8, m=2){
  t<-theta
  cn<-t/(1-t) 
  lcn<-log(cn)
  x<-qgamma(u,m)/(1-t)
  y<-qgamma(v,m) 
  tmpf<-function(k)
    pgamma(x,m+k, log.p=TRUE)+lcn*k-cn*y-lgamma(k+1)
  out<-exp(tmpf(0))
  eps<-1e-10
  del<-exp(tmpf(1))*y
  k<-1
  while(max(del, na.rm =TRUE)>eps){
    out<-out+del
    k<-k+1
    del<-exp(tmpf(k))*y^k
    #cat(max(del, na.rm =TRUE),"\n")
  }
  out
}
u<-seq(0,1,length=20)
pcop.nakagami.cond(u, v=.5, theta=.8, m=2)


#' @keywords internal
#' @noRd
qcop.nakagami.cond<-function(p, v=.5, theta=.8, m=2){
  np<-length(p)
  a<-min(p[p!=0:1])/2
  b<-(1+max(p[p!=0:1]))/2
  q<-rep(0,np)
  for(i in 1:np){
    if(p[i]==0 || p[i]==1) q[i]<-p[i]
    else q[i]<-uniroot(function(u) pcop.nakagami.cond(u, v=v, theta=theta, m=m)-p[i], c(a,b))$root
  }
  q
}
#' @keywords internal
#' @noRd
rcop.nakagami.cond<-function(n, v=.5, theta=.8, m=2){
  qcop.nakagami.cond(runif(n), v, theta, m)
}


#############################
# Plackett
#############################

#' @keywords internal
#' @noRd
dcop.plackett.cond<-function(u, v=.5, theta=1.8){
  dcop.plackett(u, v, theta)
}
#' @keywords internal
#' @noRd
pcop.plackett.cond<-function(u, v=.5, theta=1.8){
  .5-.5*(1+(theta-1)*(u+v)-2*u*theta)/sqrt((1+(theta-1)*(u+v))^2-4*theta*(theta-1)*u*v)
}
#' @keywords internal
#' @noRd
#qcop.plackett.cond<-function(p, v=.5, theta=1.8){
#  np<-length(p)
#  a<-min(p[p!=0:1])/2
#  b<-(1+max(p[p!=0:1]))/2
#  q<-rep(0,np)
#  for(i in 1:np){
#    if(p[i]==0 || p[i]==1) q[i]<-p[i]
#    else q[i]<-uniroot(function(u) pcop.plackett.cond(u, v=v, theta=theta)-p[i], c(a,b))$root
#  }
#  q
#}
qcop.plackett.cond<-function(p, v=.5, theta=1.8){
  if(theta==1) p
  else{
      y<-1+(theta-1)*v
      A<--4*(p*(1-p)*(theta-1)^2+theta)
      B<-2*((theta-1)*(1-2*p)^2*(y-2*theta*v)+(theta+1)*y)
      C<--4*p*(1-p)*y^2
      q<-apply(round((-B+sqrt(B^2-4*A*C)%o%c(-1,1))/(2*A),10), 1, function(x) x[x>=0 & x<=1])
      (p<.5)*apply(q,2,min)+(p>=.5)*apply(q,2,max)
  }
}
#p<-seq(0,1,length=50)
#v<-.3; theta<-.4
#qcop.plackett.cond(p, v, theta)
#qcop.plackett.cond1(p, v, theta)

#' @keywords internal
#' @noRd
rcop.plackett.cond<-function(n, v=.5, theta=1.8){
  qcop.plackett.cond(runif(n), v, theta)
}

#############################
#  Student's t
#############################

#' @keywords internal
#' @noRd
dcop.t.cond<-function(u, v=.5, theta=.7, df=4){
  dcop.t(u, v, theta, df)
}
#' @importFrom stats pt 
#' @keywords internal
#' @noRd
pcop.t.cond<-function(u, v=.5, theta=.7, df=4){
  tv<-qt(v,df)
  tu<-qt(u,df)
  pt(sqrt((df+1)/(1-theta^2)/(df+tv^2))*(tu-theta*tv),df=df+1)
}
#' @keywords internal
#' @noRd
qcop.t.cond<-function(p, v=.5, theta=.7, df=4){
  tv<-qt(v,df)
  pt(qt(p,df=df+1)*sqrt((1-theta^2)*(df+tv^2)/(df+1))+theta*tv,df)
}
#' @keywords internal
#' @noRd
rcop.t.cond<-function(n, v=.5, theta=.7, df=4){
  qcop.t.cond(runif(n), v, theta, df)
}
