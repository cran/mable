#############################################################
## Maximum Approximate Bernstein Likelihood Estimate
##  of Multivariate Density Function
##
## References:
##  Wang, T. and Guan, Z.,(2019) Bernstein Polynomial Model    
##     for Nonparametric Multivariate Density, Statistics,     
##              Vol. 53, no. 2, 321-338                        
#############################################################
# setwd("C:\\Users\\zguan\\Documents\\papers\\proportional odds ratio\\parametrization using bernstein polynomials\\multivariate-bernstein-poly-model")
# setwd("C:\\Documents and Settings\\zguan\\My Documents\\papers\\proportional odds ratio\\parametrization using bernstein polynomials\\multivariate-bernstein-poly-model")
# dyn.load("mable-multivar")
## Call C functions by .C
# R --arch x64 CMD SHLIB mable-multivar.c
# R CMD SHLIB mable-multivar.c
#############################################################
#' Maximum Approximate Bernstein Likelihood Estimate
#'  of Multivariate Density Function
#' @param x an \code{n x d} matrix or \code{data.frame} of multivariate sample of size \code{n}
#' @param M a positive integer or a vector of \code{d} positive integers specify  
#'    the maximum candidate or the given model degrees for the joint density.
#' @param search logical, whether to search optimal degrees using \code{M} as maximum candidate 
#'    degrees or not but use \code{M} as the given model degrees for the joint density.
#' @param interval a vector of two endpoints or a \code{d x 2} matrix, each row containing 
#'    the endpoints of support/truncation interval for each marginal density.
#'    If missing, the i-th row is assigned as \code{c(min(x[,i]), max(x[,i]))}.
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit
#' and the convergence criterion \code{eps}. Default is \code{\link{mable.ctrl}}. See Details.
#' @param progress if TRUE a text progressbar is displayed
#' @details 
#'   A \eqn{d}-variate density \eqn{f} on a hyperrectangle \eqn{[a, b]
#'   =[a_1, b_1] \times \cdots \times [a_d, b_d]} can be approximated 
#'   by a mixture of \eqn{d}-variate beta densities on \eqn{[a, b]}, 
#'   \eqn{\beta_{mj}(x) = \prod_{i=1}^d\beta_{m_i,j_i}[(x_i-a_i)/(b_i-a_i)]/(b_i-a_i)},
#'   with proportion \eqn{p(j_1, \ldots, j_d)}, \eqn{0 \le j_i \le m_i, i = 1, \ldots, d}. 
#'   If \code{search=TRUE}, we start with \eqn{M_0}\code{=rep(2,d)} and select  
#'   \eqn{M_1} which maximizes the likelohood with degrees \eqn{M_1=M_0+e_i}, where
#'   \eqn{e_i}, \eqn{i=1,\ldots,d}, form the basis of \eqn{R^n}.  The above procedure
#'   are repeated at least four times to \eqn{M_s} and loglikelihood \eqn{\ell_s}, 
#'   \eqn{s=0, 1,\ldots,k}, and stop whenever the p-value of the change point of
#'   \eqn{\ell_{s+1}-\ell_s} is small or a component of \eqn{M_k} reached its maximum 
#'   value specified by \code{M}. For each \eqn{M_s+e_i} the data are fitted using 
#'   EM algorithm and the multivariate Bernstein polynomial model with a vector of the    
#'   mixture proportions \eqn{p(j_1, \ldots, j_d)}, arranged in the 
#'   column-major order of \eqn{j = (j_1, \ldots, j_d)}, \eqn{p_0, \ldots, p_{K-1}}, 
#'   where \eqn{K=\prod_{i=1}^d (m_i+1)}, to obtain likelihood.
#' @return  A list with components
#' \itemize{
#'  \item \code{dim} the dimension \code{d} of the data
#'  \item \code{m} a vector of the selected optimal degrees by the method of change-point
#'  \item \code{p} a vector of the mixture proportions \eqn{p(j_1, \ldots, j_d)}, arranged in the 
#'   column-major order of \eqn{j = (j_1, \ldots, j_d)}, \eqn{0 \le j_i \le m_i, i = 1, \ldots, d}.
#'  \item \code{mloglik}  the maximum log-likelihood at an optimal degree \code{m}
#'  \item \code{pval}  the p-values of change-points for choosing the optimal degrees for the 
#'    marginal densities
#'  \item \code{M} the vector \code{(m1, m2, ... , md)}, where \code{mi} is the largest candidate 
#'    degree when the search stoped for the \code{i}-th marginal density
#'  \item \code{interval} support hyperrectangle \eqn{[a, b]=[a_1, b_1] \times \cdots \times [a_d, b_d]}
#'  \item \code{convergence} An integer code. 0 indicates successful completion(the EM iteration is   
#'    convergent). 1 indicates that the iteration limit \code{maxit} had been reached in the EM iteration;
#' }
#' @examples
#' ## Old Faithful Data
#' \donttest{
#'  a<-c(0, 40); b<-c(7, 110)
#'  #ans<-mable.mvar(faithful, M = c(100, 100), interval = cbind(a, b))
#'  ans<- mable.mvar(faithful, M = c(46,19), search =FALSE, interval = cbind(a,b))
#'  plot(ans, which="density") 
#'  plot(ans, which="cumulative")
#' }
#' @keywords distribution nonparametric multivariate 
#' @concept multivariate Bernstein polynomial model
#' @concept density estimation
#' @author Zhong Guan <zguan@iusb.edu>
#' @references Wang, T. and Guan, Z.,(2019) Bernstein Polynomial Model for Nonparametric Multivariate Density,    
#'    \emph{Statistics}, Vol. 53, no. 2, 321-338  
#' @export
mable.mvar<-function(x, M, search=TRUE, interval=NULL,  
        controls = mable.ctrl(), progress=TRUE){
    data.name<-deparse(substitute(x))
    n<-nrow(x)
    d<-ncol(x)
    xNames<-names(x)
    if(is.null(xNames)) 
        for(i in 1:d) xNames[i]<-paste("x",i,sep='')
    x<-as.matrix(x)   
    m<-rep(0, d)
    if(is.null(interval))  
        interval<-cbind(apply(x, 2, min), apply(x, 2, max))
    else if(is.matrix(interval)){
        if(nrow(interval)!=d || ncol(interval)!=2)
            stop("Invalid argument 'interval'.\n")
    }
    else{
        if(length(interval)!=2) 
            stop("Invalid argument 'interval'.\n")
        else if(any(apply(x, 2, min)<interval[1])||any(apply(x, 2, max)>interval[2]))
                stop("Invalid argument 'interval'.\n")
        else interval<-matrix(rep(interval, d), nrow=d, byrow=TRUE)
    }
    for(i in 1:d) x[,i]<-(x[,i]-interval[i,1])/(interval[i,2]-interval[i,1])
    if(missing(M) || length(M)==0) stop("'M' is missing.\n")
    else if(length(M)==1) M<-rep(M,d)
    else if(length(M)!=d){
        if(search){
            if(max(M)>=5) warning("Length of 'M' does not match 'ncol(x)'. Use largest 'M'.\n")
            else stop("'M' are too small for searching optimal degrees.\n")
            M<-rep(max(M),d)}
        else stop("Invalide 'M'.\n")
    }
    else M<-M[1:d]
    kmax<-sum(M-2)
    pval<-rep(0, kmax)
    K<-prod(M+1)
    lk<-rep(0, kmax+1)
    lr<-rep(0, kmax)
    chpts<-rep(0, kmax)
    p<-rep(1,K)/K
    convergence<-0
    level<-controls$sig.level
    ## Call C mable_mvar
    res<-.C("mable_mvar",
      as.integer(M), as.integer(n), as.integer(d), as.integer(search), as.double(p), 
      as.double(x), as.integer(controls$maxit), as.double(controls$eps), 
      as.double(lk), as.double(lr), as.double(pval), as.integer(chpts), 
      as.logical(progress), as.integer(convergence), as.double(level))
      m<-res[[1]]
      k<-sum(m-2)
      K<-prod(m+1)
    out<-list(p=res[[5]][1:K], mloglik=res[[9]][max(res[[12]][1:k])], lk=res[[9]][1:(k+1)], 
        lr=res[[10]][1:k], pval=res[[11]][1:k], chpts=res[[12]][1:k], interval=interval, 
        m=m, dim=d, xNames=xNames,  convergence=res[[14]])
    out$data.type<-"mvar"
    class(out)<-"mable"
    return(out)
}
#########################################################
#' Multivariate Mixture Beta Distribution
#' @description Density, distribution function,  and 
#' pseudorandom number generation for the multivariate Bernstein polynomial model, 
#' mixture of multivariate beta distributions, with given mixture proportions 
#' \eqn{p = (p_0, \ldots, p_{K-1})}, given degrees \eqn{m = (m_1, \ldots, m_d)},
#'  and support \code{interval}.
#' @param x a matrix with \code{d} columns or a vector of length \code{d} within 
#'   support hyperrectangle \eqn{[a, b] = [a_1, b_1] \times \cdots \times [a_d, b_d]}
#' @param p a vector of \code{K} values. All components of \code{p} must be 
#'  nonnegative and sum to one for the mixture multivariate beta distribution. See 'Details'. 
#' @param m a vector of degrees, \eqn{(m_1, \ldots, m_d)} 
#' @param n sample size
#' @param interval a vector of two endpoints or a \code{d x 2} matrix, each row containing 
#'    the endpoints of support/truncation interval for each marginal density.
#'    If missing, the i-th row is assigned as \code{c(min(x[,i]), max(x[,i]))}.
#' @details 
#'  \code{dmixmvbeta()} returns a linear combination \eqn{f_m} of \eqn{d}-variate beta densities 
#'  on \eqn{[a, b]}, \eqn{\beta_{mj}(x) = \prod_{i=1}^d\beta_{m_i,j_i}[(x_i-a_i)/(b_i-a_i)]/(b_i-a_i)},   
#'  with coefficients \eqn{p(j_1, \ldots, j_d)}, \eqn{0 \le j_i \le m_i, i = 1, \ldots, d}, where
#'  \eqn{[a, b] = [a_1, b_1] \times \cdots \times [a_d, b_d]} is a hyperrectangle, and the  
#'  coefficients are arranged in the column-major order of \eqn{j = (j_1, \ldots, j_d)}, 
#'  \eqn{p_0, \ldots, p_{K-1}},  where \eqn{K = \prod_{i=1}^d (m_i+1)}. 
#'  \code{pmixmvbeta()} returns a linear combination \eqn{F_m} of the distribution
#'  functions of \eqn{d}-variate beta distribution.
#'
#'  If all \eqn{p_i}'s are nonnegative and sum to one, then \code{p}
#'  are the mixture proportions of the mixture multivariate beta distribution.
#' @importFrom stats rbeta
#' @export
dmixmvbeta<-function(x, p, m, interval=NULL){
    d<-length(m)
    if(any(m-floor(m)!=0) || d==0) stop("Invalid argument 'm'.\n")
    km<-c(1,cumprod(m+1))
    K<-length(p)
    if(K==0) stop("Missing mixture proportions 'p' without default.")
    if(K!=km[d+1]) stop("Invalid argument 'p': length must equal to 'prod(m+1)'.\n")
    if(any(p<0)) warning("Argument 'p' has negative component(s).\n")
    #else if(abs(sum(p)-1)>.Machine$double.eps)
    #    warning("Sum of 'p's is not 1.\n")
    if(is.null(interval))  
        interval<-cbind(rep(0,d), rep(1,d))
    else if(is.matrix(interval)){
        if(nrow(interval)!=d || ncol(interval)!=2)
            stop("Invalid argument 'interval'.\n")
    }
    else{
        if(length(interval)!=2) 
            stop("Invalid argument 'interval'.\n")
        else interval<-matrix(rep(interval, d), nrow=d, byrow=TRUE)
    }
    int<-interval
    if(is.matrix(x)){
        nx<-nrow(x)
        if(d!=ncol(x)) stop("Wrong dim of 'x'.\n")
        if(any(x-matrix(rep(int[,1],nx), ncol=d, byrow=TRUE)<0) 
            || any(x-matrix(rep(int[,2],nx), ncol=d, byrow=TRUE)>0))
            stop("'x' must be in the hyperrectangle 'interval'.\n")
        for(i in 1:d) x[,i]<-(x[,i]-int[i,1])/(int[i,2]-int[i,1])
    }
    else {
        if(d!=length(x)) stop("'x' is a vector. Its length must be the same as 'm'.\n") 
        nx<-1}
#    obj<-list(m=m, p=p, interval=interval)
#    ans<-mvbern.poly(x, obj, density=TRUE)
    fb<-rep(0, nx)
    density<-TRUE
    res<-.C("mable_mvdf",
      as.integer(d), as.integer(m), as.integer(km), as.integer(nx), as.double(x),  
      as.double(p), as.double(fb), as.logical(density))
    fb<-res[[7]]/prod(int[,2]-int[,1])
    return(fb)
}
#' @rdname dmixmvbeta
#' @export
pmixmvbeta<-function(x, p, m, interval=NULL){
    d<-length(m)
    if(any(m-floor(m)!=0) || d==0) stop("Invalid argument 'm'.\n")
    km<-c(1,cumprod(m+1))
    K<-length(p)
    if(K==0) stop("Missing mixture proportions 'p' without default.")
    if(K!=km[d+1]) stop("Invalid argument 'p': length must equal to prod(m+1).\n")
    if(any(p<0)) warning("Argument 'p' has negative component(s).\n")
    #else if(abs(sum(p)-1)>.Machine$double.eps)
    #    warning("Sum of 'p's is not 1.\n")
    if(is.null(interval))  
        interval<-cbind(rep(0,d), rep(1,d))
    else if(is.matrix(interval)){
        if(nrow(interval)!=d || ncol(interval)!=2)
            stop("Invalid argument 'interval'.\n")
    }
    else{
        if(length(interval)!=2) 
            stop("Invalid argument 'interval'.\n")
        else interval<-matrix(rep(interval, d), nrow=d, byrow=TRUE)
    }
    int<-interval
    if(is.matrix(x)){
        nx<-nrow(x)
        if(d!=ncol(x)) stop("Wrong dim of 'x'.\n")
        if(any(x-matrix(rep(int[,1],nx), ncol=d, byrow=TRUE)<0) 
            || any(x-matrix(rep(int[,2],nx), ncol=d, byrow=TRUE)>0))
            stop("'x' must be in the hyperrectangle 'interval'.\n")
        for(i in 1:d) x[,i]<-(x[,i]-int[i,1])/(int[i,2]-int[i,1])
    }
    else {
        if(d!=length(x)) stop("'x' is a vector. Its length must be the same as 'm'.\n") 
        nx<-1}
#    obj<-list(m=m, p=p, interval=interval)
#    ans<-mvbern.poly(x, obj, density=FALSE)
    fb<-rep(0, nx)
    density<-FALSE
    res<-.C("mable_mvdf",
      as.integer(d), as.integer(m), as.integer(km), as.integer(nx), as.double(x),  
      as.double(p), as.double(fb), as.logical(density))
    fb<-res[[7]] 
    return(fb)
}
#########################################################
# generating prn from mixture of multivariate beta 
#   with degrees m=(m1,...,md)
#' @rdname dmixmvbeta
#' @export
rmixmvbeta<-function(n, p, m, interval=NULL){
    d<-length(m)
    if(any(m-floor(m)!=0) || d==0) stop("Invalid argument 'm'.\n")
    K<-length(p)
    if(K==0) stop("Missing mixture proportions 'p' without default.")
    if(K!=prod(m+1)) stop("Invalid argument 'p': length must equal to prod(m+1).\n")
    if(any(p<0)) stop("Negative component(s) of argument 'p' is not allowed.\n")
    else if(abs(sum(p)-1)>.Machine$double.eps){
        warning("Sum of 'p's is not 1. Dividing 'p's by the total.\n")
        p<-p/sum(p)
    }
    if(is.null(interval))  
        interval<-cbind(rep(0,d), rep(1,d))
    else if(is.matrix(interval)){
        if(nrow(interval)!=d || ncol(interval)!=2)
            stop("Invalid argument 'interval'.\n")
    }
    else{
        if(length(interval)!=2) 
            stop("Invalid argument 'interval'.\n")
        else interval<-matrix(rep(interval, d), nrow=d, byrow=TRUE)
    }
    if(d==1) x<-rmixbeta(n, p, interval)
    else{
        km<-cumprod(m+1)
        # column-major order
        ii<-matrix(0, nrow=K, ncol=d)
        for(k in 1:d) ii[,k]<-rep(0:m[k], each=km[k]/(m[k]+1))
        w<-sample(1:K, n, replace = TRUE, prob = p)
        x<-rep(NULL, d)
        for(i in 1:n){
            tmp<-rbeta(d, shape1 = ii[w[i],]+1, shape2 = m-ii[w[i],]+1)
            x<-rbind(x, interval[,1]+(interval[,2]-interval[,1])*tmp)
        }
    }
    return(x)
}
