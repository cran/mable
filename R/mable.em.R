mable.em <-
function(x, m, a=0, b=1, maxit = 50000, eps = 1.0e-9){
    n<-length(x)
    if(a>=b) stop("a must be smaller than b")
    x<-(x-a)/(b-a)
    llik<-0
    p<-rep(1,m+1)/(m+1)
    ## Call C mable_em
    res<-.C("mable_em",
      as.integer(m), as.integer(n), as.double(p), as.double(x), as.integer(maxit),
      as.double(eps),  as.double(llik))
    ans<-list(phat=res[[3]], llik=res[[7]], support=c(a,b))
    ans
}
