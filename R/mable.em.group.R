mable.em.group <-
function(x, breaks, m,  a=0, b=1, maxit = 50000, eps = 1.0e-9){
    N<-length(x)
    llik<-0
    p<-rep(1,m+1)/(m+1)
    t<-(breaks-a)/(b-a)
    if(any(t<0) | any(t>1)) stop("components of 'breaks' must be in (a,b)")
    ## Call C mable_em_group
    res<-.C("mable_em_group",
      as.integer(m), as.double(x), as.integer(N), as.double(p), as.double(t), as.integer(maxit),
      as.double(eps),  as.double(llik))
    ans<-list(phat=res[[4]], llik=res[[8]], support=c(a,b))
    ans
}
