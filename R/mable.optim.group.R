mable.optim.group <-
function(x, breaks, m0, m1,  a=0, b=1, maxit = 50000, eps = 1.0e-9){
    N<-length(x)
    llik<-0
    Llik<-rep(0, m1-m0+1)
    lr<-rep(0,m1-m0)
    p<-rep(1, m1+1)
    opti.m<-0
    t<-(breaks-a)/(b-a)
    if(any(t<0) | any(t>1)) stop("components of 'breaks' must be in (a,b)")
    ## Call C mable_optim_group
    res<-.C("mable_optim_group",
      as.integer(m0), as.integer(m1), as.integer(N), as.double(p), as.double(t), as.double(x),  as.integer(maxit),
      as.double(eps), as.double(Llik), as.double(lr),   as.double(llik),   as.integer(opti.m))
    opti.m<-res[[12]]
    m<-opti.m+1
    ans<-list(Llik=res[[9]],  phat=res[[4]][1:m], llik=res[[11]], opti.m=opti.m, lr=-res[[10]], support=c(a,b))
    ans
}
