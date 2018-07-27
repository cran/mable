mable.optim <-
function(x, m0, m1, a=0, b=1, maxit = 50000, eps = 1.0e-9, tini=1e-4){
    n<-length(x)
    if(a>=b) stop("a must be smaller than b")
    x<-(x-a)/(b-a)
    llik<-0;
    Llik<-rep(0, m1-m0+1)
    lr<-rep(0,m1-m0)
    p<-rep(1, m1+1)
    opti.m<-0
    ## Call C mable_optim
    res<-.C("mable_optim",
      as.integer(m0), as.integer(m1), as.integer(n), as.double(p),  as.double(x),  as.integer(maxit),
      as.double(eps), as.double(Llik), as.double(lr), as.integer(opti.m), as.double(tini))
    opti.m<-res[[10]]
    llik<-res[[8]][opti.m-m0+1]
    mhat<-opti.m+1
    ans<-list(Llik=res[[8]],  phat=res[[4]][1:mhat], opti.m=opti.m, llik=llik, lr=res[[9]], support=c(a,b))
    ans
}
