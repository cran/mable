bern.approx <-
function(x, p, a=0, b=1, density=TRUE){
    m<-length(p)-1
    n<-length(x)
    if(a>=b) stop("a must be smaller than b")
    u<-(x-a)/(b-a)
    cdf<-!density
    res<-.C("mable_approx",
      as.double(u), as.double(p), as.integer(m), as.integer(n), as.integer(cdf))
    if(density) res[[1]]/(b-a)
    else res[[1]]
}
