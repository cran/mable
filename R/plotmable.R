plotmable <-
function(mable.fit, density=TRUE, nx=512, add = TRUE, ...){
   phat<-mable.fit$phat
   support<-mable.fit$support
   m<-length(phat)-1
   a<-support[1]
   b<-support[2]
   x<-seq(a, b,len=nx)
   if(add) lines(x, bern.approx(x, phat, a, b, density), ...)
   else plot(x, bern.approx(x, phat, a, b, density), ...)
}
