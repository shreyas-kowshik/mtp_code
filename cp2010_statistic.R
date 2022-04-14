Gn<-function(x){
    n<-length(x)

    t<-array(0,dim=c(2))
    xc<-x-mean(x) 
    tnk<-(cumsum(xc)/sqrt(n))
    tnk<-tnk[1:n-1]
    vnk<-numeric(0)

    for(k in 1: (n-1)){
        ss<-1:k  # index 
        sf<-x[ss]-mean(x[ss]) # first part 
        sb<-x[-ss]-mean(x[-ss]) # end part 
        vnk[k]<-(sum((cumsum(sf))^2)+ sum((cumsum(sb))^2))/n^2
    }

    st<-tnk*tnk/vnk
    t[1]<-which(st==max(st))
    t[2]<-st[t[1]]
    t
}
