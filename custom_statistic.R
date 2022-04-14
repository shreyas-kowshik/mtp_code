source("variance_est.R")
source("utils.R")

generate_covar_matrix_ar1<-function(n, rho, sigma) {
    M<-array(sigma^2/(1-rho^2),dim=c(n,n))
    for(i in 1: (n)){
        for(j in 1:(n)) {
            M[i,j]<-M[i,j]*rho^(abs(i-j)) #*(abs(i-j)<=delta)
        }
    }
    M
}

# Exact variance of numerator for AR(1) process
variance_ar1<-function(n, rho, sigma, M=-1){
    # M is passed for caching #
    if(M==-1) {
        M<-generate_covar_matrix_ar1(n, rho, sigma)
    }

    gamma_sum=array(0,dim=c(n-1))
    for(r in 1 :(n-1)){
        v<-c(rep((1-r/n),r), rep((-r/n),(n-r)))
        gamma_sum[r]<-t(v)%*%M%*%v/n
    }
    gamma_sum
}

# New statistic #
Hn_ar1_exact<-function(x, sigma, phi) {
    n<-length(x)

    t<-array(0,dim=c(2))
    xc<-x-mean(x)
    tnr<-(cumsum(xc)/n)[1:n-1]
    
    # gamma_sum<-(sigma/(1 - phi))^2
    den<-variance_ar1(n, phi, sigma)
    
    vals<-abs(sqrt(n) * tnr)/(sqrt(den))
    t[1]<-which(vals==max(vals))
    t[2]<-vals[t[1]]
    t
}

Hn_estimated<-function(x, lambdas, kernel, bandwidth, delta=-1, epsilon=1e-6) {
    n<-length(x)

    if(delta==-1) {
        delta=n-1
    }

    t<-array(0,dim=c(2))
    xc<-x-mean(x)
    tnr<-(cumsum(xc)/n)[1:n-1]
    
    out<-transformed_variance_estimate(x, lambdas, kernel, bandwidth)
    den<-out[1:n-1] # Denominator val
    M_z<-resize(out[n:n+(n*n)-1],n,n)
    Br_tilde<-resize(out[n+(n*n):length(out)],n,n)

    # Add epsilon
    den<-den+epsilon
    
    vals<-abs(sqrt(n) * tnr)/(sqrt(den))
    t[1]<-which(vals==max(vals))
    t[2]<-vals[t[1]]
    t
}
