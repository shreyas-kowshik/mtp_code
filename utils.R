# AR(1) generation function
"
n : Number of points to sample
k : Changepoint location
m1 : Mean before changepoint
m2 : Mean after changepoint
rho : AR(1) coefficient
sigma : WN(0,sigma^2)
"
ar1_data<-function(n, k, m1, m2, rho, sigma){
    Y=arima.sim(model=list(ar=c(rho)),n=n,sd=sigma)  
    X=Y*0
    X[1:k]=Y[1:k]+m1
    X[(k+1):n]=Y[(k+1):n]+m2 
    X
}
