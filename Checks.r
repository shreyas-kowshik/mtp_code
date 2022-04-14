
setwd('/media/shreyas/Data/mtp/Sem-10/workstation_files/shreyas_changepoint')

library("matrixStats")
library("cointReg")
library("ramify")

# Include files #
source("utils.R")
source("kernel.R")
source("variance_est.R")
source("custom_statistic.R")
source("cp2010_statistic.R")

n<-200
k<-100
rho<-0.7
sigma<-0.5
m1<-1.0
k<-as.integer(n/2) + 10

lambdas<-seq(0.1,3.0,by=0.2)
kernel<-"ba"
bandwidth<-9.0

# Simulate and compute statistics #
simulate<-function(n,k,mu1,mu2,itrn1,rho,sigma){
  # One array per statistic
  g<-array(0,dim=c(itrn1,2))
  h<-array(0,dim=c(itrn1,2)) # For custom statistic

  for (i in  1 : itrn1){
    x<-gdata(n,k,mu1,mu2,rho,sigma)
    g[i,]<-gn(x,n)
    # Add new statistic here
    h[i,]<-Hn(x,n,sigma,rho)
  }
 
  list(g, h)
}


null_simul_hn<-function(n,k,mu1,itrn1,rho,sigma){
  q<-c(0.9,0.95,0.99)

  # One array per statistic
  h<-array(0,dim=c(itrn1,2)) # For custom statistic

  for (i in  1 : itrn1){
    x<-ar1_data(n, k, mu1, mu1, rho, sigma)
    # Add new statistic here
    h[i,]<-Hn_estimated(x, lambdas, kernel, bandwidth)
  }
  h
}

#------Compute power-------
comp_pow<-function(n,k,m1,m2,itrn2,u1,u2,rho,sigma){
  # One array per statistic
  v<-array(0,dim = c(itrn2,3))
  h<-array(0,dim = c(itrn2,3))

  for (i in  1 : itrn2){
    x<-ar1_data(n, k, m1, m2, rho, sigma)
    v[i,1:2]<-Gn(x)
    v[i,3]<-(v[i,2]>u1)
      
    # Add a new statistic here
    h[i,1:2]<-Hn_estimated(x, lambdas, kernel, bandwidth)
    h[i,3]<-(h[i,2]>u2)
  }
  
  list(v, h)
  
}


### Parallel computation code ###

library("parallel")

n<-200
k<-100
rho<-0.7
sigma<-0.5
m1<-1.0
k<-as.integer(n/2) + 10

lambdas<-seq(0.1,3.0,by=0.2)
kernel<-"ba"
bandwidth<-9.0
itrn2<-2

hn_sim_stats<-readRDS(file="hn_estimated_sims_500_v2.RDS")
gn_sim_stats_<-readRDS(file="sim_stats_null_G_v2.RDS")
gn_sim_stats<-gn_sim_stats_[17,,]
q<-c(0.9,0.95,0.99)
hcut<-quantile(hn_sim_stats[,2], probs=q)
gcut<-quantile(gn_sim_stats[,2], probs=q)
u1<-gcut[2]
u2<-hcut[2]

mm<-seq(0.0,2.0,by=0.10)

alt_sim_stats_gn<-array(0, dim=c(length(mm), itrn2, 3))
alt_sim_stats_hn<-array(0, dim=c(length(mm), itrn2, 3))
alt_sim_stats_gnhn<-array(0, dim=c(length(mm), 2, itrn2, 3))

parallel_simulation<-function(i) {
    # print(i)
    # print(mm[i])
    stats<-comp_pow(n,k,m1,mm[i],itrn2,gcut[2],hcut[2],rho,sigma)
    # alt_sim_stats_gn[i,,]<-stats[[1]]
    # alt_sim_stats_hn[i,,]<-stats[[2]]
    alt_sim_stats_gnhn[i,1,,]<-stats[[1]]
    alt_sim_stats_gnhn[i,2,,]<-stats[[2]]
    
    # return(alt_sim_stats_gn, alt_sim_stats_hn)
    return(alt_sim_stats_gnhn)
}

input_arr<-array(0, dim=c(length(mm)))
for (i in 1:length(mm)) {
    input_arr[i]<-i
}

begin<-Sys.time()
r <- mclapply(input_arr, parallel_simulation, mc.cores=6)      ## Split this job across 10 cores
end<-Sys.time()
print(end-begin)

for (i in 1:length(mm)) {
    alt_sim_stats_gn[i,,]<-r[[i]][i,1,,]
    alt_sim_stats_hn[i,,]<-r[[i]][i,2,,]
    alt_sim_stats_gnhn[i,1,,]<-r[[i]][i,1,,]
    alt_sim_stats_gnhn[i,2,,]<-r[[i]][i,2,,]
}

pow_gn<-array(0, dim=c(length(mm)))
pow_hn<-array(0, dim=c(length(mm)))

for (i in 1:length(mm)) {
    pow_gn[i]<-mean(alt_sim_stats_gn[i,,3])
    pow_hn[i]<-mean(alt_sim_stats_hn[i,,3])
}
# plot(mm, mean(alt_sim_stats_gn[,,3]))

saveRDS(alt_sim_stats_gnhn, file = "alt_sim_stats_gnhn_500_14April.RDS")

plot(mm, pow_gn, type='l', col='red')
lines(mm, pow_hn, col='green')



pow_gn

pow_hn

