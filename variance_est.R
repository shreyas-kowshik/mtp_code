source("kernel.R")

# Takes data as input from a time series and returns variance estimate without using transformation of series #
non_transformed_variance_estimate<-function(x, kernel, bandwidth, delta=-1){
    n<-length(x)

    if(delta==-1) {
        delta=n-1
    }

    gamma_sum=array(0,dim=c(n-1))

    w<-get_kernel_weights(n, kernel, bandwidth)

    for(r in 1 :(n-1)){
        a<-as.numeric(acf(x,type=c("covariance"),lagax=delta,plot = F,demean = T)$acf)
        a<-((a)[1:min(length(a), delta+1)]) 
        aa<-c(a,rep(0,(n-length(a)))) * w

        M<-array(0,dim=c(n,n))
        for(i in 1: (n)){
            for(j in 1:(n)) {
                if (is.na(aa[abs(i-j)+1])) {
                    M[i,j]<-0.0
                }
                else {
                    M[i,j]<-aa[abs(i-j)+1]
                }
            }
        }

        v<-c(rep((1-r/n),r), rep((-r/n),(n-r)))
        gamma_sum[r]<-t(v)%*%M%*%v/n
    }
    gamma_sum
}

# Utility functions for Br consistency regularization #

# Takes sample covariance matrices M_x and M_z without and with transformations and returns Br_tilde leading to #
# the closest quadratic form value without the transformation #
get_br_tilde_opt<-function(M_x, M_z, ar, eps, Br, lambdas) {
    best_br<-1.0
    eps<-1e-6
    min_val<-10000000.0
    
    for (i in 1:length(lambdas)) {
        k<-lambdas[i]
        flag<-1
        Br_tilde<-Br+k*diag(n)

        tryCatch(
        expr = {
            x_qcqp<-solve(Br_tilde,ar)
        },
        error = function(e){ 
            # (Optional)
            # Do this if an error is caught...
            print('caught')
            flag<-0
        },
        warning = function(w){
            # (Optional)
            # Do this if an warning is caught...
            # print('warning')
        },
        finally = {
            # (Optional)
            # Do this at the end before quitting the tryCatch structure...
            # print('finally')
        }
        )
    
        # print(flag)
        if(flag==0) {
            next
        }

        Brx<-(Br_tilde)%*%x_qcqp
        vals<-t(x_qcqp)%*%M_z%*%x_qcqp - t(ar)%*%M_x%*%ar
        vals<-abs(vals)
        constraint_vals<-sum((Brx - ar)^2)

        if (constraint_vals<eps && vals<min_val) {
            min_val<-vals
            best_br<-(Br+k*diag(n))
        }
    }
    
    best_br
}

transformed_variance_estimate<-function(x, lambdas, kernel, bandwidth, delta=-1) {
    n<-length(x)

    if(delta==-1) {
        delta=n-1
    }

    gamma_sum=array(0,dim=c(n-1))
    num_steps_skip<-0

    w<-get_kernel_weights(n, kernel, bandwidth)

    a<-as.numeric(acf(x,type=c("covariance"),lagax=delta,plot = F,demean = T)$acf)
    a<-((a)[1:min(length(a), delta+1)]) 
    aa<-c(a,rep(0,(n-length(a)))) * w # Multiply by w to enforce weights

    M<-array(0,dim=c(n,n))
    for(i in 1: (n)){
        for(j in 1:(n)) {
            if (is.na(aa[abs(i-j)+1])) {
                M[i,j]<-0.0
            }
            else {
                M[i,j]<-aa[abs(i-j)+1]
            }
        }
    }
    M_x<-M

    for(r in 1 :(n-1)){
        first_r<-r-num_steps_skip
        z<-array(0,dim=c(n))

        # Build Br matrix for transformation #
        Br<-array(0,dim=c(n,n))
        r_mat<-array(-1.0/r,dim=c(max(first_r,1),max(first_r,1)))
        Br[1:max(first_r,1),1:max(first_r,1)]<-r_mat
        nr_mat<-array(-1.0/(n-r),dim=c(n-first_r,n-first_r))
        Br[(first_r+1):n,(first_r+1):n]<-nr_mat
        Br=Br+diag(n)
            
        z<-Br%*%x
            
        a<-as.numeric(acf(z,type=c("covariance"),lagax=delta,plot = F,demean = F)$acf)
        a<-((a)[1:min(length(a), delta+1)]) 
        aa<-c(a,rep(0,(n-length(a)))) * w

        M<-array(0,dim=c(n,n))
        for(i in 1: (n)){
            for(j in 1:(n)) {
                if (is.na(aa[abs(i-j)+1])) {
                    M[i,j]<-0.0
                }
                else {
                    M[i,j]<-aa[abs(i-j)+1]
                }
            }
        }
        M_z<-M


        # Construct v to multiply with covariance matrix #
        # B_r*v = ar (v is the solution of this)
        ar<-c(rep((1-r/n),r), rep((-r/n),(n-r)))
            
        # Get Br_tilde matrix from optimization function #
        Br_tilde<-get_br_tilde_opt(M_x,M_z,ar,1e-6,Br,lambdas)
        v<-solve(Br_tilde,ar)
            
        if(abs(sum((Br_tilde%*%v)-ar))>1e-6) {
            stop("Br*v not equal to ar")
        }

        gamma_sum[r]<-t(v)%*%M_z%*%v/n
    }
    return(c(gamma_sum, M_z, Br_tilde)) 
}
