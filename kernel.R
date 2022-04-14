# Official Code #
"
n : Number of points
bandwidth : Kernel bandwidth value
kernel : One of [tr, ba, pa, bo, da, qs]
"
getLongRunWeights <- function(n, bandwidth, kernel) {
  w <- numeric(n - 1)
  bw <- bandwidth

  if (kernel == "tr") {
    w <- w + 1
    upper <- min(bw, n - 1)
  }

  else if (kernel == "ba") {
    upper <- ceiling(bw) - 1
    if (upper > 0) {
      j <- 1:upper
    } else {
      j <- 1
    }
    w[j] <- 1 - j / bw
  }

  else if (kernel == "pa") {
    upper1 <- floor(bw / 2)
    if (upper1 > 0) {
      j <- 1:upper1
    } else {
      j <- 1
    }
    jj <- j / bw
    w[j] <- 1 - 6 * jj^2 + 6 * jj^3
    j2 <- (floor(bw / 2) + 1):bw
    jj2 <- j2 / bw
    w[j2] <- 2 * (1 - jj2)^3
    upper <- ceiling(bw) - 1
  }

  else if (kernel == "bo") {
    upper <- ceiling(bw) - 1
    if (upper > 0) {
      j <- 1:upper
    } else {
      j <- 1
    }
    jj <- j / bw
    w[j] <- (1 - jj) * cos(pi * jj) + sin(pi * jj) / pi
  }

  else if (kernel == "da") {
    upper <- n - 1
    j <- 1:upper
    w[j] <- sin(pi * j / bw) / (pi * j / bw)
  }

  else if (kernel == "qs") {
    sc <- 1.2 * pi
    upper <- n - 1
    j <- 1:upper
    jj <- j / bw
    w[j] <- 25 / (12 * pi^2 * jj^2) * (sin(sc * jj) / (sc * jj) - cos(sc * jj))
  }

  if (upper <= 0) upper <- 1

  return(list(w = w, upper = upper))
}

# Set up kernel weights #
get_kernel_weights<-function(n, kernel, bandwidth, delta=-1) {
    if(delta==-1) {
        delta=n-1
    }

    weights<-getLongRunWeights(n, bandwidth, kernel=kernel)
    w <- weights[[1]]
    upper <- weights[[2]]

    w<-c(1.0, w)[1:(delta+1)]
    w<-c(w,rep(0,(n-1-delta)))

    # Normalize #
    norm_factor<-1.0
    if (kernel == "tr") {
        norm_factor<-0.5
    }
    
    else if (kernel=="ba") {
        norm_factor<-1.5
    }

    else if (kernel=="pa") {
        norm_factor<-(1.0/0.539285)
    }

    else {
        stop("Kernel not implemented")
    }

    w<-w*norm_factor
    w
}
