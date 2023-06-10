source('tropical_clustering.R')

######## Experiments
res <- matrix(rep(0, 600), 6, 100)
res2 <- matrix(rep(0, 600), 6, 100)

### R = 0.25
T1 <- ape::read.nexus("data/R025genetrees_S1.dat")
T2 <- ape::read.nexus("data/R025genetrees_S2.dat")

L <- T2[[1]]$tip.label
n <- 10 ## number of leaves
m <- n - 1 ## m is the number of leaves in subtrees m-dissimilarity map
d <- choose(n, m)
star.tree <- rep(2, d)
N1 <- length(T1) ## sample size of the data set
D1 <- matrix(rep(0, N1*choose(n, 2)), N1, choose(n, 2))
N2 <- length(T2) ## sample size of the data set
D2 <- matrix(rep(0, N2*choose(n, 2)), N2, choose(n, 2))

imputeSize1 <- round(N1*0.1)
imputeSize2 <- round(N2*0.1)

TT1 <- list()
TT2 <- list()

for(i in 1:N1){
  u <- force.ultrametric(T1[[i]], "nnls")
  u$edge.length <- u$edge.length/max(u$edge.length)
  TT1[[i]] <- u
  
  D0 <- cophenetic(u)[L, ]
  D0 <- 2*cophenetic(u)[L, ]/max(cophenetic(u)[L, ])
  D1[i, ] <- normaliz.tree(D0[lower.tri(t(D0))], 2)
}

for(i in 1:N2){
  u <- force.ultrametric(T2[[i]], "nnls")
  u$edge.length <- u$edge.length/max(u$edge.length)
  TT2[[i]] <- u
  
  D0 <- cophenetic(u)[L, ]
  D0 <- 2*cophenetic(u)[L, ]/max(cophenetic(u)[L, ])
  D2[i, ] <- normaliz.tree(D0[lower.tri(t(D0))], 2)
}

for(i in 1:100){
  index1 <- sample(1:N1, 20)
  index2 <- sample(1:N2, 20)
  D <- rbind(D1[index1,], D2[index2,])
  
  res[1, i] <- Tropical.DIANA(D, n)
  res2[1, i] <- DIANA(D, n)
}

### R = 0.5
T1 <- ape::read.nexus("data/R05genetrees_S1.dat")
T2 <- ape::read.nexus("data/R05genetrees_S2.dat")

L <- T2[[1]]$tip.label
n <- 10 ## number of leaves
m <- n - 1 ## m is the number of leaves in subtrees m-dissimilarity map
d <- choose(n, m)
star.tree <- rep(2, d)
N1 <- length(T1) ## sample size of the data set
D1 <- matrix(rep(0, N1*choose(n, 2)), N1, choose(n, 2))
N2 <- length(T2) ## sample size of the data set
D2 <- matrix(rep(0, N2*choose(n, 2)), N2, choose(n, 2))

imputeSize1 <- round(N1*0.1)
imputeSize2 <- round(N2*0.1)

TT1 <- list()
TT2 <- list()

for(i in 1:N1){
  u <- force.ultrametric(T1[[i]], "nnls")
  u$edge.length <- u$edge.length/max(u$edge.length)
  TT1[[i]] <- u
  
  D0 <- cophenetic(u)[L, ]
  D0 <- 2*cophenetic(u)[L, ]/max(cophenetic(u)[L, ])
  D1[i, ] <- normaliz.tree(D0[lower.tri(t(D0))], 2)
}

for(i in 1:N2){
  u <- force.ultrametric(T2[[i]], "nnls")
  u$edge.length <- u$edge.length/max(u$edge.length)
  TT2[[i]] <- u
  
  D0 <- cophenetic(u)[L, ]
  D0 <- 2*cophenetic(u)[L, ]/max(cophenetic(u)[L, ])
  D2[i, ] <- normaliz.tree(D0[lower.tri(t(D0))], 2)
}

for(i in 1:100){
  index1 <- sample(1:N1, 20)
  index2 <- sample(1:N2, 20)
  D <- rbind(D1[index1,], D2[index2,])
  
  res[2, i] <- Tropical.DIANA(D, n)
  res2[2, i] <- DIANA(D, n)
}

### R = 1
T1 <- ape::read.nexus("data/R1genetrees_S1.dat")
T2 <- ape::read.nexus("data/R1genetrees_S2.dat")

L <- T2[[1]]$tip.label
n <- 10 ## number of leaves
m <- n - 1 ## m is the number of leaves in subtrees m-dissimilarity map
d <- choose(n, m)
star.tree <- rep(2, d)
N1 <- length(T1) ## sample size of the data set
D1 <- matrix(rep(0, N1*choose(n, 2)), N1, choose(n, 2))
N2 <- length(T2) ## sample size of the data set
D2 <- matrix(rep(0, N2*choose(n, 2)), N2, choose(n, 2))

imputeSize1 <- round(N1*0.1)
imputeSize2 <- round(N2*0.1)

TT1 <- list()
TT2 <- list()

for(i in 1:N1){
  u <- force.ultrametric(T1[[i]], "nnls")
  u$edge.length <- u$edge.length/max(u$edge.length)
  TT1[[i]] <- u
  
  D0 <- cophenetic(u)[L, ]
  D0 <- 2*cophenetic(u)[L, ]/max(cophenetic(u)[L, ])
  D1[i, ] <- normaliz.tree(D0[lower.tri(t(D0))], 2)
}

for(i in 1:N2){
  u <- force.ultrametric(T2[[i]], "nnls")
  u$edge.length <- u$edge.length/max(u$edge.length)
  TT2[[i]] <- u
  
  D0 <- cophenetic(u)[L, ]
  D0 <- 2*cophenetic(u)[L, ]/max(cophenetic(u)[L, ])
  D2[i, ] <- normaliz.tree(D0[lower.tri(t(D0))], 2)
}

for(i in 1:100){
  index1 <- sample(1:N1, 20)
  index2 <- sample(1:N2, 20)
  D <- rbind(D1[index1,], D2[index2,])
  
  res[3, i] <- Tropical.DIANA(D, n)
  res2[3, i] <- DIANA(D, n)
}
### R = 2
T1 <- ape::read.nexus("data/R2genetrees_S1.dat")
T2 <- ape::read.nexus("data/R2genetrees_S2.dat")

L <- T2[[1]]$tip.label
n <- 10 ## number of leaves
m <- n - 1 ## m is the number of leaves in subtrees m-dissimilarity map
d <- choose(n, m)
star.tree <- rep(2, d)
N1 <- length(T1) ## sample size of the data set
D1 <- matrix(rep(0, N1*choose(n, 2)), N1, choose(n, 2))
N2 <- length(T2) ## sample size of the data set
D2 <- matrix(rep(0, N2*choose(n, 2)), N2, choose(n, 2))

imputeSize1 <- round(N1*0.1)
imputeSize2 <- round(N2*0.1)

TT1 <- list()
TT2 <- list()

for(i in 1:N1){
  u <- force.ultrametric(T1[[i]], "nnls")
  u$edge.length <- u$edge.length/max(u$edge.length)
  TT1[[i]] <- u
  
  D0 <- cophenetic(u)[L, ]
  D0 <- 2*cophenetic(u)[L, ]/max(cophenetic(u)[L, ])
  D1[i, ] <- normaliz.tree(D0[lower.tri(t(D0))], 2)
}

for(i in 1:N2){
  u <- force.ultrametric(T2[[i]], "nnls")
  u$edge.length <- u$edge.length/max(u$edge.length)
  TT2[[i]] <- u
  
  D0 <- cophenetic(u)[L, ]
  D0 <- 2*cophenetic(u)[L, ]/max(cophenetic(u)[L, ])
  D2[i, ] <- normaliz.tree(D0[lower.tri(t(D0))], 2)
}

for(i in 1:100){
  index1 <- sample(1:N1, 20)
  index2 <- sample(1:N2, 20)
  D <- rbind(D1[index1,], D2[index2,])
  
  res[4, i] <- Tropical.DIANA(D, n)
  res2[4, i] <- DIANA(D, n)
}

### R = 5
T1 <- ape::read.nexus("data/R5genetrees_S1.dat")
T2 <- ape::read.nexus("data/R5genetrees_S2.dat")

L <- T2[[1]]$tip.label
n <- 10 ## number of leaves
m <- n - 1 ## m is the number of leaves in subtrees m-dissimilarity map
d <- choose(n, m)
star.tree <- rep(2, d)
N1 <- length(T1) ## sample size of the data set
D1 <- matrix(rep(0, N1*choose(n, 2)), N1, choose(n, 2))
N2 <- length(T2) ## sample size of the data set
D2 <- matrix(rep(0, N2*choose(n, 2)), N2, choose(n, 2))

imputeSize1 <- round(N1*0.1)
imputeSize2 <- round(N2*0.1)

TT1 <- list()
TT2 <- list()

for(i in 1:N1){
  u <- force.ultrametric(T1[[i]], "nnls")
  u$edge.length <- u$edge.length/max(u$edge.length)
  TT1[[i]] <- u
  
  D0 <- cophenetic(u)[L, ]
  D0 <- 2*cophenetic(u)[L, ]/max(cophenetic(u)[L, ])
  D1[i, ] <- normaliz.tree(D0[lower.tri(t(D0))], 2)
}

for(i in 1:N2){
  u <- force.ultrametric(T2[[i]], "nnls")
  u$edge.length <- u$edge.length/max(u$edge.length)
  TT2[[i]] <- u
  
  D0 <- cophenetic(u)[L, ]
  D0 <- 2*cophenetic(u)[L, ]/max(cophenetic(u)[L, ])
  D2[i, ] <- normaliz.tree(D0[lower.tri(t(D0))], 2)
}

for(i in 1:100){
  index1 <- sample(1:N1, 20)
  index2 <- sample(1:N2, 20)
  D <- rbind(D1[index1,], D2[index2,])
  
  res[5, i] <- Tropical.DIANA(D, n)
  res2[5, i] <- DIANA(D, n)
}

### R = 10
T1 <- ape::read.nexus("data/R10genetrees_S1.dat")
T2 <- ape::read.nexus("data/R10genetrees_S2.dat")

L <- T2[[1]]$tip.label
n <- 10 ## number of leaves
m <- n - 1 ## m is the number of leaves in subtrees m-dissimilarity map
d <- choose(n, m)
star.tree <- rep(2, d)
N1 <- length(T1) ## sample size of the data set
D1 <- matrix(rep(0, N1*choose(n, 2)), N1, choose(n, 2))
N2 <- length(T2) ## sample size of the data set
D2 <- matrix(rep(0, N2*choose(n, 2)), N2, choose(n, 2))

imputeSize1 <- round(N1*0.1)
imputeSize2 <- round(N2*0.1)

TT1 <- list()
TT2 <- list()

for(i in 1:N1){
  u <- force.ultrametric(T1[[i]], "nnls")
  u$edge.length <- u$edge.length/max(u$edge.length)
  TT1[[i]] <- u
  
  D0 <- cophenetic(u)[L, ]
  D0 <- 2*cophenetic(u)[L, ]/max(cophenetic(u)[L, ])
  D1[i, ] <- normaliz.tree(D0[lower.tri(t(D0))], 2)
}

for(i in 1:N2){
  u <- force.ultrametric(T2[[i]], "nnls")
  u$edge.length <- u$edge.length/max(u$edge.length)
  TT2[[i]] <- u
  
  D0 <- cophenetic(u)[L, ]
  D0 <- 2*cophenetic(u)[L, ]/max(cophenetic(u)[L, ])
  D2[i, ] <- normaliz.tree(D0[lower.tri(t(D0))], 2)
}

for(i in 1:100){
  index1 <- sample(1:N1, 20)
  index2 <- sample(1:N2, 20)
  D <- rbind(D1[index1,], D2[index2,])
  
  res[6, i] <- Tropical.DIANA(D, n)
  res2[6, i] <- DIANA(D, n)
}


x <- c(0.25, 0.5, 1, 2, 5, 10)

plot(x,apply(res, 1, mean),type="l",col="red")
lines(x,apply(res2, 1, mean),col="blue")
legend("bottomright", legend=c("Tropical Metric", "Euclidean"),
       col=c("red", "blue"),  lty = 1:1)


plot(x,apply(res, 1, mean),type="l",col="red")
lines(x,apply(res2, 1, mean),col="blue")
lines(x, apply(res, 1, mean)-1.96*apply(res, 1, sd), pch = 18, col = "red", type = "b", lty = 2)
lines(x, apply(res, 1, mean)+1.96*apply(res, 1, sd), pch = 18, col = "red", type = "b", lty = 2)
lines(x, apply(res2, 1, mean)-1.96*apply(res2, 1, sd), pch = 18, col = "blue", type = "b", lty = 2)
lines(x, apply(res2, 1, mean)+1.96*apply(res2, 1, sd), pch = 18, col = "blue", type = "b", lty = 2)

plot(x,apply(res, 1, mean),type="l",col="red")
lines(x,apply(res2, 1, mean),col="blue")
lines(x, apply(res, 1, mean)-1*apply(res, 1, sd), pch = 18, col = "red", type = "b", lty = 2)
lines(x, apply(res, 1, mean)+1*apply(res, 1, sd), pch = 18, col = "red", type = "b", lty = 2)
lines(x, apply(res2, 1, mean)-1*apply(res2, 1, sd), pch = 18, col = "blue", type = "b", lty = 2)
lines(x, apply(res2, 1, mean)+1*apply(res2, 1, sd), pch = 18, col = "blue", type = "b", lty = 2)
