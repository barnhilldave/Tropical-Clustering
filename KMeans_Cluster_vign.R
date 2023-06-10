
source('tropical_clustering.R')
########## Experiment 1
set.seed(124)
P<-matrix(c(0,-100,-100,0,100,0,0,0,100),3,3,TRUE)
x0<-c(0,0,0)
y0<-c(0,0,0)
z0<-c(0,0,0)
mu1<-c(0,-10,-20)
sig1<-5
mu2<-c(0,20,30)
sig2<-8
mu3<-c(0,30,10)
sig3<-3
N<-50
har_norms1<-matrix(0,N,3)
har_norms2<-matrix(0,N,3)
har_norms3<-matrix(0,N,3)

for(i in 1:N){
  print(i)
  x<-TropicalPolytope.extrapolation.HAR_NORM(P, x0, I = 50,k=2,mu1,sig1)
  y<-TropicalPolytope.extrapolation.HAR_NORM(P, y0, I = 50,k=2,mu2,sig2)
  z<-TropicalPolytope.extrapolation.HAR_NORM(P, z0, I = 50,k=2,mu3,sig3)
  x0<-x
  y0<-y
  z0<-z
  har_norms1[i,]<-x
  har_norms2[i,]<-y
  har_norms3[i,]<-z
}

hars<-rbind(har_norms1,har_norms2,har_norms3)

plot(hars[,2],hars[,3],asp=1,xlab='x2',ylab='x3')

res<-TKmeans(hars,3,10)
try<-res[[1]]
cen<-res[[2]]
plot(try[,2],try[,3],col=try[,4],asp=1)
plot(try[,2],try[,3],col=try[,4],asp=1,xlab='x2',ylab='x3')
points(cen[,2],cen[,3],col=c('purple','hotpink','orange'),pch=19)

################# Experiment 2 ####################
P<-matrix(c(0,-100,-100,0,100,0,0,0,100),3,3,TRUE)
x0<-c(0,0,0)
y0<-c(0,0,0)
z0<-c(0,0,0)
mu1<-c(0,-5,-5)
sig1<-4
mu2<-c(0,5,5)
sig2<-4
mu3<-c(0,10,0)
sig3<-4
N<-50
har_norms11<-matrix(0,N,3)
har_norms22<-matrix(0,N,3)
har_norms33<-matrix(0,N,3)

for(i in 1:N){
  print(i)
  x<-TropicalPolytope.extrapolation.HAR_NORM(P, x0, I = 50,k=2,mu1,sig1)
  y<-TropicalPolytope.extrapolation.HAR_NORM(P, y0, I = 50,k=2,mu2,sig2)
  z<-TropicalPolytope.extrapolation.HAR_NORM(P, z0, I = 50,k=2,mu3,sig3)
  x0<-x
  y0<-y
  z0<-z
  har_norms11[i,]<-x
  har_norms22[i,]<-y
  har_norms33[i,]<-z
}

hars2<-rbind(har_norms11,har_norms22,har_norms33)
plot(hars2[,2],hars2[,3],asp=1,xlab='x2',ylab='x3')
points(har_norms11[,2],har_norms11[,3],col='blue')
points(har_norms22[,2],har_norms22[,3],col="red")
points(har_norms33[,2],har_norms33[,3],col="green")


res2<-TKmeans(hars2,3)
try<-res2[[1]]
cen<-res2[[2]]
plot(try[,2],try[,3],col=try[,4],asp=1)
plot(try[,2],try[,3],col=try[,4],asp=1,xlab='x2',ylab='x3')
points(cen[,2],cen[,3],col=c('purple','hotpink','orange'),pch=19)

t_clu<-c(rep(1,50),rep(2,50),rep(3,50))
hars_ag<-cbind(hars2,t_clu)
which(hars_ag[,4]!=try[,4])
plot(hars_ag[,2],hars_ag[,3],asp=1,col=hars_ag[,4],xlab='x2',ylab='x3')
###################### IRIS data
iris1<-iris[,-ncol(iris)]

iris12<-scale(iris1)
for(i in 1:ncol(iris1)){
  iris1[,i]<-scale(iris1[,i])
}

iris2<-cbind(rep(0,nrow(iris12)),iris12)
iris2<-cbind(rep(0,nrow(iris1)),iris1)
res_iris<-TKmeans(iris2,3,15)

blend_i<-cbind(res_iris[[1]][,ncol(res_iris[[1]])],iris[,ncol(iris)])

blend_i[which(blend_i[,1]==1),]
blend_i[which(blend_i[,1]==2),]
blend_i[which(blend_i[,1]==3),]

res_iris1<-kmeans(iris12,3)
blend_i_euc<-cbind(res_iris1[[1]],iris[,ncol(iris)])
blend_i_euc[which(blend_i_euc[,1]==1),]
blend_i_euc[which(blend_i_euc[,1]==2),]
blend_i_euc[which(blend_i_euc[,1]==3),]
############# Wine data
library(HDclassif)
data(wine)

wine1<-wine[,-1]

for(i in 1:ncol(wine1)){
  wine1[,i]<-scale(wine1[,i])
}

wine2<-cbind(rep(0,nrow(wine1)),wine1)

res_wine_euc<-kmeans(wine1,3)

blend_euc<-cbind(res_wine_euc[[1]],wine[,1])

blend_euc[blend_euc[,1]==1,]

res_wine<-TKmeans(wine2,3,20)

blend<-cbind(res_wine[[1]][,ncol(res_wine[[1]])],wine[,1])

############# KMeans Clustering Ultrametrics ###################

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

D <- rbind(D1[1:40,], D2[1:40,])
TEST10<-TKmeans_Ultra(D,2,10,10)
TEST10[[1]][,ncol(TEST10[[1]])]

### R = .25
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

D <- rbind(D1[1:40,], D2[1:40,])
TEST.25<-TKmeans_Ultra(D,2,10)
TEST.25[[1]][,ncol(TEST.25[[1]])]

### R = .5
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

D <- rbind(D1[1:40,], D2[1:40,])
TEST.5<-TKmeans_Ultra(D,2,10,50)
TEST.5[[1]][,ncol(TEST.5[[1]])]

TEST.5K<-TKmeans(D,2)
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

D <- rbind(D1[1:40,], D2[1:40,])
TEST1<-TKmeans_Ultra(D,2,10,10)
TEST1[[1]][,ncol(TEST1[[1]])]
kmeans(D,2)
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

D <- rbind(D1[1:40,], D2[1:40,])
TEST2<-TKmeans_Ultra(D,2,10,10)
TEST2[[1]][,ncol(TEST2[[1]])]
TEST2[[4]]
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

D <- rbind(D1[1:40,], D2[1:40,])
TEST5<-TKmeans_Ultra(D,2,10,10)
TEST5[[1]][,ncol(TEST5[[1]])]
