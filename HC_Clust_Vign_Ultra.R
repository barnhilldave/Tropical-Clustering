
source('tropical_clustering.R')
############## Hierarchical Clustering Ultrametrics #################

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
output1 <- tropical.HC(D)
#mean(lower_tri(pairws.distances(D[output1[[2]][[79]][[1]],])))
#mean(lower_tri(pairws.distances(D[output1[[2]][[79]][[2]],])))

#DD11<-D[output1[[2]][[79]][[1]],]
#DD22<-D[output1[[2]][[79]][[2]],]
#d1star<-100000



#trop_bet_dist(A[V[[1]],],A[V[[2]],])


#d2star<-100000
#for(i in 1:nrow(DD22)){
#  pD22<-project_pi(DD11,DD22[i,])
#  td<-trop.dist(DD22[i,],pD22)
#  if(td<d2star){
#    d2star<-td
#  }
#}

#output <- tropical.HC(D, method="max")
#output2 <- tropical.HC(D, method="min")

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
output1 <- tropical.HC(D)
output1 <- tropical.HC(D,method="max")
output1 <- tropical.HC(D,method="min")

output1[[2]][[80]]

### R = .5
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

### R = 1
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
### R = 2
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

### R = 5
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
######### Wine data
wines<-wine[,-1]
wines1<-cbind(rep(0,nrow(wines)),wines)
W<-TKmeans(wines1,3) #.713
WKM<-kmeans(wines1,3) # .7022



iris2<-iris[,-ncol(iris)]
iris3<-cbind(rep(0,nrow(iris2)),iris2)
iris3[,2:ncol(iris3)]<-scale(iris3[,2:ncol(iris3)])
tropical.HC(iris3)
IKM<-kmeans(iris3,3)
ITKM<-TKmeans(iris3,3,1000)

as.matrix(ITKM[[1]][,ncol(ITKM[[1]])])



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
pwD<-as.dist(pairws.distances(D))
e_D<-dist(D)
############ Average
output110 <- tropical.HC(D)
output110[[2]]

hd110<-hclust(pwD,method='average')
dd110<-diana(pwD)
plot(dd110)
plot(hd110)

h110_ed<-hclust(e_D,method='average')
dd110_ed<-diana(e_D)
plot(dd110_ed)
plot(h110_ed)
################ Complete
output210 <- tropical.HC(D,method="max")
output210[[2]]

hd110c<-hclust(pwD,method='complete')
plot(hd110c)

h110_edc<-hclust(e_D,method='complete')
plot(h110_edc)
############### Single
output310 <- tropical.HC(D,method="min")
output310[[2]]

hd110s<-hclust(pwD,method='single')
plot(hd110s)

h110_eds<-hclust(e_D,method='single')
plot(h110_eds)


########### Average
output110 <- tropical.HC(D)
output110[[2]]
########### Complete
output210 <- tropical.HC(D,method="max")
output210[[2]]
########### Single
output310 <- tropical.HC(D,method="min")
output310[[2]]

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
pwD<-as.dist(pairws.distances(D))
e_D<-dist(D)
############ Average
output15 <- tropical.HC(D)
output15[[2]]

hd15<-hclust(pwD,method='average')
dd15<-diana(pwD)
plot(hd15)
plot(dd15)

h15_ed<-hclust(e_D,method='average')
dd15_ed<-diana(e_D)
plot(dd15_ed)
plot(h15_ed)
################ Complete
output25 <- tropical.HC(D,method="max")
output25[[2]]

hd15c<-hclust(pwD,method='complete')
plot(hd15c)

h15_edc<-hclust(e_D,method='complete')
plot(h15_edc)
############### Single
output35 <- tropical.HC(D,method="min")
output35[[2]]

hd15s<-hclust(pwD,method='single')
plot(hd15s)

h15_eds<-hclust(e_D,method='single')
plot(h15_eds)



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
pwD<-as.dist(pairws.distances(D))
e_D1.25<-dist(D)
################# Average Linkage
output1.25 <- tropical.HC(D)
output1.25[[2]]

hd1.25<-hclust(pwD,method='average')
dd1.25<-diana(pwD)
plot(hd1.25)
plot(dd1.25)

h1.25_ed<-hclust(e_D1.25,method='average')
dd1.25_ed<-diana(e_D1.25)
plot(dd1.25_ed)
plot(h1.25_ed)
################# Complete Linkage
output2.25 <- tropical.HC(D,method="max")
output2.25[[2]]
hd2.25c<-hclust(pwD,method='complete')
plot(hd2.25c)
e_D2.25<-dist(D)
h2.25_edc<-hclust(e_D2.25,method='complete')
plot(h2.25_edc)
############# Single Linkage
output3.25 <- tropical.HC(D,method="min")
output3.25[[2]]
hd3.25c<-hclust(pwD,method='single')
plot(hd3.25c)
e_D3.25<-dist(D)
h3.25_edc<-hclust(e_D3.25,method='single')
plot(h3.25_edc)
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
pwD<-as.dist(pairws.distances(D))
e_D<-dist(D)
############ Average
output1.5 <- tropical.HC(D)
output1.5[[2]]

hd1.5<-hclust(pwD,method='average')
dd1.5<-diana(pwD)
plot(dd1.5)
plot(hd1.5)

h1.5_ed<-hclust(e_D,method='average')
dd1.5_ed<-diana(e_D)
plot(dd1.5_ed)
plot(h1.5_ed)
################ Complete
output2.5 <- tropical.HC(D,method="max")
output2.5[[2]]

hd1.5c<-hclust(pwD,method='complete')
plot(hd1.5c)

h1.5_edc<-hclust(e_D,method='complete')
plot(h1.5_edc)
############### Single
output3.5 <- tropical.HC(D,method="min")
output3.5[[2]]

hd1.5s<-hclust(pwD,method='single')
plot(hd1.5s)

h1.5_eds<-hclust(e_D,method='single')
plot(h1.5_eds)
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
pwD<-as.dist(pairws.distances(D))
e_D<-dist(D)
############ Average
output11 <- tropical.HC(D)
output11[[2]]

hd11<-hclust(pwD,method='average')
dd11<-diana(pwD)
plot(dd11)
plot(hd11)

h11_ed<-hclust(e_D,method='average')
dd11_ed<-diana(e_D)
plot(dd11_ed)
plot(h11_ed)
################ Complete
output21 <- tropical.HC(D,method="max")
output21[[2]]

hd11c<-hclust(pwD,method='complete')
plot(hd11c)

h11_edc<-hclust(e_D,method='complete')
plot(h11_edc)
############### Single
output31 <- tropical.HC(D,method="min")
output31[[2]]

hd11s<-hclust(pwD,method='single')
plot(hd11s)

h11_eds<-hclust(e_D,method='single')
plot(h11_eds)
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
pwD<-as.dist(pairws.distances(D))
e_D<-dist(D)
############ Average
output12 <- tropical.HC(D)
output12[[2]]

hd12<-hclust(pwD,method='average')
dd12<-diana(pwD)
plot(dd12)
plot(hd12)

h12_ed<-hclust(e_D,method='average')
dd12_ed<-diana(e_D)
plot(dd12_ed)
plot(h12_ed)
################ Complete
output22 <- tropical.HC(D,method="max")
output22[[2]]

hd12c<-hclust(pwD,method='complete')
plot(hd12c)

h12_edc<-hclust(e_D,method='complete')
plot(h12_edc)
############### Single
output32 <- tropical.HC(D,method="min")
output32[[2]]

hd12s<-hclust(pwD,method='single')
plot(hd12s)

h12_eds<-hclust(e_D,method='single')
plot(h12_eds)



