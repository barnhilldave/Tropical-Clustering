###########################################
## Author :  Ruriko Yoshida and David Barnhill
## Date   :  May 26th 2023
## Update :  June 10th 2023
## Program:  Tropical clustering via tropical polytopes
## Input  :  
## Output :  
## Execute: type in R as 
##
#############################################


library(ape)
library(phangorn)
library(phytools)
# New library for calculating ROC curves
library(ROCR)
library(TreeTools)
library(gtools)
library(mvtnorm)
library(caper)
library(metaSEM)
library(lpSolve)

set.seed(213)

project_pi<-function(D_s,D){
    d <- dim(D_s)
    lambda <- rep(0, d[1])
    for(i in 1:d[1])
        lambda[i] <- min(D - D_s[i,])

    x <- rep(0, d[2])
    for(i in 1:d[2])
        x[i] <- max((lambda+D_s)[, i])
    return(x)
}

normaliz.tree <- function(D, h = 1){
    d <- length(D)
    a <- max(D) - h
    x <- D - rep(a, d)
    for(i in 1:d)
        if(x[i] < 0)
            x[i] <- 0
    return(x)
}

trop.dist <- function(D1, D2){
    return(max(D2 - D1) - min(D2 - D1))
}

normaliz.vector <- function(D){
    return(D - rep(D[1], length(D)))
}

pairws.distances <- function(D1){
## D1 is a matrix.  the rows are observations.
    d1 <- dim(D1)
    D <- matrix(rep(0, d1[1]*d1[1]), d1[1], d1[1])
    for(i in 1:d1[1])
        for(j in 1:d1[1])
            D[i, j] <- trop.dist(D1[i, ], D1[j, ])
    return(D)
}

trop_wi_dist<-function(D1,method='avg'){
  if(method=='avg'){
  md<-mean(lower_tri(pairws.distances(D1)))}
  if(method=='max'){
    md<-max(lower_tri(pairws.distances(D1)))
  }
  return(md)
}

trop_bet_dist<-function(D1,D2){
  dstar<-0
  for(i in 1:nrow(D2)){
    pD2<-project_pi(D1,D2[i,])
    td<-trop.dist(D2[i,],pD2)
    dstar<-dstar+td
  }
  return(dstar/nrow(D2))
}

over_bet_KM<-function(A,C){
  bets<-c()
  for (i in 1:C){
    bet<-c()
    for(j in 1:C){
      if(i!=j){
        bet1<-mean(c(trop_bet_dist(A[which(A[,ncol(A)]==i),1:(ncol(A)-1)],A[which(A[,ncol(A)]==j),1:(ncol(A)-1)]),trop_bet_dist(A[which(A[,4]==j),1:(ncol(A)-1)],A[which(A[,4]==i),1:(ncol(A)-1)])))
        bet<-append(bet,bet1)
      }
    }
    mbet<-mean(bet)
    bets<-append(bets,mbet)
  }
  return(bets)  
}

over_bet_HC<-function(A,V){
  bets<-c()
  for (i in 1:length(V)){
    bet<-c()
    for(j in 1:length(V)){
      if(i!=j){
        bet1<-mean(c(trop_bet_dist(A[V[[i]],],A[V[[j]],]),trop_bet_dist(A[V[[j]],],A[V[[i]],])))
        bet<-append(bet,bet1)
      }
    }
    mbet<-mean(bet)
    bets<-append(bets,mbet)
  }
return(bets)  
}

cluster.ratio_KM<-function(A,C,method='avg'){
  wits<-c()
  for (i in 1:C){
    wit<-trop_wi_dist(A[which(A[,ncol(A)]==i),1:(ncol(A)-1)],method='avg')
    wits<-append(wits,wit)
  }
  bets<-over_bet_KM(A,C)
  rats<-c()
  for(i in 1:C){
    bet_p<-bets[-i]
    rat<-wits[i]/mean(bet_p)
    rats<-append(rats,rat)
  }
  return(rats)
}

cluster.ratio_HC<-function(A,V,method='avg'){
  wits<-c()
  for (i in 1:length(V)){
    wit<-trop_wi_dist(A[V[[i]],],method='avg')
    wits<-append(wits,wit)
  }
  bets<-over_bet(A,V)
  rats<-c()
  for(i in 1:length(V)){
    bet_p<-bets[-i]
    rat<-wits[i]/mean(bet_p)
    rats<-append(rats,rat)
  }
  return(rats)
}


# A is a matrix with rows as observations
Trop_FW<-function(A){
  perms<-permutations(ncol(A),2)
  mins<-apply(A,2,min)
  if(any(mins<0)){
    mins_ind<-which(mins<0)
    A[,mins_ind]<--rep(mins[mins_ind],each=nrow(A))+A[,mins_ind]
  }
  
  yys<-matrix(0,nrow(perms),ncol(A),TRUE)
  for(i in 1:nrow(yys)){
    yys[i,perms[i,1]]=1
    yys[i,perms[i,2]]=-1
  }
  L<-matrix(rep(t(yys),nrow(A)),ncol=ncol(yys),byrow=TRUE)
  
  R<-matrix(0,0,nrow(A),TRUE)
  for(j in 1:nrow(A)){
    rrs<-matrix(0,nrow(perms),nrow(A),TRUE)
    rrs[,j]=-1
    R<-rbind(R,rrs)
  }
  f.con<-cbind(L,R)
  zer<-c(1,rep(0,(ncol(f.con)-1)))
  f.con<-rbind(f.con,zer)
  
  f.rhs<-c()
  for(j in 1:nrow(A)){
    for(i in 1:nrow(perms)){
      diff<--A[j,perms[i,1]]+A[j,perms[i,2]]
      f.rhs<-append(f.rhs,-diff)
    }
  }
  f.rhs<-append(f.rhs,0)
  f.dir <- c(rep("<=",nrow(f.con)))
  f.obj<-c(rep(0,(ncol(f.con)-ncol(R))),rep(1,ncol(R)))
  res<-lp ("min", f.obj, f.con, f.dir, f.rhs)
  FW<-res$solution[1:ncol(A)]
  if(any(mins<0)){
    FW[mins_ind]<-FW[mins_ind]+mins[mins_ind]
  }
  return(FW)
}

########### Tropical K-Means Clustering #################
TKmeans<-function(A,C,M){
  or_cents<-matrix(0,C,dim(A)[2],TRUE)
  clus<-sample(1:C,nrow(A),TRUE)
  try<-cbind(A,clus)
  CL=FALSE
  ct<-0
  while(CL==FALSE|ct<=M){
    for (i in 1:C){
      c<-Trop_FW(try[which(try[,ncol(try)]==i),1:(ncol(try)-1)])
      or_cents[i,]<-c
    }
    try1<-try
    for (j in 1:nrow(try)){
      try1[j,ncol(try1)]<-which.min(apply(or_cents,1,function(x) trop.dist(x,A[j,])))
    }
    
    if(all(try[,ncol(try)]==try1[,ncol(try)])){
      try=try1
      CL=TRUE
      ct<-ct+1
    }
    else{
      try=try1
      ct<-ct+1
    }
    print(ct)
  }
  RES<-list(try,or_cents,ct)
  return(RES)
}


# Tropical K-means for Ultrametrics using FW

############# KMeans Clustering Ultrametrics ###################

## Project a point onto the space of ultrametrics
pro.ultrametric <- function(x0, n,h = 1){
  d <- length(x0)
  x <- normaliz.tree(x0, h)
  a <- h
  D0 <- vec2symMat(x,diag=FALSE)
  for(k in 1:n)
    D0[k, k] <- 0
  tree <- upgma(D0)
  tree$edge.length <- tree$edge.length/max(tree$edge.length)
  DD <- cophenetic(tree)
  v <- DD[lower.tri(t(D))] #rep(0, d)
  return(normaliz.tree(v, h))
}

TKmeans_Ultra<-function(A,C,n,st){
  h<-max(A)
  or_cents<-matrix(0,C,dim(A)[2],TRUE)
  clus<-sample(1:C,nrow(A),TRUE)
  try<-cbind(A,clus)
  CL=FALSE
  ct<-0
  tries<-matrix(0,0,nrow(A),TRUE)
  while(CL==FALSE & ct<st){
    for (i in 1:C){
      if(nrow(try[which(try[,ncol(try)]==i),1:(ncol(try)-1)])>1){
        print('if')
        c<-Trop_FW(try[which(try[,ncol(try)]==i),1:(ncol(try)-1)])
        or_cents[i,]<-pro.ultrametric(c,n,h)
      }
      else if(length(try[which(try[,ncol(try)]==i),1:(ncol(try)-1)])>0){
        print('elseif')
        c<-try[which(try[,ncol(try)]==i),1:(ncol(try)-1)]
        or_cents[i,]<-pro.ultrametric(c,n,h)
      }
      else{
        print('break')
        break
      }
    }
    try1<-try
    for (j in 1:nrow(try)){
      try1[j,ncol(try1)]<-which.min(apply(or_cents,1,function(x) trop.dist(x,A[j,])))
    }
    
    if(all(try[,ncol(try)]==try1[,ncol(try)])){
      try=try1
      tries<-rbind(tries,try[,ncol(try)])
      CL=TRUE
      ct<-ct+1
    }
    else{
      try=try1
      tries<-rbind(tries,try[,ncol(try)])
      ct<-ct+1
    }
    print(ct)
  }
  RES<-list(try,or_cents,ct,tries)
  return(RES)
}


################# Begin Tropical HC #############

####### Agglomerative
tropical.complete.linkage <- function(D1, D2){
## D1 is a set of vertices for a tropical polytope (raws are obs)
## D2 is a set of observations of one cluster
    x.star <- 0 ## max distance
    d1 <- dim(D1)
    d2 <- dim(D2)
    for(i in 1:d2[1]){
        x <- trop.dist(project_pi(D1, D2[i, ]), D2[i, ])
        if(x > x.star)
            x.star <- x
    }
    for(i in 1:d1[1]){
        x <- trop.dist(project_pi(D2, D1[i, ]), D1[i, ])
        if(x > x.star)
            x.star <- x
    }
    return(x.star)
}

tropical.minimum.linkage <- function(D1, D2){
## D1 is a set of vertices for a tropical polytope (raws are obs)
## D2 is a set of observations of one cluster
    x.star <-trop.dist(project_pi(D1, D2[1, ]), D2[1, ]) ## max distance
    d1 <- dim(D1)
    d2 <- dim(D2)
    for(i in 1:d2[1]){
        x <- trop.dist(project_pi(D1, D2[i, ]), D2[i, ])
        if(x < x.star)
            x.star <- x
    }
    for(i in 1:d1[1]){
        x <- trop.dist(project_pi(D2, D1[i, ]), D1[i, ])
        if(x < x.star)
            x.star <- x
    }
    return(x.star)
}

tropical.average.linkage <- function(D1, D2){
## D1 is a set of vertices for a tropical polytope (raws are obs)
## D2 is a set of observations of one cluster
    x.star <- 0 ## max distance
    d1 <- dim(D1)
    d2 <- dim(D2)
    x <- rep(0, (d1[1] + d2[1]))
    for(i in 1:d2[1]){
        x[i] <- trop.dist(project_pi(D1, D2[i, ]), D2[i, ])
    }
    for(i in 1:d1[1]){
        x[i + d2[1]] <- trop.dist(project_pi(D2, D1[i, ]), D1[i, ])
    }
    return(mean(x))
}

finding.pair <- function(D1, method="average"){
## D1 is a list of matrices.  
    ### method is one of "average", "min", "max"
    d1 <- length(D1)
    best.pair.index <- c(1, 2)
    d.star <- tropical.average.linkage(D1[[1]], D1[[2]]) + 100000
    for(i in 1:d1)
        for(j in 1:d1){
            d <- d.star
            if(method == "average")
                if(i != j)
                    d <- tropical.average.linkage(D1[[i]], D1[[j]])
            if(method == "min")
                if(i != j)
                    d <- tropical.minimum.linkage(D1[[i]], D1[[j]])            
            if(method == "max")
                if(i != j)
                    d <- tropical.complete.linkage(D1[[i]], D1[[j]])
            if(d < d.star){
                d.star <- d
                best.pair.index[1] <- i
                best.pair.index[2] <- j
            }
                
        }
    return(list(best.pair.index, d.star))
}

make.list.matrices <- function(D, index){
    D1 <- list()
    for(j in 1:length(index)){
        D2 <- as.matrix(t(D[index[[j]][1],]))
        if(length(index[[j]]) > 1)
            for(k in 2:length(index[[j]]))
                D2 <- rbind(D2, D[index[[j]][k],])
        D1[[j]] <- D2
    }
    return(D1)
}

tropical.HC <- function(D, method="average"){
 ### method is one of "average", "min", "max"
    d <- dim(D)
    distance <- rep(0, d[1])
    index.list <- list()
    index <- list()
    for(i in 1:d[1])
        index[[i]] <- c(i)
    index.list[[1]] <- index
    for(i in 2:d[1]){
        D1 <- make.list.matrices(D, index)
        best.pair <- finding.pair(D1, method)
        #print(best.pair[[1]])
        index[[best.pair[[1]][1]]] <- cbind(index[[best.pair[[1]][1]]], index[[best.pair[[1]][2]]])
        index <- index[-best.pair[[1]][2]]
        #print(index[best.pair[[1]][2]])
        index.list[[i]] <- index
        distance[i] <- best.pair[[2]]
    }
    return(list(distance, index.list))
}

########### Divisive

Tropical.DIANA <- function(D, n, k=2){
  ## n is the number of leaves
  ## D is a dataset (matrix) whose row are ultrametrics
  ## k is the number of clusters
  d <- dim(D)
  x <- 0
  M <- matrix(rep(0, d[1]*d[1]), d[1], d[1])
  for(i in 1:d[1])
    for(j in 1:d[1])
      M[i, j] <- trop.dist(D[i, ], D[j, ])
  cl <- diana(M, diss=TRUE)
  grp <- cutree(cl, k = k)
  x <- (sum(grp[1:20]==grp[1])+sum(grp[21:40]==grp[21]))/40
  if(x < 0.5)
    x <- 1 - x
  return(x)
}

DIANA <- function(D, n, k=2){
  ## n is the number of leaves
  ## D is a dataset (matrix) whose row are ultrametrics
  ## k is the number of clusters
  cl <- diana(D)
  grp <- cutree(cl, k = k)
  x <- (sum(grp[1:20]==grp[1])+sum(grp[21:40]==grp[21]))/40
  if(x < 0.5)
    x <- 1 - x
  return(x)
}


