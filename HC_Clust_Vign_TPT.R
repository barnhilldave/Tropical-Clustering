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

XX<-tropical.HC(hars,method='average')
tropical.HC(as.matrix(iris3[,-1]),method='average')

plot(hars[,2],hars[,3],asp=1,xlab='x2',ylab='x3')
points(har_norms1[,2],har_norms1[,3],col='blue')
points(har_norms2[,2],har_norms2[,3],col="red")
points(har_norms3[,2],har_norms3[,3],col="green")

plot(XX[[1]],type='l')
abline(v=139)

length(XX[[2]][[148]])
plot(hars[,2],hars[,3],xlab='x2',ylab='x3',asp=1)
for(i in 1:length(XX[[2]][[148]])){
  points(hars[XX[[2]][[148]][[i]],2],hars[XX[[2]][[148]][[i]],3],col=(i+4))
}

YY<-tropical.HC(hars,method='max')



plot(hars[,2],hars[,3],xlab='x2',ylab='x3',asp=1)
for(i in 1:length(YY[[2]][[148]])){
  points(hars[YY[[2]][[148]][[i]],2],hars[YY[[2]][[148]][[i]],3],col=(i+4))
}
length(YY[[2]][[148]])
over_bet(hars,YY[[2]][[148]])

plot(YY[[1]],type='l')

points(hars[YY[[2]][[148]][[1]],2],hars[YY[[2]][[148]][[1]],3],col='blue')

points(hars[YY[[2]][[148]][[2]],2],hars[YY[[2]][[148]][[2]],3],col='red')
points(hars[YY[[2]][[148]][[3]],2],hars[YY[[2]][[148]][[3]],3],col='green')

ZZ<-tropical.HC(hars,method='min')
plot(hars[,2],hars[,3],xlab='x2',ylab='x3',asp=1)
for(i in 1:length(ZZ[[2]][[148]])){
  points(hars[ZZ[[2]][[148]][[i]],2],hars[ZZ[[2]][[148]][[i]],3],col=(i+4))
}
plot(ZZ[[1]],type='l')
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

XX<-tropical.HC(hars,method='average')
for(i in 1:length(XX[[2]][[148]])){
  points(hars[XX[[2]][[148]][[i]],2],hars[XX[[2]][[148]][[i]],3],col=(i+4))
}
YY<-tropical.HC(hars,method='max')
for(i in 1:length(YY[[2]][[148]])){
  points(hars[YY[[2]][[148]][[i]],2],hars[YY[[2]][[148]][[i]],3],col=(i+4))
}
ZZ<-tropical.HC(hars,method='min')
plot(hars[,2],hars[,3],xlab='x2',ylab='x3',asp=1)
for(i in 1:length(ZZ[[2]][[148]])){
  points(hars[ZZ[[2]][[148]][[i]],2],hars[ZZ[[2]][[148]][[i]],3],col=(i+4))
}