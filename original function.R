############################the cumulative distribution of H(y) 

#Arguments:
#hy: the most updated H
#PI: the most updated mixing probabilities
#alphabeta: the most updated basis coefficient of mean functions
#sigma: the most updated variances of errors
#Lamda: the most updated eigenvalues
#Alpha: the most updated basis coefficient of eigenfunctions
#Group_NUM: the most updated number of clusters

#return: the cumulative distribution of H(y) 

original<-function(hy,PI,alphabeta,sigma,Lamda,Alpha,Group_NUM){
  out=0
  for(i in 1:n){
    XXX<-splinex1(TT[[i]])
    for(j in 1:length(TT[[i]])){
      for(k in 1:Group_NUM){
        Phi<-XXX%*%t(Alpha[,,k])
        sigmak<-sqrt(Phi[j,]%*%Lamda[,,k]%*%t(t(Phi[j,]))+sigma[k])
        out=out+PI[k]*pnorm((hy-XXX[j,]%*%alphabeta[,k])/sigmak)
      }
    }
  }
  return(out)
}