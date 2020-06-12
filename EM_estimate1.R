############################the function for estimating Omega_n given H. It uses function defined in the file identify.R.
#Arguments:
#TT: the observed time points for n curves
#YY: the observed functional data
#Group_NUM: the most updated number of clusters
#PI_int: the most updated mixing probabilities
#Beta: the most updated basis coefficient of mean functions
#gamma_int: the most updated variances of errors
#Lamda: the most updated eigenvalues
#Alpha: the most updated basis coefficient of eigenfunctions
#HYY: the most updated H
#lambda: the penalty parameter


#return:
#Group_number:  the estimated number of clusters
#PI: the estimated mixing probabilities
#beta: the estimated basis coefficient of mean functions
#sigma: the estimated variances of errors
#lamda: the estimated eigenvalues
#eigencoefficient: the estimated basis coefficient of eigenfunctions
#Sstep: the number of Iterations
#del: the estimated probabilities for individuals belong to each group

EM_estimate<-function(TT,YY,Group_NUM,PI_int,Beta,gamma_int,Lamda,Alpha,HYY,lambda)
{
  delta<-matrix(0,nrow=length(YY),ncol=Group_NUM) 
  Beta_XY<-matrix(0,nrow=q_n,ncol=Group_NUM) 
  Beta_matrix<-array(0,dim=c(Group_NUM,q_n,q_n)) 
  Sigma_XY<-matrix(0,nrow=1,ncol=Group_NUM)
  Sigma_matrix<-matrix(0,nrow=1,ncol=Group_NUM)
  lamda_XY<-matrix(0,nrow=K_g,ncol=Group_NUM)
  lamda_matrix<-matrix(0,nrow=K_g,ncol=Group_NUM)
  Alpha_XY<-array(0,dim=c(q_n,K_g,Group_NUM))
  Alpha_matrix<-array(0,dim=c(q_n,q_n,K_g,Group_NUM))
  f_den<-matrix(0,nrow=length(YY),ncol=Group_NUM)
  
  SStep<-0
  Epsilon<-1
  G_1<-sum(PI_int!=0)
  Beta_1<-Beta
  gamma_1<-gamma_int
  Lamda_1<-Lamda
 
  PI_1=PI_int
    for (k in 1:Group_NUM)
    {
      Alpha[,,k]<-t(qr.Q(qr(t(Alpha[,,k]))))
    }
    Alpha_1<-Alpha
 ############Start a loop to find the estimates
 while((SStep<50)&(Epsilon>10^(-5))){
    SStep<-SStep+1
    print(SStep)
   likelihood_em=0
  
    for (i in 1:length(YY)){
      XXX<-splinex1(TT[[i]])
      Phi<-matrix(NA,length(TT[[i]]),K_g)
MU<-matrix(NA,length(TT[[i]]),1)
      HYYY<-HYY[identify(YY[[i]],allpoint)]
      for (k in 1:Group_NUM){
          Phi<-XXX%*%t(Alpha[,,k])
        Delta<-matrix(nrow=length(TT[[i]]),ncol=length(TT[[i]]))
        Delta<-Phi%*%Lamda[,,k]%*%t(Phi)+gamma_int[k]*diag(1,ni[i])##
        
        MU<-XXX%*%Beta[,k]
        f_den[i,k]<-(2*pi)^(-length(TT[[i]])/2)/sqrt(det(Delta))*exp(-1/2*t(HYYY-MU)%*%solve(Delta)%*%(HYYY-MU))
      }
       likelihood_em=likelihood_em+log(f_den[i,]%*%PI_int) 
      }
    
    for (i in 1:length(YY)){
      XXX<-splinex1(TT[[i]])
   Phi<-matrix(NA,length(TT[[i]]),K_g)
MU<-matrix(NA,length(TT[[i]]),1)
      HHYY<-HYY[identify(YY[[i]],allpoint)]
      for (k in 1:Group_NUM){
          Phi<-XXX%*%t(Alpha[,,k])
        Delta<-matrix(nrow=length(TT[[i]]),ncol=length(TT[[i]]))
        Delta<-Phi%*%Lamda[,,k]%*%t(Phi)+gamma_int[k]*diag(1,ni[i])##
        delta[i,k]<-PI_int[k]*f_den[i,k]/((t(as.matrix(PI_int))%*%as.matrix(f_den[i,])))
          
        
        MU<-XXX%*%Beta[,k]
        times<-matrix(HHYY-MU,ni[i],1)
        Beta_XY[,k]<-Beta_XY[,k]+delta[i,k]*t(XXX)%*%solve(Delta)%*%HHYY##beta_k
        Beta_matrix[k,,]<-Beta_matrix[k,,]+delta[i,k]*t(XXX)%*%solve(Delta)%*%XXX##beta_K
        Sigma_XY[k]<-Sigma_XY[k]+delta[i,k]*sum(diag(solve(Delta)%*%solve(Delta)%*%(Phi%*%Lamda[,,k]%*%t(Phi)-times%*%t(times))))
        Sigma_matrix[k]<-Sigma_matrix[k]+delta[i,k]*sum(diag(solve(Delta)%*%solve(Delta)))
        for (g in 1:K_g){
          
          lamda_XY[g,k]<-lamda_XY[g,k]+delta[i,k]*t(Phi[,g])%*%solve(Delta)%*%(Delta-Lamda[g,g,k]*Phi[,g]%*%t(Phi[,g])-times%*%t(times))%*%solve(Delta)%*%Phi[,g]
          lamda_matrix[g,k]<-lamda_matrix[g,k]+delta[i,k]*t(Phi[,g])%*%solve(Delta)%*%Phi[,g]%*%t(Phi[,g])%*%solve(Delta)%*%Phi[,g]
         
         Alpha_XY[,g,k]<-Alpha_XY[,g,k]+delta[i,k]*t(XXX)%*%solve(Delta)%*%times%*%t(times)%*%solve(Delta)%*%Phi[,g]
          Alpha_matrix[,,g,k]<-Alpha_matrix[,,g,k]+delta[i,k]*t(XXX)%*%solve(Delta)%*%XXX
        }
      }
    }

    PII=rep(NA,Group_NUM)
    for (k in 1:Group_NUM){
      PII[k]<-max(0,1/(1-lambda*Group_NUM)*(1/length(YY)*sum(delta[,k])-lambda))
    }
    PII<-PII/sum(PII)
 
    for (k in 1:Group_NUM){
      Beta[,k]=solve(Beta_matrix[k,,])%*%Beta_XY[,k]
      gamma_int[k]=-Sigma_XY[k]/(Sigma_matrix[k])
      for (g in 1:K_g){
        Lamda[g,g,k]<-lamda_XY[g,k]/(-lamda_matrix[g,k])
        Alpha[g,,k]<-solve(Alpha_matrix[,,g,k])%*%Alpha_XY[,g,k]
      }
    }
    ############identification for groups
    if(Group_NUM==3){
      WZ=order(c(splinex1(0)%*%Beta[1:q_n,1],splinex1(0)%*%Beta[1:q_n,2],splinex1(0)%*%Beta[1:q_n,3]))
      Beta=Beta[,WZ]
      gamma_int=gamma_int[WZ]
      Lamda=Lamda[,,WZ]
      Alpha=Alpha[,,WZ]
      PII=PII[WZ]
    }  

if((is.nan(sum(delta)))&(SStep==1))  {  Updated_Est<-list(Group_number=G_1,PI=PI_int,beta=Beta_1,sigma=gamma_1,lamda=Lamda_1,eigencoefficient=Alpha_1,Sstep=100,del=delta)
    break}
if(is.nan(sum(delta))) {  Updated_Est<-list(Group_number=G_1,PI=PI_int,beta=Beta_1,sigma=gamma_1,lamda=Lamda_1,eigencoefficient=Alpha_1,Sstep=100,del=delta_1)
    break}

 ############identification for eigenfunctions
    for (k in 1:Group_NUM)
    {
      Alpha[,,k]<-t(qr.Q(qr(t(Alpha[,,k]))))
    
    }
   eigen11<-function(t) splinex1(t)%*%Alpha[1,1:q_n,1]
    eigen12<-function(t) splinex1(t)%*%Alpha[2,1:q_n,1]
    eigen21<-function(t) splinex1(t)%*%Alpha[1,1:q_n,2]
    eigen22<-function(t) splinex1(t)%*%Alpha[2,1:q_n,2]
    eigen31<-function(t) splinex1(t)%*%Alpha[1,1:q_n,3]
    eigen32<-function(t) splinex1(t)%*%Alpha[2,1:q_n,3]
if(is.nan(fderiv(eigen11,0.6))){break}
if(fderiv(eigen11,0.6)>0) {Alpha[1,,1]<--Alpha[1,,1]}
 if(fderiv(eigen12,0.6)>0) {Alpha[2,,1]<--Alpha[2,,1]} 

if(fderiv(eigen21,0.6)<0) {Alpha[1,,2]<--Alpha[1,,2]}
 if(fderiv(eigen22,0.6)>0) {Alpha[2,,2]<--Alpha[2,,2]} 

if(fderiv(eigen31,0.6)<0) {Alpha[1,,3]<--Alpha[1,,3]}
 if(fderiv(eigen32,0.6)>0) {Alpha[2,,3]<--Alpha[2,,3]} 

    yucha<-c(max(abs(Beta-Beta_1)),max(abs(gamma_int-gamma_1)),max(abs(PII-PI_1)),max(abs(Lamda-Lamda_1)),max(abs(Alpha-Alpha_1)))
    Epsilon<-max(max(abs(Beta-Beta_1)),max(abs(gamma_int-gamma_1)),max(abs(Lamda-Lamda_1)),max(abs(PII-PI_1)),max(abs(Alpha-Alpha_1)))

############keep the groups with nonzero probabilities
    G_1<-sum(PII!=0) 
    Group_NUM<-sum(PII!=0)
    Beta<-as.matrix(Beta[,PII!=0])
    Beta_1<-Beta
    gamma_int<-gamma_int[PII!=0]
    gamma_1<-gamma_int
    Lamda<-Lamda[,,PII!=0]
    Lamda_1<-Lamda
    Alpha<-Alpha[,,PII!=0]
    Alpha_1<-Alpha
    PI_1=PII[PII!=0]
    PI_int<-PII[PII!=0]
    delta_1<-delta
    df<-sum(PI_int!=0)
    Updated_Est<-list(Group_number=G_1,PI=PI_int,beta=Beta_1,sigma=gamma_1,lamda=Lamda_1,eigencoefficient=Alpha_1,Sstep=SStep,del=delta)
    
     if(G_1<3){break}
    f_den<-matrix(0,nrow=length(YY),ncol=Group_NUM)
    delta<-matrix(nrow=length(YY),ncol=Group_NUM)
    Beta_XY<-matrix(0,nrow=q_n,ncol=Group_NUM)
    Beta_matrix<-array(0,dim=c(Group_NUM,q_n,q_n))
    Sigma_XY<-matrix(0,nrow=1,ncol=Group_NUM)
    Sigma_matrix<-matrix(0,nrow=1,ncol=Group_NUM)
    lamda_XY<-matrix(0,nrow=K_g,ncol=Group_NUM)
    lamda_matrix<-matrix(0,nrow=K_g,ncol=Group_NUM)
    Alpha_XY<-array(0,dim=c(q_n,K_g,Group_NUM))
    Alpha_matrix<-array(0,dim=c(q_n,q_n,K_g,Group_NUM))
  }
  return(Updated_Est)
}
