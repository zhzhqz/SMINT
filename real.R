rm(list=ls())
############This is the main code  for ADNI study data. It uses functions defined in the file original function.R, first_der function.R, H_est.R, identify.R.
########the true clutser id: trueclass
########the estimated clutser id: categorial
########the estimated number of clusters: G_1
########the estimated mixing probabilities: PI_int
########the estimated basis coefficient of mean functions: Beta_1
########the estimated varainces of errors: gamma_1
########the estimated eigenvalue: Lamda_1
########the estimated basis coefficient of eigenfunctions: Alpha_1
########the estimated H: HY[[l+1]]

library(MASS)
library(splines)
library(pracma)
library(fda)
library(rJava)
library(xlsx)
source("original function.R")
source("first_der function.R")
source("H_est.R")
source("identify.R")
YY<-read.xlsx2("ADNIdata.xlsx",sheetIndex=1) 
n<-nrow(YY)######sample size
nii<-ncol(YY)######the number of observed points
data<-matrix(NA,n,nii)
for (i in 1:nrow(YY))
{ data[i,]=as.numeric(t(YY[i,1:nii])) }#####the functional observation
ni<-rep(nii,nrow(YY))
tt<-seq(0,1,length.out=nii)
TT<-list()
for (i in 1:n)
  TT[[i]]<-tt
ymm<-data
ntrain<-n
lambda<-0.002######penalty parameter
q_n<-6###the number of spline basis functions
num=6###########initial number of clusters
K_g<-2######the number of eigengunctions
sstep<-rep(0,n)

###############the orthonormal spline basis functions
kontx=c(0.4,0.75)
splinex<-function(t) bs(t,df=NULL,kontx,degree=3,intercept=T,Boundary.knots=c(0,1))
work<-matrix(0,q_n,q_n)
pos<-matrix(0,q_n,1)
for (i in 1:10000)
{work=work+splinex(i/10000)[1,]%o%splinex(i/10000)[1,]
}
working=work/10000
VV<-solve(sqrtm(working)$B)
splinex1<- function(t) splinex(t)%*%VV

categorial<-rep(NA,ntrain)

#######################initialization BBeta, SSigma, LLamda, eigencoef, PI_int for the basis coefficient of mean functions, the varainces of errors, eigenvalues, the basis coefficient of eigenfunctions, the mixing probabilities
gamma_int<-rep(0.1,num)####
BBeta_XY<-matrix(0,nrow=q_n,ncol=num) ###
BBeta_matrix<-array(0,dim=c(num,q_n,q_n))  
BBeta<-matrix(0,nrow=q_n,ncol=num) 

Sigma_XY<-matrix(0,nrow=1,ncol=num)
Sigma_matrix<-matrix(0,nrow=1,ncol=num)
SSigma<-matrix(0,nrow=1,ncol=num)

lamda_XY<-matrix(0,nrow=K_g,ncol=num)
lamda_matrix<-matrix(0,nrow=K_g,ncol=num)

eigencoef<-array(0,dim=c(K_g,q_n,num))
SSigma<-gamma_int

Group_N<-rep(NA,500)
PI_MATRIX<-list()
ALPHABETA<-list()
ERRORSIGMA<-list()
EIGENLAMDA<-list()
EIGENCOEFFI<-list()
HY<-list()


SSigma<-rep(0.1,num)
LLamda<-array(0,dim=c(K_g,K_g,num))
eigenfun11<-function(t) sqrt(2)*cos(1*pi*t)
eigenfun12<-function(t) sqrt(2)*sin(1*pi*t)
H<-function(t) 80*t
for (i in 1:num) 
{
  LLamda[1,1,i]<-1
  LLamda[2,2,i]<-0.5
} 
ytrain<-ymm
kgroup<-kmeans(H(ytrain),num)
PI_int<-kgroup$size/sum(kgroup$size)

a<-splinex1(tt)
b1<-eigenfun11(tt)
b2<-eigenfun12(tt)
Hessian<-solve(t(a)%*%a)%*%t(a)
for (k in 1:num)
{
  IND<-kgroup$cluster==k 
  Ygroup<-ytrain[IND,] 
  
  b<-colMeans(Ygroup)
  BBeta[,k]<-Hessian%*%b
  eigencoef[1,,k]<-Hessian%*%b1
  eigencoef[2,,k]<-Hessian%*%b2
}
Group_N[[1]]=num
PI_MATRIX[[1]]=PI_int
ALPHABETA[[1]]=BBeta
ERRORSIGMA[[1]]=SSigma
EIGENLAMDA[[1]]=LLamda
EIGENCOEFFI[[1]]=eigencoef

############initialization H
geshu=50
lowbound=max(sort(ymm,decreasing=F)[1:20]) 
upbbound=max(sort(ymm,decreasing=T)[1:20])
point=seq(lowbound,upbbound,length=geshu)
DESIMATRX_low=cbind(rep(1,5),point[1:5])##
DESIMATRX_up=cbind(rep(1,10),point[(geshu-9):geshu])
textpoint1=seq(upbbound,max(ymm),length=50)
textpoint2=seq(min(ymm),lowbound,length=50)
allpoint=c(textpoint2,point,textpoint1)
HY[[1]]=H(allpoint)

Alpha<-eigencoef
Group_NUM<-num
Lamda<-LLamda
Beta<-BBeta
delta<-matrix(0,nrow=n,ncol=Group_NUM) 
gamma_int<-SSigma
G_1<-Group_NUM
Beta_1<-Beta
gamma_1<-gamma_int
Lamda_1<-Lamda
PI_1=PI_int
for (k in 1:Group_NUM)
{
  Alpha[,,k]<-t(qr.Q(qr(t(Alpha[,,k]))))
}

Alpha_1<-Alpha
XXX<-splinex1(tt)


############################begin iteration 
for(l in 1:17)
{
  
  HYY<-HY[[l]]
  SStep<-0
  Epsilon<-1
  ############################estimate Omega_n given H
  while((SStep<300)&(Epsilon>10^(-3))){
    SStep<-SStep+1
    print(SStep)
    likelihood=0
    f_den<-matrix(0,nrow=n,ncol=Group_NUM)
    delta<-matrix(nrow=n,ncol=Group_NUM)
    Beta_XY<-matrix(0,nrow=q_n,ncol=Group_NUM)
    Beta_matrix<-array(0,dim=c(Group_NUM,q_n,q_n))
    Sigma_XY<-matrix(0,nrow=1,ncol=Group_NUM)
    Sigma_matrix<-matrix(0,nrow=1,ncol=Group_NUM)
    lamda_XY<-matrix(0,nrow=K_g,ncol=Group_NUM)
    lamda_matrix<-matrix(0,nrow=K_g,ncol=Group_NUM)
    Alpha_XY<-array(0,dim=c(q_n,K_g,Group_NUM))
    Alpha_matrix<-array(0,dim=c(q_n,q_n,K_g,Group_NUM))
    
    for (i in 1:n){
      HYYY<-HYY[identify(ytrain[i,],allpoint)]
      for (k in 1:Group_NUM){
        
        Phi<-XXX%*%t(Alpha[,,k])
        Delta<-matrix(nrow=nii,ncol=nii)
        Delta<-Phi%*%Lamda[,,k]%*%t(Phi)+gamma_int[k]*diag(1,ni[1])##
        
        MU<-XXX%*%Beta[,k]
        f_den[i,k]<-(2*pi)^(-nii/2)/sqrt(det(Delta))*exp(-1/2*t(HYYY-MU)%*%solve(Delta)%*%(HYYY-MU))
      }
      likelihood=likelihood+log(f_den[i,]%*%PI_int)
    }
    
    for (i in 1:n){
      HHYY<-HYY[identify(ytrain[i,],allpoint)]
      for (k in 1:Group_NUM){
        Phi<-XXX%*%t(Alpha[,,k])
        Delta<-matrix(nrow=nii,ncol=nii)
        Delta<-Phi%*%Lamda[,,k]%*%t(Phi)+gamma_int[k]*diag(1,nii)##
        delta[i,k]<-PI_int[k]*f_den[i,k]/((t(as.matrix(PI_int))%*%as.matrix(f_den[i,])))
        
        MU<-XXX%*%Beta[,k]
        times<-matrix(HHYY-MU,nii,1)
        Beta_XY[,k]<-Beta_XY[,k]+delta[i,k]*t(XXX)%*%solve(Delta)%*%HHYY
        Beta_matrix[k,,]<-Beta_matrix[k,,]+delta[i,k]*t(XXX)%*%solve(Delta)%*%XXX
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
      PII[k]<-max(0,1/(1-lambda*Group_NUM)*(1/length(ymm)*sum(delta[,k])-lambda))
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
    
    for (k in 1:Group_NUM)
    {
      Alpha[,,k]<-t(qr.Q(qr(t(Alpha[,,k]))))
      Alpha[1,,k]<-Alpha[1,,k]*sign(Alpha[1,1,k])
      Alpha[2,,k]<-Alpha[2,,k]*sign(Alpha[2,1,k])
    }
    yucha<-c(max(abs(Beta-Beta_1)),max(abs(gamma_int-gamma_1)),max(abs(PII-PI_1)),max(abs(Lamda-Lamda_1)),max(abs(Alpha-Alpha_1)))
    Epsilon<-max(max(abs(Beta-Beta_1)),max(abs(gamma_int-gamma_1)),max(abs(Lamda-Lamda_1)),max(abs(PII-PI_1)),max(abs(Alpha-Alpha_1)))
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
    df<-sum(PI_int!=0)
  }  
  
  ##############the estimaed results for Omega_n 
  Group_N[l+1]<-G_1
  PI_MATRIX[[l+1]]<-PI_int
  ALPHABETA[[l+1]]<-Beta_1
  ERRORSIGMA[[l+1]]<-gamma_1
  EIGENLAMDA[[l+1]]<-Lamda_1
  EIGENCOEFFI[[l+1]]<-Alpha_1
  
  ############################estimate H given Omega_n 
  HY[[l+1]]<-H_est(ytrain,PI_MATRIX[[l+1]],ALPHABETA[[l+1]],ERRORSIGMA[[l+1]],EIGENLAMDA[[l+1]],EIGENCOEFFI[[l+1]],Group_NUM,HY[[l]][51:100],100)
  H_diff=mean(abs(HY[[l+1]]-HY[[l]]))
  
  dif=10000
  if(Group_N[l+1]-Group_N[l]==0)
  {
    dif=mean(c(mean(abs(sort(PI_MATRIX[[l+1]])-sort(PI_MATRIX[[l]]))),mean(abs(ALPHABETA[[l+1]]-ALPHABETA[[l]])),mean(abs(ERRORSIGMA[[l+1]]-ERRORSIGMA[[l]])),mean(abs(EIGENLAMDA[[l+1]]-EIGENLAMDA[[l]])),mean(abs(EIGENCOEFFI[[l+1]]-EIGENCOEFFI[[l]]))))
  }
  if(dif<1e-4&H_diff<1e-3){ break }
}
#########################the end of iteration

#############the Bayes' optimal allocation rule

categorial<-rep(NA,768)
for (i in 1:n)
{
  categorial[i]<-which.max(delta[i,])
}
trueclass<-c(rep(1,172),rep(2,218),rep(3,378)) 
sum(trueclass!=categorial) 