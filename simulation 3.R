rm(list=ls())
############This is the main code  for simulation 3. It uses functions defined in the files EM_estimate1.R, original function.R, first_der function.R, H_est.R, identify.R.
library(MASS)
library(splines)
library(pracma)
library(copula)
source("EM_estimate1.R")
source("original function.R")
source("first_der function.R")
source("H_est.R")
source("identify.R")

############This function generated simulated data with errors from gamma distribution
sim_data = function(n1,PI,gamma,LLamda,gammashape)
{
  T1<-sort(runif(n1))
  Phi<-matrix(0,n1,K_g)
  index<-runif(1)
  INDEX=1 
  G<-length(PI)-1
  for (l in 1:G){
    if ( index>=sum(PI[1:l]) )
    { INDEX <- INDEX+1 }
  }
  xi<-mvrnorm(1,rep(0,K_g),LLamda[,,INDEX])
  paramMargins<-list()
  for (si in 1:n1)
  {
    paramMargins[[si]]<-list(shape=gammaeshape, scale=1)
  }
  myCop<-normalCopula(param=c(0.1),dim=n1,dispstr="ar1")
  myMvd<-mvdc(copula=myCop, margins=rep("gamma",n1),
              paramMargins)
  Z2<-rMvdc(1,myMvd)
  Z3=Z2-gammaeshape
  error1=as.numeric(Z3%*%diag(rep(1/sqrt(gammaeshape),n1))*sqrt(gamma[INDEX]))
  
  if(INDEX==1){
    Phi[,1]<-eigenfun11(T1)
    Phi[,2]<-eigenfun12(T1)
    Yt=SQ(g1f(T1)+xi[1]*Phi[,1]+xi[2]*Phi[,2]+error1)##
  }
  else{
    if(INDEX==2){
      Phi[,1]<-eigenfun21(T1)
      Phi[,2]<-eigenfun22(T1)
      Yt=SQ(g2f(T1)+xi[1]*Phi[,1]+xi[2]*Phi[,2]+error1)###
    }
    else{
      Phi[,1]<-eigenfun31(T1)
      Phi[,2]<-eigenfun32(T1)
      Yt=SQ(g3f(T1)+xi[1]*Phi[,1]+xi[2]*Phi[,2]+error1)
    }
  }
  data<-list(Y=Yt,TT=T1,IND=INDEX)
}

sim_data_all<-function(n,ni,PI,gamma,LLamda,gammashape){
  data_all<-list()
  for (i in 1:n){
    sample<-sim_data(ni[i],PI,gamma,LLamda,gammashape)
    data_all[[i]]<-sample
  }
  return(data_all)
}

n=400######sample size
lambda<-0.06######penalty parameter
PI<-c(1/3,1/3,1/3) #####the mixing probabilities
q_n<-6###the number of spline basis functions
num=7##########initial number of clusters
K_g<-2#######the number of eigengunctions
gammaeshape=1######the shape parameter of gamma distribution
sstep<-rep(0,300)
ggamma<-c(0.1,0.15,0.2) ######## the varainces of errors
g1f<-function(t) t+sin(pi*t)########mean function for cluster 1
g2f<-function(t) exp(t)########mean function for cluster 2
g3f<-function(t) 2*t^2+2########mean function for cluster 3
LLamda<-array(0,dim=c(K_g,K_g,num))

##############the eigenfunctions
eigenfun11<-function(t) sqrt(2)*cos(1*pi*t)
eigenfun12<-function(t) sqrt(2)*sin(1*pi*t)
eigenfun21<-function(t) sqrt(2)*cos(2*pi*t)
eigenfun22<-function(t) sqrt(2)*cos(1*pi*t)
eigenfun31<-function(t) sqrt(2)*cos(2*pi*t)
eigenfun32<-function(t) sqrt(2)*sin(1*pi*t)

H<-function(t) 3*log(t)####the transformation function in Case 1
SQ<-function(t) exp(t/3)  #inverse function of H


###############the orthonormal spline basis functions
kontx=c(0.4,0.75)###
splinex<-function(t) bs(t,df=NULL,kontx,degree=3,intercept=T,Boundary.knots=c(0,1))
work<-matrix(0,q_n,q_n)
pos<-matrix(0,q_n,1)
for (i in 1:10000)
{work=work+splinex(i/10000)[1,]%o%splinex(i/10000)[1,]
}
working=work/10000
VV<-solve(sqrtm(working)$B)
splinex1<- function(t) splinex(t)%*%VV

num_group=rep(NA,300)
prob=matrix(NA,300,length(PI))
parmet=array(NA,dim=c(300,q_n,length(PI)))##
stderror=matrix(NA,300,length(PI))
eigenvalue=array(NA,dim=c(300,K_g,K_g,length(PI)))
eigencoe=array(NA,dim=c(300,K_g,q_n,length(PI)))
H_fun=matrix(NA,300,170)
cluster_result<-matrix(NA,300,800)

###########simulate 300 datasets
for(iter in 1:300){
  ####################the eigenvalues
  LLamda[1,1,1]<-1
  LLamda[2,2,1]<-0.25
  
  LLamda[1,1,2]<-1.1
  LLamda[2,2,2]<-0.2
  
  LLamda[1,1,3]<-0.9
  LLamda[2,2,3]<-0.15
  
  for (i in 4:num) {
    LLamda[1,1,i]<-1
    LLamda[2,2,i]<-0.2
  } 
  ####################generating data
  set.seed(iter)
  ni<-ceiling(runif(n,8,12))####the number of observed time points 
  data<-sim_data_all(n,ni,PI,ggamma,LLamda,gammaeshape)
  YY<-list()
  TT<-list()
  IND<-list()
  n<-length(data)
  for (i in 1:n) {
    YY[[i]]<-data[[i]]$Y######the functional observation
    TT[[i]]<-data[[i]]$TT####obvserved time
    IND[[i]]<-data[[i]]$IND####true clutser id of the data with the individuals in the same order
  }
  true_index<-rep(NA,n)
  
  #######################initialization BBeta, SSigma, LLamda, eigencoef, PI_int for the basis coefficient of mean functions, the varainces of errors, eigenvalues, the basis coefficient of eigenfunctions, the mixing probabilities
  gamma_int<-rep(0.15,num)
  BBeta_XY<-matrix(0,nrow=q_n,ncol=num) 
  BBeta_matrix<-array(0,dim=c(num,q_n,q_n))  
  BBeta<-matrix(0,nrow=q_n,ncol=num) 
  
  Sigma_XY<-matrix(0,nrow=1,ncol=num)
  Sigma_matrix<-matrix(0,nrow=1,ncol=num)
  SSigma<-matrix(0,nrow=1,ncol=num)
  
  lamda_XY<-matrix(0,nrow=K_g,ncol=num)
  lamda_matrix<-matrix(0,nrow=K_g,ncol=num)
  
  eigencoef<-array(0,dim=c(K_g,q_n,num))
  
  YYlen<-rep(NA,400)
  for (i in 1:length(YY))
    YYlen[i] <-length(unlist(YY[[i]]))
  lenmin<-YYlen[which.min(YYlen[1:400])]
  lenmax<-YYlen[which.max(YYlen[1:400])]
  
  ymm<-matrix(0,nrow=length(YY),ncol=lenmin)
  for (i in 1:length(YY))
  { ymm[i,]=YY[[i]][1:lenmin] 
  true_index[i]=IND[[i]]
  }
  
  kgroup<-kmeans(H(ymm),num)
  PI_int<-kgroup$size/sum(kgroup$size)
  
  for (k in 1:num){
    IND<-kgroup$cluster==k 
    Ygroup<-YY[IND] 
    Tgroup<-TT[IND]
    for (i in 1:length(Ygroup)){
      Phi<-matrix(0,length(Tgroup[[i]]),K_g)
      Phi[,1]<-eigenfun11(Tgroup[[i]])
      Phi[,2]<-eigenfun12(Tgroup[[i]])
      Delta<-matrix(nrow=length(Tgroup[[i]]),ncol=length(Tgroup[[i]]))
      Delta<-Phi%*%LLamda[,,k]%*%t(Phi)+gamma_int[k]*diag(1,length(Tgroup[[i]]))
      
      BX<-splinex1(Tgroup[[i]])###
      BBeta_XY[,k]<-BBeta_XY[,k]+t(BX)%*%solve(Delta)%*%(H(Ygroup[[i]]))
      BBeta_matrix[k,,]<-BBeta_matrix[k,,]+t(BX)%*%solve(Delta)%*%BX
    }
    BBeta[,k]=solve(BBeta_matrix[k,,]+diag(0.000001,q_n))%*%BBeta_XY[,k]
  }
  samp=seq(0.001,0.999,length=100)
  desig<-splinex1(samp)
  project<-solve(t(desig)%*%desig)%*%t(desig)
  
  for (k in 1:num){
    eigencoef[1,,k]<-project%*%eigenfun11(samp)
    eigencoef[2,,k]<-project%*%eigenfun12(samp)
  }
  
  for (k in 1:num){
    SSigma[k]<-ggamma[1]
  }
  Group_NUM<-rep(NA,500)
  PI_MATRIX<-list()
  ALPHABETA<-list()
  ERRORSIGMA<-list()
  EIGENLAMDA<-list()
  EIGENCOEFFI<-list()
  HY<-list()
  Group_NUM[[1]]=num
  PI_MATRIX[[1]]=PI_int
  ALPHABETA[[1]]=BBeta
  ERRORSIGMA[[1]]=SSigma
  EIGENLAMDA[[1]]=LLamda
  EIGENCOEFFI[[1]]=eigencoef
  
  ############initialization H
  geshu=70
  lowbound=max(sort(unlist(YY),decreasing=F)[1:20]) 
  upbbound=min(sort(unlist(YY),decreasing=T)[1:10])
  point=seq(lowbound,upbbound,length=geshu)
  DESIMATRX_low=cbind(rep(1,5),point[1:5])##
  DESIMATRX_up=cbind(rep(1,10),point[(geshu-9):geshu])
  textpoint1=seq(upbbound,max(unlist(YY)),length=50)
  textpoint2=seq(min(unlist(YY)),lowbound,length=50)
  allpoint=c(textpoint2,point,textpoint1)
  HY[[1]]=H(allpoint)
  dif=10000
  likelibefore=100
  
  ############################begin iteration 
  for(i in 1:10)####
  {
    ############################estimate Omega_n given H
    EM_est<-EM_estimate(TT,YY,Group_NUM[i],PI_MATRIX[[i]],ALPHABETA[[i]],ERRORSIGMA[[i]],EIGENLAMDA[[i]],EIGENCOEFFI[[i]],HY[[i]],lambda)
    Group_NUM[i+1]<-EM_est$Group_number
    PI_MATRIX[[i+1]]<-EM_est$PI
    ALPHABETA[[i+1]]<-EM_est$beta
    ERRORSIGMA[[i+1]]<-EM_est$sigma
    EIGENLAMDA[[i+1]]<-EM_est$lamda
    EIGENCOEFFI[[i+1]]<-EM_est$eigencoefficient
    delta<-EM_est$del
    sstep[iter]=i
   if(EM_est$Sstep==100) {  HY[[i+1]]<-HY[[i]]
    break}
    if(Group_NUM[i+1]<length(PI)) {  HY[[i+1]]<-HY[[i]]
    break}
    
    ############################estimate H given Omega_n 
    HY[[i+1]]<-H_est(YY,PI_MATRIX[[i+1]],ALPHABETA[[i+1]],ERRORSIGMA[[i+1]],EIGENLAMDA[[i+1]],EIGENCOEFFI[[i+1]],Group_NUM[i+1],HY[[i]][51:120],10)
    
    #############the computation of probabilities for individuals belong to each group   
    delta<-matrix(0,nrow=length(YY),ncol=EM_est$Group_number)
    f_den<-matrix(0,nrow=length(YY),ncol=EM_est$Group_number) 
    lHYY<-unlist(HY[[i+1]])
    lAlpha<-EM_est$eigencoefficient
    lLamda<-EM_est$lamda
    lgamma_int<-EM_est$sigma
    lBeta<-EM_est$beta
    lPI_int<-EM_est$PI
    likeliafter=0
    
    for (j in 1:length(YY)){
      XXX<-splinex1(TT[[j]])
      Phi<-matrix(NA,length(TT[[j]]),K_g)
      MU<-matrix(NA,length(TT[[j]]),1)
      lHYYY<-lHYY[identify(YY[[j]],allpoint)]
      for (k in 1:EM_est$Group_number){
        lPhi<-XXX%*%t(lAlpha[,,k])
        Delta<-matrix(nrow=length(TT[[j]]),ncol=length(TT[[j]]))
        Delta<-lPhi%*%lLamda[,,k]%*%t(lPhi)+lgamma_int[k]*diag(1,ni[j])##
        lMU<-XXX%*%lBeta[,k]
        f_den[j,k]<-(2*pi)^(-length(TT[[j]])/2)/sqrt(det(Delta))*exp(-1/2*t(lHYYY-lMU)%*%solve(Delta)%*%(lHYYY-lMU))
      }
    }
    for (j in 1:length(YY)){
      for (k in 1:EM_est$Group_number){
        delta[j,k]<-lPI_int[k]*f_den[j,k]/((t(lPI_int)%*%f_den[j,]))
        likeliafter=likeliafter+delta[j,k]*log(f_den[j,k]%*%lPI_int[k])
      }
    }
    
    cat("iter=",iter,"\n")
    H_diff=mean(abs(HY[[i+1]][51:120]-HY[[i]][51:120]))/mean(abs(HY[[i]][51:120]))
    if(Group_NUM[i+1]-Group_NUM[i]==0){dif=mean(c(mean(abs(sort(PI_MATRIX[[i+1]])-sort(PI_MATRIX[[i]]))),mean(abs(ALPHABETA[[i+1]]-ALPHABETA[[i]])),mean(abs(ERRORSIGMA[[i+1]]-ERRORSIGMA[[i]])),mean(abs(EIGENLAMDA[[i+1]]-EIGENLAMDA[[i]])),mean(abs(EIGENCOEFFI[[i+1]]-EIGENCOEFFI[[i]]))))}
    if(dif<1e-3&H_diff<1e-3){ break }
    
  }
  #########################the end of iteration
  
  num_group[iter]=Group_NUM[i+1]#######the estimated number of clusters
  if(Group_NUM[i+1]!=3){
    prob[iter,]=NA
    parmet[iter,,]=NA
    stderror[iter,]=NA
    eigenvalue[iter,,,]=NA
    eigencoe[iter,,,]=NA
    H_fun[iter,]=NA
  }
  if(Group_NUM[i+1]!=3){next}
  
 #############the Bayes' optimal allocation rule
  categorial<-rep(NA,n)
  for (j in 1:n)
  {
    categorial[j]<-which.max(delta[j,])
  }
  cluster_result[iter,1:400]<-true_index#######the true clutser id
  cluster_result[iter,401:800]<-categorial#######the estimated clutser id
  ##############the estimaed results
  prob[iter,]=PI_MATRIX[[i+1]]#######the probabilities
  parmet[iter,,]=ALPHABETA[[i+1]]#######the basis coefficient of mean functions
  stderror[iter,]=ERRORSIGMA[[i+1]]#######the varainces of errors
  eigenvalue[iter,,,]=EIGENLAMDA[[i+1]]#######the diagonal matrix of eigenvalue
  eigencoe[iter,,,]=EIGENCOEFFI[[i+1]]#######the basis coefficient of eigenfunctions
  H_fun[iter,]=HY[[i+1]]#######the value of H
}
