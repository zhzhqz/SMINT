rm(list=ls())
############This is the main code  for simulation 2. It uses functions defined in the files EM_estimate2.R, original function.R, first_der function.R, H_est.R, identify.R.
########the true clutser id: cluster_result[,1:n]
########the estimated clutser id: cluster_result[,(n+1):(2*n)]
library(MASS)
library(splines)
library(pracma)
source("EM_estimate2.R")
source("original function.R")
source("first_der function.R")
source("H_est.R")
source("identify.R")

############This function generated simulated data
sim_data = function(n1,PI,gamma,LLamda)
{
  T1<-sort(runif(n1))
  Phi<-matrix(NA,n1,K_g)
  index<-runif(1)
  INDEX=1 
  G<-length(PI)-1
  for (l in 1:G){
    if ( index>=sum(PI[1:l]) )
    { INDEX <- INDEX+1 }
  }
  xi<-mvrnorm(1,rep(0,K_g),LLamda[,,INDEX])
  Sigma_r=matrix(nrow=n1,ncol=n1)
  Sigma_r=gamma[INDEX]*diag(1,n1)
  error1<-mvrnorm(1,rep(0,n1),Sigma_r)
  if(INDEX==1){
    Phi[,1]<-eigenfun11(T1)
    Phi[,2]<-eigenfun12(T1)
    Yt=SQ(g1f(T1)+xi[1]*Phi[,1]+xi[2]*Phi[,2]+error1)
  }
  else{
    if(INDEX==2){
      Phi[,1]<-eigenfun21(T1)
      Phi[,2]<-eigenfun22(T1)
      Yt=SQ(g2f(T1)+xi[1]*Phi[,1]+xi[2]*Phi[,2]+error1)
    }
    else{
      Phi[,1]<-eigenfun31(T1)
      Phi[,2]<-eigenfun32(T1)
      Yt=SQ(g3f(T1)+xi[1]*Phi[,1]+xi[2]*Phi[,2]+error1)
    }
  }
  data<-list(Y=Yt,TT=T1,IND=INDEX)
}

sim_data_all<-function(n,ni,PI,gamma,LLamda){
  data_all<-list()
  for (i in 1:n){
    sample<-sim_data(ni[i],PI,gamma,LLamda)
    data_all[[i]]<-sample
  }
  return(data_all)
}
n=400######sample size
lambda<-0.03#######penalty parameter
PI<-c(0.229,0.277,0.494)###the mixing probabilities
q_n<-6###the number of spline basis functions
num=7########initial number of clusters
K_g<-2#######the number of eigengunctions
sstep<-rep(0,300)
ggamma<-c(0.306,0.302,0.018) ######## the varainces of errors

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

###########mean functions
BBeta_esco<-matrix(c(-0.17069210,0.02013216,-0.08453613
                ,-0.05221607, 0.24048968,0.10558663
                ,0.24767664, 0.81013221,0.46133860
                ,0.16861082, 0.71978989, 0.33228626
                ,-0.11147241, 0.20341094,-0.02377541
                ,-0.12787621, 0.02899900, -0.07583228),6,3,byrow=T)

g1f<-function(t) splinex1(t)%*%BBeta_esco[,1]########mean function for cluster 1
g2f<-function(t) splinex1(t)%*%BBeta_esco[,2]########mean function for cluster 2
g3f<-function(t) splinex1(t)%*%BBeta_esco[,3]########mean function for cluster 3
LLamda<-array(0,dim=c(K_g,K_g,num))
dd<-seq(0,1,length.out=100)
##############the eigenfunctions
eigen_esco<-matrix(c(0.3042409,0.58428511,0.4206789,-0.2937568, -0.47816660, -0.27228903,
                     0.2555993, 0.04120785, -0.6963424, -0.6681271, -0.01160255,  0.03936685,
                     0.2798120, 0.55781747, 0.4003899,-0.4618208, -0.4521572, -0.1803006,
                     0.2674949, 0.09281317, -0.7061496, -0.5941956, 0.1904865, 0.1784132,
                     0.36565672,0.5364813,0.1853044,-0.4806000,-0.49674627,-0.2577055,
                     0.08910582,-0.2257959,-0.7727863,-0.5725718,0.02278733,0.1245789),6,6,byrow=T)
eigenfun11<-function(t) splinex1(t)%*%eigen_esco[1,]
eigenfun12<-function(t) splinex1(t)%*%eigen_esco[2,]
eigenfun21<-function(t) splinex1(t)%*%eigen_esco[3,]
eigenfun22<-function(t) splinex1(t)%*%eigen_esco[4,]
eigenfun31<-function(t) splinex1(t)%*%eigen_esco[5,]
eigenfun32<-function(t) splinex1(t)%*%eigen_esco[6,]

####the transformation function
kontx_H=(seq(-0.03,0.04,length=8))[2:7]
splinex_H<-function(t) bs(t,df=NULL,kontx_H,degree=3,intercept=T,Boundary.knots=c(-0.035,0.05))
H<-function(t) splinex_H(t)%*%c(-4.0793174,-3.5300838,-1.7940349,-0.9425109,-0.1370304,0.3993313,0.8920254,1.6257655, 5.0469606,-2.1223741)

#inverse function of H
kontx_SQ=(seq(-3.1,3,length=8))[2:7]
splinex_SQ<-function(t) bs(t,df=NULL,kontx_SQ,degree=3,intercept=T,Boundary.knots=c(-3.8,3.8))
SQ<-function(t) splinex_SQ(t)%*%c(-0.036439888,-0.027476215,-0.025819601,-0.015612284,-0.004704981,0.008755589,0.028535617,0.032098119,0.040827433,0.056836708)

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
  LLamda[1,1,1]<-0.163
  LLamda[2,2,1]<-0
  
  LLamda[1,1,2]<-0
  LLamda[2,2,2]<-0
  
  LLamda[1,1,3]<-0.308
  LLamda[2,2,3]<-0.16
  
  for (i in 4:num) {
    LLamda[1,1,i]<-0.2
    LLamda[2,2,i]<-0.16
  } 
  
  ####################generating data
  set.seed(iter)
  ni<-ceiling(runif(n,8,12))####the number of observed time points 
  data<-sim_data_all(n,ni,PI,ggamma,LLamda)
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
  gamma_int<-rep(0.2,num)
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
  
  zhojia<-matrix(H(c(ymm)),n,lenmin)
  kgroup<-kmeans(ymm,num)
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
    BBeta[,k]=solve(BBeta_matrix[k,,])%*%BBeta_XY[,k]
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
  lowbound=max(sort(unlist(YY),decreasing=F)[1:30])
  upbbound=min(sort(unlist(YY),decreasing=T)[1:10])
  point=seq(lowbound,upbbound,length=geshu)
  DESIMATRX_low=cbind(rep(1,5),point[1:5])##
  DESIMATRX_up=cbind(rep(1,10),point[(geshu-9):geshu])
  textpoint1=seq(upbbound,max(unlist(YY)),length=50)
  textpoint2=seq(min(unlist(YY)),lowbound,length=50)
  par(mfrow=c(3,3))
  allpoint=c(textpoint2,point,textpoint1)
  HY[[1]]=H(allpoint)
  dif=10000
  likelibefore=100
  delta<-matrix(NA,n,3)
  ############################begin iteration 
  for(i in 1:10)
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
    HY[[i+1]]<-H_est(YY,PI_MATRIX[[i+1]],ALPHABETA[[i+1]],ERRORSIGMA[[i+1]],EIGENLAMDA[[i+1]],EIGENCOEFFI[[i+1]],Group_NUM[i+1],HY[[i]][51:120],100)
    
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

