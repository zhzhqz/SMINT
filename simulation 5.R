rm(list=ls())
############This is the main code  for simulation 5. It uses functions defined in the files EM_estimate4.R, original function.R, first_der function.R, H_est.R, identify.R.
########the true clutser id: cluster_result[,1:n]
########the estimated clutser id: cluster_result[,(n+1):(2*n)]
library(MASS)
library(splines)
library(pracma)
source("EM_estimate4.R")
source("original function.R")
source("first_der function.R")
source("H_est.R")
source("identify.R")

############This function generated simulated data
sim_data = function(i,n1,PI,gamma)
{
  T1<-seq(1,21,length.out=n1)
  Phi<-matrix(0,n1,K_g)
  
  index<-runif(1)
  INDEX=1
  G<-length(PI)-1
  for (l in 1:G){
    if ( index>=sum(PI[1:l]) )
    { INDEX <- INDEX+1 }
  }
  xi<-runif(2)
  
  Sigma_r=gamma[1]*diag(1,n1)
  error1<-mvrnorm(1,rep(0,n1),Sigma_r)
  if(INDEX==1){
    Phi[,1]<-eigenfun1(T1)
    Yt=SQ(g1f(T1)+Phi%*%xi+error1)
    
  }
  else{
    Phi[,1]<-eigenfun2(T1)
    Yt=SQ(g2f(T1)+Phi%*%xi+error1)
  }
  
  data<-list(Y=Yt,TT=T1,IND=INDEX)
}

sim_data_all<-function(n,ni,PI,gamma){
  data_all<-list()
  for (i in 1:n){
    sample<-sim_data(i,ni[i],PI,gamma)
    data_all[[i]]<-sample
  }
  return(data_all)
}

sstep<-rep(0,400)
q_n<-10###the number of spline basis functions
################the orthonormal spline basis functions
kontx=c(0.1,0.3,0.45,0.60,0.75,0.9)*21
splinex<-function(t) bs(t,df=NULL,kontx,degree=3,intercept=T,Boundary.knots=c(1,21))
work<-matrix(0,q_n,q_n)
pos<-matrix(0,q_n,1)
dian<-seq(1,21,length.out=10000)
for (i in 1:10000)
{work=work+splinex(dian[i])[1,]%o%splinex(dian[i])[1,]
}
working=work/10000
VV<-solve(sqrtm(working)$B)
splinex1<- function(t) splinex(t)%*%VV

###########simulate 100 datasets
N=100
n=100
true_num<-2
K_g<-2
cluster_result<-matrix(NA,N,(2*n))
num_group=rep(NA,N)
prob=matrix(NA,N,true_num)
parmet=array(NA,dim=c(N,q_n,true_num))##
stderror=matrix(NA,N,true_num)
eigenvalue=array(NA,dim=c(N,K_g,K_g,true_num))
eigencoe=array(NA,dim=c(N,K_g,q_n,true_num))
H_fun=matrix(NA,N,170)

for (iter in 1:N)
{
  n=100##########sample size
  lambda<-0.15#########penalty parameter
  
  PI=c(1/2,1/2) ###the mixing probabilities
  K_g<-2#######the number of eigengunctions
  true_num<-2#######true number of clusters
  ggamma<-rep(0.5,4) ######## the varainces of errors
  h1<-function(t) 6-abs(t-7)######## the function h1(t) in Bouveyron et al.(2015)
  h2<-function(t) 6-abs(t-15)######## the function h2(t) in Bouveyron et al.(2015)
  
  H<-function(t) t ########the transformation function in Bouveyron et al.(2015)
  SQ<-function(t) t ########the inverse function of H in Bouveyron et al.(2015)
  
  eigenfun1<-function(t) 1-h1(t)
  eigenfun2<-function(t) 1-h2(t)
  
  g1f<-function(t) h1(t)
  g2f<-function(t) h2(t)
  
  num=7####initial number of clusters
  LLamda<-array(0,dim=c(K_g,K_g,num))
  
  ####################generating data
  set.seed(iter)
  ni<-ceiling(runif(n,101,101))####the number of observed time points 
  data<-sim_data_all(n,ni,PI,ggamma)
  YY<-list()
  TT<-list()
  IND<-list()
  for (i in 1:n) {
    YY[[i]]<-data[[i]]$Y######the functional observation
    TT[[i]]<-data[[i]]$TT####obvserved time
    IND[[i]]<-data[[i]]$IND####true clutser id of the data with the individuals in the same order
  }
  true_index<-rep(NA,n)
  
  #######################initialization BBeta, SSigma, LLamda, eigencoef, PI_int for the basis coefficient of mean functions, the varainces of errors, eigenvalues, the basis coefficient of eigenfunctions, the mixing probabilities
  gamma_int<-rep(0.5,num)
  BBeta_XY<-matrix(0,nrow=q_n,ncol=num) 
  BBeta_matrix<-array(0,dim=c(num,q_n,q_n)) 
  BBeta<-matrix(0,nrow=q_n,ncol=num) 
  
  Sigma_XY<-matrix(0,nrow=1,ncol=num)
  Sigma_matrix<-matrix(0,nrow=1,ncol=num)
  SSigma<-matrix(0,nrow=1,ncol=num)
  
  lamda_XY<-matrix(0,nrow=K_g,ncol=num)
  lamda_matrix<-matrix(0,nrow=K_g,ncol=num)
  eigencoef<-array(0,dim=c(K_g,q_n,num))
  
  ymm<-matrix(0,nrow=n,ncol=ni[1])###
  for (i in 1:n)
  { ymm[i,]=YY[[i]][1:ni[1]] 
  true_index[i]=IND[[i]]
  }
  
  kgroup<-kmeans(H(ymm),num)
  PI_int<-kgroup$size/sum(kgroup$size)
  
  Phi<-matrix(0,length(TT[[1]]),K_g)
  Phi[,1]<-eigenfun1(TT[[1]])
  
  for (i in 1:num) {
    LLamda[1,1,i]<-1/10
    LLamda[2,2,i]<-1/10
  }
  
  for (k in 1:num){
    IND<-kgroup$cluster==k 
    Ygroup<-YY[IND] ####
    Tgroup<-TT[IND]
    for (i in 1:length(Ygroup)){
      Delta<-matrix(nrow=length(Tgroup[[i]]),ncol=length(Tgroup[[i]]))
      Delta<-Phi%*%LLamda[,,k]%*%t(Phi)+gamma_int[k]*diag(1,length(Tgroup[[i]]))##
      
      BX<-splinex1(Tgroup[[i]])###
      BBeta_XY[,k]<-BBeta_XY[,k]+t(BX)%*%solve(Delta)%*%(H(Ygroup[[i]]))###
      BBeta_matrix[k,,]<-BBeta_matrix[k,,]+t(BX)%*%solve(Delta)%*%BX
    }
    BBeta[,k]=solve(BBeta_matrix[k,,]+diag(0.00001,q_n))%*%BBeta_XY[,k]###
  }
  
  samp=seq(1,21,length=100)
  desig<-splinex1(samp)
  project<-solve(t(desig)%*%desig)%*%t(desig)
  
  for (k in 1:num){
    
    eigencoef[1,,k]<-project%*%eigenfun1(samp)
    eigencoef[2,,k]<-project%*%eigenfun1(samp)
  }
  
  
  for (k in 1:num){
    SSigma[k]<-ggamma[1]
  }
  index_data<-seq(1,ni[1],by=6)
  index_data<-c(index_data,ni[1])
  YY<-list()
  TT<-list()
  
  for (i in 1:n) {
    YY[[i]]<-(data[[i]]$Y)[index_data]
    TT[[i]]<-(data[[i]]$TT)[index_data]
    
  }
  ni<-rep(18,n)
  SSigma<-gamma_int
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
  delta<-matrix(NA,n,true_num)
  xu=12
  ############################begin iteration 
  for(i in 1:10)
  {
    ############################estimate Omega_n given H
    EM_est<-EM_estimate(TT,YY,Group_NUM[i],PI_MATRIX[[i]],ALPHABETA[[i]],ERRORSIGMA[[i]],EIGENLAMDA[[i]],EIGENCOEFFI[[i]],HY[[i]],lambda,xu)
    Group_NUM[i+1]<-EM_est$Group_number
    PI_MATRIX[[i+1]]<-EM_est$PI
    ALPHABETA[[i+1]]<-EM_est$beta
    ERRORSIGMA[[i+1]]<-EM_est$sigma
    EIGENLAMDA[[i+1]]<-EM_est$lamda
    EIGENCOEFFI[[i+1]]<-EM_est$eigencoefficient
    delta<-EM_est$del
    sstep[iter]=i
    if(EM_est$Sstep==100) {break}
    if(Group_NUM[i+1]<true_num){break}
    
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
  
  #############the Bayes' optimal allocation rule
  categorial<-rep(NA,n)
  for (j in 1:n)
  {
    categorial[j]<-which.max(delta[j,])
  }
  cluster_result[iter,1:100]<-true_index########the true clutser id
  cluster_result[iter,101:200]<-categorial########the estimated clutser id
  num_group[iter]=Group_NUM[i+1] ########the estimated number of clusters
  
}

