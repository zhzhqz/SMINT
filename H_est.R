############################the function for estimating H  given Omega_n. It uses functions defined in the files original function.R and first_der function.R .

#Arguments:
#YY: the observed functional data
#PI_MATRIX: the most updated mixing probabilities
#ALPHABETA: the most updated basis coefficient of mean functions
#gamma_int: the most updated variances of errors
#EIGENLAMDA: the most updated eigenvalues
#Alpha: the most updated basis coefficient of eigenfunctions
#Group_NUM: the most updated number of clusters
#Hini: the most updated H
#iter_shu: the number of iteration

#return: the estimates of H(y) 

H_est<-function(YY,PI_MATRIX,ALPHABETA,gamma_int,EIGENLAMDA,Alpha,Group_NUM,Hini,iter_shu)
{
  HY_MATRIX<-rep(NA,geshu)
  for(j in 1:geshu){
 initial=Hini[j]
    before=0
    after=1
    s=0
    s=sum(unlist(YY)<point[j])

    StepH<-0
############Start a loop to find the estimates
     while((StepH<iter_shu)&(abs(after-before)>1e-5)){
StepH<-StepH+1
      orig<-original(initial,PI_MATRIX,ALPHABETA,gamma_int,EIGENLAMDA,Alpha,Group_NUM)
      after=initial+(s-orig)/(first_der(initial,PI_MATRIX,ALPHABETA,gamma_int,EIGENLAMDA,Alpha,Group_NUM))
      before=initial
      initial=after
    }
    HY_MATRIX[j]<-after
    if(j>1){if(HY_MATRIX[j]<HY_MATRIX[j-1]){HY_MATRIX[j]=HY_MATRIX[j-1]}}
    if(j %in% seq(10,geshu,by=10)){cat("j=",j,"\n")}
  }
  
  coe_up<-solve(t(DESIMATRX_up)%*%DESIMATRX_up)%*%t(DESIMATRX_up)%*%HY_MATRIX[(geshu-9):geshu]
  coe_low<-solve(t(DESIMATRX_low)%*%DESIMATRX_low)%*%t(DESIMATRX_low)%*%HY_MATRIX[1:5]
  tranfun_up<-function(t) c(1,t)%*%coe_up
  tranfun_low<-function(t) c(1,t)%*%coe_low
  
  jjk1=sapply(textpoint1,tranfun_up)
  jjk2=sapply(textpoint2,tranfun_low)
  transfunc<-c(jjk2,HY_MATRIX,jjk1)

  return(transfunc)
}