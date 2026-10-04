## Description: Simulation codes for eigenvalue ratio estimators and variance ratio tests.
## The code can produce the results given in Tables 1 and 2 of the paper.

## Loading packages
library(fda);library(tseries);library(sandwich);library(sde);library(variables);library(basefun);library(polynom)

## Loading precalculated critical values of the variance-ratio test
setwd("Path to the Working Directory")
inner = dget("inprod.R")
operator = dget("operator.R")
set.seed(99999)

pred_operator=function(basis_temp,eigen,TTTT,f,grid){
  ttemp=0
  for (tindex in 1:TTTT){
  temp = 0 
  for (k in 1:length(eigen)){
    temp = temp + (1/eigen[k])*inner(basis_temp[,k],f,grid)*inner(basis_temp[,k],x_mat[,tindex],grid)*y_mat[,tindex]
  }
  ttemp=ttemp+temp
  }
  return(ttemp/TTTT)
}

BBridge2=function(x0=0,y=0,t00=0,T0=1,N0=nt-1)
{
vvtem=BM(x=x0,t0=t00,T=T0,N=N0) * 1
vv=vvtem*(1/sqrt(3))
vv
}

FAR = 0; reg_iden= 0  # Notes: FAR=0, reg_iden= 1 = > Model 1;	FAR=0, reg_iden= 0 = > Model 2;	FAR=1, reg_iden= 1 = > Model 3;  FAR=1, reg_iden= 0 = > Model 4    

THRE=c(0.01,0.02)
EXPO=c(0.475,0.450,0.425)

####################################################
## Other key simulation paramters and related settings ###
####################################################
subTT=c(50,100,200,400,800); T=max(subTT)

lbnumber2=100
nt = 200 ; t = (0:(nt-1))/(nt-1); 

LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}
LBF=cbind(rep(1,length(t)),LBF)
lb=LBF
for(i in 2:lbnumber2){  
  for(j in 1:i)  { 
    if (j != i) {lb[,i] = lb[,i]-(inner(lb[,i],lb[,j],t)/inner(lb[,j],lb[,j],t))*lb[,j]  }}}

for(i in 1:lbnumber2){
  LBF[,i] = lb[,i]/(sqrt(inner(lb[,i],lb[,i],t)))
} 

nnbasis=100
inidiscard=100
kpoint = t
basis_fn = create.bspline.basis(rangeval = c(0,1), nbasis = nnbasis )

EMSE1=NULL; EMSE2=NULL; EMSE3=NULL ; KKKK1=NULL;  KKKK2=NULL;  KKKK3=NULL 
EMSE4=NULL; EMSE5=NULL; EMSE6=NULL ; KKKK4=NULL;  KKKK5=NULL;  KKKK6=NULL

#start_time <- Sys.time()
for (iteration in 1:1000)
{
#######################################################
####### Data generation ###############################
#######################################################
eigenfactors=(1:(lbnumber2+1))^(-1)
eta = matrix(rep(NA,((T+inidiscard)*nt)),ncol=T+inidiscard) 
eta[,1] = BBridge2(x=0,y=0,t0=0,T=1,N=nt-1)
eta2 = eta; eta3 = eta; x_mat = eta;  y_mat = eta;
if(FAR==1){dcfmax=0.75}
if(FAR==0){dcfmax=0.25}
dcoefficients1=runif(lbnumber2+1,-0.1,dcfmax)
dcoefficients2=runif(lbnumber2+1,-2.5,2.5)
dcoefficients2iden=runif(1,-2.5,2.5)



  for (jiter in 2:((T+inidiscard))){
 	eta2[,jiter] = sqrt(1)*BBridge2(x=0,y=0,t0=0,T=1,N=nt-1)
	eta3[,jiter] = sqrt(1)*BBridge2(x=0,y=0,t0=0,T=1,N=nt-1)
    {x_mat[,jiter] = eta2[,jiter]  + operator(LBF[,1:length(dcoefficients1)],dcoefficients1, x_mat[,jiter-1] ,t)} #else{x_mat[,jiter] = eta2[,jiter] } 
	  if (reg_iden==1){y_mat[,jiter] = dcoefficients2iden*x_mat[,jiter] +eta3[,jiter] }else{y_mat[,jiter] = x_mat[,jiter] +eta3[,jiter] + operator(LBF[,1:length(dcoefficients2)],dcoefficients2-1, x_mat[,jiter] ,t)}  	
    }

x_mat=x_mat[,(inidiscard+1):(T+inidiscard)]
y_mat=y_mat[,(inidiscard+1):(T+inidiscard)]


emse1=NULL;emse2=NULL;emse3=NULL;emse4=NULL;emse5=NULL;emse6=NULL
KKK1=NULL;KKK2=NULL;KKK3=NULL;KKK4=NULL;KKK5=NULL;KKK6=NULL;
for (subT in subTT)
{
f_mat = Data2fd(x_mat[,1:subT], argvals = kpoint, basis_fn)
pca_f=pca.fd(f_mat,nharm=trunc(subT^{1/2}),centerfns=FALSE)
fy_mat = Data2fd(y_mat[,1:subT], argvals = kpoint, basis_fn)
pca_fy=pca.fd(fy_mat,nharm=trunc(subT^{1/2}),centerfns=FALSE)
y_evector =eval.fd(pca_fy$harmonics,kpoint)
yd_mat=t(t(y_mat)%*%(y_evector)%*%solve(t(y_evector)%*%(y_evector))%*%t(y_evector))

KKtem=NULL; KKKtem=NULL
for (thre in THRE)
{
KKtem=NULL
for (expo in EXPO)
{
KK=min(max(which(pca_f$values[1:(nnbasis-1)]-pca_f$values[2:(nnbasis)] > thre*sum(pca_f$values)*subT^(-expo))),0.5*subT^(expo))
KK=trunc(KK)
KKtem=append(KKtem,KK)
}
KKKtem=rbind(KKKtem,KKtem)
}
maxKKK=max(KKKtem)
KKtem=KKKtem[1,]
KKtemA=KKKtem[2,]

KKK1=append(KKK1,KKtem[1])
KKK2=append(KKK2,KKtem[2])
KKK3=append(KKK3,KKtem[3])
KKK4=append(KKK4,KKtemA[1])
KKK5=append(KKK5,KKtemA[2])
KKK6=append(KKK6,KKtemA[3])


evalue=pca_f$values[1:maxKKK]
evector =eval.fd(pca_f$harmonics,kpoint)
evector = evector[,1:maxKKK]
gaps=NULL
preds=NULL
etas=NULL

preds=NULL

uval=sort(unique(as.vector(KKKtem)))

  gaps=array(NA,dim = c(nt,subT, length(uval)))
  
  for (ttt in 1:subT){
  ttemp=matrix(0,nrow=nt,ncol=length(uval))
  for (tindex in 1:subT){
  temp = NULL 
  {
  for (k in 1:length(evalue)){
    temp = cbind(temp, (1/evalue[k])*inner(evector[,k],x_mat[,ttt],kpoint)*inner(evector[,k],x_mat[,tindex],kpoint)*yd_mat[,tindex])
  }
  for (jjin in 1:length(uval)){ttemp[,jjin]=ttemp[,jjin]+(rowSums(temp[,1:uval[jjin]]))}
  }
  }
  ggg=(y_mat[,ttt] - ttemp/subT)
  for (jjin in 1:length(uval)){gaps[,ttt,jjin]=ggg[,jjin]}
  }
 
 gg=NULL 
 for (jjin in 1:length(uval))
 {
 gaps1_mat= Data2fd(gaps[,,jjin], argvals = kpoint, basis_fn)
 gg1=sum(pca.fd(gaps1_mat,nharm=trunc(T^{1/2}),centerfns=FALSE)$values)
 gg=append(gg,gg1)
 }
 
Kvec=c(KKKtem[1,],KKKtem[2,])
emseval=rep(0,length=6)
for(jjin in 1:length(uval))
{
emseval[which(Kvec==uval[jjin])]=gg[jjin]
}
if (subT==subTT[1]){EMSE1=rbind(EMSE1,emseval)}
if (subT==subTT[2]){EMSE2=rbind(EMSE2,emseval)}
if (subT==subTT[3]){EMSE3=rbind(EMSE3,emseval)}
if (subT==subTT[4]){EMSE4=rbind(EMSE4,emseval)}
if (subT==subTT[5]){EMSE5=rbind(EMSE5,emseval)}
}

set1=c(mean(EMSE1[,1]),mean(EMSE2[,1]),mean(EMSE3[,1]),mean(EMSE4[,1]),mean(EMSE5[,1]))
set2=c(mean(EMSE1[,2]),mean(EMSE2[,2]),mean(EMSE3[,2]),mean(EMSE4[,2]),mean(EMSE5[,2]))
set3=c(mean(EMSE1[,3]),mean(EMSE2[,3]),mean(EMSE3[,3]),mean(EMSE4[,3]),mean(EMSE5[,3]))
set4=c(mean(EMSE1[,4]),mean(EMSE2[,4]),mean(EMSE3[,4]),mean(EMSE4[,4]),mean(EMSE5[,4]))
set5=c(mean(EMSE1[,5]),mean(EMSE2[,5]),mean(EMSE3[,5]),mean(EMSE4[,5]),mean(EMSE5[,5]))
set6=c(mean(EMSE1[,6]),mean(EMSE2[,6]),mean(EMSE3[,6]),mean(EMSE4[,6]),mean(EMSE5[,6]))

if(iteration%%6==1){print(round(c(iteration,FAR,reg_iden,set1),digits=3))}
if(iteration%%6==2){print(round(c(iteration,FAR,reg_iden,set2),digits=3))}
if(iteration%%6==3){print(round(c(iteration,FAR,reg_iden,set3),digits=3))}
if(iteration%%6==4){print(round(c(iteration,FAR,reg_iden,set4),digits=3))}
if(iteration%%6==5){print(round(c(iteration,FAR,reg_iden,set5),digits=3))}
if(iteration%%6==0){print(round(c(iteration,FAR,reg_iden,set6),digits=3))}
}

result1=c(mean(EMSE1[,1]),mean(EMSE2[,1]),mean(EMSE3[,1]),mean(EMSE4[,1]),mean(EMSE5[,1]))
result1 












