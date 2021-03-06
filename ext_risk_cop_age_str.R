library(copula)
library(VineCopula)
library(mvtnorm)

# This is the simulator for age structured 2 variable pop model with fecundity and survival drawn from a specific copula family

# Input :
#       p0 : a vector of two numeric values : E0, A0 initial values for Egg and Adult
#       ext_thrs : extinction threshold, generally = 1
#       cop : a bivariate copula matrix
#       par_dist : a vector containing 4 numbers for gamma and beta distribution parameters: gshape, gscale, bshape1, bshape2
#       omg : a number : if infinity then linear model, if finite nonlinear stage structured model
#       numsims : number of simulations
#       numsteps : number of time steps
#       ploton : logical to get optional plot
#       resloc : location path to save optional plot

# Output :
#       A list of three : resE, resA : each is a numsims by numsteps+1 matrix for Egg and Adult stage
#                         extdf : a dataframe (nrow= numsteps+1, ncol=3) : extinction risk for adult and its 95% CI

sim_age_str<-function(p0,ext_thrs,cop,par_dist,omg,numsims,numsteps,ploton,resloc){
  
  gshape<-par_dist[1]
  gscale<-par_dist[2]
  bshape1<-par_dist[3]
  bshape2<-par_dist[4]
  
  #initialization for Egg and Adult pop variables
  resE<-matrix(NA,nrow=numsims, ncol=numsteps+1)
  resA<-matrix(NA,nrow=numsims, ncol=numsteps+1)
  
  resE[,1]<-p0[1]
  resA[,1]<-p0[2]
  
  # make the transformation on cop
  
  # 1st column : fecundity of adult : gamma distributed
  x <- qgamma(cop[,1],shape=gshape,scale=gscale)
  
  # 2nd column : survival of eggs : beta distributed
  y <- qbeta(cop[,2],shape1=bshape1,shape2=bshape2) 
  
  transcop<-cbind(x,y)
  
  if(ploton==ploton){
    pdf(paste(resloc,"fA_sE_gshape_",gshape,"_gscale_",gscale,"_bshape1_",bshape1,"_bshape2_",bshape2,".pdf",sep=""),height=5,width=10)
    op<-par(mfrow=c(1,2))
    hist(x,100,main="",col="grey",border=F,xlab="Fecundity of adults")
    #hist(x,100,main="fecundity of adult : gamma dist.",col="grey",border=F,xlab="")
    #mtext(paste0("gshape = ",gshape," , gscale = ", gscale),line=-1) 
    
    hist(y,100,main="",col="grey",border=F,xlab="Survival of eggs")
    #mtext(paste0("bshape1 = ",bshape1," , bshape2 = ", bshape2),line=-1) 
    par(op)
    dev.off()
  }
  
  Ls<-split(as.data.frame(transcop), rep(1:numsteps, each = numsims)) 
  
  #step the populations forward
  for (tct in 1:numsteps){
    
    L<-as.matrix(Ls[[tct]]) # this is the Leslie matrix (dim = numsims by numsims)
    fA<-L[,1]
    sE<-L[,2]
    
    resE[,tct+1]<-fA * resA[,tct] * exp(-(resA[,tct]/omg))
    resA[,tct+1]<-sE * resE[,tct] 
    
    #see if there is extinction on total Adult and Egg population
    h<-resA[,tct+1]+resE[,tct+1]
    id_ext<-which(h<ext_thrs)
    resA[id_ext,tct+1]<-0
    resE[id_ext,tct+1]<-0
  }
  
  erA<-get_er(mat=resA)
  
  CI0.025<-qbinom(p=0.025,size=numsims,prob=erA)/numsims
  CI0.975<-qbinom(p=0.975,size=numsims,prob=erA)/numsims
  
  if(ploton==ploton){
    pdf(paste(resloc,"ext_riskA_vs_time_omg_",omg,".pdf",sep=""),height=5,width=5)
    op<-par(mar=c(4,4,2,2),mgp=c(2.5,0.5,0))
    plot(c(0:numsteps),erA,xlab="time",ylab="Extinction risk",cex.lab=1.5,cex=0.3,pch=16,col="black",ylim=c(0,1),panel.first = grid())
    arrows(x0=c(1:numsteps),y0=CI0.025[2:(numsteps+1)],x1=c(1:numsteps),y1=CI0.975[2:(numsteps+1)],
           angle=90,length=0.01,col="grey",code=3)
    par(op)
    dev.off()
  }
  
  extdf<-as.data.frame(cbind(CI0.025,erA,CI0.975))
  
  return(list(resE=resE,
              resA=resA,
              extdf=extdf))
}
#----------------------------------------------------------------------------
# function to calculate extinction risk
get_er<-function(mat){
  er<-(apply(FUN=sum,X=(mat==0),MARGIN=2))/nrow(mat)
  return(er)
}

#--------------------------------

#copC<-claytonCopula(par)
#tC<-tau(copC)
#tC

#rC  <- rotCopula(claytonCopula(2), flip = T)
#tau(rC)
#iTau(rC,tau=0.6) # This should be the same parameter as of clayton
#-------------------------------------------------------------------------

