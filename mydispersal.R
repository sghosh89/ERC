#--------------------------------------------------------------------------
source("ExtremeTailDep.R")
#--------------------------------------------------------------------------
library(copula)
#-------------------------------------------------------------------------------
# A function to generate dispersal matrix D (a square matrix : numlocs by numlocs)
# Input :
#       1)numlocs : number of location : an integer
#       2) d : degree of dispersal : a number in [0,1]
#       3) disp_everywhere :logical:
#             if T : gives D for linear chain model with equal dispersal everywhere
#             if F : gives D for linear chain model with equal dispersal only to nearest neighbor location
# Output : dispersal matrix D (a square matrix : numlocs by numlocs)

Dispersal_mat<-function(numlocs,d,disp_everywhere){
  D<-matrix(0,numlocs,numlocs)
  # A companion matrix (delta) that indicates how "off" a diagonal is:
  mat <- matrix(seq_len(numlocs ^ 2), ncol = numlocs)
  delta <- row(mat) - col(mat)
  
  diag(D)<-(1-d)
  
  if(disp_everywhere==T){
    D[delta!=0]<-d/(numlocs-1)
  }else{
    D[delta ==1 | delta==-1] <- d/2
    D[numlocs-1,numlocs]<-d
    D[2,1]<-d
  } 
  
  stopifnot(all(colSums(D)==1)) # this is a check
  return(D)
}
#-----------------------------------------------------------------------------------
#The simulator of the multi-location model with a dispersal matrix.
#-------Args--------------------------------------------------
#p0       A length numlocs vector holding initial populations
#ns       An numsims by numlocs by numsteps array of the epsilons, 
#           where numsteps is the number of time steps you want
#D        Dispersal matrix
#params   a vector with model parameters in order as written in equation (e.g. r,a,b... or r,K,L)
#ext_thrs  a threshold below which populations go extinct
#model    a character specifying model name
#--------Output------------------------------------------------------
#A numsims by numlocs by numsteps+1 array of populations
#----------------------------------------------------------------------
popsim_ml_D<-function(p0,ns,D,params,ext_thrs,model){
  numsims<-dim(ns)[1]
  numlocs<-dim(ns)[2]
  numsteps<-dim(ns)[3]
  
  #initialization
  res<-array(NA,c(numsims,numlocs,numsteps+1))
  res[,,1]<-rep(p0,each=numsims)
  
  #step the populations forward
  for (tct in 1:numsteps){
    
    #stochastic ricker model : if r=0 then back to LC model
    
    #growth rates based on the noise
    lam_sto<-exp(ns[,,tct])  #numsims by numlocs matrix
    
    #growth prior to dispersal
    if(model=="ricker"){
      r<-params[1]
      K<-params[2]
      lam_ricker<-exp(r*(1-(res[,,tct]/K)))*res[,,tct]
      res[,,tct+1]<-lam_ricker*lam_sto #numsims by numlocs matrix
    }else if(model=="verhulst"){
      r<-params[1]
      K<-params[2]
      lam_verhulst<-(1+(r*(1-(res[,,tct]/K))))*res[,,tct]
      res[,,tct+1]<-lam_verhulst*lam_sto #numsims by numlocs matrix
    }else if(model=="hassell"){
      r<-params[1]
      a<-params[2]
      b<-params[3]
      lam_hassell<-(r*res[,,tct])/((1+(a*res[,,tct]))^b)
      res[,,tct+1]<-lam_hassell*lam_sto #numsims by numlocs matrix
    }else if(model=="msmith"){
      r<-params[1]
      a<-params[2]
      b<-params[3]
      lam_msmith<-(r*res[,,tct])/(1+((a*res[,,tct])^b))
      res[,,tct+1]<-lam_msmith*lam_sto #numsims by numlocs matrix
    }else if(model=="pennycuick"){
      r<-params[1]
      a<-params[2]
      b<-params[3]
      lam_pennycuick<-(r*res[,,tct])/(1+exp(-a*(1-(res[,,tct]/b))))
      res[,,tct+1]<-lam_pennycuick*lam_sto #numsims by numlocs matrix
    }else if(model=="malthus"){
      r<-params[1]
      K<-params[2]
      L<-params[3]
      
      x1<-res[,,tct]
      id0<-which(x1==0,arr.ind = T)
      idp<-which(x1!=0,arr.ind = T)
      
      x2<-x1
      
      x2[idp]<-r*x1[idp]*(K-x1[idp]+(L*log(x1[idp])))
      x2[id0]<-0 # just to ensure extinct value
      
      #lam_malthus<-r*res[,,tct]*(K-res[,,tct]+(L*log(res[,,tct])))
      lam_malthus<-x2
      
      res[,,tct+1]<-lam_malthus*lam_sto #numsims by numlocs matrix
      
    }else if(model=="austinbrewer"){
      r<-params[1]
      K<-params[2]
      s<-params[3]
      lam_austin<-(1+(r*(K-res[,,tct])*(1-exp(-(s*res[,,tct])))))*res[,,tct]
      res[,,tct+1]<-lam_austin*lam_sto #numsims by numlocs matrix
    }else if(model=="varley"){
      r<-params[1]
      b<-params[2]
      C<-params[3]
      x1<-res[,,tct]
      x2<-C*matrix(1,nrow(x1),ncol(x1))
      idl<-which(x1<=x2,arr.ind = T)
      x1[idl]<-r*x1[idl]
      idg<-which(x1>x2,arr.ind = T)
      x1[idg]<-r*(x1[idg])^(1-b)
      res[,,tct+1]<-x1*lam_sto #numsims by numlocs matrix
    }else{ 
      warning("model not specified",immediate.=T,call.=T)
    }
    
    #after dispersal
    tres<-t(res[,,tct+1]) # numlocs by numsims matrix after taking transpose
    temp<-D%*%tres  # output from a matrix multiplication [(numlocs by numlocs) * (numlocs by numsims)]
    res[,,tct+1]<-t(temp) # a matrix (numsims by numlocs) 
    
    #see if there is extinction
    h<-res[,,tct+1]
    h[h<ext_thrs]<-0
    res[,,tct+1]<-h
  }
  
  return(res)
}

# function to calculate extinction risk
#
#Args
#sims     A numeric array, assumed to be numsims by numlocs by numsteps
#
#Output
#Produces a vector of length numsteps with the entry in location i being
#the fraction of simulations which have gone extinct by time step i. 
#"Extinct" here means all locations have population 0.
#
extrisk<-function(sims){
  totpop<-apply(FUN=sum,X=sims,MARGIN=c(1,3))
  return(apply(FUN=sum,X=(totpop==0),MARGIN=2)/dim(sims)[1])
}

# function to calculate average ACF (lag=1) from all numlocs after specified time steps
# input : 
#        1) sims : an array (numsims by numlocs by numsteps)

get_avg_acf<-function(sims){
  #x<-sims[,,-1] # to get rid off the initial p0 values
  x<-sims # dim(x) = numsims x numlocs x numsteps
  acfs<-c() # to store acf from timeseries of each location for each simulation
  
  for(i in c(1:dim(x)[1])){ # loop for each numsim
    m<-x[i,,] # matrix of dim = numlocs by numsteps
    macf<-apply(X=m,MARGIN = 1,FUN=acf,lag.max = 1,type="correlation",plot=F) # along row : time series from each location
    s<-unlist(macf)
    id<-which(names(s)=="acf2")
    sl1<-as.numeric(s[id])
    acfs<-c(acfs,sl1)
  }
  # which(is.na(acfs))
  ans<-mean(acfs,na.rm=T) # This line to remove NaN which is produced due to pop = 0 throughout the time at any patch
  # you can see that in the row of matrix m (for a given i)
  return(ans)
}
# I checked : ACF of noise (right tail) = ACF of noise (left tail)

#-------------------------------------------------------------------------------------------------
# function to optionally plot extinction risk agsinst time and give you back extinction risk 
#                                                            after specified numsteps 
# Input :
#       1) numsims : an integer : number of simulations
#       2) numsteps : an integer : number of time steps
#       3) numlocs : an integer : number of locations
#       4) D : dispersal matrix which is the output of Dispersal_mat function
#       5) p0 : initial pop value to start with
#       6) params : a vector containing model parameters
#       7) ext_thrs : extinction threshold below which populations go extinct
#       8) scl : a scaling factor which is multiplied with noise generated
#       9) srho : given Spearman's rho to generate moderate tail dep. copula
#       10) model : a character specifying model name
#       11) ploton : logical(T or F) to get optional plot
#       12) resloc : location to save plots
#       13) getacf : these are tags

plotter_ext_risk<-function(numsims,numsteps,numlocs,D,p0,params,ext_thrs,scl,srho,model,ploton,resloc,getacf){
  
  # ========================== when noises are extreme tail dep. ========================================
  ns1<-retd(n=numsteps*numsims,d=numlocs,rl=1,mn=0,sdev=1)# a righttail dep matrix(numpoints by numlocs,     
  # numpoints=numsteps*numsims)
  
  ns1<-array(ns1,c(numsteps,numsims,numlocs))# convert to an array (numsteps by numsims by numlocs)
  ns1<-aperm(ns1,c(2,3,1)) # convert to an array (numsims by numlocs by numsteps)
  
  ns1<-scl*ns1
  pops1<-popsim_ml_D(p0=rep(p0,numlocs),ns=ns1,D=D,params=params,ext_thrs=ext_thrs,model=model)
  risk_xright<-extrisk(pops1) # a vector
  
  ns2<-(-ns1)
  pops2<-popsim_ml_D(p0=rep(p0,numlocs),ns=ns2,D=D,params=params,ext_thrs=ext_thrs,model=model)
  risk_xleft<-extrisk(pops2) # a vector
  
  #------------------------------------------------------
  # This two number must be same and I checked that
  #noise_acf_xright<-get_avg_acf(sims=ns1) # a number
  #noise_acf_xleft<-get_avg_acf(sims=ns2) # a number
  #--------------------------------------------------------
  if(getacf==T){
    pop_acf_xright<-get_avg_acf(sims=pops1) # a number
    pop_acf_xleft<-get_avg_acf(sims=pops2) # a number
  }else{
    pop_acf_xright<-NA # a number
    pop_acf_xleft<-NA
  }
  
  # ========================== when noises have moderate tail dep. ========================================
  
  # for clayton : left tail dep.
  copC<-claytonCopula(3)
  parC<-iRho(copC,rho=srho)
  
  # first generate Clayton copula with dimensions=numlocs
  cc<-claytonCopula(par=parC,dim=numlocs)
  ns3<-rCopula(numsteps*numsims,cc)  # numpoints=(numsteps*numsims)
  
  # flip Clayton copula with dimensions=numlocs
  ns4<-1-ns3
  
  # make marginal normally distributed
  ns3<-qnorm(ns3)
  ns4<-qnorm(ns4)
  
  ns3<-array(ns3,c(numsteps,numsims,numlocs))# convert to an array (numsteps by numsims by numlocs)
  ns3<-aperm(ns3,c(2,3,1)) # convert to an array (numsims by numlocs by numsteps)
  
  ns3<-scl*ns3
  pops3<-popsim_ml_D(p0=rep(p0,numlocs),ns=ns3,D=D,params=params,ext_thrs=ext_thrs,model=model)
  risk_mleft<-extrisk(pops3) # a vector
  
  ns4<-array(ns4,c(numsteps,numsims,numlocs))# convert to an array (numsteps by numsims by numlocs)
  ns4<-aperm(ns4,c(2,3,1)) # convert to an array (numsims by numlocs by numsteps)
  
  ns4<-scl*ns4
  pops4<-popsim_ml_D(p0=rep(p0,numlocs),ns=ns4,D=D,params=params,ext_thrs=ext_thrs,model=model)
  risk_mright<-extrisk(pops4) # a vector
  
  # =========================================================================================================
  if(ploton==T){
    
    # noise time series plot for last simulation (numsim=numsims) from each patches
    pdf(paste0(resloc,"r_",params[1],"_noise_lastsim_timeseries_all_locs.pdf",sep=""),height=10,width=10)
    op<-par(mfrow=c(2,2),mar=c(3.5,4.5,2,3.5),mgp=c(1.9,0.5,0))
    
    ns_sim1<-ns2[numsims,,]
    ylm<-max(c(abs(min(ns_sim1)),abs(max(ns_sim1))))
    plot(0,0,xlim = c(0,numsteps),ylim = c(-ylm,ylm),type = "n",xlab="Time",ylab="Noise : extreme LTdep.")
    cl <- rainbow(numlocs)
    for(i in 1:numlocs){
      lines(c(0:numsteps),c(0,ns_sim1[i,]),col=cl[i],type="l") # at t=0 noise was also 0, for t>0 , there was finite noise in each time step
    }
    
    ns_sim1<-ns1[numsims,,]
    ylm<-max(c(abs(min(ns_sim1)),abs(max(ns_sim1))))
    plot(0,0,xlim = c(0,numsteps),ylim = c(-ylm,ylm),type = "n",xlab="Time",ylab="Noise : extreme RTdep.")
    cl <- rainbow(numlocs)
    for(i in 1:numlocs){
      lines(c(0:numsteps),c(0,ns_sim1[i,]),col=cl[i],type="l") # at t=0 noise was also 0, for t>0 , there was finite noise in each time step
    }
    
    ns_sim3<-ns3[numsims,,]
    ylm<-max(c(abs(min(ns_sim3)),abs(max(ns_sim3))))
    plot(0,0,xlim = c(0,numsteps),ylim = c(-ylm,ylm),type = "n",xlab="Time",ylab="Noise : moderate LTdep.")
    cl <- rainbow(numlocs)
    for(i in 1:numlocs){
      lines(c(0:numsteps),c(0,ns_sim3[i,]),col=cl[i],type="l") # at t=0 noise was also 0, for t>0 , there was finite noise in each time step
    }
    
    ns_sim4<-ns4[numsims,,]
    ylm<-max(c(abs(min(ns_sim4)),abs(max(ns_sim4))))
    plot(0,0,xlim = c(0,numsteps),ylim = c(-ylm,ylm),type = "n",xlab="Time",ylab="Noise : moderate RTdep.")
    cl <- rainbow(numlocs)
    for(i in 1:numlocs){
      lines(c(0:numsteps),c(0,ns_sim4[i,]),col=cl[i],type="l") # at t=0 noise was also 0, for t>0 , there was finite noise in each time step
    }
    
    par(op)
    dev.off()
    
    #--------------------------------------------------------------------------
    
    # population time series plot for last simulation (numsim=numsims) from each patches
    pdf(paste0(resloc,"r_",params[1],"_pops_lastsim_timeseries_all_locs.pdf",sep=""),height=10,width=10)
    op<-par(mfrow=c(2,2),mar=c(3.5,4.5,2,3.5),mgp=c(1.9,0.5,0))
    
    pop_sim2<-pops2[numsims,,]
    ylm<-max(c(abs(min(pop_sim2)),abs(max(pop_sim2))))
    plot(0,0,xlim = c(0,numsteps),ylim = c(0,ylm),type = "n",xlab="Time",ylab="Population for extreme LTdep. noise")
    cl <- rainbow(numlocs)
    for(i in 1:numlocs){
      lines(c(0:numsteps),pop_sim2[i,],col=cl[i],type="l")
    }
    
    pop_sim1<-pops1[numsims,,]
    ylm<-max(c(abs(min(pop_sim1)),abs(max(pop_sim1))))
    plot(0,0,xlim = c(0,numsteps),ylim = c(0,ylm),type = "n",xlab="Time",ylab="Population for extreme RTdep. noise")
    cl <- rainbow(numlocs)
    for(i in 1:numlocs){
      lines(c(0:numsteps),pop_sim1[i,],col=cl[i],type="l")
    }
    
    pop_sim3<-pops3[numsims,,]
    ylm<-max(c(abs(min(pop_sim3)),abs(max(pop_sim3))))
    plot(0,0,xlim = c(0,numsteps),ylim = c(0,ylm),type = "n",xlab="Time",ylab="Population for moderate LTdep. noise")
    cl <- rainbow(numlocs)
    for(i in 1:numlocs){
      lines(c(0:numsteps),pop_sim3[i,],col=cl[i],type="l")
    }
    
    pop_sim4<-pops4[numsims,,]
    ylm<-max(c(abs(min(pop_sim4)),abs(max(pop_sim4))))
    plot(0,0,xlim = c(0,numsteps),ylim = c(0,ylm),type = "n",xlab="Time",ylab="Population for moderate RTdep. noise")
    cl <- rainbow(numlocs)
    for(i in 1:numlocs){
      lines(c(0:numsteps),pop_sim4[i,],col=cl[i],type="l")
    }
    
    par(op)
    dev.off()
    
    #-------------------------------------------------
    
    # extinction risk time series plot for last simulation (numsim=numsims) from each patches
    pdf(paste0(resloc,"r_",params[1],"_extrisk_vs_time.pdf",sep=""),height=10,width=10)
    op<-par(mfrow=c(2,2),mar=c(3.5,4.5,2,3.5),mgp=c(1.9,0.5,0))
    
    plot(c(0:numsteps),risk_xleft,type='b',xlab='Time',ylab='Extinction risk',col='red',panel.first = grid())
    mtext("Noise : extreme LTdep.")
    
    plot(c(0:numsteps),risk_xright,type='b',
         xlab='Time',ylab='Extinction risk',col='blue',panel.first = grid())
    mtext("Noise : extreme RTdep.")
    
    plot(c(0:numsteps),risk_mleft,type='b',
         xlab='Time',ylab='Extinction risk',col='red',panel.first = grid())
    mtext("Noise : moderate LTdep.")
    
    plot(c(0:numsteps),risk_mright,type='b',
         xlab='Time',ylab='Extinction risk',col='blue',panel.first = grid())
    mtext("Noise : moderate RTdep.")
    
    par(op)
    dev.off()  
    
  }
  
  return(list(risk_xright_after_numsteps=risk_xright[numsteps+1],
              risk_xleft_after_numsteps=risk_xleft[numsteps+1],
              pop_avg_acf_xright_after_numsteps=pop_acf_xright,
              pop_avg_acf_xleft_after_numsteps=pop_acf_xleft,
              risk_mleft_after_numsteps=risk_mleft[numsteps+1],
              risk_mright_after_numsteps=risk_mright[numsteps+1]
  ))
  
}
#----------------------------------------------------------------------------------------------------

# function to plot extinction risk agsinst degree of dispersal (d) after a certain numstep (time)

# Input :
#       1) numsims : an integer : number of simulations
#       2) numsteps : an integer : number of time steps
#       3) numlocs : an integer : number of locations
#       4) p0 : initial population to start with
#       5) params : a vector with model parameters
#       6) ext_thrs extinction threshold below which populations go extinct
#       7) scl : a scaling factor which is multiplied with noise generated
#       8) srho : given Spearman's rho to generate moderate tail dep. copula
#       9) model : a character specifying model name
#       10) disp_everywhere :logical:
#             if T : gives D for linear chain model with equal dispersal everywhere
#             if F : gives D for linear chain model with equal dispersal only to nearest neighbor location
#       11) plotteron : logical to get optional plot
#       12) resloc : location to save plots
#       13) getacf : these are tags

varying_d<-function(numsims,numsteps,numlocs,p0,params,ext_thrs=ext_thrs,scl,srho,model,disp_everywhere,plotteron,resloc,getacf){
  risk_xright<-c()
  risk_xleft<-c()
  risk_mright<-c()
  risk_mleft<-c()
  pop_avg_acf_xright<-c()
  pop_avg_acf_xleft<-c()
  noise_avg_acf_xright<-c()
  noise_avg_acf_xleft<-c()
  
  if(disp_everywhere==F){
    tempo2<-paste(resloc,"local_disp_d_",sep="")
  }else{
    tempo2<-paste(resloc,"global_disp_d_",sep="")
  }
  
  d_seq<-seq(from=0,to=1,by=0.1)
  
  for(d in d_seq){
    D_mat<-Dispersal_mat(numlocs=numlocs,d=d,disp_everywhere=disp_everywhere)
    
    tempo3<-paste(tempo2,d,sep="")
    
    if (!dir.exists(tempo3)){
      dir.create(tempo3)
    }
    resloc2<-paste(tempo3,"/",sep="")
    
    riskrl<- plotter_ext_risk(numsims=numsims,numsteps=numsteps,numlocs=numlocs,D=D_mat,p0=p0,params=params,
                              ext_thrs=ext_thrs,scl=scl,srho=srho,model=model,ploton=T,resloc=resloc2,getacf=getacf)
    
    risk_xr<-riskrl$risk_xright_after_numsteps
    risk_xl<-riskrl$risk_xleft_after_numsteps
    pop_avg_acf_xr<-riskrl$pop_avg_acf_xright_after_numsteps
    pop_avg_acf_xl<-riskrl$pop_avg_acf_xleft_after_numsteps
    #noise_avg_acf_r<-riskrl$noise_avg_acf_xright_after_numsteps
    #noise_avg_acf_l<-riskrl$noise_avg_acf_xleft_after_numsteps
    risk_mr<-riskrl$risk_mright_after_numsteps
    risk_ml<-riskrl$risk_mleft_after_numsteps
    
    risk_xright<-c(risk_xright,risk_xr)
    risk_xleft<-c(risk_xleft,risk_xl)
    pop_avg_acf_xright<-c(pop_avg_acf_xright,pop_avg_acf_xr)
    pop_avg_acf_xleft<-c(pop_avg_acf_xleft,pop_avg_acf_xl)
    #noise_avg_acf_xright<-c(noise_avg_acf_xright,noise_avg_acf_xr)
    #noise_avg_acf_xleft<-c(noise_avg_acf_xleft,noise_avg_acf_xl)
    
    risk_mright<-c(risk_mright,risk_mr)
    risk_mleft<-c(risk_mleft,risk_ml)
  }
  
  if(plotteron==T){
    
    if(getacf==T){
      op<-par(mfrow=c(2,1), mar=c(5.2,4.2,1,1.2))
      plot(d_seq,risk_xleft,xlim=c(0,1),ylim=c(0,1),type="b",col="red",panel.first = grid(),
           xlab="d",ylab="ext_risk",cex.lab=1.5,cex.axis=1.5)
      lines(d_seq,risk_xright,type="b",col="blue")
      legend("bottomleft", c("noise: extreme LT","noise: extreme RT"), lty=c(1,1), pch=c(1,1), 
             col=c('red', 'blue'), horiz=T, bty='n', cex=1.2)
      mtext(paste0("r = ", r))
      
      plot(d_seq,pop_avg_acf_xleft,xlim=c(0,1),ylim=c(-0.1,0.1),type="b",col="red",panel.first = grid(),pch=16,
           xlab="d",ylab="avg. ACF",cex.lab=1.5,cex.axis=1.5)
      lines(d_seq,pop_avg_acf_xright,type="b",col="blue",pch=16)
      legend("bottomleft", c("noise: extreme LT","noise: extreme RT"), lty=c(1,1), pch=c(16,16), 
             col=c('red', 'blue'), horiz=T, bty='n', cex=1.2)
      #mtext(paste0("r = ", r))
      
      par(op)
    }else{
      op<-par(mfrow=c(2,1), mar=c(5.2,4.2,1,1.2))
      plot(d_seq,risk_xleft,xlim=c(0,1),ylim=c(0,1),type="b",col="red",panel.first = grid(),
           xlab="d",ylab="ext_risk",cex.lab=1.5,cex.axis=1.5)
      lines(d_seq,risk_xright,type="b",col="blue")
      legend("bottomleft", c("noise : extreme LT","noise : extreme RT"), lty=c(1,1), pch=c(1,1), 
             col=c('red', 'blue'), horiz=T, bty='n', cex=1.2)
      mtext(paste0("r = ", r))
      
      plot(d_seq,risk_mleft,xlim=c(0,1),ylim=c(0,1),type="b",col="red",panel.first = grid(),
           xlab="d",ylab="ext_risk",cex.lab=1.5,cex.axis=1.5)
      lines(d_seq,risk_mright,type="b",col="blue")
      legend("bottomleft", c("noise : moderate LT","noise : moderate RT"), lty=c(1,1), pch=c(1,1), 
             col=c('red', 'blue'), horiz=T, bty='n', cex=1.2)
      mtext(paste0("r = ", r))
      
      par(op)
    }
  }
  
  return(data.frame(d_seq=d_seq,
                    risk_xleft=risk_xleft,
                    risk_xright=risk_xright,
                    pop_avg_acf_xleft=pop_avg_acf_xleft,
                    pop_avg_acf_xright=pop_avg_acf_xright,
                    risk_mleft=risk_mleft,
                    risk_mright=risk_mright
  ))
  
  
}
#------------------------------------------------------------------------------------------------------------
