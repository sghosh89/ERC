#--------------------------------------------------------------------------
source("ExtremeTailDep.R")
#--------------------------------------------------------------------------
#The simulator of the multi-location model.
#
#Args
#p0       A length numlocs vector holding initial populations
#ns       An numsims by numlocs by numsteps array of the epsilons, 
#           where numsteps is the number of time steps you want
#
#Output
#A numsims by numlocs by numsteps+1 array of populations
#
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
      lam_verhulst<-1+(r*(1-res[,,tct]))*res[,,tct]
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
      lam_malthus<-r*res[,,tct]*(K-res[,,tct]+(L*log(res[,,tct])))
      res[,,tct+1]<-lam_malthus*lam_sto #numsims by numlocs matrix
    }else if(model=="austinbrewer"){
      r<-params[1]
      K<-params[2]
      s<-params[3]
      lam_austin<-(1+(r*(K-res[,,tct])*(1-exp(-s*res[,,tct]))))*res[,,tct] 
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
extrisk<-function(sims){
  totpop<-apply(FUN=sum,X=sims,MARGIN=c(1,3))
  return(apply(FUN=sum,X=(totpop==0),MARGIN=2)/dim(sims)[1])
}

# function to calculate average ACF (lag=1) from all numlocs after specified time steps
# input : 
#        1) sims : an array (numsims by numlocs by numsteps+1)
#        2) nts : numsteps after which you want the acf; generally equal to numsteps+1 as to compare with extinction risk
get_avg_acf<-function(sims,nts){
  m<-sims[,,nts]
  macf<-apply(X=m,MARGIN = 2,FUN=acf,lag.max = 1,type="correlation",plot=F)
  acfs<-c() # initialize to store acf with lag1 for each col of matrix m
  for(i in c(1:ncol(m))){
    acfs<-c(acfs,macf[[i]]$acf[2])
  }
  ans<-mean(acfs)
  return(ans)
}

#get_avg_acf(sims=sims,nts=numsteps+1)


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
#       9) model : a character specifying model name
#       10) ploton : logical(T or F) to get optional plot
#       11) resloc : location to save plots

plotter_ext_risk<-function(numsims,numsteps,numlocs,D,p0,params,ext_thrs,scl,model,ploton,resloc){
 
  ns1<-retd(n=numsteps*numsims,d=numlocs,rl=1)# a righttail dep matrix(numpoints by numlocs,     
  #                                                      numpoints=numsteps*numsims)
  
  ns1<-array(ns1,c(numsteps,numsims,numlocs))# convert to an array (numsteps by numsims by numlocs)
  ns1<-aperm(ns1,c(2,3,1)) # convert to an array (numsims by numlocs by numsteps)
  
  ns1<-scl*ns1
  pops1<-popsim_ml_D(p0=rep(p0,numlocs),ns=ns1,D=D,params=params,ext_thrs=ext_thrs,model=model)
  risk_right<-extrisk(pops1) # a vector
  acf_right<-get_avg_acf(sims=pops1,nts=numsteps+1) # a number
  
  ns2<-(-ns1)
  pops2<-popsim_ml_D(p0=rep(p0,numlocs),ns=ns2,D=D,params=params,ext_thrs=ext_thrs,model=model)
  risk_left<-extrisk(pops2) # a vector
  acf_left<-get_avg_acf(sims=pops2,nts=numsteps+1) # a number
  
  if(ploton==T){
    # noise time series plot for last simulation (numsim=numsims) from each patches
    pdf(paste0(resloc,"params_",params,"_noise_lastsim_timeseries_all_locs.pdf",sep=""),height=5,width=10)
    op<-par(mfrow=c(1,2),mar=c(3.5,4.5,2,3.5),mgp=c(1.9,0.5,0))
    
    ns_sim1<-ns1[numsims,,]
    ylm<-max(c(abs(min(ns_sim1)),abs(max(ns_sim1))))
    plot(0,0,xlim = c(1,numsteps),ylim = c(-ylm,ylm),type = "n",xlab="time",ylab="noise_right_tail")
    cl <- rainbow(numlocs)
    for(i in 1:numlocs){
      lines(c(1:numsteps),ns_sim1[i,],col=cl[i],type="l")
    }
    
    ns_sim1<-ns2[numsims,,]
    ylm<-max(c(abs(min(ns_sim1)),abs(max(ns_sim1))))
    plot(0,0,xlim = c(1,numsteps),ylim = c(-ylm,ylm),type = "n",xlab="time",ylab="noise_left_tail")
    cl <- rainbow(numlocs)
    for(i in 1:numlocs){
      lines(c(1:numsteps),ns_sim1[i,],col=cl[i],type="l")
    }
    
    par(op)
    dev.off()
    
    pdf(paste0(resloc,"params_",params,"_pops_lastsim_timeseries_all_locs.pdf",sep=""),height=5,width=10)
    op<-par(mfrow=c(1,2),mar=c(3.5,4.5,2,3.5),mgp=c(1.9,0.5,0))
    
    pop_sim1<-pops1[numsims,,]
    ylm<-max(c(abs(min(pop_sim1)),abs(max(pop_sim1))))
    plot(0,0,xlim = c(1,numsteps+1),ylim = c(-ylm,ylm),type = "n",xlab="time",ylab="pops_with_noise_right_tail")
    cl <- rainbow(numlocs)
    for(i in 1:numlocs){
      lines(c(1:c(numsteps+1)),pop_sim1[i,],col=cl[i],type="l")
    }
    
    pop_sim2<-pops2[numsims,,]
    ylm<-max(c(abs(min(pop_sim2)),abs(max(pop_sim2))))
    plot(0,0,xlim = c(1,numsteps+1),ylim = c(-ylm,ylm),type = "n",xlab="time",ylab="pops_with_noise_left_tail")
    cl <- rainbow(numlocs)
    for(i in 1:numlocs){
      lines(c(1:c(numsteps+1)),pop_sim2[i,],col=cl[i],type="l")
    }
    
    par(op)
    dev.off()
    
    pdf(paste0(resloc,"params_",params,"_extrisk_vs_time.pdf",sep=""),height=5,width=10)
    op<-par(mfrow=c(1,2),mar=c(3.5,4.5,2,3.5),mgp=c(1.9,0.5,0))
    plot(1:length(risk_right),risk_right,type='b',
         xlab='Time',ylab='Risk',col='blue',panel.first = grid())
    mtext("Right tail-dep")
    
    plot(1:length(risk_left),risk_left,type='b',xlab='Time',ylab='Risk',col='red',panel.first = grid())
    mtext("Left tail-dep")
    
    par(op)
    dev.off()  
    
  }
  
  return(list(risk_right_after_numsteps=risk_right[numsteps+1],
              risk_left_after_numsteps=risk_left[numsteps+1],
              avg_acf_right_after_numsteps=acf_right,
              avg_acf_left_after_numsteps=acf_left))
  
}
#----------------------------------------------------------------------------------------------------

# function to plot extinction risk agsinst degree of dispersal (d) after a certain numstep (time)

# Input :
#       1) numsims : an integer : number of simulations
#       2) numsteps : an integer : number of time steps
#       3) numlocs : an integer : number of locations
#       4) p0 : initial population to start with
#       5) params : a vector with model parameters
#       6) scl : a scaling factor which is multiplied with noise generated
#       7) ext_thrs extinction threshold below which populations go extinct
#       8) model : a character specifying model name
#       9) disp_everywhere :logical:
#             if T : gives D for linear chain model with equal dispersal everywhere
#             if F : gives D for linear chain model with equal dispersal only to nearest neighbor location
#       10) plotteron : logical to get optional plot

varying_d<-function(numsims,numsteps,numlocs,p0,params,ext_thrs=ext_thrs,scl,model,disp_everywhere,plotteron,resloc){
  risk_right<-c()
  risk_left<-c()
  avg_acf_right<-c()
  avg_acf_left<-c()
  
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
                              ext_thrs=ext_thrs,scl=scl,model=model,ploton=T,resloc=resloc2)
    
    risk_r<-riskrl$risk_right_after_numsteps
    risk_l<-riskrl$risk_left_after_numsteps
    avg_acf_r<-riskrl$avg_acf_right_after_numsteps
    avg_acf_l<-riskrl$avg_acf_left_after_numsteps
    risk_right<-c(risk_right,risk_r)
    risk_left<-c(risk_left,risk_l)
    avg_acf_right<-c(avg_acf_right,avg_acf_r)
    avg_acf_left<-c(avg_acf_left,avg_acf_l)
  }
  
  if(plotteron==T){
    op<-par(mfrow=c(2,1), mar=c(5.2,4.2,1,1.2))
    plot(d_seq,risk_left,xlim=c(0,1),ylim=c(0,1),type="b",col="red",panel.first = grid(),
         xlab="d",ylab="ext_risk",cex.lab=1.5,cex.axis=1.5)
    lines(d_seq,risk_right,type="b",col="blue")
    legend("bottomleft", c("noise : left","noise : right"), lty=c(1,1), pch=c(1,1), 
           col=c('red', 'blue'), horiz=T, bty='n', cex=1.2)
    mtext(paste0("r = ", r))
    
    plot(d_seq,avg_acf_left,xlim=c(0,1),ylim=c(-0.1,0.1),type="b",col="red",panel.first = grid(),pch=16,
         xlab="d",ylab="avg. ACF",cex.lab=1.5,cex.axis=1.5)
    lines(d_seq,avg_acf_right,type="b",col="blue",pch=16)
    legend("bottomleft", c("noise : left","noise : right"), lty=c(1,1), pch=c(16,16), 
           col=c('red', 'blue'), horiz=T, bty='n', cex=1.2)
    #mtext(paste0("r = ", r))
    
    par(op)
  }
  
  return(data.frame(d_seq=d_seq,
                    risk_left=risk_left,
                    risk_right=risk_right,
                    avg_acf_left=avg_acf_left,
                    avg_acf_right=avg_acf_right))
  

}
#------------------------------------------------------------------------------------------------------------
