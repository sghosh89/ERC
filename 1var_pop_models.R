# Exploring Ricker model
Ricker<-function(r,K,p0,lensim){
  pop<-c(p0)
  time<-c(0)
  # po<-rep(K,r)
  # print(po)
  
  for(it in c(1:lensim)){
    pt<-p0*exp(r*(1-(p0/K)))
    time<-c(time,it)
    pop<-c(pop,pt)
    p0<-pt
  }
  
  plot(time,pop,type="b")
  
}

#call the function
Ricker(r=1.1,K=50,p0=5,lensim=1000) # bifurcation point r=1
#-------------------------------------------------------

# Exploring Hassell model
# help: https://jmahaffy.sdsu.edu/courses/s00a/math121/lectures/qual_discrete/qualdiscrete.html#Hassell
#       https://www.jstor.org/stable/pdf/3863.pdf?refreqid=excelsior%3A348907abb1bca7f039109d94db5667b5
#       http://lab.rockefeller.edu/cohenje/PDFs/231CohenUnexpectedDominanceHighFreqChaoticNonlinPpnModelsLet.pdf

Hassell<-function(r,a,b,p0,lensim){
  pop<-c(p0)
  time<-c(0)
  
  for(it in c(1:lensim)){
    pt<-(r*p0)/((1+(a*p0))^b)
    time<-c(time,it)
    pop<-c(pop,pt)
    p0<-pt
  }
  
  r1<-((b/(b-1)))^b
  r2<-((b/(b-2)))^b
  plot(time,pop,type="b")
  return(pop)
  
  
  #provided b>2 (to get finite r2), these are the stability conditions:
  cat("For stable, monotonic approach towards eqm. upper limit of r should be",r1,"\n")

  cat("For stable, oscillatory approach towards eqm. upper limit of r should be within ",r1,"to",r2,"\n")
  
}

#call the function
r<-2.72
a<-0.5
b<-100
K_e<-((r^(1/b))-1)/a
K_e
x<-Hassell(r=r, a=a, b=b, p0=K_e-0.01, lensim=50) #b~1.5 : contest, b~100 : scramble
diff(x)
Hassell(r=5, a=0.01, b=2.5, p0=0.1, lensim=200) # see change : r=5,a=0.5 and vary b 1.5 to 5.5
#----------------------------------------------------------------------------------------------------
# Exploring Maynard Smith model
Msmith<-function(r,a,b,p0,lensim){
  pop<-c(p0)
  time<-c(0)
  
  for(it in c(1:lensim)){
    pt<-(r*p0)/(1+((a*p0)^b))
    time<-c(time,it)
    pop<-c(pop,pt)
    p0<-pt
  }
  
  plot(time,pop,type="b")
  
}

#call the function
r<-1.2
a<-0.5
b<-4
K_e<-((r-1)^(1/b))/a
K_e
Msmith(r=r,a=a,b=b,p0=K_e-0.1,lensim=200) #see changes as r=5,a=0.5 and vary b : 1,1.5,2,4
                                          #see changes as a=0.5,b=4 and vary r : 1.2, 1.5, 2
#--------------------------------------------------------------------------------------------------

# Exploring Verhulst model
Verhulst<-function(r,p0,lensim){
  pop<-c(p0)
  time<-c(0)
  # po<-rep(K,r)
  # print(po)
  
  for(it in c(1:lensim)){
    pt<-p0*(1+r*(1-p0))
    time<-c(time,it)
    pop<-c(pop,pt)
    p0<-pt
  }
  
  plot(time,pop,type="b")
}

#call the function
Verhulst(r=1.4,p0=0.4,lensim=400) # eqm. point =1

#-------------------------------------------------------------
# Exploring pennycuick model
Pennycuick<-function(r,a,b,p0,lensim){
  pop<-c(p0)
  time<-c(0)
  
  for(it in c(1:lensim)){
    pt<-(r*p0)/(1+exp(-a*(1-(p0/b))))
    time<-c(time,it)
    pop<-c(pop,pt)
    p0<-pt
  }
  
  plot(time,pop,type="b")
  
}

#call the function
a<-0.1
b<-0.1
r<-2
K_e<-(b/a)*(a+log(r-1))
K_e
Pennycuick(r=r,a=a,b=b,p0=K_e-0.01,lensim=200) 

#-----------------------------------------
# Exploring Malthus model
Malthus<-function(r,K,L,p0,lensim){
  pop<-c(p0)
  time<-c(0)
  
  for(it in c(1:lensim)){
    pt<-r*p0*(K-p0+(L*log(p0)))
    time<-c(time,it)
    pop<-c(pop,pt)
    p0<-pt
  }
  
  plot(time,pop,type="b")
  
}

#call the function
r<-0.0032
K<-600
L<-4

Malthus(r=r,K=K,L=L,p0=0.1,lensim=200) 

#-----------------------------------------
# Exploring Austin-Brewer model
AustinBrewer<-function(r,K,s,p0,lensim){
  pop<-c(p0)
  time<-c(0)
  
  for(it in c(1:lensim)){
    pt<-p0*(1+(r*(K-p0)*(1-exp(-s*p0)))) 
    time<-c(time,it)
    pop<-c(pop,pt)
    p0<-pt
  }
  
  plot(time,pop,type="b")
  
}

#call the function
s<-0.13
K<-100
rc<-1/(K*(1-exp(-s*K)))
rc
r<-0.014
AustinBrewer(r=r,K=K,s=s,p0=1,lensim=200) 

#----------------------------------------------------
# Exploring Varley model

Varley<-function(r,b,C,p0,lensim){
  pop<-c(p0)
  time<-c(0)
  
  for(it in c(1:lensim)){
    if(p0<=C){
      pt<-r*p0
    }else{
      pt<-r*(p0^(1-b))
    }
    time<-c(time,it)
    pop<-c(pop,pt)
    p0<-pt
  }
  
  plot(time,pop,type="b")
  
}

#call the function
b<-1.4
r<-4
C<-1
p0<-r^(1/b)
p0
Varley(r=r,b=b,C=C,p0=p0-1,lensim=200) 


































