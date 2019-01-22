# Bifurcation diagram for 1variable model
#                    Ricker model : p(t+1)=p(t)*exp(r(1-(p(t)/K))) # rc=1
#                    Hassell model : p(t+1)=r*p(t)/(1+(a*p(t)))^b  # rc=(b/(b-1))^b
#                    Maynard-smith model : p(t+1)=(r*p(t))/(1+((a*p(t))^b)) # rc=b/(b-1)
#                    Pennycuick model : p(t+1)=(r*p(t))/(1+exp(-a*(1-(p(t)/b)))) # rc calculate numerically


fn_d1_eqm<-function(params,model){
  if(model=="ricker"){
    fn<-expression(p*exp(r*(1-(p/K))))
    fn_d1<-D(fn,"p")
    r<-params[1]
    K<-params[2]
    p_e<-K
    res<-eval(fn_d1,list(r=r,K=K,p=p_e))
  }
  
  if(model=="hassell"){
    fn<-expression(r*p/(1+(a*p))^b)
    fn_d1<-D(fn,"p")
    r<-params[1]
    a<-params[2]
    b<-params[3]
    p_e<-((r^(1/b))-1)/a
    res<-eval(fn_d1,list(r=r,a=a,b=b,p=p_e))
  }
  
  if(model=="msmith"){
    fn<-expression((r*p)/(1+((a*p)^b)))
    fn_d1<-D(fn,"p")
    r<-params[1]
    a<-params[2]
    b<-params[3]
    p_e<-((r-1)^(1/b))/a
    res<-eval(fn_d1,list(r=r,a=a,b=b,p=p_e))
  }
  
  if(model=="pennycuick"){
    fn<-expression((r*p)/(1+exp(-a*(1-(p/b)))))
    fn_d1<-D(fn,"p")
    r<-params[1]
    a<-params[2]
    b<-params[3]
    p_e<-(b/a)*(a+log(r-1))
    res<-eval(fn_d1,list(r=r,a=a,b=b,p=p_e))
  }
  
  return(res)
}

#----------------------------------------------------------------
# calling for ricker
#rlist<-sort(c(0,seq(from=0.1,to=1.9,by=0.2),1))
#K<-50 
#res_ricker<-c() 
#for(r in rlist){
#  ans<-fn_d1_eqm(params=c(r,K),model="ricker")
#  res_ricker<-c(res_ricker,ans)
#}

#plot(rlist,res_ricker,col="red",type="b",xlab="r",ylab="ricker_1st_der_at_eqm")
#abline(h=0)

#formatC(res_ricker[7],digits = 60, format="f") #formatting upto 60th decimal : eqm point
#---------------------------------------------------------------------------------------------------
# calling for hassell
#rlist<-c(seq(from=1.5,to=2.5,by=0.5),seq(from=2.731,to=2.732,by=0.00000001),seq(from=3,to=7.5,by=0.5))
#a<-0.5
#b<-100
#res_hassell<-c() 
#for(r in rlist){
#  ans<-fn_d1_eqm(params=c(r,a,b),model="hassell")
#  res_hassell<-c(res_hassell,ans)
#}

#plot(rlist,res_hassell,col="red",type="l",xlab="r",ylab="hassell_1st_der_at_eqm")
#abline(h=0)
#rc_numerics<-rlist[which.min(abs(res_hassell-0))]
#formatC(rc_numerics,digits = 10, format="f") # compare with rc_analytical

#rc_analytical<-(b/(b-1))^b
#formatC(rc_analytical,digits = 10, format="f")
#--------------------------------------------------------------------------------------------------------

# calling for msmith 
#a<-1
#b<-3
#rc_analytical<-(b/(b-1))
#rlist<-sort(c(seq(from=1.2,to=2.8,by=0.2),rc_analytical))

#res_msmith<-c() 
#for(r in rlist){
#  ans<-fn_d1_eqm(params=c(r,a,b),model="msmith")
#  res_msmith<-c(res_msmith,ans)
#}

#plot(rlist,res_msmith,col="red",type="b",xlab="r",ylab="msmith_1st_der_at_eqm")
#abline(h=0)
#---------------------------------------------------------------------------------------
# calling for pennycuick

#a<-0.1
#b<-0.1
#eps<-1e-10
#rlist<-seq(from=4.3231657,to=4.3231658,by=eps)
#res_pennycuick<-c() 
#for(r in rlist){
#  ans<-fn_d1_eqm(params=c(r,a,b),model="pennycuick")
#  res_pennycuick<-c(res_pennycuick,ans)
#}

#plot(rlist,res_pennycuick,col="red",type="l",xlab="r",ylab="pennycuick_1st_der_at_eqm")
#abline(h=0)

#rc_numerics<-rlist[which.min(abs(res_pennycuick-0))]
#formatC(rc_numerics,digits = 10, format="f") # compare with rc_analytical







