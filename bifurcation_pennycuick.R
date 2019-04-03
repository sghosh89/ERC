# This function evaluates 1st order derrivative at equilibrium for pennycuick model
# Args
#     params : a vector containing model parameters
#     model : a character with model name (default = "pennycuick")
fn_d1_eqm<-function(params,model="pennycuick"){
  fn<-expression((r*p)/(1+exp(-a*(1-(p/b)))))
  fn_d1<-D(fn,"p")
  r<-params[1]
  a<-params[2]
  b<-params[3]
  p_e<-(b/a)*(a+log(r-1))
  res<-eval(fn_d1,list(r=r,a=a,b=b,p=p_e))
  return(res)
}
