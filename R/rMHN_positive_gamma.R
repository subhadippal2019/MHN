
#' @export
rExtendedGamma_positive_b_alpha<-function(alpha,a,b,ifPlot=F){
  #alpha>, a>0, b>0
  if(ifPlot){
    log_f<-function(xx){  return(   log(xx)*(alpha-1)-a*xx^2+b*xx  )};
    plot(log_f,0,100,xlab="log x")
  }

  # test the theory of the function

  if(1-(a>0)*(b>0)*(alpha>0)){  print("The parameters a or b and alpha have to be positive."); return(NULL); }


  fraction_part_alpha=alpha-floor(alpha);

  if(fraction_part_alpha==0){ return( rExtendedGamma_positive_b_integer_alpha(alpha,a,b)) }

  if(fraction_part_alpha>0){
    rt=b^2/(4*a);
    Integer_part_alpha=floor(alpha); N_upper=Integer_part_alpha+1;
    log_P=generateP_alt(rt,N_upper,iflog = TRUE)
    K=sum(exp(log_P))

    Scaled_prob=exp(log_P-max(log_P));P=(Scaled_prob)/sum(Scaled_prob)
    ACCEPT1=F


      while(!ACCEPT1){

        #browser()
        if(K==Inf){bin=0}
        if(K!=Inf){bin<-rbinom(1,1,1/(alpha*K+1))}


        if(!bin){
          i=sample(1:N_upper,size=1,prob=P)
          x_candidate=1+sqrt(   rgamma(1,shape=i/2, rate=rt)   )
          if(log(runif(1))<= -(1-fraction_part_alpha)*log(x_candidate)){
            x=x_candidate; ACCEPT1=T
          }
        }

        if(bin){
          x_candidate=rbeta(1,alpha,1)
          if(log(runif(1))<= -rt*(x_candidate-1)^2){
            x=x_candidate;ACCEPT1=T
          }
        }
      }


      return(x*b/(2*a))
  }#End fractional_part>0
}



###########################################################################################################
###########################################################################################################
###########################################################################################################

rExtendedGamma_positive_b_integer_alpha<-function(alpha,a,b){
  #alpha= positive intger, a>0, b>0

  if( (alpha-floor(alpha))>0){ print("alpha must be integer to use this function. use the function rExtendedGamma_positive_b_alpha instad"); return(NULL);}
  if(1-(a>0)*(b>0)*(alpha>0)){  print("The parameters a or b and alpha have to be positive."); return(NULL); }


  rt=b^2/(4*a)
  log_P=generateP_alt(rt,alpha,iflog = TRUE)

  K=sum(exp(log_P))
  Scaled_prob=exp(log_P-max(log_P));P=(Scaled_prob)/sum(Scaled_prob)
  ACCEPT=F

     while(!ACCEPT){
      #browser()
      if(K==Inf){bin=0}
      if(K!=Inf){bin<-rbinom(1,1,1/(alpha*K+1))}


      if(!bin){
        i=sample(1:alpha,size=1,prob=P)
        x=1+sqrt(   rgamma(1,shape=i/2, rate=rt)   )
        ACCEPT=T

      }

      if(bin){
        x_candidate=rbeta(1,alpha,1)
        if(log(runif(1))<= -rt*(x_candidate-1)^2){
          x=x_candidate;
          ACCEPT=T
        }
      }
    }
  return(x*b/(2*a))
}








###########################################################################################################
###########################################################################################################
###########################################################################################################
####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
###########################################################################################################
################################## Additional required functions ##########################################
###########################################################################################################



generateP_alt<-function(rt,n,iflog=FALSE){
  #generateP_alt function can manage higher values of the arguments. try generateP_alt(100,1000)
  log_p=0
  for(i in 1: (n) ){
    log_p[i]=log(.5)+log_comb((n-1),(i-1))  +   lgamma(i/2) -(i/2)*log(rt)
  }

  if(iflog){ return(log_p)  }
  if(!iflog){return(exp(log_p))}
}

generateP<-function(rt,n){
  p=0
  for(i in 1: (n) ){
    p[i]=.5*comb((n-1),(i-1))   *   gamma(i/2) / (rt)^(i/2)
  }
  return(p)
}



#########################

comb = function(n, x) {
  return(factorial(n) / (factorial(x) * factorial(n-x)))
}

log_comb = function(n, x) {
  return(lgamma(n+1) -(lgamma(x+1)+ lgamma(n-x+1)))
}


f<-function(x){  f=x^(alpha-1)*exp(-(a*x^2+b*x));  return(f);}

log_f<-function(x){   log_f=(alpha-1)*log(x) -(a*x^2+b*x);}





###########Test Functions
#
# n=100;a=.0001;b=.1
# samlpe_g= GenerateG(alpha,a,-b)
# samlpe_g=scaled_f_sampler(n,a,b)
# # for the case b>0, the sampler is extremely robust. try GenerateG(1000000,100000000,1000000000). to s the result.
# #GenerateG(1000000,100000000,100)
# #GenerateG(1000000,1,1000)
#
# for(iii in 1:1000){
#
#   samlpe_g[iii]=GenerateG(alpha,a,-b)
# }
#
# par(mfrow=c(1,2))
#
# plot(density(samlpe_g))
# plot(function(x){x^(n-1)*exp(-a*x^2-b*x)},xlim=c(300,700))
#



#
#
#
#
# log_p_i<-function(i,rt,n){
#   log(.5)+log_comb((n-1),(i-1))  +   lgamma(i/2) -(i/2)*log(rt)
# }
# generateP_alt1<-function(rt,n,iflog=FALSE){
#   #generateP_alt function can manage higher values of the arguments. try generateP_alt(100,1000)
#   #log_p=0
#
#   seq_i=1:n
#
#   log_p=apply(as.array(seq_i),1,log_p_i,rt=rt,n=n)
#
#   #for(i in 1: (n) ){
#   # log_p[i]=log(.5)+log_comb((n-1),(i-1))  +   lgamma(i/2) -(i/2)*log(rt)
#   #}
#
#   if(iflog){ return(log_p)  }
#   if(!iflog){return(exp(log_p))}
# }
#
#
#
#


