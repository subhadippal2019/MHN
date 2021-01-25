

####################################################
# FoxWright Method 1
####################################################

#' @export
FoxWright1<-function(alpha,beta,gamma,eps=.00001,log=FALSE){
  u=alpha;v=gamma;w=beta;
  j=0
  x0=0
  x1=lgamma(u/2+j/2)+j*log(abs(v)/sqrt(w))-lgamma(1)-lfactorial(j)
  while(abs(x1-x0)>eps && exp(x1)!=0){
    j=j+1
    x0=x1
    x1=lgamma(u/2+j/2)+j*log(abs(v)/sqrt(w))-lgamma(1)-lfactorial(j)
  }

  SumTimes=j
  ######################

  i=seq(0,SumTimes,by=1)
  if(v>0){
    logitem=lgamma(u/2+i/2)+i*log(v/sqrt(w))-lgamma(1)-lfactorial(i)
    Max_of_terms=max(logitem)
    logitem_Minus_Maxlog=logitem-Max_of_terms
    item_star=exp(logitem_Minus_Maxlog)
    fox=exp(Max_of_terms)*sum(item_star)
    logfox=log(fox)}
  else if (v<0) {
    logitem=lgamma(u/2+i/2)+i*log(abs(v)/sqrt(w))-lgamma(1)-lfactorial(i)
    Max_of_terms=max(logitem)
    logitem_Minus_Maxlog=logitem-Max_of_terms
    item_star=exp(logitem_Minus_Maxlog)
    fox=exp(Max_of_terms)*sum(item_star*(-1)^i)
    logfox=log(fox)
  } else {
    fox=0
    logfox=-Inf
  }

  result=ifelse(log==FALSE,fox,logfox)
  return(result)

}

####################################################
# FoxWright Method 2
####################################################

#' @export
FoxWright2<-function(alpha,beta,gamma,eps=.00001,log=FALSE){
    u=alpha;v=gamma;w=beta;

  RobustFox<-function(u,v,w,log,SumTimes){
    i=seq(0,SumTimes,by=1)
    if(v>0){
      logitem=lgamma(u/2+i/2)+i*log(v/sqrt(w))-lgamma(1)-lfactorial(i)
      item=exp(logitem)
      fox=sum(item)
      logfox=log(fox)}
    else if (v<0) {
      logitem=lgamma(u/2+i/2)+i*log(abs(v)/sqrt(w))-lgamma(1)-lfactorial(i)
      item=exp(logitem)*(-1)^i
      fox=sum(item)
      logfox=log(fox)
    } else {
      fox=0
      logfox=-Inf
    }

    result=ifelse(log==F,fox,logfox)
    return(result)
  }


  # eps: converge level (difference of two numbers)
  ## find the SumTimes which makes RobustFox function converge
  i=1
  x0=0
  x1=RobustFox(u,v,w,SumTimes=i,log)
  while(abs(x1-x0)>eps){
    i=i+1
    x0=x1
    x1=RobustFox(u,v,w,SumTimes=i,log)
  }

  return(FoxWright=x1)
}

####################################################
# FoxWright Method 3
####################################################

#install.packages('Rmpfr')
#' @export
FoxWright3<-function(alpha,beta,gamma,eps=.00001,log=FALSE){
  u=alpha;v=gamma;w=beta;
  j=1
  x0=0
  x1=lgamma(u/2+j/2)+j*log(abs(v)/sqrt(w))-lgamma(1)-lfactorial(j)
  while(abs(x1-x0)>eps && round(exp(x1),2)!=0){
    j=j+1
    x0=x1
    x1=lgamma(u/2+j/2)+j*log(abs(v)/sqrt(w))-lgamma(1)-lfactorial(j)
  }

  SumTimes=j
  ########################################
  i=seq(0,SumTimes,by=1)
  if(v>0){
    logitem=lgamma(u/2+i/2)+i*log(v/sqrt(w))-lgamma(1)-lfactorial(i)

    logitem <- mpfr(logitem, precBits = 106)
    item=exp(logitem)
    fox=sum(item)
    logfox=log(fox)}
  else if (v<0) {
    logitem=lgamma(u/2+i/2)+i*log(abs(v)/sqrt(w))-lgamma(1)-lfactorial(i)
    logitem <- mpfr(logitem, precBits = 106)
    item=exp(logitem)*(-1)^i
    fox=sum(item)   # this is a list type
    logfox=log(fox)
  } else {
    fox=0
    logfox=-Inf
  }


  fox <- capture.output(fox)[2]
  fox <- substr(fox,5,nchar(fox))

  logfox <- capture.output(logfox)[2]
  logfox <- substr(logfox,5,nchar(logfox))

  result=ifelse(log==FALSE,fox,logfox)
  return(noquote(result))
}

####################################################
# FoxWright Method 4
####################################################

#' @export
FoxWright<-function(alpha,beta,gamma,eps=.00001,log=FALSE){
  u=alpha;v=gamma;w=beta;
  ###### part 1 #######
  k=0
  x1=lgamma(u/2+k)-lgamma(u/2)+k*log(v^2)-k*log(w)-lfactorial(2*k)
  x0=x1+1
  while(exp(lgamma(u/2))*exp(x1)!=0 |
        abs(exp(x1)-exp(x0))>0 |
        k<20){
    k=k+1
    x0=x1
    x1=lgamma(u/2+k)-lgamma(u/2)+k*log(v^2)-k*log(w)-lfactorial(2*k)
  }
  T1=k
  ###### part 2 #######
  k=0
  x1=lgamma((u+1)/2+k)-lgamma((u+1)/2)+k*log(v^2)-k*log(w)-lfactorial(2*k+1)
  x0=x1+1
  while(exp(lgamma((u+1)/2))*(v/sqrt(w))*exp(x1)!=0 |
        abs(exp(x1)-exp(x0))>0 |
        k<20){
    k=k+1
    x0=x1
    x1=lgamma((u+1)/2+k)-lgamma((u+1)/2)+k*log(v^2)-k*log(w)-lfactorial(2*k+1)
  }
  T2=k
  ###### SUM #######
  SumTimes=max(T1,T2)
  #  print(paste("SumTimes=",SumTimes))

  #k=seq(0,SumTimes,by=1)
  if(SumTimes<180){k=seq(0,SumTimes,by=1)} else {k=seq(0,SumTimes+300,by=1)}
  #  print(paste("k=",k[length(k)]))

  logM1=lgamma(u/2+k)-lgamma(u/2)+k*log(v^2)-k*log(w)-lfactorial(2*k)

  logM2=lgamma((u+1)/2+k)-lgamma((u+1)/2)+k*log(v^2)-k*log(w)-lfactorial(2*k+1)

  fox=exp(lgamma(u/2))*sum(exp(logM1))+exp(lgamma((u+1)/2))*(v/sqrt(w))*sum(exp(logM2))

  #fox=sum(c(exp(lgamma(u/2))*exp(logM1),exp(lgamma((u+1)/2))*(v/sqrt(w))*exp(logM2)))
  result=ifelse(log==F,fox,log(fox))
  return(result)
}


