

#' @examples
#' ss=Gamma_Positive_G(N=1000, alpha=20, beta=1, gamma=1)
#' hist(ss$sample,prob=TRUE, col="grey")
#' lines(density(ss$sample))
#' @export
Gamma_Positive_G<-function(N,alpha,beta,gamma){ # N: sample size
  delta_h=beta+(2*gamma^2-sqrt(4*gamma^4+32*alpha*beta*gamma^2))/(8*alpha)

  t <- rgamma(N,shape=alpha,rate=delta_h)
  U <- runif(N,0,1)
  Rejection_Rate <- exp(-(beta-delta_h)+gamma*sqrt(t)+gamma^2/(4*(beta-delta_h)))
  sample_t <- t[U<Rejection_Rate]  # in this step we can't gurantee the sample size =N
  sample_x=sqrt(sample_t)


  size=N-length(sample_x)
  i=1
  while(size!=0){
    t <- rgamma(size,shape=alpha,rate=delta_h)
    U <- runif(size,0,1)
    Rejection_Rate <- exp(-(beta-delta_h)+gamma*sqrt(t)+gamma^2/(4*(beta-delta_h)))
    sample_t <- t[U<Rejection_Rate]
    sample_x=c(sample_x,sqrt(sample_t))

    size=N-length(sample_x)

    i=i+1
    print(i)
  }

  sample=as.vector(sample_x)
  return(sample)
}





#' @examples
#' ss=Gamma_Positive_N(N=1000, alpha=20, beta=1, gamma=1)
#' hist(ss$sample,prob=TRUE, col="grey")
#' lines(density(ss$sample))
#' @export
Gamma_Positive_N<-function(N,alpha,beta,gamma){ # N: sample size
  m=(gamma+sqrt(gamma^2+8*beta*(alpha-1)))/(4*beta);m # mode

  x <- rnorm(N,mean=m, sd=1/(2*beta))
  U <- runif(N,0,1)
  Rejection_Rate <- (x/m)^(alpha-1)*exp((2*beta*m-gamma)*(m-x))
  sample <- x[U<Rejection_Rate]  # in this step we can't gurantee the sample size =N


  size=N-length(sample)
  i=1
  while(size!=0){
    x <- rnorm(size,mean=m, sd=1/(2*beta))
    U <- runif(size,0,1)
    Rejection_Rate <- (x/m)^(alpha-1)*exp((2*beta*m-gamma)*(m-x))
    sample0 <- x[U<Rejection_Rate]
    sample=c(sample,sample0)
    size=N-length(sample)

    i=i+1
    # print(i)
  }
  sample=as.vector(sample)
  return(sample)
}





Acceptence_Rate_G<-function(alpha,beta,gamma){

  delta_h=beta+(2*gamma^2-sqrt(4*gamma^4+32*alpha*beta*gamma^2))/(8*alpha)

  set.seed(100)
  t=runif(50000,0,1000)

  fx=t^(alpha/2-1)*exp(-beta*t+gamma*sqrt(t))
  Int_fx=mean(fx)

  gx=t^(alpha/2-1)*exp(-delta_h*t)
  Int_gx=mean(gx)


  Acceptance_Rate=Int_fx/(exp(gamma^2/(4*(beta-delta_h)))*Int_gx)

  return(Acceptance_Rate)
}




Acceptence_Rate_N<-function(alpha,beta,gamma){
  m=(gamma+sqrt(gamma^2+8*beta*(alpha-1)))/(4*beta);m # mode
  set.seed(100)
  x=runif(10000,0,100)
  gx=exp(-beta*(x-m)^2)
  Int_gx=mean(gx)

  fx=x^(alpha-1)*exp(-beta*x^2+gamma*x)
  Int_fx=mean(fx)

  Acceptance_Rate=Int_fx/(m^(alpha-1)*exp(-beta*m^2+gamma*m)*Int_gx)

  return(Acceptance_Rate)
}






#' @examples
#' ss=rMHN(N=1000, alpha=3, beta=1, gamma=1)
#' ss=rMHN(N=1000, alpha=3, beta=1, gamma=-1)
#' hist(ss$sample,prob=TRUE, col="grey")
#' lines(density(ss$sample))
#' @export
rMHN_old<-function(N,alpha,beta,gamma){
  if(gamma<=0){
    sample=replicate(N,rExtendedGamma(alpha,a,b=abs(gamma)))
    AC_Rate=Acceptence_Rate_Neg(alpha,beta,gamma)
  }

  else if(gamma>0){

    Acceptence_Rate_G=Acceptence_Rate_G(alpha,beta,gamma)
    Acceptence_Rate_N=Acceptence_Rate_N(alpha,beta,gamma)

    if(Acceptence_Rate_G>Acceptence_Rate_N){
      sample= Gamma_Positive_G(N,alpha,beta,gamma)
      AC_Rate= Acceptence_Rate_G }
    else if((Acceptence_Rate_G<=Acceptence_Rate_N)) {
      sample=Gamma_Positive_N(N,alpha,beta,gamma)
      AC_Rate=Acceptence_Rate_N }

  }

  return(list(sample=sample,Acceptence_Rate=AC_Rate))
}

#ss=rMHN(N=10000,alpha=32,beta=20,gamma=-3)
#ss$sample
#ss$Acceptence_Rate
