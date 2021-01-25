


Psi_1<-function(alpha, x){
  val=0*seq(1:alpha)
  #browser()

  if(alpha==1){
    val= 2*sqrt(pi)*exp(x^2/4)*(pnorm(x/sqrt(2)))
  }
  if(alpha==2){
    val[1]= 2*sqrt(pi)*exp(x^2/4)*(pnorm(x/sqrt(2)))
    val[2]= x*val[1]/2+1
  }

  if(alpha>2){
    val[1]= 2*sqrt(pi)*exp(x^2/4)*(pnorm(x/sqrt(2)))
    val[2]= x*val[1]/2+1
    for(i in 3:alpha){
      val[i]=x*val[i-1]/2+ (i-2)*val[i-2]/2
    }

  }
  return(val[alpha])

}


Psi_ratio_a_plus1_by_a<-function(alpha, x){

  term_0=exp(x^2/4)*(   2*sqrt(pi)*pnorm(x/sqrt(2))   )
  psi_1=x/2+ exp(-x^2/4)/(   2*sqrt(pi)*(1-pnorm(-x/sqrt(2)))   )
  #term=0*(1:(alpha-1))
  term=psi_1
  if(alpha>1){
      for(i in 2: as.integer(alpha)){
        term=x/2+.5*(i-1)/term
      }
  }


  return(term)

}



FoxWright_MC<-function(alpha, x){


}




