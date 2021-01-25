
#' rExtendedGamma
#'
#' @export
rExtendedGamma<-function(alpha, a, b){
  # alpha>0, a>0 and b = any real

  if(1-(a>0)*(alpha>0)){  print("The parameters a and alpha required to be positive."); return(NULL) }


  if(b<0){
    return(rExtendedGamma_negative_b(alpha=alpha, a=a, b=abs(b)))
  }
  if(b>0){
   return( rExtendedGamma_positive_b_alpha(alpha=alpha, a=a,b =b))
  }

  if(b==0){
    shape_par=alpha/2; rate_par=a;
    Exp_g<-  rgamma(1,shape=shape_par,rate=rate_par)
    g=sqrt(Exp_g)
    return(g)
  }

}


