# This part is new part of the analysis 
# We can use reads instead of hard calls
# For this to work we need to estimate genotype frequencies.

# let's use EM algorithm

# Data consists of matrix n by 3
# p = P(G=0)
# q = P(G=1)

EM_calc = function(M){
  p_0 = 0.15
  q_0 = 0.15
  
  q_n = 1
  p_n = 0
  d_n = 0
  k = 0
  while ((p_n - p_0)^2 + (q_n - q_0)^2>0.000000001){
    d_0 = 1-p_0 - q_0
    v = c(p_0,q_0,d_0)
    p_D = M%*%(v)
    E_p = M[,1]*p_0/p_D
    E_q = M[,2]*q_0/p_D
    p_n= p_0
    q_n = q_0
    d_n = 1-q_0-p_0
    p_0 = sum(E_p)/length(E_p)
    q_0 = sum(E_q)/length(E_q)
    k = k+1
    if (k==5000){
      cat('hi','\n')
      return (c(p_0,q_0,1-p_0-q_0))
    }
  }
  
  return (c(p_0,q_0,1-p_0-q_0))
  
}




EM_calc_HWE = function(M){
  
  p_t<-0.5
  
  p0_t<-(1-p_t)^2
  p1_t<-2*p_t*(1-p_t)
  p2_t<-p_t^2
  
  p_n<-0.7
  
  k<-0
  
  while ( (p_t-p_n)^2 >0.000000001 ){
    
    p_t<-p_n
    
    v<-c(p0_t,p1_t,p2_t)
    
    p_D<-M%*%v
    E_0<-M[,1]*p0_t/p_D
    E_1<-M[,2]*p1_t/p_D
    E_2<-M[,3]*p2_t/p_D
    
    # solve for new p
    eq_solve<-function(p){
      val<-sum(E_0)/(p-1) + sum(E_1)*(1-2*p)/(2*p*(1-p)) + sum(E_2)/p
      return(val)
    }
    
    p_n<-uniroot(eq_solve,c(0.000001,0.999999))$root
    
    # update p0, p1, p2
    p0_t<-(1-p_n)^2
    p1_t<-2*p_n*(1-p_n)
    p2_t<-p_n^2
    
    k<-k+1
    if (k==5000){
      cat('hi','\n')
      return (c(p_n))
    }
    
  }
  
  return (c(p_n))
  
}