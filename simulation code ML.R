# simulation code with ML
expit<-function(x){return(exp(x)/(1+exp(x)))}
logit<-function(x){return(log(x/(1-x)))}

library(truncnorm)

# use some of Derbach et al.'s code
source('em.r')

get_L = function(vect_row_read,A1,A2,error){
  if (vect_row_read == A1){
    p1 = 1-error
  }else{
    p1 = error/3
  }
  
  if (vect_row_read == A2){
    p2 = 1-error
  }else{
    p2 = error/3
  }
  p = 1/2*p1 + 1/2*p2
  
  return (p)
}

get_Mr = function(Error,vect_row_reads,rdv){
  L = length(rdv)
  M = NULL
  LG  =c('TT','CT','CC')
  S = 0
  for (i in 1:L){
    m = NULL
    for  (j in 1:3){
      LL = 1
      for (kk in 1:rdv[i]){
        G = LG[j]
        A1 = substring(G,1,1)
        A2 = substring(G,2,2)
        LL = LL*get_L(vect_row_reads[S+kk],A1,A2,Error[S+kk])
      }
      m =c(m,LL)
    }
    S =S+rdv[i]
    M = rbind(M,m)
  }
  return (M)
}

# function for ML simulation
full_MLE_optim<-function(n_case,n_control,rd_case_mean,rd_case_sd,rd_con_mean,rd_con_sd,p,e_case_mean,e_case_sd,e_con_mean,e_con_sd,beta0,beta1,R,alpha=0.05){
  
  # arguments:
  # n_case          number of cases
  # n_control       number of controls
  # rd_case_mean    mean of read depth for cases
  # rd_case_sd      sd of read depth for cases
  # rd_con_mean     mean of read depth for controls
  # rd_con_sd       sd of read depth for controls
  # p               minor allele frequency
  # e_case_mean     mean of error rate for cases
  # e_case_sd       sd of error rate for cases
  # e_con_mean      mean of error rate for controls
  # e_con_sd        sd of error rate for controls
  # beta0           beta0 of logistic model for phenotype
  # beta1           beta1 of logistic model for phenotype
  # R               number of replications
  # alpha           alpha level for number of rejection in function print out (p-values are still returned)
  
  gca<-cbind(rep(0,n_case),rep(1,n_case),rep(2,n_case))
  gcon<-cbind(rep(0,n_control),rep(1,n_control),rep(2,n_control))
  
  pos_g<-c("A","G","T","C")
  
  # get beta0_cc, the intercept within case-control first
  P_Y_1<-sum(expit(beta0+beta1*c(0,1,2))*c((1-p)^2,2*p*(1-p),p^2))
  beta0_cc<-beta0+log( n_case/n_control * (1-P_Y_1)/P_Y_1 )
  
  P_G_Y_0<-(1-expit(beta0+beta1*c(0,1,2)))*c((1-p)^2,2*p*(1-p),p^2)/(1-P_Y_1)
  P_G_Y_1<-expit(beta0+beta1*c(0,1,2))*c((1-p)^2,2*p*(1-p),p^2)/P_Y_1
  
  P_G_S_1<-P_G_Y_0*n_control/(n_case+n_control)+P_G_Y_1*n_case/(n_case+n_control)
  
  print(P_G_Y_0)
  print(P_G_Y_1)
  
  est<-NULL
  
  orig_est<-NULL
  
  EM_est<-NULL
  
  info_Total<-matrix(rep(0,9),3,3)
  
  cov_est_info_Total<-matrix(rep(0,4),2,2)
  
  score<-NULL
  
  log_L<-NULL
  
  bad_counter<-0
  
  p_val_true<-NULL
  p_val<-NULL
  
  # count number of rejections less than alpha
  rej_T<-0
  rej_info<-0
  
  plot(0,0,xlim=c(0,R),ylim=c(0,R))
  
  for(i in 1:R){
    if( (i %% 10) == 0 ){points(i,i)}
    
    a<-"B"
    
    while(typeof(a)=="character"){
      
      # generate genotype of controls
      g_control<-as.vector(t(rmultinom(n_control,1,P_G_Y_0))%*%c(0,1,2))
      
      # generate genotype of cases
      g_case<-as.vector(t(rmultinom(n_case,1,P_G_Y_1))%*%c(0,1,2))
      
      # get reads for controls
      P_D_G_control<-NULL
      for(j in 1:n_control){
        # get genotype
        if(g_control[j]==0){
          cur_g<-c("T","T")
        } else if(g_control[j]==1){
          cur_g<-c("C","T")
        } else{
          cur_g<-c("C","C")
        }
        
        # generate read depth
        cur_rd<-round(rtruncnorm(1,a=1,b=Inf,mean=rd_con_mean,sd=rd_con_sd))
        
        # generate errors
        cur_error<-rtruncnorm(cur_rd,a=0,b=1,mean=e_con_mean,sd=e_con_sd)
        
        # generate reads
        cur_read<-NULL
        for(k in 1:cur_rd){
          look_at<-sample(1:2,1)
          if(runif(1)>cur_error[k]){
            cur_read<-c(cur_read,cur_g[look_at])
          } else{
            pos_read<-pos_g[which(pos_g!=cur_g[look_at])]
            cur_read<-c(cur_read,pos_read[sample(1:3,1)])
          }
        }
        
        P_D_G_control<-rbind(P_D_G_control,as.vector(get_Mr(cur_error,cur_read,cur_rd)))
      }
      
      # get reads for cases
      P_D_G_case<-NULL
      for(j in 1:n_case){
        # get genotype
        if(g_case[j]==0){
          cur_g<-c("T","T")
        } else if(g_case[j]==1){
          cur_g<-c("C","T")
        } else{
          cur_g<-c("C","C")
        }
        
        # generate read depth
        cur_rd<-round(rtruncnorm(1,a=1,b=Inf,mean=rd_case_mean,sd=rd_case_sd))
        
        # generate errors
        cur_error<-rtruncnorm(cur_rd,a=0,b=1,mean=e_case_mean,sd=e_case_sd)
        
        # generate reads
        cur_read<-NULL
        for(k in 1:cur_rd){
          look_at<-sample(1:2,1)
          if(runif(1)>cur_error[k]){
            cur_read<-c(cur_read,cur_g[look_at])
          } else{
            pos_read<-pos_g[which(pos_g!=cur_g[look_at])]
            cur_read<-c(cur_read,pos_read[sample(1:3,1)])
          }
        }
        
        P_D_G_case<-rbind(P_D_G_case,as.vector(get_Mr(cur_error,cur_read,cur_rd)))
      }
      
      neg_log_Likelihood<-function(t){
        fc<-function(x) expit(x)*expit(t[2])+expit(x+t[1])*(expit(t[2]+exp(t[3]))-expit(t[2]))+expit(x+2*t[1])*(1-expit(t[2]+exp(t[3])))-n_case/(n_case+n_control)
        
        b0cc<-uniroot(f=fc,interval=c(-100,100))$root
        
        a<-c( F1 = (-1)*sum( log( as.vector( ( P_D_G_case*expit(b0cc+t[1]*gca) )%*%c( expit(t[2]) , expit(t[2]+exp(t[3]))-expit(t[2]) , 1-expit(t[2]+exp(t[3])) ) ) ) ) 
              - sum( log( as.vector( ( P_D_G_control*(1-expit(b0cc+t[1]*gcon)) )%*%c( expit(t[2]) , expit(t[2]+exp(t[3]))-expit(t[2]) , 1-expit(t[2]+exp(t[3])) ) ) ) ) )
        
        return(a)
      }
      
      a<-tryCatch({
        sqrt(solve(optim(c(0,-1,0),neg_log_Likelihood,hessian=T)$hessian)[1,1])
        m<-optim(c(0,-1,0),neg_log_Likelihood,hessian=T)$hessian
        ma<-m[1:2,1:2]
        mb<-matrix(m[1:2,3])
        mc<-t(matrix(m[3,1:2]))
        md<-matrix(m[3,3])
        bwi<-solve(ma) + solve(ma)%*%mb%*%solve(md-mc%*%solve(ma)%*%mb)%*%mc%*%solve(ma)
        sqrt(bwi[1,1])
      }, warning = function(w){
        return("B")
      }, error = function(e){
        return("B")
      }
      )
      
      if(typeof(a) == "character"){
        bad_counter<-bad_counter+1
      }
    }
    
    mod<-optim(c(0,-1,0),neg_log_Likelihood,hessian=T)
    
    moda<-m[1:2,1:2]
    modb<-matrix(m[1:2,3])
    modc<-t(matrix(m[3,1:2]))
    modd<-matrix(m[3,3])
    
    cov_est<-solve(moda) + solve(moda)%*%modb%*%solve(modd-modc%*%solve(moda)%*%modb)%*%modc%*%solve(moda)
    
    cur_est<-mod$par
    
    fb0cc<-function(x) expit(x)*expit(cur_est[2])+expit(x+cur_est[1])*(expit(cur_est[2]+exp(cur_est[3]))-expit(cur_est[2]))+expit(x+2*cur_est[1])*(1-expit(cur_est[2]+exp(cur_est[3])))-n_case/(n_case+n_control)
    
    root<-uniroot(f=fb0cc,interval=c(-100,100))$root
    
    est<-rbind(est,c(root,cur_est[1],expit(cur_est[2]),expit(cur_est[2]+exp(cur_est[3]))-expit(cur_est[2])))
    
    orig_est<-rbind(orig_est,mod$par)
    
    # fit true G regression
    y<-c(rep(1,n_case),rep(0,n_control))
    g_all<-c(g_case,g_control)
    reg<-glm(y~g_all,family=binomial)
    
    
    if(summary(reg)$coefficients[2,4]<alpha){
      rej_T<-rej_T+1
    }
    
    p_val_true<-c(p_val_true,summary(reg)$coefficients[2,4])
    
    if(abs(cur_est[1]/sqrt(cov_est[1,1]))>qnorm(1-alpha/2)){
      rej_info<-rej_info+1
    }
    
    cov_est_info_Total<-cov_est_info_Total+cov_est
    info_Total<-info_Total+mod$hessian
    
    #print(mod$hessian)
    
    p_val<-c(p_val,2*(1-pnorm(abs(cur_est[1]/sqrt(cov_est[1,1])))))
  }
  
  # Outputs
  print(c(beta0_cc,beta1,P_G_S_1[1:2]))
  print("Average of transformed estimates:")
  print(apply(est,2,mean))
  print("Average of original estimates:")
  print(apply(orig_est,2,mean))
  print("Covariance of transformed estimates:")
  print((n_case+n_control)*cov(est))
  print("Covariance of original estimates:")
  print((n_case+n_control)*cov(orig_est))
  #print((n_case+n_control)*cov_est_delta_Total/R)
  print("Covariance of original estimates:")
  print(info_Total/(R*(n_case+n_control)))
  print("Average over inverse information (blockwise):")
  print((n_case+n_control)*cov_est_info_Total/R)
  print("True G MLE:")
  print(rej_T/R)
  print("Reads MLE:")
  print(rej_info/R)
  print("END OF THIS RUN")
  print("END OF THIS RUN")
  print("END OF THIS RUN")
  
  par(mfrow=c(1,2))
  qqplot(-log10(1:R/R),-log10(p_val_true),pch=1,main="MLE with True G",xlab="Theoretical",ylab="Empirical")
  abline(0,1,col="red")
  qqplot(-log10(1:R/R),-log10(p_val),pch=1,main="MLE with Reads",xlab="Theoretical",ylab="Empirical")
  abline(0,1,col="red")
  par(mfrow=c(1,1))
  
  par(mfrow=c(1,3))
  hist(orig_est[,1])
  hist(orig_est[,2])
  hist(orig_est[,3])
  par(mfrow=c(1,1))
  
  # return values:
  # est           all ML estimates, and estimates of genotype distribution (order of columns are beta0, beta1, P(G=0), P(G=1))
  # p_val_true    all p-values with true genotype
  # p_val         all ML p-values
  
  return(list(est,p_val_true,p_val))
}