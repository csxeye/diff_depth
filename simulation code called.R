# simulation code with called genotype
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

analysis_w_called<-function(n_case,n_control,rd_case_mean,rd_case_sd,rd_con_mean,rd_con_sd,p,e_case_mean,e_case_sd,e_con_mean,e_con_sd,beta0,beta1,R){
  
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
  
  pos_g<-c("A","G","T","C")
  
  # get beta0_cc, the intercept within case-control first
  P_Y_1<-sum(expit(beta0+beta1*c(0,1,2))*c((1-p)^2,2*p*(1-p),p^2))
  
  all_est<-NULL
  
  bad_counter<-0
  
  # save p values
  p_val_called<-NULL
  
  for(i in 1:R){
    if( (i %% 1000) == 0 ){print(i)}
    
    a<-"B"
    
    while( (typeof(a)=="character") ){
      
      # generate genotype of controls
      P_G_Y_0<-(1-expit(beta0+beta1*c(0,1,2)))*c((1-p)^2,2*p*(1-p),p^2)/(1-P_Y_1)
      g_control<-as.vector(t(rmultinom(n_control,1,P_G_Y_0))%*%c(0,1,2))
      
      # generate genotype of cases
      P_G_Y_1<-expit(beta0+beta1*c(0,1,2))*c((1-p)^2,2*p*(1-p),p^2)/P_Y_1
      g_case<-as.vector(t(rmultinom(n_case,1,P_G_Y_1))%*%c(0,1,2))
      
      # get reads controls
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
      
      # get called genotypes
      P_G_EM_case<-EM_calc(P_D_G_case)
      P_G_EM_control<-EM_calc(P_D_G_control)
      
      called_control<-NULL
      for(j in 1:n_control){
        cur_P_D_and_G<-P_D_G_control[j,]*P_G_EM_control
        
        if(max(cur_P_D_and_G)==cur_P_D_and_G[1]){
          called_control<-c(called_control,0)
        } else if(max(cur_P_D_and_G)==cur_P_D_and_G[2]){
          called_control<-c(called_control,1)
        } else{
          called_control<-c(called_control,2)
        }
      }
      
      called_case<-NULL
      for(j in 1:n_case){
        cur_P_D_and_G<-P_D_G_case[j,]*P_G_EM_case
        
        if(max(cur_P_D_and_G)==cur_P_D_and_G[1]){
          called_case<-c(called_case,0)
        } else if(max(cur_P_D_and_G)==cur_P_D_and_G[2]){
          called_case<-c(called_case,1)
        } else{
          called_case<-c(called_case,2)
        }
      }
      
      # fit regression
      y<-c(rep(1,n_case),rep(0,n_control))
      called_all<-c(called_case,called_control)
      
      #################
      #################
      a<-tryCatch({
        glm(y~called_all,family=binomial)
      }, warning = function(w){
        return("B")
      }, error = function(e){
        return("B")
      }
      )
      
      
      if( (typeof(a)=="character") ){
        bad_counter<-bad_counter+1
      }
      #################
      #################
      
    }
    
    reg_called<-glm(y~called_all,family=binomial)
    
    all_est<-rbind(all_est,c(reg_called$coefficients))
    
    p_val_called<-c(p_val_called,summary(reg_called)$coefficients[2,4])
    
  }
  
  # Outputs
  print("END OF THIS RUN")
  print("END OF THIS RUN")
  print("END OF THIS RUN")
  
  # return values:
  # all_est       all naive estimates with called genotype
  # p_val_called  all p-values with called genotype
  
  return(list(all_est,p_val_called))
}