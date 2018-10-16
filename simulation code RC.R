# simulation code with RC
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

# function for RC simulation
regression_w_E_G_D<-function(n_case,n_control,rd_case_mean,rd_case_sd,rd_con_mean,rd_con_sd,p,e_case_mean,e_case_sd,e_con_mean,e_con_sd,beta0,beta1,R,alpha=0.05){
  
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
  
  pos_g<-c("A","G","T","C")
  
  # get beta0_cc, the intercept within case-control first
  P_Y_1<-sum(expit(beta0+beta1*c(0,1,2))*c((1-p)^2,2*p*(1-p),p^2))
  beta0_cc<-beta0+log( n_case/n_control * (1-P_Y_1)/P_Y_1 )
  
  all_est<-NULL
  est<-NULL
  est_naive<-NULL
  mod_based_cov<-cbind(c(0,0),c(0,0))
  obs_info<-cbind(c(0,0),c(0,0))
  sandwich_cov<-matrix(rep(0,16),4,4)
  
  bun_total<-matrix(rep(0,16),4,4)
  meat_total<-matrix(rep(0,16),4,4)
  
  # count number of rejections less than alpha
  rej_T<-0
  rej_EGD_N<-0
  rej_EGD_S<-0
  
  bad_counter<-0
  
  # save p values
  p_val_true<-NULL
  p_val_stan<-NULL
  p_val_sand<-NULL
  
  #plot(0,0,xlim=c(0,R),ylim=c(0,R))
  
  for(i in 1:R){
    if( (i %% 1000) == 0 ){print(i)}
    
    a<-"B"
    b<-"B"
    
    while( (typeof(a)=="character") | (typeof(b)=="character") ){
      
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
      
      # get E(G|D)
      P_G_EM<-EM_calc(rbind(P_D_G_control,P_D_G_case))
      
      E_G_D_control<-NULL
      for(j in 1:n_control){
        P_D_and_G<-P_D_G_control[j,]*P_G_EM
        E_G_D<-sum(c(0,1,2)*P_D_and_G/sum(P_D_and_G))
        E_G_D_control<-c(E_G_D_control,E_G_D)
      }
      
      E_G_D_case<-NULL
      for(j in 1:n_case){
        P_D_and_G<-P_D_G_case[j,]*P_G_EM
        E_G_D<-sum(c(0,1,2)*P_D_and_G/sum(P_D_and_G))
        E_G_D_case<-c(E_G_D_case,E_G_D)
      }
      
      # fit regression
      y<-c(rep(1,n_case),rep(0,n_control))
      E_G_D_all<-c(E_G_D_case,E_G_D_control)
      
      g_all<-c(g_case,g_control)
      
      
      #################
      #################
      a<-tryCatch({
        glm(y~E_G_D_all,family=binomial)
      }, warning = function(w){
        return("B")
      }, error = function(e){
        return("B")
      }
      )
      
      
      
      b<-tryCatch({
        glm(y~g_all,family=binomial)
      }, warning = function(w){
        return("B")
      }, error = function(e){
        return("B")
      }
      )
      
      if( (typeof(a)=="character") | (typeof(b)=="character") ){
        bad_counter<-bad_counter+1
      }
      #################
      #################
      
    }
    
    reg_naive<-glm(y~E_G_D_all,family=binomial)
    est_naive<-rbind(est_naive,reg_naive$coefficients)
    
    all_est<-rbind(all_est,c(reg_naive$coefficients,P_G_EM[-3]))
    
    est_naive_cur<-reg_naive$coefficients
    
    reg<-glm(y~g_all,family=binomial)
    est<-rbind(est,reg$coefficients)
    
    # get model based covariance estimate, for regression with E(G|D)
    mod_based_cov<-mod_based_cov+summary(reg_naive)$cov.scaled
    obs_info<-obs_info+solve(summary(reg_naive)$cov.scaled)/(n_case+n_control)
    
    # get sandwich estimator for variance
    P_D_0_EM<-as.vector(P_D_G_control%*%P_G_EM)
    P_D_1_EM<-as.vector(P_D_G_case%*%P_G_EM)
    
    bun<-matrix(rep(0,16),4,4)
    meat<-matrix(rep(0,16),4,4)
    
    for(j in 1:n_control){
      bun[1:2,1:2]<-bun[1:2,1:2]-expit( (est_naive_cur%*%c(1,E_G_D_control[j]))[1,1] )*( 1-expit( (est_naive_cur%*%c(1,E_G_D_control[j]))[1,1] ) )*cbind(c(1,E_G_D_control[j]),c(E_G_D_control[j],E_G_D_control[j]^2))
      bun[3:4,3:4]<-bun[3:4,3:4]-(1/P_D_0_EM[j])^2*cbind(c( (P_D_G_control[j,1]-P_D_G_control[j,3])^2,(P_D_G_control[j,1]-P_D_G_control[j,3])*(P_D_G_control[j,2]-P_D_G_control[j,3]) ),c( (P_D_G_control[j,1]-P_D_G_control[j,3])*(P_D_G_control[j,2]-P_D_G_control[j,3]),(P_D_G_control[j,2]-P_D_G_control[j,3])^2 ))
      
      d_E_G_D_control_p0<-(-1)*P_D_G_control[j,2]*P_G_EM[2]*(P_D_G_control[j,1]-P_D_G_control[j,3])/P_D_0_EM[j]^2 - 2*P_D_G_control[j,3]*( P_D_0_EM[j]+P_G_EM[3]*(P_D_G_control[j,1]-P_D_G_control[j,3]) )/P_D_0_EM[j]^2
      d_E_G_D_control_p1<-( P_D_0_EM[j]*P_D_G_control[j,2] - P_D_G_control[j,2]*P_G_EM[2]*(P_D_G_control[j,2]-P_D_G_control[j,3]) )/P_D_0_EM[j]^2 - 2*P_D_G_control[j,3]*( P_D_0_EM[j] + P_G_EM[3]*(P_D_G_control[j,2]-P_D_G_control[j,3]) )/P_D_0_EM[j]^2
      
      d_E_G_D_control_p<-c(d_E_G_D_control_p0,d_E_G_D_control_p1)
      
      bun[1,3:4]<-bun[1,3:4] - est_naive_cur[2]*expit( (est_naive_cur%*%c(1,E_G_D_control[j]))[1,1] )*( 1-expit( (est_naive_cur%*%c(1,E_G_D_control[j]))[1,1] ) )*d_E_G_D_control_p
      bun[2,3:4]<-bun[2,3:4] + ( ( 0 - expit( (est_naive_cur%*%c(1,E_G_D_control[j]))[1,1] ) ) - est_naive_cur[2]*E_G_D_control[j]*expit( (est_naive_cur%*%c(1,E_G_D_control[j]))[1,1] )*( 1-expit( (est_naive_cur%*%c(1,E_G_D_control[j]))[1,1] ) ) )*d_E_G_D_control_p
      
      psi<-c( 0 - expit( (est_naive_cur%*%c(1,E_G_D_control[j]))[1,1] ), E_G_D_control[j]*( 0 - expit( (est_naive_cur%*%c(1,E_G_D_control[j]))[1,1] ) ) , (P_D_G_control[j,1]-P_D_G_control[j,3])/P_D_0_EM[j] , (P_D_G_control[j,2]-P_D_G_control[j,3])/P_D_0_EM[j] )
      
      meat<-meat + psi%*%t(psi)
    }
    
    for(j in 1:n_case){
      bun[1:2,1:2]<-bun[1:2,1:2]-expit( (est_naive_cur%*%c(1,E_G_D_case[j]))[1,1] )*( 1-expit( (est_naive_cur%*%c(1,E_G_D_case[j]))[1,1] ) )*cbind(c(1,E_G_D_case[j]),c(E_G_D_case[j],E_G_D_case[j]^2))
      bun[3:4,3:4]<-bun[3:4,3:4]-(1/P_D_1_EM[j])^2*cbind(c( (P_D_G_case[j,1]-P_D_G_case[j,3])^2,(P_D_G_case[j,1]-P_D_G_case[j,3])*(P_D_G_case[j,2]-P_D_G_case[j,3]) ),c( (P_D_G_case[j,1]-P_D_G_case[j,3])*(P_D_G_case[j,2]-P_D_G_case[j,3]),(P_D_G_case[j,2]-P_D_G_case[j,3])^2 ))
      
      d_E_G_D_case_p0<-(-1)*P_D_G_case[j,2]*P_G_EM[2]*(P_D_G_case[j,1]-P_D_G_case[j,3])/P_D_1_EM[j]^2 - 2*P_D_G_case[j,3]*( P_D_1_EM[j]+P_G_EM[3]*(P_D_G_case[j,1]-P_D_G_case[j,3]) )/P_D_1_EM[j]^2
      d_E_G_D_case_p1<-( P_D_1_EM[j]*P_D_G_case[j,2] - P_D_G_case[j,2]*P_G_EM[2]*(P_D_G_case[j,2]-P_D_G_case[j,3]) )/P_D_1_EM[j]^2 - 2*P_D_G_case[j,3]*( P_D_1_EM[j] + P_G_EM[3]*(P_D_G_case[j,2]-P_D_G_case[j,3]) )/P_D_1_EM[j]^2
      
      d_E_G_D_case_p<-c(d_E_G_D_case_p0,d_E_G_D_case_p1)
      
      bun[1,3:4]<-bun[1,3:4] - est_naive_cur[2]*expit( (est_naive_cur%*%c(1,E_G_D_case[j]))[1,1] )*( 1-expit( (est_naive_cur%*%c(1,E_G_D_case[j]))[1,1] ) )*d_E_G_D_case_p
      bun[2,3:4]<-bun[2,3:4] + ( ( 1 - expit( (est_naive_cur%*%c(1,E_G_D_case[j]))[1,1] ) ) - est_naive_cur[2]*E_G_D_case[j]*expit( (est_naive_cur%*%c(1,E_G_D_case[j]))[1,1] )*( 1-expit( (est_naive_cur%*%c(1,E_G_D_case[j]))[1,1] ) ) )*d_E_G_D_case_p
      
      psi<-c( 1 - expit( (est_naive_cur%*%c(1,E_G_D_case[j]))[1,1] ), E_G_D_case[j]*( 1 - expit( (est_naive_cur%*%c(1,E_G_D_case[j]))[1,1] ) ) , (P_D_G_case[j,1]-P_D_G_case[j,3])/P_D_1_EM[j] , (P_D_G_case[j,2]-P_D_G_case[j,3])/P_D_1_EM[j] )
      
      meat<-meat + psi%*%t(psi)
    }
    
    bun_total<-bun_total+bun/(n_case+n_control)
    meat_total<-meat_total+meat/(n_case+n_control)
    
    sand_cov_est<-(solve(bun)%*%meat%*%t(solve(bun)))
    
    sandwich_cov<-sandwich_cov+sand_cov_est
    
    if(summary(reg)$coefficients[2,4]<alpha){
      rej_T<-rej_T+1
    }
    
    p_val_true<-c(p_val_true,summary(reg)$coefficients[2,4])
    
    if(summary(reg_naive)$coefficients[2,4]<alpha){
      rej_EGD_N<-rej_EGD_N+1
    }
    
    p_val_stan<-c(p_val_stan,summary(reg_naive)$coefficients[2,4])
    
    if( abs(est_naive_cur[2]/sqrt(sand_cov_est[2,2])) > qnorm(1-alpha/2) ){
      rej_EGD_S<-rej_EGD_S+1
    }
    
    p_val_sand<-c( p_val_sand , 2*(1-pnorm(abs(est_naive_cur[2]/sqrt(sand_cov_est[2,2])))) )
  }
  
  # get expect information matrix for the estimator with true G
  i1<-sum(expit(beta0+beta1*c(0,1,2))*(1-expit(beta0+beta1*c(0,1,2)))*P_G_Y_0)*n_control/(n_control+n_case)+sum(expit(beta0+beta1*c(0,1,2))*(1-expit(beta0+beta1*c(0,1,2)))*P_G_Y_1)*n_case/(n_control+n_case)
  i2<-sum(c(0,1,2)*expit(beta0+beta1*c(0,1,2))*(1-expit(beta0+beta1*c(0,1,2)))*P_G_Y_0)*n_control/(n_control+n_case)+sum(c(0,1,2)*expit(beta0+beta1*c(0,1,2))*(1-expit(beta0+beta1*c(0,1,2)))*P_G_Y_1)*n_case/(n_control+n_case)
  i3<-sum(c(0,1,2)^2*expit(beta0+beta1*c(0,1,2))*(1-expit(beta0+beta1*c(0,1,2)))*P_G_Y_0)*n_control/(n_control+n_case)+sum(c(0,1,2)^2*expit(beta0+beta1*c(0,1,2))*(1-expit(beta0+beta1*c(0,1,2)))*P_G_Y_1)*n_case/(n_control+n_case)
  i<-cbind(c(i1,i2),c(i2,i3))
  
  # Outputs
  print("Asymptotic limit for estimator w/G:______________________________________________________________________")
  print(c(beta0_cc,beta1))
  print("Mean of estimator w/G:___________________________________________________________________________________")
  print(apply(est,2,mean))
  print("Mean of estimator w/E(G|D):______________________________________________________________________________")
  print(apply(est_naive,2,mean))
  print("#########################################################################################################")
  print("#########################################################################################################")
  print("Normalized covariance of estimator w/E(G|D) over all replications:_______________________________________")
  print(cov(sqrt(n_case+n_control)*(all_est)))
  print("Average over model-based covariance estimate:____________________________________________________________")
  print((n_case+n_control)*mod_based_cov/R)
  print("Average over model-based observed information:___________________________________________________________")
  print(obs_info/R)
  print("Average over sandwich covariance estimate:_______________________________________________________________")
  print((n_case+n_control)*sandwich_cov/R)
  print("Bun average:_____________________________________________________________________________________________")
  print(bun_total/R)
  print("Meat Average:____________________________________________________________________________________________")
  print(meat_total/R)
  print("#########################################################################################################")
  print("#########################################################################################################")
  print("Normlized covariance of estimator w/G over all replications (check stability):___________________________")
  print(cov(sqrt(n_case+n_control)*(est-cbind(rep(beta0_cc,R),rep(beta1,R)))))
  print("Inverse information matrix:______________________________________________________________________________")
  print(solve(i))
  print("#########################################################################################################")
  print("#########################################################################################################")
  print("Proportion of rejection with true G:_____________________________________________________________________")
  print(rej_T/R)
  print("Proportion of rejection with RC (naive var est):_________________________________________________________")
  print(rej_EGD_N/R)
  print("Proportion of rejection with RC (sand var est):__________________________________________________________")
  print(rej_EGD_S/R)
  print("END OF THIS RUN")
  print("END OF THIS RUN")
  print("END OF THIS RUN")
  
  par(mfrow=c(1,3))
  qqplot(-log10(1:R/R),-log10(p_val_true),pch=1,main="True G",xlab="Theoretical",ylab="Empirical")
  abline(0,1,col="red")
  qqplot(-log10(1:R/R),-log10(p_val_stan),pch=1,main="Standard",xlab="Theoretical",ylab="Empirical")
  abline(0,1,col="red")
  qqplot(-log10(1:R/R),-log10(p_val_sand),pch=1,main="Sandwich",xlab="Theoretical",ylab="Empirical")
  abline(0,1,col="red")
  par(mfrow=c(1,1))
  
  # return values:
  # all_est       all RC estimates, and estimates of genotype distribution (order of columns are beta0, beta1, P(G=0), P(G=1))
  # p_val_true    all p-values with true genotype
  # p_val_stan    all RC p-values with naive variance
  # p_val_sand    all RC p-values with sandwich variance
  
  return(list(all_est,p_val_true,p_val_stan,p_val_sand))
}