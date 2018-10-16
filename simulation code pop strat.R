# code for simulation with population stratification
library(truncnorm)

source('em.r')

###
###
###
###
###
###

# generate data
expit<-function(x){return(exp(x)/(1+exp(x)))}
logit<-function(x){return(log(x/(1-x)))}

# Fst value
Fst<-0.01

# number of cases
n_case<-300

# number of controls
n_control<-600

# number of SNPs
n_SNP<-20000

# param for generating phenotype, note that beta1 is effect size of population membership here
beta0<-logit(0.2)
beta1<-1

# generate phenotype with case-control sampling
# just have g_control and g_case by the end of the case-control sampling
pop_anc<-c(rep(0,2000),rep(1,2000))

pop_y<-rbinom(4000,1,expit(beta0+beta1*pop_anc))

while( !( ( sum(pop_y)>=n_case ) & ( (length(pop_y)-sum(pop_y))>=n_control ) ) ){
  anc_extra<-c(rep(0,1000),rep(1,1000))
  y_extra<-rbinom(2000,1,expit(beta0+beta1*anc_extra))
  
  pop_anc<-c(pop_anc,anc_extra)
  pop_y<-c(pop_y,y_extra)
}

cases<-which(pop_y==1)
controls<-which(pop_y==0)

case_sampled<-sample(cases,n_case)
control_sampled<-sample(controls,n_control)

anc_case<-pop_anc[case_sampled]
anc_control<-pop_anc[control_sampled]

y<-c(rep(1,n_case),rep(0,n_control))

pheno_anc<-cbind(c(anc_case,anc_control),y)

pheno_anc<-pheno_anc[order(pheno_anc[,1]),]

pop1_num<-length(pheno_anc[,1])-sum(pheno_anc[,1])
pop2_num<-sum(pheno_anc[,1])

# generate genotypes
pop1_genotype<-NULL
pop2_genotype<-NULL

for(i in 1:n_SNP){
  
  if((i %% 50)==0){print(i)}
  
  # ancestral allele frequency
  cur_anc_maf<-0.2
  
  # generate allele frequency of two populations with Balding-Nichols model
  pop1_maf<-rbeta(1,cur_anc_maf*(1-Fst)/Fst,(1-cur_anc_maf)*(1-Fst)/Fst)
  pop2_maf<-rbeta(1,cur_anc_maf*(1-Fst)/Fst,(1-cur_anc_maf)*(1-Fst)/Fst)
  
  pop1_genotype<-rbind(pop1_genotype,rbinom(pop1_num,2,pop1_maf))
  pop2_genotype<-rbind(pop2_genotype,rbinom(pop2_num,2,pop2_maf))
}





###
###
###

###
###
###
# now generate the reads and get genotype likelihood and conditional expectation here

# use some of Derbach et al.'s code
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

# mean and sd of read depth for cases
rd_case_mean<-20
rd_case_sd<-1

# mean and sd of read depth for controls
rd_con_mean<-3
rd_con_sd<-1

# mean and sd of error rate for cases
e_case_mean<-0.005
e_case_sd<-0.0001

# mean and sd of error rate for controls
e_con_mean<-0.05
e_con_sd<-0.001

# first sort the genotypes by case and control instead
colnames(pheno_anc)[1]<-"anc"

combined_pheno_anc_genotype<-rbind(t(pheno_anc),cbind(pop1_genotype,pop2_genotype))

combined_pheno_anc_genotype<-combined_pheno_anc_genotype[,order(combined_pheno_anc_genotype[2,])]

control_genotype<-combined_pheno_anc_genotype[3:(n_SNP+2),1:n_control]
case_genotype<-combined_pheno_anc_genotype[3:(n_SNP+2),(1+n_control):(n_control+n_case)]

# now start the generation of reads
pos_g<-c("A","G","T","C")

case_E_G_D<-NULL
control_E_G_D<-NULL

for(i in 1:50){
  
  print(i)
  
  cur_case_GL<-NULL
  cur_control_GL<-NULL
  
  for(l in ((n_SNP/50)*(i-1)+1):((n_SNP/50)*i)){
    
    g_case<-case_genotype[l,]
    g_control<-control_genotype[l,]
    
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
    
    case_E_G_D<-rbind(case_E_G_D,E_G_D_case)
    control_E_G_D<-rbind(control_E_G_D,E_G_D_control)
    
    cur_case_GL<-rbind(cur_case_GL,P_D_G_case)
    cur_control_GL<-rbind(cur_control_GL,P_D_G_control)
    
  }
  
  assign(paste0("case_GL_",i),cur_case_GL)
  assign(paste0("control_GL_",i),cur_control_GL)
  
}

###
###
###
###
###
###

# analysis without PCs first

p_val_stan<-NULL # RC with naive variance p-values
p_val_sand<-NULL # RC with sandwich variance p-values

for(i in 1:50){
  
  cur_ca_GL<-get(paste0("case_GL_",i))
  cur_co_GL<-get(paste0("control_GL_",i))
  
  for(ii in 1:(n_SNP/50)){
    
    start_time<-proc.time()
    
    E_G_D_control<-control_E_G_D[(i-1)*(n_SNP/50)+ii,]
    E_G_D_case<-case_E_G_D[(i-1)*(n_SNP/50)+ii,]
    
    P_D_G_control<-cur_co_GL[(ii-1)*n_control+(1:n_control),]
    P_D_G_case<-cur_ca_GL[(ii-1)*n_case+(1:n_case),]
    
    
    E_G_D_all<-c(E_G_D_case,E_G_D_control)
    y<-c(rep(1,n_case),rep(0,n_control))
    
    reg_naive<-glm(y~E_G_D_all,family=binomial)
    
    est_naive_cur<-reg_naive$coefficients
    
    # get sandwich estimator for variance
    P_G_EM<-EM_calc(rbind(P_D_G_control,P_D_G_case))
    
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
    
    sand_cov_est<-(solve(bun)%*%meat%*%t(solve(bun)))
    
    p_val_stan<-c(p_val_stan,summary(reg_naive)$coefficients[2,4])
    
    p_val_sand<-c( p_val_sand , 2*(1-pnorm(abs(est_naive_cur[2]/sqrt(sand_cov_est[2,2])))) )
    
    if((ii %% 50)==0){
      print(proc.time()-start_time)
    }
    
  }
}

p_val_ML<-NULL # ML p-values

gca<-cbind(rep(0,n_case),rep(1,n_case),rep(2,n_case))
gcon<-cbind(rep(0,n_control),rep(1,n_control),rep(2,n_control))

for(i in 1:50){
  
  cur_ca_GL<-get(paste0("case_GL_",i))
  cur_co_GL<-get(paste0("control_GL_",i))
  
  for(ii in 1:(n_SNP/50)){
    
    start_time<-proc.time()
    
    P_D_G_control<-cur_co_GL[(ii-1)*n_control+(1:n_control),]
    P_D_G_case<-cur_ca_GL[(ii-1)*n_case+(1:n_case),]
    
    neg_log_Likelihood<-function(t){
      fc<-function(x) expit(x)*expit(t[2])+expit(x+t[1])*(expit(t[2]+exp(t[3]))-expit(t[2]))+expit(x+2*t[1])*(1-expit(t[2]+exp(t[3])))-n_case/(n_case+n_control)
      
      b0cc<-uniroot(f=fc,interval=c(-100,100))$root
      
      a<-c( F1 = (-1)*sum( log( as.vector( ( P_D_G_case*expit(b0cc+t[1]*gca) )%*%c( expit(t[2]) , expit(t[2]+exp(t[3]))-expit(t[2]) , 1-expit(t[2]+exp(t[3])) ) ) ) ) 
            - sum( log( as.vector( ( P_D_G_control*(1-expit(b0cc+t[1]*gcon)) )%*%c( expit(t[2]) , expit(t[2]+exp(t[3]))-expit(t[2]) , 1-expit(t[2]+exp(t[3])) ) ) ) ) )
      
      return(a)
    }
    
    a<-tryCatch({
      mod<-optim(c(0,-1,0),neg_log_Likelihood,hessian=T)
      m<-mod$hessian
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
    
    if(typeof(a) != "character"){
      cov_est<-bwi
      
      cur_est<-mod$par
      
      p_val_ML<-c(p_val_ML,2*(1-pnorm(abs(cur_est[1]/sqrt(cov_est[1,1])))))
    }
    
    if((ii %% 50)==0){
      print(proc.time()-start_time)
    }
    
  }
}

p_val<-NULL # p-value with true genotype

for(i in 1:n_SNP){
  p_val<-c(p_val,summary(glm(combined_pheno_anc_genotype[2,]~combined_pheno_anc_genotype[2+i,],family=binomial))$coefficients[2,4])
}


###
###
###
###
###
###
# now w PCs

# first normalize E(G|D)
EGD_normalized_combined<-NULL

for(i in 1:n_SNP){
  #print(i)
  cur_normalized<-c(control_E_G_D[i,],case_E_G_D[i,])
  cur_normalized<-cur_normalized-mean(cur_normalized)
  cur_normalized<-cur_normalized/sd(cur_normalized)
  EGD_normalized_combined<-rbind(EGD_normalized_combined,cur_normalized)
}

# compute PCs
svd_g_mat<-svd(t(EGD_normalized_combined))

# compute p values again

p_val_stan_pc<-NULL # RC with naive variance p-values with PCs
p_val_sand_pc<-NULL # RC with sandwich variance p-values with PCs

for(i in 1:50){
  
  cur_ca_GL<-get(paste0("case_GL_",i))
  cur_co_GL<-get(paste0("control_GL_",i))
  
  for(ii in 1:(n_SNP/50)){
    
    start_time<-proc.time()
    
    E_G_D_control<-control_E_G_D[(i-1)*(n_SNP/50)+ii,]
    E_G_D_case<-case_E_G_D[(i-1)*(n_SNP/50)+ii,]
    
    P_D_G_control<-cur_co_GL[(ii-1)*n_control+(1:n_control),]
    P_D_G_case<-cur_ca_GL[(ii-1)*n_case+(1:n_case),]
    
    
    E_G_D_all<-c(E_G_D_control,E_G_D_case)
    y<-c(rep(0,n_control),rep(1,n_case))
    
    reg_naive<-glm(y~E_G_D_all+svd_g_mat$u[,1]+svd_g_mat$u[,2],family=binomial)
    
    est_naive_cur<-reg_naive$coefficients
    
    # get sandwich estimator for variance
    P_G_EM<-EM_calc(rbind(P_D_G_control,P_D_G_case))
    
    P_D_0_EM<-as.vector(P_D_G_control%*%P_G_EM)
    P_D_1_EM<-as.vector(P_D_G_case%*%P_G_EM)
    
    bun<-matrix(rep(0,36),6,6)
    meat<-matrix(rep(0,36),6,6)
    
    for(j in 1:n_control){
      
      con_lin_pred<-(est_naive_cur%*%c(1,E_G_D_control[j],svd_g_mat$u[j,1],svd_g_mat$u[j,2]))[1,1]
      
      bun[1:4,1:4]<-bun[1:4,1:4]-expit( con_lin_pred )*( 1-expit( con_lin_pred ) )*(c(1,E_G_D_control[j],svd_g_mat$u[j,1],svd_g_mat$u[j,2])%*%t(c(1,E_G_D_control[j],svd_g_mat$u[j,1],svd_g_mat$u[j,2])))
      bun[5:6,5:6]<-bun[5:6,5:6]-(1/P_D_0_EM[j])^2*cbind(c( (P_D_G_control[j,1]-P_D_G_control[j,3])^2,(P_D_G_control[j,1]-P_D_G_control[j,3])*(P_D_G_control[j,2]-P_D_G_control[j,3]) ),c( (P_D_G_control[j,1]-P_D_G_control[j,3])*(P_D_G_control[j,2]-P_D_G_control[j,3]),(P_D_G_control[j,2]-P_D_G_control[j,3])^2 ))
      
      d_E_G_D_control_p0<-(-1)*P_D_G_control[j,2]*P_G_EM[2]*(P_D_G_control[j,1]-P_D_G_control[j,3])/P_D_0_EM[j]^2 - 2*P_D_G_control[j,3]*( P_D_0_EM[j]+P_G_EM[3]*(P_D_G_control[j,1]-P_D_G_control[j,3]) )/P_D_0_EM[j]^2
      d_E_G_D_control_p1<-( P_D_0_EM[j]*P_D_G_control[j,2] - P_D_G_control[j,2]*P_G_EM[2]*(P_D_G_control[j,2]-P_D_G_control[j,3]) )/P_D_0_EM[j]^2 - 2*P_D_G_control[j,3]*( P_D_0_EM[j] + P_G_EM[3]*(P_D_G_control[j,2]-P_D_G_control[j,3]) )/P_D_0_EM[j]^2
      
      d_E_G_D_control_p<-c(d_E_G_D_control_p0,d_E_G_D_control_p1)
      
      bun[1,5:6]<-bun[1,5:6] - est_naive_cur[2]*expit( con_lin_pred )*( 1-expit( con_lin_pred ) )*d_E_G_D_control_p
      bun[2,5:6]<-bun[2,5:6] + ( ( 0 - expit( con_lin_pred ) ) - est_naive_cur[2]*E_G_D_control[j]*expit( con_lin_pred )*( 1-expit( con_lin_pred ) ) )*d_E_G_D_control_p
      bun[3,5:6]<-bun[3,5:6] - svd_g_mat$u[j,1]*est_naive_cur[2]*expit( con_lin_pred )*( 1-expit( con_lin_pred ) )*d_E_G_D_control_p
      bun[4,5:6]<-bun[4,5:6] - svd_g_mat$u[j,2]*est_naive_cur[2]*expit( con_lin_pred )*( 1-expit( con_lin_pred ) )*d_E_G_D_control_p
      
      psi<-c( 0 - expit( con_lin_pred ) , E_G_D_control[j]*( 0 - expit( con_lin_pred ) ) , svd_g_mat$u[j,1]*( 0 - expit( con_lin_pred ) ) , svd_g_mat$u[j,2]*( 0 - expit( con_lin_pred ) ) , (P_D_G_control[j,1]-P_D_G_control[j,3])/P_D_0_EM[j] , (P_D_G_control[j,2]-P_D_G_control[j,3])/P_D_0_EM[j] )
      
      meat<-meat + psi%*%t(psi)
    }
    
    for(j in 1:n_case){
      
      case_lin_pred<-(est_naive_cur%*%c(1,E_G_D_case[j],svd_g_mat$u[j+n_control,1],svd_g_mat$u[j+n_control,2]))[1,1]
      
      bun[1:4,1:4]<-bun[1:4,1:4]-expit( case_lin_pred )*( 1-expit( case_lin_pred ) )*(c(1,E_G_D_case[j],svd_g_mat$u[j+n_control,1],svd_g_mat$u[j+n_control,2])%*%t(c(1,E_G_D_case[j],svd_g_mat$u[j+n_control,1],svd_g_mat$u[j+n_control,2])))
      bun[5:6,5:6]<-bun[5:6,5:6]-(1/P_D_1_EM[j])^2*cbind(c( (P_D_G_case[j,1]-P_D_G_case[j,3])^2,(P_D_G_case[j,1]-P_D_G_case[j,3])*(P_D_G_case[j,2]-P_D_G_case[j,3]) ),c( (P_D_G_case[j,1]-P_D_G_case[j,3])*(P_D_G_case[j,2]-P_D_G_case[j,3]),(P_D_G_case[j,2]-P_D_G_case[j,3])^2 ))
      
      d_E_G_D_case_p0<-(-1)*P_D_G_case[j,2]*P_G_EM[2]*(P_D_G_case[j,1]-P_D_G_case[j,3])/P_D_1_EM[j]^2 - 2*P_D_G_case[j,3]*( P_D_1_EM[j]+P_G_EM[3]*(P_D_G_case[j,1]-P_D_G_case[j,3]) )/P_D_1_EM[j]^2
      d_E_G_D_case_p1<-( P_D_1_EM[j]*P_D_G_case[j,2] - P_D_G_case[j,2]*P_G_EM[2]*(P_D_G_case[j,2]-P_D_G_case[j,3]) )/P_D_1_EM[j]^2 - 2*P_D_G_case[j,3]*( P_D_1_EM[j] + P_G_EM[3]*(P_D_G_case[j,2]-P_D_G_case[j,3]) )/P_D_1_EM[j]^2
      
      d_E_G_D_case_p<-c(d_E_G_D_case_p0,d_E_G_D_case_p1)
      
      bun[1,5:6]<-bun[1,5:6] - est_naive_cur[2]*expit( case_lin_pred )*( 1-expit( case_lin_pred ) )*d_E_G_D_case_p
      bun[2,5:6]<-bun[2,5:6] + ( ( 1 - expit( case_lin_pred ) ) - est_naive_cur[2]*E_G_D_case[j]*expit( case_lin_pred )*( 1-expit( case_lin_pred ) ) )*d_E_G_D_case_p
      bun[3,5:6]<-bun[3,5:6] - svd_g_mat$u[j+n_control,1]*est_naive_cur[2]*expit( case_lin_pred )*( 1-expit( case_lin_pred ) )*d_E_G_D_case_p
      bun[4,5:6]<-bun[4,5:6] - svd_g_mat$u[j+n_control,2]*est_naive_cur[2]*expit( case_lin_pred )*( 1-expit( case_lin_pred ) )*d_E_G_D_case_p
      
      psi<-c( 1 - expit( case_lin_pred ) , E_G_D_case[j]*( 1 - expit( case_lin_pred ) ) , svd_g_mat$u[j+n_control,1]*( 1 - expit( case_lin_pred ) ) , svd_g_mat$u[j+n_control,2]*( 1 - expit( case_lin_pred ) ) , (P_D_G_case[j,1]-P_D_G_case[j,3])/P_D_1_EM[j] , (P_D_G_case[j,2]-P_D_G_case[j,3])/P_D_1_EM[j] )
      
      meat<-meat + psi%*%t(psi)
    }
    
    sand_cov_est<-(solve(bun)%*%meat%*%t(solve(bun)))
    
    p_val_stan_pc<-c(p_val_stan_pc,summary(reg_naive)$coefficients[2,4])
    
    p_val_sand_pc<-c(p_val_sand_pc , 2*(1-pnorm(abs(est_naive_cur[2]/sqrt(sand_cov_est[2,2])))) )
    
    if((ii %% 50)==0){
      print(proc.time()-start_time)
    }
    
  }
}

p_val_ML_pc<-NULL # ML p-values with PCs

gca<-cbind(rep(0,n_case),rep(1,n_case),rep(2,n_case))
gcon<-cbind(rep(0,n_control),rep(1,n_control),rep(2,n_control))

control_PC<-svd_g_mat$u[1:n_control,1:2]
case_PC<-svd_g_mat$u[n_control+(1:n_case),1:2]

for(i in 1:50){
  
  cur_ca_GL<-get(paste0("case_GL_",i))
  cur_co_GL<-get(paste0("control_GL_",i))
  
  for(ii in 1:(n_SNP/50)){
    
    start_time<-proc.time()
    
    P_D_G_control<-cur_co_GL[(ii-1)*n_control+(1:n_control),]
    P_D_G_case<-cur_ca_GL[(ii-1)*n_case+(1:n_case),]
    
    control_PC<-svd_g_mat$u[1:n_control,1:2]
    case_PC<-svd_g_mat$u[n_control+(1:n_case),1:2]
    
    neg_log_Likelihood_simple<-function(t){
      
      control_P_G_g_X_S1<-c( expit(t[5]) , expit(t[5]+exp(t[6]))-expit(t[5]) , 1-expit(t[5]+exp(t[6])) )
      
      case_P_G_g_X_S1<-c( expit(t[5]) , expit(t[5]+exp(t[6]))-expit(t[5]) , 1-expit(t[5]+exp(t[6])) )
      
      a<-c( F1 = (-1)*sum( log( as.vector( ( P_D_G_case*expit(t[1]+t[2]*gca+t[3]*case_PC[,1]+t[4]*case_PC[,2]) )%*%case_P_G_g_X_S1 ) ) ) 
            - sum( log( as.vector( ( P_D_G_control*(1-expit(t[1]+t[2]*gcon+t[3]*control_PC[,1]+t[4]*control_PC[,2])) )%*%control_P_G_g_X_S1 ) ) ) )
      
      return(a)
    }
    
    a<-tryCatch({
      mod<-optim(c(-1,0,0,0,0,0),neg_log_Likelihood_simple,control=list(abstol=10^(-12),reltol=10^(-12),maxit=20000),hessian=T)
      m<-mod$hessian
      sqrt(solve(m)[2,2])
    }, warning = function(w){
      return("B")
    }, error = function(e){
      return("B")
    }
    )
    
    if(typeof(a) != "character"){
      cov_est<-solve(m)
      
      cur_est<-mod$par
      
      p_val_ML_pc<-c(p_val_ML_pc,2*(1-pnorm(abs(cur_est[2]/sqrt(cov_est[2,2])))))
    }
    
    if((ii %% 50)==0){
      print(proc.time()-start_time)
    }
    
  }
}


p_val_pc<-NULL # p-value with true genotype with PCs

for(i in 1:n_SNP){
  p_val_pc<-c(p_val_pc,summary(glm(combined_pheno_anc_genotype[2,]~combined_pheno_anc_genotype[2+i,]+svd_g_mat$u[,1]+svd_g_mat$u[,2],family=binomial))$coefficients[2,4])
}