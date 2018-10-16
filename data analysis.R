# script for data analysis

# use some of Derbach et al.'s code
source('em.r')

expit<-function(x){return(exp(x)/(1+exp(x)))}
logit<-function(x){return(log(x/(1-x)))}

# set directory to the results of process_vcf_files_filter.py

# get panel data from 1000 genomes
otg_panel<-read.table("phase1_integrated_calls.20101123.ALL.panel.txt",header=F,fill=T)
CEU_GBR_names<-as.character(otg_panel[which(otg_panel[,2]=="CEU" | otg_panel[,2]=="GBR"),1])

# get E(G|reads) for PC computation
incon_rs_num<-0
incon_a_def<-0
in_bad_snp_list<-0

num_pass_filter<-0

stored_cn_af_EGD<-list()

for(cn in 1:22){
  print(paste("Chromosome ",cn,sep=""))
  
  ali_data<-read.csv(paste0(cn,"_ALI_data.csv"),header=T)
  
  otg_data<-read.csv(paste0(cn,"_OTG_data.csv"),header=T)
  
  snp_info<-read.csv(paste0(cn,"_snp_info.csv"),header=T)
  
  cur_cn_af_EGD<-list(NULL,NULL,NULL,NULL,NULL,NULL)
  
  for(i in 1:dim(snp_info)[1]){
    
    if( nchar(toString(snp_info[i,2]))>0 & nchar(toString(snp_info[i,4]))>0 & toString(snp_info[i,2])!=toString(snp_info[i,4]) ){
      
      incon_rs_num<-incon_rs_num+1
      
    } else if( toString(snp_info[i,6])!=toString(snp_info[i,8]) | toString(snp_info[i,7])!=toString(snp_info[i,9]) ){
      
      incon_a_def<-incon_a_def+1
      
    } else{
      
      
      P_D_G_case<-ali_data[(sum(snp_info[0:(i-1),3])+1):sum(snp_info[1:i,3]),1:4]
      P_D_G_control<-otg_data[(sum(snp_info[0:(i-1),5])+1):sum(snp_info[1:i,5]),1:4]
      
      all_zeros<-NULL
      
      for(j in 1:dim(P_D_G_control)[1]){
        if( all(P_D_G_control[j,2:4]==c(0,0,0)) ){
          all_zeros<-c(all_zeros,j)
        }
      }
      
      if(length(all_zeros)!=0){
        P_D_G_control<-P_D_G_control[-all_zeros,]
      }
      
      cur_control_names<-P_D_G_control[,1]
      cur_case_names<-P_D_G_case[,1]
      
      P_D_G_control<-P_D_G_control[,2:4]
      P_D_G_control<-P_D_G_control[which(cur_control_names %in% CEU_GBR_names),]
      P_D_G_case<-P_D_G_case[,2:4]
      
      P_D_G_case<-10^(-P_D_G_case/10)
      P_D_G_control<-10^P_D_G_control
      
      n_case<-dim(P_D_G_case)[1]
      n_control<-dim(P_D_G_control)[1]
      
      P_G_EM<-EM_calc(rbind(P_D_G_control,P_D_G_case))
      
      P_G_EM_case<-EM_calc(P_D_G_case)
      P_G_EM_control<-EM_calc(P_D_G_control)
      
      if( (P_G_EM_control[1]>0.0025) & (P_G_EM_control[2]>0.095) & (P_G_EM_control[2]<=0.5) & (P_G_EM_control[3]>0.0025) & (P_G_EM_case[1]>0.0025) & (P_G_EM_case[2]>0.095) & (P_G_EM_case[2]<=0.5) & (P_G_EM_case[3]>0.0025) & (n_case==89) & (n_control==174) ){
        
        E_G_D_control<-NULL
        for(j in 1:n_control){
          P_D_and_G<-P_D_G_control[j,]*P_G_EM
          E_G_D<-sum(c(0,1,2)*P_D_and_G/sum(P_D_and_G))
          E_G_D_control<-c(E_G_D_control,E_G_D)
        }
        
        E_G_D_control<-E_G_D_control[order(as.character(cur_control_names[which(cur_control_names %in% CEU_GBR_names)]))]
        
        E_G_D_case<-NULL
        for(j in 1:n_case){
          P_D_and_G<-P_D_G_case[j,]*P_G_EM
          E_G_D<-sum(c(0,1,2)*P_D_and_G/sum(P_D_and_G))
          E_G_D_case<-c(E_G_D_case,E_G_D)
        }
        
        E_G_D_case<-E_G_D_case[order(cur_case_names)]
        
        cur_cn_af_EGD[[1]]<-c(cur_cn_af_EGD[[1]],toString(snp_info[i,4]))
        cur_cn_af_EGD[[2]]<-rbind(cur_cn_af_EGD[[2]],P_G_EM)
        cur_cn_af_EGD[[3]]<-rbind(cur_cn_af_EGD[[3]],P_G_EM_case)
        cur_cn_af_EGD[[4]]<-rbind(cur_cn_af_EGD[[4]],P_G_EM_control)
        cur_cn_af_EGD[[5]]<-rbind(cur_cn_af_EGD[[5]],E_G_D_control)
        cur_cn_af_EGD[[6]]<-rbind(cur_cn_af_EGD[[6]],E_G_D_case)
        
      }
      
      if( (P_G_EM_control[1]>0.0025) & (P_G_EM_control[2]>0.095) & (P_G_EM_control[2]<=0.5) & (P_G_EM_control[3]>0.0025) & (P_G_EM_case[1]>0.0025) & (P_G_EM_case[2]>0.095) & (P_G_EM_case[2]<=0.5) & (P_G_EM_case[3]>0.0025) ){
        num_pass_filter<-num_pass_filter+1
      }
      
      
    }
  }
  
  stored_cn_af_EGD<-c(stored_cn_af_EGD,list(cur_cn_af_EGD))
  
}


###
###
###
# Now do svd for the PCs
control_E_G_D_for_PC<-NULL
case_E_G_D_for_PC<-NULL

for(i in 1:22){
  
  control_E_G_D_for_PC<-rbind(control_E_G_D_for_PC,stored_cn_af_EGD[[i]][[5]])
  case_E_G_D_for_PC<-rbind(case_E_G_D_for_PC,stored_cn_af_EGD[[i]][[6]])
  
}

EGD_normalized_combined<-NULL

for(i in 1:dim(control_E_G_D_for_PC)[1]){
  print(i)
  cur_normalized<-c(control_E_G_D_for_PC[i,],case_E_G_D_for_PC[i,])
  cur_normalized<-cur_normalized-mean(cur_normalized)
  cur_normalized<-cur_normalized/sd(cur_normalized)
  EGD_normalized_combined<-rbind(EGD_normalized_combined,cur_normalized)
}

svd_g_mat<-svd(t(EGD_normalized_combined))

###
###
###
###
###
###
# Now do analysis with PCs

control_PC<-svd_g_mat$u[1:n_control,1:2]
case_PC<-svd_g_mat$u[(n_control+1):(n_control+n_case),1:2]

# get sample names onto PCs
ali_data<-read.csv(paste0("1_ALI_data.csv"),header=T)
otg_data<-read.csv(paste0("1_OTG_data.csv"),header=T)
snp_info<-read.csv(paste0("1_snp_info.csv"),header=T)

for(i in 1:dim(snp_info)[1]){
  
  cur_case_names<-as.character(ali_data[(sum(snp_info[0:(i-1),3])+1):sum(snp_info[1:i,3]),1])
  
  if(length(cur_case_names)==89){
    case_names<-sort(cur_case_names)
    break
  }
  
}

for(i in 1:dim(snp_info)[1]){
  
  cur_control_names<-as.character(otg_data[(sum(snp_info[0:(i-1),5])+1):sum(snp_info[1:i,5]),1])
  
  if( length(which(cur_control_names %in% CEU_GBR_names))==174 ){
    control_names<-sort(cur_control_names[which(cur_control_names %in% CEU_GBR_names)])
    break
  }
  
}

rownames(control_PC)<-control_names
rownames(case_PC)<-case_names

# now actual computation
res_w_PC<-list() # list for all results with PC

RC_bad_snp_PC<-list() # names of SNPs that run into error in code

for(cn in 1:22){
  print(paste("Chromosome ",cn,":______________________________________________________________________________________________________",sep=""))
  
  ali_data<-read.csv(paste0(cn,"_ALI_data.csv"),header=T)
  
  otg_data<-read.csv(paste0(cn,"_OTG_data.csv"),header=T)
  
  snp_info<-read.csv(paste0(cn,"_snp_info.csv"),header=T)
  
  cur_RC_bad_snp_PC<-NULL
  
  cur_res_w_PC<-list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)
  
  for(i in 1:dim(snp_info)[1]){
    
    if( nchar(toString(snp_info[i,2]))>0 & nchar(toString(snp_info[i,4]))>0 & toString(snp_info[i,2])!=toString(snp_info[i,4]) ){
      
      print("Different RS numbers:")
      print(snp_info[i,])
      
    } else if( toString(snp_info[i,6])!=toString(snp_info[i,8]) | toString(snp_info[i,7])!=toString(snp_info[i,9]) ){
      
      print("Inconsistent definition of major and/or minor alleles:")
      print(snp_info[i,])
      
    } else{
      
      P_D_G_case<-ali_data[(sum(snp_info[0:(i-1),3])+1):sum(snp_info[1:i,3]),1:4]
      P_D_G_control<-otg_data[(sum(snp_info[0:(i-1),5])+1):sum(snp_info[1:i,5]),1:4]
      
      all_zeros<-NULL
      
      for(j in 1:dim(P_D_G_control)[1]){
        if( all(P_D_G_control[j,2:4]==c(0,0,0)) ){
          all_zeros<-c(all_zeros,j)
        }
      }
      
      if(length(all_zeros)!=0){
        P_D_G_control<-P_D_G_control[-all_zeros,]
      }
      
      cur_control_names<-P_D_G_control[,1]
      cur_case_names<-P_D_G_case[,1]
      
      P_D_G_control<-P_D_G_control[,2:4]
      P_D_G_control<-P_D_G_control[which(cur_control_names %in% CEU_GBR_names),]
      cur_control_names<-cur_control_names[which(cur_control_names %in% CEU_GBR_names)]
      P_D_G_case<-P_D_G_case[,2:4]
      
      P_D_G_case<-10^(-P_D_G_case/10)
      P_D_G_control<-10^P_D_G_control
      
      n_case<-dim(P_D_G_case)[1]
      n_control<-dim(P_D_G_control)[1]
      
      # get RC p-value
      P_G_EM<-EM_calc(rbind(P_D_G_control,P_D_G_case))
      
      P_G_EM_case<-EM_calc(P_D_G_case)
      P_G_EM_control<-EM_calc(P_D_G_control)
      
      if( (P_G_EM_control[1]>0.0025) & (P_G_EM_control[2]>0.095) & (P_G_EM_control[2]<=0.5) & (P_G_EM_control[3]>0.0025) & (P_G_EM_case[1]>0.0025) & (P_G_EM_case[2]>0.095) & (P_G_EM_case[2]<=0.5) & (P_G_EM_case[3]>0.0025) ){
        
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
        
        cur_PC1<-c( case_PC[paste0(cur_case_names,""),1] , control_PC[paste0(cur_control_names,""),1] )
        cur_PC2<-c( case_PC[paste0(cur_case_names,""),2] , control_PC[paste0(cur_control_names,""),2] )
        
        y<-c(rep(1,n_case),rep(0,n_control))
        E_G_D_all<-c(E_G_D_case,E_G_D_control)
        
        a<-tryCatch({
          reg_naive<-glm(y~E_G_D_all+cur_PC1+cur_PC2,family=binomial)
          est_naive_cur<-reg_naive$coefficients
          
          # get sandwich estimator for variance
          P_D_0_EM<-as.vector(P_D_G_control%*%P_G_EM)
          P_D_1_EM<-as.vector(P_D_G_case%*%P_G_EM)
          
          bun<-matrix(rep(0,36),6,6)
          meat<-matrix(rep(0,36),6,6)
          
          for(j in 1:n_control){
            
            con_lin_pred<-(est_naive_cur%*%c(1,E_G_D_control[j],control_PC[paste0(cur_control_names,"")[j],1],control_PC[paste0(cur_control_names,"")[j],2]))[1,1]
            
            bun[1:4,1:4]<-bun[1:4,1:4]-expit( con_lin_pred )*( 1-expit( con_lin_pred ) )*(c(1,E_G_D_control[j],control_PC[paste0(cur_control_names,"")[j],1],control_PC[paste0(cur_control_names,"")[j],2])%*%t(c(1,E_G_D_control[j],control_PC[paste0(cur_control_names,"")[j],1],control_PC[paste0(cur_control_names,"")[j],2])))
            bun[5:6,5:6]<-bun[5:6,5:6]-(1/P_D_0_EM[j])^2*cbind(c( (P_D_G_control[j,1]-P_D_G_control[j,3])^2,(P_D_G_control[j,1]-P_D_G_control[j,3])*(P_D_G_control[j,2]-P_D_G_control[j,3]) ),c( (P_D_G_control[j,1]-P_D_G_control[j,3])*(P_D_G_control[j,2]-P_D_G_control[j,3]),(P_D_G_control[j,2]-P_D_G_control[j,3])^2 ))
            
            d_E_G_D_control_p0<-(-1)*P_D_G_control[j,2]*P_G_EM[2]*(P_D_G_control[j,1]-P_D_G_control[j,3])/P_D_0_EM[j]^2 - 2*P_D_G_control[j,3]*( P_D_0_EM[j]+P_G_EM[3]*(P_D_G_control[j,1]-P_D_G_control[j,3]) )/P_D_0_EM[j]^2
            d_E_G_D_control_p1<-( P_D_0_EM[j]*P_D_G_control[j,2] - P_D_G_control[j,2]*P_G_EM[2]*(P_D_G_control[j,2]-P_D_G_control[j,3]) )/P_D_0_EM[j]^2 - 2*P_D_G_control[j,3]*( P_D_0_EM[j] + P_G_EM[3]*(P_D_G_control[j,2]-P_D_G_control[j,3]) )/P_D_0_EM[j]^2
            
            d_E_G_D_control_p<-c(d_E_G_D_control_p0,d_E_G_D_control_p1)
            
            bun[1,5:6]<-bun[1,5:6] - est_naive_cur[2]*expit( con_lin_pred )*( 1-expit( con_lin_pred ) )*d_E_G_D_control_p
            bun[2,5:6]<-bun[2,5:6] + ( ( 0 - expit( con_lin_pred ) ) - est_naive_cur[2]*E_G_D_control[j]*expit( con_lin_pred )*( 1-expit( con_lin_pred ) ) )*d_E_G_D_control_p
            bun[3,5:6]<-bun[3,5:6] - control_PC[paste0(cur_control_names,"")[j],1]*est_naive_cur[2]*expit( con_lin_pred )*( 1-expit( con_lin_pred ) )*d_E_G_D_control_p
            bun[4,5:6]<-bun[4,5:6] - control_PC[paste0(cur_control_names,"")[j],2]*est_naive_cur[2]*expit( con_lin_pred )*( 1-expit( con_lin_pred ) )*d_E_G_D_control_p
            
            psi<-c( 0 - expit( con_lin_pred ) , E_G_D_control[j]*( 0 - expit( con_lin_pred ) ) , control_PC[paste0(cur_control_names,"")[j],1]*( 0 - expit( con_lin_pred ) ) , control_PC[paste0(cur_control_names,"")[j],2]*( 0 - expit( con_lin_pred ) ) , (P_D_G_control[j,1]-P_D_G_control[j,3])/P_D_0_EM[j] , (P_D_G_control[j,2]-P_D_G_control[j,3])/P_D_0_EM[j] )
            
            meat<-meat + psi%*%t(psi)
          }
          
          for(j in 1:n_case){
            
            case_lin_pred<-(est_naive_cur%*%c(1,E_G_D_case[j],case_PC[paste0(cur_case_names,"")[j],1],case_PC[paste0(cur_case_names,"")[j],2]))[1,1]
            
            bun[1:4,1:4]<-bun[1:4,1:4]-expit( case_lin_pred )*( 1-expit( case_lin_pred ) )*(c(1,E_G_D_case[j],case_PC[paste0(cur_case_names,"")[j],1],case_PC[paste0(cur_case_names,"")[j],2])%*%t(c(1,E_G_D_case[j],case_PC[paste0(cur_case_names,"")[j],1],case_PC[paste0(cur_case_names,"")[j],2])))
            bun[5:6,5:6]<-bun[5:6,5:6]-(1/P_D_1_EM[j])^2*cbind(c( (P_D_G_case[j,1]-P_D_G_case[j,3])^2,(P_D_G_case[j,1]-P_D_G_case[j,3])*(P_D_G_case[j,2]-P_D_G_case[j,3]) ),c( (P_D_G_case[j,1]-P_D_G_case[j,3])*(P_D_G_case[j,2]-P_D_G_case[j,3]),(P_D_G_case[j,2]-P_D_G_case[j,3])^2 ))
            
            d_E_G_D_case_p0<-(-1)*P_D_G_case[j,2]*P_G_EM[2]*(P_D_G_case[j,1]-P_D_G_case[j,3])/P_D_1_EM[j]^2 - 2*P_D_G_case[j,3]*( P_D_1_EM[j]+P_G_EM[3]*(P_D_G_case[j,1]-P_D_G_case[j,3]) )/P_D_1_EM[j]^2
            d_E_G_D_case_p1<-( P_D_1_EM[j]*P_D_G_case[j,2] - P_D_G_case[j,2]*P_G_EM[2]*(P_D_G_case[j,2]-P_D_G_case[j,3]) )/P_D_1_EM[j]^2 - 2*P_D_G_case[j,3]*( P_D_1_EM[j] + P_G_EM[3]*(P_D_G_case[j,2]-P_D_G_case[j,3]) )/P_D_1_EM[j]^2
            
            d_E_G_D_case_p<-c(d_E_G_D_case_p0,d_E_G_D_case_p1)
            
            bun[1,5:6]<-bun[1,5:6] - est_naive_cur[2]*expit( case_lin_pred )*( 1-expit( case_lin_pred ) )*d_E_G_D_case_p
            bun[2,5:6]<-bun[2,5:6] + ( ( 1 - expit( case_lin_pred ) ) - est_naive_cur[2]*E_G_D_case[j]*expit( case_lin_pred )*( 1-expit( case_lin_pred ) ) )*d_E_G_D_case_p
            bun[3,5:6]<-bun[3,5:6] - case_PC[paste0(cur_case_names,"")[j],1]*est_naive_cur[2]*expit( case_lin_pred )*( 1-expit( case_lin_pred ) )*d_E_G_D_case_p
            bun[4,5:6]<-bun[4,5:6] - case_PC[paste0(cur_case_names,"")[j],2]*est_naive_cur[2]*expit( case_lin_pred )*( 1-expit( case_lin_pred ) )*d_E_G_D_case_p
            
            psi<-c( 1 - expit( case_lin_pred ) , E_G_D_case[j]*( 1 - expit( case_lin_pred ) ) , case_PC[paste0(cur_case_names,"")[j],1]*( 1 - expit( case_lin_pred ) ) , case_PC[paste0(cur_case_names,"")[j],2]*( 1 - expit( case_lin_pred ) ) , (P_D_G_case[j,1]-P_D_G_case[j,3])/P_D_1_EM[j] , (P_D_G_case[j,2]-P_D_G_case[j,3])/P_D_1_EM[j] )
            
            meat<-meat + psi%*%t(psi)
          }
          
          sand_cov_est<-(solve(bun)%*%meat%*%t(solve(bun)))
          
        }, warning = function(w){
          return("B")
        }, error = function(e){
          return("B")
        }
        )
        
        if(typeof(a)!="character"){
          
          #
          # do ML here
          #
          ML_p_val<-999
          
          b<-tryCatch({
            
            gca<-cbind(rep(0,n_case),rep(1,n_case),rep(2,n_case))
            gcon<-cbind(rep(0,n_control),rep(1,n_control),rep(2,n_control))
            
            neg_log_Likelihood_2d_int<-function(t){
              
              control_P_G_g_X_S1<-cbind( expit(t[5]+t[6]*control_PC[paste0(cur_control_names,""),1]+t[7]*control_PC[paste0(cur_control_names,""),2]+t[8]*control_PC[paste0(cur_control_names,""),1]^2+t[9]*control_PC[paste0(cur_control_names,""),2]^2+t[10]*control_PC[paste0(cur_control_names,""),1]*control_PC[paste0(cur_control_names,""),2]),
                                         expit(t[5]+t[6]*control_PC[paste0(cur_control_names,""),1]+t[7]*control_PC[paste0(cur_control_names,""),2]+t[8]*control_PC[paste0(cur_control_names,""),1]^2+t[9]*control_PC[paste0(cur_control_names,""),2]^2+t[10]*control_PC[paste0(cur_control_names,""),1]*control_PC[paste0(cur_control_names,""),2]+exp(t[11]+t[12]*control_PC[paste0(cur_control_names,""),1]+t[13]*control_PC[paste0(cur_control_names,""),2]+t[14]*control_PC[paste0(cur_control_names,""),1]^2+t[15]*control_PC[paste0(cur_control_names,""),2]^2+t[16]*control_PC[paste0(cur_control_names,""),1]*control_PC[paste0(cur_control_names,""),2]))-expit(t[5]+t[6]*control_PC[paste0(cur_control_names,""),1]+t[7]*control_PC[paste0(cur_control_names,""),2]+t[8]*control_PC[paste0(cur_control_names,""),1]^2+t[9]*control_PC[paste0(cur_control_names,""),2]^2+t[10]*control_PC[paste0(cur_control_names,""),1]*control_PC[paste0(cur_control_names,""),2]),
                                         1-expit(t[5]+t[6]*control_PC[paste0(cur_control_names,""),1]+t[7]*control_PC[paste0(cur_control_names,""),2]+t[8]*control_PC[paste0(cur_control_names,""),1]^2+t[9]*control_PC[paste0(cur_control_names,""),2]^2+t[10]*control_PC[paste0(cur_control_names,""),1]*control_PC[paste0(cur_control_names,""),2]+exp(t[11]+t[12]*control_PC[paste0(cur_control_names,""),1]+t[13]*control_PC[paste0(cur_control_names,""),2]+t[14]*control_PC[paste0(cur_control_names,""),1]^2+t[15]*control_PC[paste0(cur_control_names,""),2]^2+t[16]*control_PC[paste0(cur_control_names,""),1]*control_PC[paste0(cur_control_names,""),2])) )
              
              case_P_G_g_X_S1<-cbind( expit(t[5]+t[6]*case_PC[paste0(cur_case_names,""),1]+t[7]*case_PC[paste0(cur_case_names,""),2]+t[8]*case_PC[paste0(cur_case_names,""),1]^2+t[9]*case_PC[paste0(cur_case_names,""),2]^2+t[10]*case_PC[paste0(cur_case_names,""),1]*case_PC[paste0(cur_case_names,""),2]),
                                      expit(t[5]+t[6]*case_PC[paste0(cur_case_names,""),1]+t[7]*case_PC[paste0(cur_case_names,""),2]+t[8]*case_PC[paste0(cur_case_names,""),1]^2+t[9]*case_PC[paste0(cur_case_names,""),2]^2+t[10]*case_PC[paste0(cur_case_names,""),1]*case_PC[paste0(cur_case_names,""),2]+exp(t[11]+t[12]*case_PC[paste0(cur_case_names,""),1]+t[13]*case_PC[paste0(cur_case_names,""),2]+t[14]*case_PC[paste0(cur_case_names,""),1]^2+t[15]*case_PC[paste0(cur_case_names,""),2]^2+t[16]*case_PC[paste0(cur_case_names,""),1]*case_PC[paste0(cur_case_names,""),2]))-expit(t[5]+t[6]*case_PC[paste0(cur_case_names,""),1]+t[7]*case_PC[paste0(cur_case_names,""),2]+t[8]*case_PC[paste0(cur_case_names,""),1]^2+t[9]*case_PC[paste0(cur_case_names,""),2]^2+t[10]*case_PC[paste0(cur_case_names,""),1]*case_PC[paste0(cur_case_names,""),2]),
                                      1-expit(t[5]+t[6]*case_PC[paste0(cur_case_names,""),1]+t[7]*case_PC[paste0(cur_case_names,""),2]+t[8]*case_PC[paste0(cur_case_names,""),1]^2+t[9]*case_PC[paste0(cur_case_names,""),2]^2+t[10]*case_PC[paste0(cur_case_names,""),1]*case_PC[paste0(cur_case_names,""),2]+exp(t[11]+t[12]*case_PC[paste0(cur_case_names,""),1]+t[13]*case_PC[paste0(cur_case_names,""),2]+t[14]*case_PC[paste0(cur_case_names,""),1]^2+t[15]*case_PC[paste0(cur_case_names,""),2]^2+t[16]*case_PC[paste0(cur_case_names,""),1]*case_PC[paste0(cur_case_names,""),2])) )
              
              ngl<-c( F1 = (-1)*sum( log( as.vector( ( P_D_G_case*expit(t[1]+t[2]*gca+t[3]*cbind(case_PC[paste0(cur_case_names,""),1],case_PC[paste0(cur_case_names,""),1],case_PC[paste0(cur_case_names,""),1])+t[4]*cbind(case_PC[paste0(cur_case_names,""),2],case_PC[paste0(cur_case_names,""),2],case_PC[paste0(cur_case_names,""),2]))*case_P_G_g_X_S1 )%*%c(1,1,1) ) ) ) 
                      - sum( log( as.vector( ( P_D_G_control*(1-expit(t[1]+t[2]*gcon+t[3]*cbind(control_PC[paste0(cur_control_names,""),1],control_PC[paste0(cur_control_names,""),1],control_PC[paste0(cur_control_names,""),1])+t[4]*cbind(control_PC[paste0(cur_control_names,""),2],control_PC[paste0(cur_control_names,""),2],control_PC[paste0(cur_control_names,""),2])))*control_P_G_g_X_S1 )%*%c(1,1,1) ) ) ) )
              
              return(ngl)
            }
            
            mod<-optim(c(-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),neg_log_Likelihood_2d_int,control=list(abstol=10^(-12),reltol=10^(-12),maxit=20000),hessian=T)
            m<-mod$hessian
            cur_est<-mod$par
            cov_est<-solve(m)
            
            ML_p_val<-2*(1-pnorm(abs(cur_est[2]/sqrt(cov_est[2,2]))))
            
          }, warning = function(w){
            return("B")
          }, error = function(e){
            return("B")
          }
          )
          ###
          ###
          ###
          
          #
          # do naive here
          #
          naive_p_val<-999
          
          c<-tryCatch({
            
            r_g_case<-NULL
            
            for(iii in 1:n_case){
              
              cur_P_D_and_G<-P_D_G_case[iii,]*P_G_EM_case
              
              if(max(cur_P_D_and_G)==cur_P_D_and_G[1]){
                r_g_case<-c(r_g_case,0)
              } else if(max(cur_P_D_and_G)==cur_P_D_and_G[2]){
                r_g_case<-c(r_g_case,1)
              } else{
                r_g_case<-c(r_g_case,2)
              }
            }
            
            r_g_control<-NULL
            
            for(iii in 1:n_control){
              
              cur_P_D_and_G<-P_D_G_control[iii,]*P_G_EM_control
              
              if(max(cur_P_D_and_G)==cur_P_D_and_G[1]){
                r_g_control<-c(r_g_control,0)
              } else if(max(cur_P_D_and_G)==cur_P_D_and_G[2]){
                r_g_control<-c(r_g_control,1)
              } else{
                r_g_control<-c(r_g_control,2)
              }
            }
            
            y<-c(rep(1,n_case),rep(0,n_control))
            g_called<-c(r_g_case,r_g_control)
            
            naive_p_val<-summary(glm(y~g_called,family=binomial()))$coef[2,4]
            
          }, warning = function(w){
            return("B")
          }, error = function(e){
            return("B")
          }
          )
          ###
          ###
          ###
          
          p_val_sand<-2*(1-pnorm(abs(est_naive_cur[2]/sqrt(sand_cov_est[2,2]))))
          
          cur_res_w_PC[[1]]<-c(cur_res_w_PC[[1]],toString(snp_info[i,4])) # name of SNP
          
          cur_res_w_PC[[2]]<-c(cur_res_w_PC[[2]],p_val_sand) # RC with sandwich variance p-value
          cur_res_w_PC[[3]]<-rbind(cur_res_w_PC[[3]],est_naive_cur) # RC estimates
          
          cur_res_w_PC[[4]]<-c(cur_res_w_PC[[4]],ML_p_val) # ML p-value
          cur_res_w_PC[[5]]<-c(cur_res_w_PC[[5]],naive_p_val) # p-value with called genotype
          
          cur_res_w_PC[[6]]<-c(cur_res_w_PC[[6]],n_case) # number of cases
          cur_res_w_PC[[7]]<-c(cur_res_w_PC[[7]],n_control) # number of controls
          
          cur_res_w_PC[[8]]<-rbind(cur_res_w_PC[[8]],P_G_EM) # estimate of P(G) overall
          cur_res_w_PC[[9]]<-rbind(cur_res_w_PC[[9]],P_G_EM_case) # estimate of P(G) in cases
          cur_res_w_PC[[10]]<-rbind(cur_res_w_PC[[10]],P_G_EM_control) # estimate of P(G) in controls
        } else{
          cur_RC_bad_snp_PC<-c(cur_RC_bad_snp_PC,toString(snp_info[i,4]))
        }
        
      }
    }
  }
  
  RC_bad_snp_PC<-c(RC_bad_snp_PC,list(cur_RC_bad_snp_PC))
  
  res_w_PC<-c(res_w_PC,list(cur_res_w_PC))
  
}