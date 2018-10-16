# script for data analysis with called genotype

expit<-function(x){return(exp(x)/(1+exp(x)))}
logit<-function(x){return(log(x/(1-x)))}

# set directory to the results of process_vcf_files_filter_called.py

# get panel data from 1000 genomes
otg_panel<-read.table("phase1_integrated_calls.20101123.ALL.panel.txt",header=F,fill=T)
CEU_GBR_names<-as.character(otg_panel[which(otg_panel[,2]=="CEU" | otg_panel[,2]=="GBR"),1])

# get called genotype for PC computation
incon_rs_num<-0
incon_a_def<-0
in_bad_snp_list<-0

stored_cn_af_G<-list()

for(cn in 1:22){
  print(paste("Chromosome ",cn,sep=""))
  
  ali_data<-read.csv(paste0("called_",cn,"_ALI_data.csv"),header=T)
  
  otg_data<-read.csv(paste0("called_",cn,"_OTG_data.csv"),header=T)
  
  snp_info<-read.csv(paste0("called_",cn,"_snp_info.csv"),header=T)
  
  cur_cn_af_G<-list(NULL,NULL,NULL,NULL,NULL)
  
  for(i in 1:dim(snp_info)[1]){
    
    if( nchar(toString(snp_info[i,2]))>0 & nchar(toString(snp_info[i,4]))>0 & toString(snp_info[i,2])!=toString(snp_info[i,4]) ){
      
      incon_rs_num<-incon_rs_num+1
      
    } else if( toString(snp_info[i,6])!=toString(snp_info[i,8]) | toString(snp_info[i,7])!=toString(snp_info[i,9]) ){
      
      incon_a_def<-incon_a_def+1
      
    } else{
      
      
      G_case<-ali_data[(sum(snp_info[0:(i-1),3])+1):sum(snp_info[1:i,3]),1:2]
      G_control<-otg_data[(sum(snp_info[0:(i-1),5])+1):sum(snp_info[1:i,5]),1:2]
      
      cur_control_names<-as.character(G_control[,1])
      cur_case_names<-G_case[,1]
      
      G_control<-G_control[which(cur_control_names %in% CEU_GBR_names),]
      cur_control_names<-cur_control_names[which(cur_control_names %in% CEU_GBR_names)]
      
      n_case<-dim(G_case)[1]
      n_control<-dim(G_control)[1]
      
      temp_G_control<-NULL
      
      for(ct in 1:n_control){
        
        if(G_control[ct,2] == "0/0" | G_control[ct,2] == "0|0"){
          temp_G_control<-c(temp_G_control,0)
        } else if(G_control[ct,2] == "1/1" | G_control[ct,2] == "1|1"){
          temp_G_control<-c(temp_G_control,2)
        } else{
          temp_G_control<-c(temp_G_control,1)
        }
        
      }
      
      temp_G_case<-NULL
      
      for(cs in 1:n_case){
        
        if(G_case[cs,2] == "0/0" | G_case[cs,2] == "0|0"){
          temp_G_case<-c(temp_G_case,0)
        } else if(G_case[cs,2] == "1/1" | G_case[cs,2] == "1|1"){
          temp_G_case<-c(temp_G_case,2)
        } else{
          temp_G_case<-c(temp_G_case,1)
        }
        
      }
      
      P_G_control<-c( length(which(temp_G_control==0)),length(which(temp_G_control==1)),length(which(temp_G_control==2)) )/n_control
      
      P_G_case<-c( length(which(temp_G_case==0)),length(which(temp_G_case==1)),length(which(temp_G_case==2)) )/n_case
      
      if( (P_G_control[1]>0.0025) & (P_G_control[2]>0.095) & (P_G_control[2]<=0.5) & (P_G_control[3]>0.0025) & (P_G_case[1]>0.0025) & (P_G_case[2]>0.095) & (P_G_case[2]<=0.5) & (P_G_case[3]>0.0025) & (n_case==89) & (n_control==174) ){
        
        # order according to name
        temp_G_control<-temp_G_control[order(cur_control_names)]
        temp_G_case<-temp_G_case[order(cur_case_names)]
        
        cur_cn_af_G[[1]]<-c(cur_cn_af_G[[1]],toString(snp_info[i,4]))
        cur_cn_af_G[[2]]<-rbind(cur_cn_af_G[[2]],P_G_control)
        cur_cn_af_G[[3]]<-rbind(cur_cn_af_G[[3]],P_G_case)
        cur_cn_af_G[[4]]<-rbind(cur_cn_af_G[[4]],temp_G_control)
        cur_cn_af_G[[5]]<-rbind(cur_cn_af_G[[5]],temp_G_case)
        
      }
      
      
    }
  }
  
  stored_cn_af_G<-c(stored_cn_af_G,list(cur_cn_af_G))
  
}


###
###
###
# Now do svd for the PCs
control_G_for_PC<-NULL
case_G_for_PC<-NULL

for(i in 1:22){
  
  control_G_for_PC<-rbind(control_G_for_PC,stored_cn_af_G[[i]][[4]])
  case_G_for_PC<-rbind(case_G_for_PC,stored_cn_af_G[[i]][[5]])
  
}

G_normalized_combined<-NULL

for(i in 1:dim(control_G_for_PC)[1]){
  print(i)
  cur_normalized<-c(control_G_for_PC[i,],case_G_for_PC[i,])
  cur_normalized<-cur_normalized-mean(cur_normalized)
  cur_normalized<-cur_normalized/sd(cur_normalized)
  G_normalized_combined<-rbind(G_normalized_combined,cur_normalized)
}

svd_g_mat<-svd(t(G_normalized_combined))

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

bad_snp_PC<-list() # names of SNPs that run into error in code

for(cn in 1:22){
  print(paste("Chromosome ",cn,":______________________________________________________________________________________________________",sep=""))
  
  ali_data<-read.csv(paste0("called_",cn,"_ALI_data.csv"),header=T)
  
  otg_data<-read.csv(paste0("called_",cn,"_OTG_data.csv"),header=T)
  
  snp_info<-read.csv(paste0("called_",cn,"_snp_info.csv"),header=T)
  
  cur_bad_snp_PC<-NULL
  
  cur_res_w_PC<-list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)
  
  for(i in 1:dim(snp_info)[1]){
    
    if( nchar(toString(snp_info[i,2]))>0 & nchar(toString(snp_info[i,4]))>0 & toString(snp_info[i,2])!=toString(snp_info[i,4]) ){
      
      print("Different RS numbers:")
      print(snp_info[i,])
      
    } else if( toString(snp_info[i,6])!=toString(snp_info[i,8]) | toString(snp_info[i,7])!=toString(snp_info[i,9]) ){
      
      print("Inconsistent definition of major and/or minor alleles:")
      print(snp_info[i,])
      
    } else{
      
      G_case<-ali_data[(sum(snp_info[0:(i-1),3])+1):sum(snp_info[1:i,3]),1:2]
      G_control<-otg_data[(sum(snp_info[0:(i-1),5])+1):sum(snp_info[1:i,5]),1:2]
      
      cur_control_names<-as.character(G_control[,1])
      cur_case_names<-G_case[,1]
      
      G_control<-G_control[which(cur_control_names %in% CEU_GBR_names),]
      cur_control_names<-cur_control_names[which(cur_control_names %in% CEU_GBR_names)]
      
      n_case<-dim(G_case)[1]
      n_control<-dim(G_control)[1]
      
      temp_G_control<-NULL
      
      for(ct in 1:n_control){
        
        if(G_control[ct,2] == "0/0" | G_control[ct,2] == "0|0"){
          temp_G_control<-c(temp_G_control,0)
        } else if(G_control[ct,2] == "1/1" | G_control[ct,2] == "1|1"){
          temp_G_control<-c(temp_G_control,2)
        } else{
          temp_G_control<-c(temp_G_control,1)
        }
        
      }
      
      temp_G_case<-NULL
      
      for(cs in 1:n_case){
        
        if(G_case[cs,2] == "0/0" | G_case[cs,2] == "0|0"){
          temp_G_case<-c(temp_G_case,0)
        } else if(G_case[cs,2] == "1/1" | G_case[cs,2] == "1|1"){
          temp_G_case<-c(temp_G_case,2)
        } else{
          temp_G_case<-c(temp_G_case,1)
        }
        
      }
      
      P_G_control<-c( length(which(temp_G_control==0)),length(which(temp_G_control==1)),length(which(temp_G_control==2)) )/n_control
      
      P_G_case<-c( length(which(temp_G_case==0)),length(which(temp_G_case==1)),length(which(temp_G_case==2)) )/n_case
      
      # pass filter for analysis
      if( (P_G_control[1]>0.0025) & (P_G_control[2]>0.095) & (P_G_control[2]<=0.5) & (P_G_control[3]>0.0025) & (P_G_case[1]>0.0025) & (P_G_case[2]>0.095) & (P_G_case[2]<=0.5) & (P_G_case[3]>0.0025) ){
        
        cur_PC1<-c( case_PC[paste0(cur_case_names,""),1] , control_PC[paste0(cur_control_names,""),1] )
        cur_PC2<-c( case_PC[paste0(cur_case_names,""),2] , control_PC[paste0(cur_control_names,""),2] )
        
        y<-c(rep(1,n_case),rep(0,n_control))
        G_all<-c(temp_G_case,temp_G_control)
        
        a<-tryCatch({
          # do analysis here
          reg<-glm(y~G_all+cur_PC1+cur_PC2,family=binomial)
          
        }, warning = function(w){
          return("B")
        }, error = function(e){
          return("B")
        }
        )
        
        if(typeof(a)!="character"){
          
          called_pv<-summary(reg)$coef[2,4]
          cur_est<-reg$coef
          
          cur_res_w_PC[[1]]<-c(cur_res_w_PC[[1]],toString(snp_info[i,4])) # name of SNP
          
          cur_res_w_PC[[2]]<-c(cur_res_w_PC[[2]],called_pv) # p-value with called genotype
          cur_res_w_PC[[3]]<-rbind(cur_res_w_PC[[3]],cur_est) # estimates with called genotype
          
          cur_res_w_PC[[4]]<-c(cur_res_w_PC[[4]],n_case) # number of cases
          cur_res_w_PC[[5]]<-c(cur_res_w_PC[[5]],n_control) # number of controls
          
          cur_res_w_PC[[6]]<-rbind(cur_res_w_PC[[6]],P_G_case) # estimate of P(G) in cases
          cur_res_w_PC[[7]]<-rbind(cur_res_w_PC[[7]],P_G_control) # estimate of P(G) in controls
        } else{
          cur_bad_snp_PC<-c(cur_bad_snp_PC,toString(snp_info[i,4]))
        }
        
      }
    }
  }
  
  bad_snp_PC<-c(bad_snp_PC,list(cur_bad_snp_PC))
  
  res_w_PC<-c(res_w_PC,list(cur_res_w_PC))
  
}