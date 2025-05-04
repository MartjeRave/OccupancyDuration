######################################################
######### Comparing Bias Correction ##################
######################################################
library(abind)
library(sqldf)
library(lubridate)
library(tidyverse)
library(plyr)
library(zoo)
library(glmnet)
library(readxl)
library(glmmLasso)
library(shapefiles)
library(sp)
library(psych)
library(ggplot2)
#library(rgeos)
library(skellam)
#library(maptools)
library(RColorBrewer)
library(gridExtra)
#library(rgeos)
library(sf)
library(terra)
install.packages("rgdal")
library(leaflet)
library(fitdistrplus)
library(lme4)
library(skellam)
library(mgcv)
library(splines)
library(dplyr)
library(plotmo)
library(readxl)
library(corrplot)
library(forecast)
library(RColorBrewer)
library(MLmetrics)
library(Hmisc)
library(psych)
library(ggdendro)
library(cluster)
library(ISLR)
library(gam)
library(VGAM)
library(cowplot)
library(miceadds)
library(gtable)
library(ggpubr)
library(scales)
library(colorspace)
library(readr)
library(taskscheduleR)
library(Rmisc)
library(sf)
library(readr)
library(readxl)
library(dplyr)
library(lubridate)
library(lazyWeave)
library(dplyr)
library(kableExtra)
library(gtools)
library(mgcv)
library(cowplot)
library(mefa)
library(quadprog)
library(data.table)
library(abind)
library(mefa)
library(quadprog)
library(data.table)
library(abind)

######################################################
##### Offset calculation through lagged incoming #####
######################################################

Off_cal<-function(data_off=datalist[[j]], w_hat=w_true){
  
  data_tot<-data_off[, !grepl("Incoming_", names(data_off))]
  data_temp<-data_off[, c("date", "districtId", "Incoming")]
  for(i in 1:length(w_hat)){
    data_temp<-data_off[, c("date", "districtId", "Incoming")]
    data_temp$date<-data_off$date+i
    names(data_temp)<-c("date", "districtId", paste0("Incoming_", i, "_"))
    data_tot<-data_tot%>%
      full_join(data_temp, by=c("date", "districtId"))
  }
  
  Offset_data<-data_tot%>%
    drop_na()
  return(Offset_data)
  
  
}

######################################################
############### Simulating data  #####################
######################################################

sim_data_offset<-function(w=w_true,
                          districtId_sim=1:200, 
                          date_sim=1:(200+ncol(w_true_Grid)),
                          ss=123, deterministic=TRUE,
                          with_coef=TRUE){
  # set seed   
  set.seed(ss)
  
  # generate M and X  
  
  if(with_coef){
    N<-rgamma(length(districtId_sim)*
                length(date_sim),0.1,0.5)
    M<-rep(rgamma(length(districtId_sim),1,3), each=length(date_sim))
    Incoming <- rpois(length(M), lambda=exp(0.5+1*M+0.2*N))    
    
    data_sim<-as.data.frame(cbind("date"=rep(date_sim, 
                                             times=length(districtId_sim)), 
                                  "districtId"=rep(districtId_sim, 
                                                   each=length(date_sim)), 
                                  M, N, Incoming))
  }else{
    Incoming <- rpois(length(date_sim)*length(districtId_sim), lambda=10)
    data_sim<-as.data.frame(cbind("date"=rep(date_sim, 
                                             times=length(districtId_sim)), 
                                  "districtId"=rep(districtId_sim, 
                                                   each=length(date_sim)), 
                                  Incoming))
  }
  
  if(deterministic){
    
    all_length<-suppressMessages(as.data.frame(cbind("total_sample"=
                                                       sample(1:length(w), 
                                                              size=sum(data_sim$Incoming), 
                                                              prob=w, replace=TRUE),
                                                     "date"=rep(data_sim$date, 
                                                                times=data_sim$Incoming),
                                                     "districtId"=rep(data_sim$districtId, 
                                                                      times=data_sim$Incoming)))%>%
                                   dplyr::group_by(date, districtId, total_sample))
    all_length<-sqldf("SELECT COUNT(*) AS count, date, districtId, total_sample FROM all_length GROUP BY total_sample, date, districtId")
    
    data_tot<-data_sim
    for(i in 1:length(w)){
      lag<-all_length[all_length$total_sample==i,]%>%
        dplyr::mutate(date=date+i)%>%
        dplyr::select(date, count, districtId)
      names(lag)<-c("date", paste0("Incom_lag",i), "districtId")
      data_tot<- suppressMessages( data_tot%>%
                                     full_join(lag, by=c("date", "districtId"))%>%
                                     na.replace(replace=0))
    }
    data_tot<-data_tot%>%
      filter(date>=min(data_sim$date)+length(w)&
               date<=max(data_sim$date))%>%
      arrange(districtId, date)
    data_tot$Outgoing<-rowSums(data_tot[,grepl("lag", names(data_tot))])
    data_tot$Diff<-data_tot$Incoming-data_tot$Outgoing
  }else{
    # data_tot<-Off_cal(data_off=data_sim, w_hat=w)
    # Offset<-t(w%*%t(data_tot[,grepl("Incoming_", names(data_tot))]))
    # Outgoing <- rpois(nrow(data_tot), lambda=Offset)
    # data_tot<-cbind(data_tot, Outgoing)
    # data_tot$Diff<-data_tot$Incoming-data_tot$Outgoing
  }
  
  return(data_tot)
}

######################################################
###############  Score Function ######################
######################################################

Score_func_w<-function(w_score_fun=c(w_true,0), 
                       data_score_fun=datalist[[j]]){
  
  Sum_score<-as.vector(t(w_score_fun%*%t(data_score_fun[,grepl("Incoming_",
                                                               names(data_score_fun))])))
  
  model_score<-mgcv::gam(Outgoing~1,
                         offset = log(Sum_score+0.1),
                         data=data_score_fun,
                         family=poisson())
  
  terms_score<-as.vector(predict.gam(model_score, 
                                     type="response"))
  score_vec<-vector(length=length(w_score_fun)-1)
  
  R_td_score<-as.vector(data_score_fun$Outgoing)
  
  I_tN_score<-as.vector(data_score_fun[,grepl(paste0("Incoming_", 
                                                     length(w_score_fun)),
                                              names(data_score_fun))])
  
  for(s in 1:length(score_vec)){
    I_ts_score<-as.vector(data_score_fun[,grepl(paste0("Incoming_", 
                                                       s, "_"),
                                                names(data_score_fun))])
    score_vec[s]<-sum((R_td_score*(I_ts_score-I_tN_score)/(Sum_score+0.00001))-
                        (I_ts_score-I_tN_score)*terms_score)
    
  }
  
  return(score_vec)
}





######################################################
###############  Negative Second Derivative ##########
######################################################

Neg2_der_w<-function(w_derivative=start_w, 
                     data_derivative=Off_data_temp){
  
  Sum_der<-as.vector(t(w_derivative%*%t(data_derivative[,grepl("Incoming_",
                                                               names(data_derivative))])))
  
  
  der_mat<-matrix(nrow=length(w_derivative)-1, ncol=length(w_derivative)-1)
  
  R_td_der<-as.vector(data_derivative$Outgoing)
  I_tN_der<-as.vector(data_derivative[,grepl(paste0("Incoming_", 
                                                    length(w_derivative)),
                                             names(data_derivative))])
  
  for(k in 1:ncol(der_mat)){
    for(s in 1:nrow(der_mat)){
      I_ts_der<-as.vector(data_derivative[,grepl(paste0("Incoming_", 
                                                        s, "_"),
                                                 names(data_derivative))])
      I_tk_der<-as.vector(data_derivative[,grepl(paste0("Incoming_", 
                                                        k, "_"),
                                                 names(data_derivative))])
      der_mat[s,k]<--sum((R_td_der*(I_ts_der-I_tN_der)*(I_tk_der-I_tN_der)/
                            ((Sum_der)^2+0.00001)))
      
    }
  }
  
  return(-der_mat)
}

######################################################
############# Optimizer Function  ####################
######################################################

optimize_Like<-function(w_opt_prev=w_try4_new, 
                        data_like_opp=datalist[[l]]){
  
  Off_data_temp<-Off_cal(w_hat=w_opt_prev, 
                         data_off=data_like_opp[, 
                                                !grepl("Incom_", names(data_like_opp))])
  
  score_w<-Score_func_w(w_score_fun=w_opt_prev,
                        data_score_fun=Off_data_temp)
  
  FisherInf<-Neg2_der_w(w_derivative=w_opt_prev,
                        data_derivative=Off_data_temp)
  
  Dmat <- FisherInf
  d<- (score_w+ FisherInf%*%w_opt_prev[-length(w_opt_prev)])
  b<-c(rep(0, length(w_opt_prev)-1), rep(-1, length(w_opt_prev)-1))
  
  
  Amat<- cbind(diag(1, length(w_opt_prev)-1), diag(-1, length(w_opt_prev)-1))
  
  Opt<-solve.QP(Dmat=Dmat,
                dvec=d,
                Amat=Amat,
                bvec=b)
  
  return(list("w_sol"=Opt$solution, "Infor"=Dmat, "w_sol_uncon"=Opt$unconstrained.solution))
}

######################################################
############# Optimizer run  #########################
######################################################

SeeSaw<-function(length_start=5, 
                 data_seesaw=datalist[[l]][,c("date", "districtId", "M", "N",
                                              "Incoming", "Outgoing", "Diff")], 
                 ss_seesaw=123,
                 no_start=NULL){
  set.seed(ss_seesaw)
  if(is.null(no_start)){
    start_w<-runif(n=length_start-1, 0, max=1/length_start)
    start_w[length(start_w)+1]<-1-sum(start_w)
    w_try4_new<-start_w
  }else{
    start_w<-no_start
    if(length_start>length(start_w)){
      for(i in 1:(length_start-length(start_w))){
        start_w<-c(start_w,1-sum(start_w))
      }  
    }
  }
  w_try4_new<-start_w
  
  for(i in 1:100){
    new_sol<-optimize_Like(w_opt_prev=w_try4_new, 
                           data_like_opp=data_seesaw)
    
    if(sum(((c(new_sol$w_sol, 1-sum(new_sol$w_sol))-w_try4_new)^2))<0.00001){
      Var<-solve(new_sol$Infor)
      break
    }else{
      w_try4_new<-c(new_sol$w_sol, 1-sum(new_sol$w_sol))      
    }}
  return(list(w=w_try4_new, VCOV_w=Var, runs=i))
}


######################################################
############# EM algorithm  ##########################
######################################################


########## Variance estimation ###########################
RubinVar<-function(Varianceruns=Run1$Coeff_var_in[200:500,,], 
                   Coefficientruns=Coeff_in[200:500,], runs=500){
  Rub1<-apply(Varianceruns,c(2,3),mean)
  betaminusbetaline<-(Coefficientruns-
                        matrix(rep(apply(Coefficientruns, 2, mean),
                                   each=nrow(Varianceruns)),
                               nrow=nrow(Varianceruns), ncol=ncol(Coefficientruns)))  
  variance_step<-array(NA,dim=c(nrow(betaminusbetaline), ncol(Coefficientruns),ncol(Coefficientruns)))
  for(i in 1:nrow(betaminusbetaline)){
    variance_step[i,,]<-t(as.matrix(betaminusbetaline[i,]))%*%
      as.matrix(betaminusbetaline[i,])
  }
  Rub2<-apply(variance_step,c(2,3),mean)
  Rubinvars<-Rub1+(1+1/nrow(Varianceruns))*(1/(nrow(Varianceruns)-1))*Rub2
  return(Rubinvars)
}





rdiffpois_array_paper =function(data_rdiff=preptrain_try,
                                IMSTEP=I_Mstep,
                                w_sim=w_start,
                                max.k=300){
  
  
  set.seed(123)
  
  
  Lagged_Data<-rep(data_rdiff[data_rdiff$date==min(data_rdiff$date),], (length(w_sim)+1))%>%
    mutate(date=date-rep(1:(length(w_sim)+1), each=length(unique(data_rdiff$districtId))))%>%
    rbind(data_rdiff)%>%
    arrange(date, districtId)
  
  
  Lagged_Data$lambda_in<-predict.gam(IMSTEP, newdata=Lagged_Data, type="response")
  
  lambda_in_array<-(Lagged_Data%>%
                      dplyr::select(date, districtId, lambda_in)%>%
                      mutate(lambda_in=ifelse(lambda_in==0, 
                                              0.000001, lambda_in))%>%
                      arrange(date, districtId)%>%
                      spread(key="date", value="lambda_in"))
  
  In_sim<-(Lagged_Data%>%
             dplyr::select(date, districtId)%>%
             mutate(in_sim=NA)%>%
             arrange(date, districtId)%>%
             spread(key="date", value="in_sim"))
  
  
  In_sim[, 2:(length(w_sim)+2)]<-apply(lambda_in_array[,2:(length(w_sim)+2)], 2, 
                                       rpois, n=nrow(In_sim))
  
  lambda_out_array<-In_sim[, -c(2:(length(w_sim)+2))]
  
  lambda_out_array[,2]<-t(w_sim%*%t(In_sim[, (length(w_sim)+2):3]))
  
  Diff_sim<-(data_rdiff%>%
               dplyr::select(date, districtId, Diff)%>%
               arrange(date, districtId)%>%
               spread(key="date", value="Diff"))
  
  Out_sim<-(data_rdiff%>%
              dplyr::select(date, districtId)%>%
              mutate(out_sim=NA)%>%
              arrange(date, districtId)%>%
              spread(key="date", value="out_sim"))
  datum<-unique(data_rdiff$date)[1]
  
  for(datum in unique(data_rdiff$date)){
    IN_poss<-matrix(rep(0:max.k, nrow(lambda_out_array)), 
                    nrow=nrow(lambda_out_array), 
                    byrow=TRUE)
    Out_poss<-IN_poss- matrix(rep(Diff_sim[, names(Diff_sim)==as.character(datum)], length(0:max.k)),
                              byrow=FALSE, ncol=length(0:max.k))
    
    lambda_out_array[, names(lambda_out_array)==as.character(datum)]<-
      t(w_sim%*%t(In_sim[,as.character((datum-1):(datum-length(w_sim)))]))
    lambda_out_array[lambda_out_array[, names(lambda_out_array)==as.character(datum)]==0,
                     names(lambda_out_array)==as.character(datum)]<-0.000001
    Prob<-apply(IN_poss, 2, dpois, lambda=lambda_in_array[, names(lambda_in_array)==as.character(datum)])*
      apply(Out_poss, 2, dpois, lambda=lambda_out_array[, names(lambda_out_array)==as.character(datum)])
    
    
    Prob<-Prob/ifelse(rowSums(Prob, na.rm=TRUE)==0, 1, rowSums(Prob, na.rm=TRUE))
    In_sim[, names(In_sim)==as.character(datum)]<- apply(Prob, 1, sample, x=0:max.k, size=1, replace=TRUE)
    
    Out_sim[, names(Out_sim)==as.character(datum)]<-In_sim[, names(In_sim)==as.character(datum)]-
      Diff_sim[, names(Diff_sim)==as.character(datum)]
    
  }
  
  
  RDIFF_Data<-suppressMessages(gather(In_sim, key="date",
                                      value = "Incoming", -1)%>%
                                 full_join(gather(Out_sim, key="date",
                                                  value = "Outgoing", -1))%>%
                                 filter(date%in% unique(data_rdiff$date))%>%
                                 mutate(date=as.numeric(date))%>%
                                 full_join(data_rdiff))
  
  RDIFF_Data_all<-suppressMessages(gather(In_sim, key="date",
                                          value = "Incoming", -1)%>%
                                     full_join(gather(Out_sim, key="date",
                                                      value = "Outgoing", -1))%>%
                                     mutate(date=as.numeric(date))%>%
                                     full_join(data_rdiff))
  
  
  return(list("data"=RDIFF_Data, "data_all"=RDIFF_Data_all))
}





######### Jack Knife Simulation #######  

sim_data_jk<-function(w=median,
                      ss=123, 
                      deterministic=TRUE,
                      coef=Coef_in,
                      covariates=datalist_prep[,c("date", 
                                                  "districtId",
                                                  "M", "N")]){
  # set seed   
  set.seed(ss)
  covariates<-covariates%>%
    rbind(rep(covariates[covariates$date==min(covariates$date),], length(w))%>%
            mutate(date=date-rep(1:length(w), each=length(unique(covariates$districtId)))))%>%
    arrange(date, districtId)
  lambda.in<-exp(t(coef%*%t(cbind(1, covariates%>%dplyr::select(-date, -districtId)))))
  Incoming<-rpois(length(lambda.in), lambda=lambda.in)
  
  data_sim<-cbind(covariates, "Incoming"=Incoming)
  
  all_length<-suppressMessages(as.data.frame(cbind("total_sample"=
                                                     sample(1:length(w), 
                                                            size=sum(data_sim$Incoming), 
                                                            prob=w, replace=TRUE),
                                                   "date"=rep(data_sim$date, 
                                                              times=data_sim$Incoming),
                                                   "districtId"=rep(data_sim$districtId, 
                                                                    times=data_sim$Incoming)))%>%
                                 dplyr::group_by(date, districtId, total_sample))
  
  all_length<-sqldf("SELECT COUNT(*) AS count, date, districtId, total_sample FROM all_length GROUP BY total_sample, date, districtId")
  
  data_tot<-data_sim
  for(i in 1:length(w)){
    lag<-all_length[all_length$total_sample==i,]%>%
      dplyr::mutate(date=date+i)%>%
      dplyr::select(date, count, districtId)
    names(lag)<-c("date", paste0("Incom_lag",i), "districtId")
    data_tot<- suppressMessages( data_tot%>%
                                   full_join(lag, by=c("date", "districtId"))%>%
                                   na.replace(replace=0))
  }
  data_tot<-data_tot%>%
    filter(date>=min(data_sim$date)+length(w)&
             date<=max(data_sim$date))%>%
    arrange(districtId, date)
  data_tot$Outgoing<-rowSums(data_tot[,grepl("lag", names(data_tot))])
  data_tot$Diff<-data_tot$Incoming-data_tot$Outgoing
  
  
  return(data_tot)
}


bias_adj<-function(w_og=omega_org, w_jk=omega_jk, max_l=max_lag, func=c("none", "sq", "lin", "exp", "abs")){
  if(!func%in%c("none", "sq", "lin", "exp", "abs")){
    print("Please define valid bias corr function")
  }
  ml<-max(length(w_og), length(w_jk), max_l)
  w_og<-c(w_og, rep(0,length=ml-length(w_og)))
  w_jk<-c(w_jk, rep(0,length=ml-length(w_jk)))
  
  if(func==("none")){
    c_hat<-0
    w_adj=w_og
  }
  if(func==("sq")){
    c_hat<-coef(lm(I((w_og-1/ml)^2)~I((w_jk-1/ml)^2)-1))
    w_adj<-ifelse(w_og<1/ml, -1, 1)*sqrt(c_hat*(w_og-1/ml)^2)+1/ml
  }
  if(func==("lin")){
    c_hat<-coef(lm(I((w_og-1/ml))~I((w_jk-1/ml))-1))
    w_adj<-c_hat*(w_og-1/ml)+1/ml
  }
  if(func==("exp")){
    c_hat<-coef(lm(I(log(w_og*ml+0.0001))~I(log(w_jk*ml+0.0001))-1))
    w_adj<-exp(c_hat*log(w_og*ml+0.0001))/ml
  }
  if(func==("abs")){
    c_hat<-coef(lm(I(abs(w_og-1/ml))~I(abs(w_jk-1/ml))-1))
    w_adj<-ifelse(w_og<1/ml, -1, 1)*sqrt((c_hat*abs(w_og-1/ml))^2)+1/ml
  }
  
  if(any(w_adj<0)|sum(w_adj)!=1){
    w_adj[w_adj<0]<-0
    w_adj<-w_adj/sum(w_adj)
  }
  return(list("w_adj"=w_adj, "c_hat"=c_hat))
}

w_true_Grid<-matrix(0, 1, 10)
w_true_Grid[1,]<-exp(0.4*seq(10,1, by=-1))/sum(exp(0.4*seq(10,1, by=-1)))



######################################################
############### Simulating data list #################
######################################################

datalist<-list()
for(j in 311:311){
  for(i in 1:nrow(w_true_Grid)){
    datalist[[length(datalist)+1]]<-sim_data_offset(w=w_true_Grid[i,],
                                                    districtId_sim=1:200, 
                                                    date_sim=1:(200+ncol(w_true_Grid)),
                                                    ss=j, 
                                                    deterministic=TRUE,
                                                    with_coef=TRUE)
    print(paste("Iteration", length(datalist), "done"))
    
  }
  
}

######################################################
####### Simulating close to real data for EM #########
######################################################

datalist_real<-list()
for(j in 1:length(datalist)){
  datalist_real[[j]]<-datalist[[j]][,c("date", "districtId", "Diff", "M", "N")]
}



######################################################
############## EM Algorithm ##########################
######################################################

l<-"sq"
omegasList_og<-list()
omegasList_jk<-list()
omegasList_adj<-list()

train_data_l=datalist_real[[1]] 
train_data_l_real=datalist[[1]]
runs1=400

formula_in=as.formula(Incoming ~ M+N) 

max_lag=12
j_k_data<-list()

train_data_l=datalist_real[[1]]
train_data_l_real=datalist[[1]]

omegasList_og<-list()
omegasList_jk<-list()
omegasList_adj<-list()

t<-1
train_data_l=datalist_real[[t]] 
train_data_l_real=datalist[[t]]

# save(train_data_l, 
#      file=paste0("C:/Users/ra98jiq/Documents/Papers/Skellam/Skellam/UpdateGoran/Jack Knife Inner/DataSim_Real_",
#                  Sys.Date(), "_", t, ".RData"))
# 
# save(train_data_l_real, 
#      file=paste0("C:/Users/ra98jiq/Documents/Papers/Skellam/Skellam/UpdateGoran/Jack Knife Inner/DataSim",
#                  Sys.Date(), "_", t, ".RData"))

train_data_l<-train_data_l%>%
  arrange(districtId, date)%>%
  mutate(Incoming=sample(10:15,1),
         Outgoing=sample(16:20,1))

I_Mstep<-mgcv::gam(formula_in, data=train_data_l, 
                   family="poisson")

w_start<-c(exp(-seq(0.1, 1.2, by=0.1))/sum(exp(-seq(0.1, 1.2, by=0.1))))


## Step 1 ######## Outer E Step ######################################### 
train_data_l_start<-rdiffpois_array_paper(data_rdiff=train_data_l%>%
                                            dplyr::select(-c("Incoming", "Outgoing")),
                                          IMSTEP=I_Mstep,
                                          w_sim=w_start)

######### Objects in which we store our results ############################

fitlist<-list()
modellist<-list()

Coeff_in<-NULL
Coeff_out<-NULL

loglike<-c()

omegas_org<-matrix(0, nrow=runs1, ncol=max_lag)
omegas_jk<-matrix(0, nrow=runs1, ncol=max_lag)
omegas_adj<-matrix(0, nrow=runs1, ncol=max_lag)

checkplot<-list()

VCOV<-array(dim=c(runs1, length(w_start)-1, length(w_start)-1))

C_hat<-vector(length=runs1)


########## EM-Loop #########################################################

for(z in 1:runs1){
  
  ## Step 1 ######## Outer E Step ######################################### 
  train_data_l_start<-rdiffpois_array_paper(data_rdiff=train_data_l%>%
                                              dplyr::select(-c("Incoming", "Outgoing")),
                                            IMSTEP=I_Mstep,
                                            w_sim=w_start)
  
  
  ## Step 2 ######## Outgoing M Step #########################################   
  
  ##### Seesaw to find best omega ####
  for(i in 1:100){
    new_sol<-optimize_Like(w_opt_prev=w_start, 
                           data_like_opp=train_data_l_start$data_all)
    
    if(sum(((c(new_sol$w_sol, 1-sum(new_sol$w_sol))-w_start)^2))<0.000000001){
      VCOV_w<-solve(new_sol$Infor)
      break
    }else{
      w_start<-c(new_sol$w_sol, 1-sum(new_sol$w_sol))      
    }
  }
  if(i==100){print("We have maxed out the runs to find omega")}
  
  # w_start<-omega_org<-round(w_start, 10)
  
  if(length(omega)>max_lag){
    print("The omega is longer than the maximum length 
                  of stay you have chosen")}
  
  
  
  ##### Incoming intensity parameter ####
  
  sim_start<-train_data_l_start$data[,c("districtId", "date", 
                                        "Incoming", 
                                        "Outgoing", 
                                        "Diff", "M", "N")]
  
  
  # ## Step 3 ######## Jack Knife Adjustment  ##################################
  # 
  # Jack_Knife_data<-sim_data_jk(w=omega_org,
  #                              ss=123, 
  #                              deterministic=TRUE,
  #                              coef=coef(I_Mstep),
  #                              covariates=sim_start[,
  #                                                   c("date", "districtId","M", "N")])
  # 
  # train_JK_temp<-Jack_Knife_data[,
  #                                c("date", "districtId","M", "N", "Diff")]
  # 
  # ## Step 4 ######## Jack Knife Adjustment Simulation  #######################
  # 
  # start_sim_jk<-rdiffpois_array_paper(data_rdiff=train_JK_temp,
  #                                     IMSTEP=I_Mstep,
  #                                     w_sim=omega_org)
  # 
  # ## Step 5 ######## Jack Knife Adjustment Omega  ############################
  # w_start_jk<-omega_org
  # ##### Seesaw to find best omega for JK ####
  # for(i in 1:100){
  #   new_sol_jk<-optimize_Like(w_opt_prev=w_start_jk,
  #                             data_like_opp=start_sim_jk$data_all)
  #   
  #   if(sum(((round(c(new_sol_jk$w_sol, 
  #                    1-sum(new_sol_jk$w_sol)), 10)-w_start_jk)^2))<0.0001){
  #     VCOV_w<-solve(new_sol_jk$Infor)
  #     break
  #   }else{
  #     w_start_jk<-round(c(new_sol_jk$w_sol, 
  #                         1-sum(new_sol_jk$w_sol)), 10)
  #   }
  # }
  # if(i==100){print("We have maxed out the runs to find omega")}
  # 
  # omega_jk<-c(w_start_jk, rep(0, max_lag-length(w_start_jk)))
  # 
  # ## Step 6 ######## Jack Knife Adjustment Omega  ############################
  # 
  # biasadj<-bias_adj(w_og=omega_org,
  #                   w_jk=omega_jk,
  #                   max_l=max_lag,
  #                   func="sq")
  # 
  # 
  # C_hat[z]<-coef(lm(I((omega_org-1/max_lag)^2)~I((omega_jk-1/max_lag)^2)-1))
  # omega_adj<-1/max_lag+ifelse(omega_org<1/max_lag, -1, 1)*sqrt(C_hat[z]*(omega_org-1/max_lag)^2)
  # omega_adj[omega_adj<0]=0
  # omega_adj<-omega_adj/sum(omega_adj)
  # ##########################
  # 
  # w_start<-omega_adj   # AENDERN ZU OMEGA ADJ MIT JK ADJUSTMENT 
  # ##########################
  ## Step 7 ######## Simulation Inner E-Step Omega  ############################
  # 
  # Inner_E_Step<-rdiffpois_array_paper(data_rdiff=sim_start%>%
  #                                       dplyr::select(-c("Incoming", "Outgoing")),
  #                                     IMSTEP=I_Mstep,
  #                                     w_sim=w_start)
  # 
  # ## Step 8 ######## M-Step for adjusted Simulated Data  #####################
  # 
  # sim_start<-cbind(Inner_E_Step$data[,c("districtId", "date", "Incoming",
  #                                       "Outgoing", "Diff", "M", "N")])


  I_Mstep<-mgcv::gam(formula_in, data=sim_start,
                     family="poisson")
  # 
  # 
  # ############################################################################
  
  fitlist[[z]]<-sim_start
  
  
  ##### Data set ####  
  
  train_data_l<-train_data_l_start<-fitlist[[z]]
  
  
  ############################################################################
  
  ##### Calculate Likelihood ####  
  loglike[z]<-sum(log(pskellam(q=train_data_l_start$Diff, 
                               lambda1=train_data_l_start$Incoming, 
                               lambda2=train_data_l_start$Outgoing)))
  
  ##### Save incoming model and omega ####  
  
  modellist[[z]]<-I_Mstep#list("I_mod"=I_Mstep,"R_mod"= Msteprun$R_Mstep)
  
  #omegas_adj[z,]<-omega_adj
  omegas_org[z,]<-c(w_start)
  #omegas_jk[z,]<-omega_jk
  #C_hat[z]<-c_hat
  VCOV[z,,]<-new_sol$Infor
  
  if(z!=1){
    Coeff_in<-rbind(Coeff_in, as.data.frame(cbind(z, t(coef(I_Mstep)))))
    
  }else{
    Coeff_in<-as.data.frame(cbind(z, t(coef(I_Mstep))))
    
  }
  print(paste("Inner iteration", z, "out of", runs1, "is done"))
  
}







varianceRuns<-200:400

coefratemed<-apply(Coeff_in[varianceRuns, 2:4], 2, median)
VCOV_Coef_Runs_b<-array(0, dim=c(length(varianceRuns), 
                                 ncol=ncol(vcov(modellist[[1]])), 
                                 nrow=nrow(vcov(modellist[[1]]))))
VCOV_Coef_Runs_s<-array(0, dim=c(length(varianceRuns), 
                                 ncol=ncol(vcov(modellist[[1]])), 
                                 nrow=nrow(vcov(modellist[[1]]))))

VCOV_Coef_Runs_w_b<-array(0, dim=c(length(varianceRuns), 11, 11))
VCOV_Coef_Runs_w_s<-array(0, dim=c(length(varianceRuns), 11, 11))

exitratemed<-apply(omegas_org[varianceRuns,], 2, median)


for(i in varianceRuns){
  VCOV_Coef_Runs_b[(i)-min(varianceRuns)+1,,]<-(as.numeric(Coeff_in[i, 2:4]-coefratemed))%*%
    t(as.numeric((Coeff_in[i, 2:4]-coefratemed)))
  VCOV_Coef_Runs_s[(i)-min(varianceRuns)+1,,]<-
    vcov(modellist[[i]])
  VCOV_Coef_Runs_w_b[(i)-min(varianceRuns)+1,,]<-((omegas_org[i,]-exitratemed)%*%t(omegas_org[i,]-exitratemed))[1:11, 1:11]
  VCOV_Coef_Runs_w_s[(i)-min(varianceRuns)+1,,]<-solve(VCOV[i,,])
}



stddev<-sqrt(diag(apply(VCOV_Coef_Runs_w_s, c(2,3), mean)+(1+(1/200))*apply(VCOV_Coef_Runs_w_b, c(2,3), sum)/(200-1)))

exitrate_med<-apply(omegas_org, 2, median)






Prediction<-fitlist[[z-1]]%>%
  mutate(Incoming_pred=Incoming, Outgoing_pred=Outgoing)%>%
  dplyr::select(-c("Incoming", "Outgoing"))%>%
  full_join(datalist[[1]]%>%
              dplyr::select(names(fitlist[[z-1]])))
ggplot()+
  geom_point(aes(x=Prediction$Incoming_pred,y=Prediction$Incoming, col="Incoming"))+
  geom_line(aes(x=Prediction$Incoming_pred,y=Prediction$Incoming_pred , col="Perf Fit"))+
  xlab("Estimation")+ylab("True")
ggplot()+
  geom_point(aes(x=Prediction$Outgoing_pred,y=Prediction$Outgoing, col="Outgoing"))+
  geom_line(aes(x=Prediction$Outgoing_pred,y=Prediction$Outgoing_pred , col="Perf Fit"))+
  xlab("Estimation")+ylab("True")



varianceRuns<-200:400

exitratemed<-apply(omegas_adj[varianceRuns,], 2, median)
VCOV_Exit_Runs<-array(0, dim=c(length(varianceRuns), 
                               ncol=ncol(VCOV[1,,]), 
                               nrow=nrow(VCOV[1,,])))
VCOV_Exit_Runs_sum<-array(0, dim=c(length(varianceRuns), 
                                   ncol=ncol(VCOV[1,,]), 
                                   nrow=nrow(VCOV[1,,])))
for(i in varianceRuns){
  VCOV_Exit_Runs[(i)-min(varianceRuns)+1,,]<-solve(VCOV[i,,])
  VCOV_Exit_Runs_sum[(i)-min(varianceRuns)+1,,]<-((omegas_adj[i,]-exitratemed)%*%t(
    omegas_adj[i,]-exitratemed))[1:11, 1:11]
}



# stddevadj<-sqrt(diag(apply(VCOV_Exit_Runs, c(2,3), mean)+((1+1/length(varianceRuns))/(length(varianceRuns)-1))*apply(VCOV_Exit_Runs_sum, c(2,3), sum)))
exitrate_med<-apply(omegas_adj, 2, median)

Exit_rates<-ggplot()+
  geom_segment(aes(x = 1, y = 1/12, xend = 1, yend = w_true_Grid[1,1], col="Pull"))+
  geom_segment(aes(x = 2, y = 1/12, xend = 2, yend = w_true_Grid[1,2], col="Pull"))+
  geom_segment(aes(x = 3, y = 1/12, xend = 3, yend = w_true_Grid[1,3], col="Pull"))+
  geom_segment(aes(x = 4, y = 1/12, xend = 4, yend = w_true_Grid[1,4], col="Pull"))+
  geom_segment(aes(x = 5, y = 1/12, xend = 5, yend = w_true_Grid[1,5], col="Pull"))+
  geom_segment(aes(x = 6, y = 1/12, xend = 6, yend = w_true_Grid[1,6], col="Pull"))+
  geom_segment(aes(x = 7, y = 1/12, xend = 7, yend = w_true_Grid[1,7], col="Pull"))+
  geom_segment(aes(x = 8, y = 1/12, xend = 8, yend = w_true_Grid[1,8], col="Pull"))+
  geom_segment(aes(x = 9, y = 1/12, xend = 9, yend = w_true_Grid[1,9], col="Pull"))+
  geom_segment(aes(x = 10, y = 1/12, xend = 10, yend = w_true_Grid[1,10], col="Pull"))+
  geom_segment(aes(x = 11, y = 1/12, xend = 11, yend = 0, col="Pull"))+
 geom_line(aes(x=1:11, y=1/12, col="1/12"), linetype=2)+
  # geom_hline(aes(yintercept=1/(max_lag), col=paste0("1/", max_lag)))+
  # geom_ribbon(aes(x=1:(max_lag-1), ymin=apply(omegas_org[ varianceRuns,], 
  #                                             2, median)[-12]-1.96*stddev[-12], 
  #                 ymax=apply(omegas_org[ varianceRuns,], 
  #                            2, median)[-12]+1.96*stddev[-12], 
  #                 fill="95% Confidence interval"))+
  # geom_line(aes(x=1:(max_lag-1), 
  #                 y=apply(omegas_adj[ varianceRuns,], 
  #                            2, median)[-12]-1.96*stddevadj[-12]))+ 
  # geom_line(aes(x=1:(max_lag-1),
  #               y=apply(omegas_adj[ varianceRuns,], 
  #                            2, median)[-12]+1.96*stddevadj[-12]))+
  geom_point(aes(x=1:(max_lag-1), y=(apply(omegas_org[ varianceRuns,], 
                                          2, median))[-12],
                col="Exit Rate"), shape=15)+
  geom_line(aes(x=1:(max_lag-1), y=(apply(omegas_org[ varianceRuns,], 
                                           2, median))[-12],
                 col="Exit Rate"))+
  theme_pubr()+
  labs(color=" ")+
  ggtitle("Exit rate")+
  geom_line(aes(x=1:11, y=c(w_true_Grid[1,], 0), col="Ground Truth"))+
  geom_point(aes(x=1:11, y=c(w_true_Grid[1,], 0), col="Ground Truth"))+
  # paste0("Median adjustment scale ", round(median(C_hat[360:400]), 
  #                                          2)))+
  ylab(expression(hat(omega)))+ 
  scale_colour_manual(values = brewer.pal(n = 5, name = "Blues")[3:6])+
  scale_fill_manual(values = brewer.pal(n = 5, name = "Blues")[3:6])+
  xlab("Lag")+
  scale_x_continuous(breaks = seq(1, max_lag, 1))+
  theme(legend.title = element_blank())



Exit_rates

ggsave(Exit_rates, file="C:/Users/ra98jiq/Downloads/exitrateNotAdjNoStd.pdf", width =5, height=3 )


True_rates<-ggplot()+
  geom_hline(aes(yintercept=1/(max_lag), col=paste0("1/", max_lag)))+
  # geom_line(aes(x=1:(max_lag-1), 
  #                 y=apply(omegas_adj[ varianceRuns,], 
  #                            2, median)[-12]-1.96*stddevadj[-12]))+ 
  # geom_line(aes(x=1:(max_lag-1),
  #               y=apply(omegas_adj[ varianceRuns,], 
  #                            2, median)[-12]+1.96*stddevadj[-12]))+
  geom_line(aes(x=1:(max_lag-1), y=c(w_true_Grid[1,], 0),
                col="Ground Truth"))+
  geom_point(aes(x=1:(max_lag-1), y=c(w_true_Grid[1,], 0),
                 col="Ground Truth"))+
  theme_pubr()+
  labs(color=" ")+
  ggtitle("Bias adjusted exit rate")+
  # paste0("Median adjustment scale ", round(median(C_hat[360:400]), 
  #                                          2)))+
  ylab("Value of exit rate")+ 
  scale_colour_manual(values = brewer.pal(n = 5, name = "Blues")[3:5])+
  scale_fill_manual(values = brewer.pal(n = 5, name = "Blues")[2:5])+
  xlab("Lag")+
  scale_x_continuous(breaks = seq(1, max_lag, 1))+
  theme(legend.title = element_blank())


# stddev<-sqrt(diag(solve(apply(VCOV[varianceRuns,,], c(2,3), mean))))
# Exit_rates<-ggplot()+
#   geom_hline(aes(yintercept=1/(max_lag), col=paste0("1/", max_lag)))+
#   geom_ribbon(aes(x=1:(max_lag-1), ymin=apply(omegas_adj[ varianceRuns,], 
#                                               2, median)[-12]-1.96*stddev[-12], 
#                   ymax=apply(omegas_adj[ varianceRuns,], 
#                              2, median)[-12]+1.96*stddev[-12], 
#                   fill="95% Confidence interval"))+
#   geom_line(aes(x=1:(max_lag-1), 
#                 y=apply(omegas_adj[ varianceRuns,], 
#                         2, median)[-12]-1.96*stddev[-12]))+ 
#   geom_line(aes(x=1:(max_lag-1),
#                 y=apply(omegas_adj[ varianceRuns,], 
#                         2, median)[-12]+1.96*stddev[-12]))+
#   geom_line(aes(x=1:(max_lag-1), y=(apply(omegas_adj[ varianceRuns,], 
#                                           2, median))[-12],
#                 col="Exit Rate"))+
#   theme_pubr()+
#   labs(color=" ")+
#   ggtitle("Bias adjusted exit rate")+
#   geom_point(aes(x=1:11, y=c(w_true_Grid[1,], 0), col="Ground Truth"))+
#   # paste0("Median adjustment scale ", round(median(C_hat[360:400]), 
#   #                                          2)))+
#   ylab("Value of exit rate")+ 
#   scale_colour_manual(values = brewer.pal(n = 5, name = "Blues")[3:5])+
#   scale_fill_manual(values = brewer.pal(n = 5, name = "Blues")[2:5])+
#   xlab("Lag")+
#   scale_x_continuous(breaks = seq(1, max_lag, 1))+
#   theme(legend.title = element_blank())
# 
# Exit_rates


chat_plot<-ggplot()+
  geom_line(aes(x= 1:400, y=sqrt(C_hat)[1:400], col="Estimated bias"), size=1.1)+
  geom_hline(aes(yintercept=1, col="No bias"))+
  coord_cartesian(ylim = c(0,1.3), xlim = c(200,400))+
  theme_pubr()+
  labs(color=" ")+
  ggtitle("b) Estimate of pull towards uniform distribution", "Over runs 200 to 400")+
  xlab("Runs")+ylab(expression(sqrt(hat(c))))+
  scale_colour_manual(values = brewer.pal(n = 5, name = "Blues")[c(5,2)])
ggsave(plot_grid(Exit_rates, chat_plot), file="C:/Users/ra98jiq/Downloads/EstimatedBias.pdf",
       width=13, height=4)



Incoming_mat<-matrix(0, 
                     nrow=nrow(fitlist[[1]]), 
                     ncol=length(round(0.5*length(fitlist)):length(fitlist)))
Outgoing_mat<-matrix(0, nrow=nrow(fitlist[[1]]), 
                     ncol=length(round(0.5*length(fitlist)):length(fitlist)))

for(i in round(0.5*length(fitlist)):length(fitlist)){
  Incoming_mat[,i-round(0.5*length(fitlist))+1]<-fitlist[[i]]$Incoming
  Outgoing_mat[,i-round(0.5*length(fitlist))+1]<-fitlist[[i]]$Outgoing
}


Inc_Out_Data<-as.data.frame(cbind("date"=fitlist[[i]]$date,
                                  "districtId"=fitlist[[i]]$districtId, 
                                  "Incoming_Est"=apply(Incoming_mat, 1, median), 
                                  "Outgoing_Est"=apply(Outgoing_mat, 1, median)))%>%
  full_join(datalist[[t]][, c("date", "districtId", "Incoming", "Outgoing")])

Max_tot<-max(Inc_Out_Data[, 3:6])
IncomTrue<-ggplot()+
  geom_line(aes(x=1:Max_tot, y=1:Max_tot, col="Perfect fit"), size=1.1)+
  geom_point(aes(y=Inc_Out_Data$Incoming_Est, 
                 x=Inc_Out_Data$Incoming, col="Incoming"), size=1.1, alpha=0.5)+
  theme_pubr()+
  labs(color=" ")+
  scale_colour_manual(values = brewer.pal(n = 5, name = "Blues")[4:6])+ 
  ggtitle("Incoming")+
  ylab("Estimated Incoming")+ 
  xlab("True Incoming") 

OutgoTrue<-ggplot()+
  geom_line(aes(x=1:Max_tot, y=1:Max_tot, col="Perfect fit"), 
            size=1.1)+
  geom_point(aes(y=Inc_Out_Data$Outgoing_Est, 
                 x=Inc_Out_Data$Outgoing, col="Outgoing"), 
             size=1.1, alpha=0.5)+
  # # geom_point(aes(y=Inc_Out_Data$Outgoing_Est[Inc_Out_Data$Incoming>20], 
  # #                x=Inc_Out_Data$Outgoing[Inc_Out_Data$Incoming>20]), col="grey", 
  # #            size=1.5)+
  # geom_point(aes(y=Inc_Out_Data$Outgoing_Est[Inc_Out_Data$Incoming>30], 
  #                x=Inc_Out_Data$Outgoing[Inc_Out_Data$Incoming>30]), col="black", 
  #            size=1.5)+
  theme_pubr()+
  scale_colour_manual(values = brewer.pal(n = 5, name = "Blues")[4:6])+ 
  labs(color=" ")+
  ggtitle("Outgoing")+
  ylab("Estimated Outgoing")+ 
  xlab("True Outgoing") 


Est_True<-plot_grid(IncomTrue, OutgoTrue, nrow=1)
ggsave(Est_True, file="C:/Users/ra98jiq/Downloads/InOutEst.pdf",
       width=8, height=4)

Exit_rates_gt<-ggplot()+
  geom_hline(aes(yintercept=1/(max_lag), col=paste0("1/", max_lag)))+
  # geom_ribbon(aes(x=1:(max_lag-1), ymin=apply(omegas_adj[ varianceRuns,], 
  #                                             2, median)[-12]-1.96*stddevadj[-12], 
  #                 ymax=apply(omegas_adj[ varianceRuns,], 
  #                            2, median)[-12]+1.96*stddevadj[-12], 
  #                 fill="95% Confidence interval"))+
  # geom_line(aes(x=1:(max_lag-1), 
  #                 y=apply(omegas_adj[ varianceRuns,], 
  #                            2, median)[-12]-1.96*stddevadj[-12]))+ 
  # geom_line(aes(x=1:(max_lag-1),
  #               y=apply(omegas_adj[ varianceRuns,], 
  #                            2, median)[-12]+1.96*stddevadj[-12]))+
# geom_line(aes(x=1:(max_lag-1), y=(apply(omegas_adj[ varianceRuns,], 
#                                         2, median))[-12],
#               col="Exit Rate"))+
theme_pubr()+
  labs(color=" ")+
  ggtitle("Ground truth exit rate")+
  geom_point(aes(x=1:11, y=c(w_true_Grid[1,], 0), col="Ground Truth"))+
  # paste0("Median adjustment scale ", round(median(C_hat[360:400]), 
  #                                          2)))+
  ylab("Value of exit rate")+ 
  scale_colour_manual(values = brewer.pal(n = 5, name = "Blues")[3:5])+
  scale_fill_manual(values = brewer.pal(n = 5, name = "Blues")[2:5])+
  xlab("Lag")+
  scale_x_continuous(breaks = seq(1, max_lag, 1))+
  theme(legend.title = element_blank())

ggsave(Exit_rates_gt, file="C:/Users/ra98jiq/Downloads/ExitRateGroundTruth.pdf",
       width=5, height=3)
ggsave(Exit_rates, file="C:/Users/ra98jiq/Downloads/ExitRateEst.pdf",
       width=5, height=3)



