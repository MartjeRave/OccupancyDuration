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
#install.packages("rgdal")
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
      dplyr::full_join(data_temp, by=c("date", "districtId"))
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
                                   dplyr::group_by(date, districtId, total_sample)%>%
                                   dplyr::summarize(count=n()))
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
  
  # model_score<-mgcv::gam(Outgoing~1,
  #                        offset = log(Sum_score+0.1),
  #                        data=data_score_fun,
  #                        family=poisson())
  # 
  # terms_score<-as.vector(predict.gam(model_score, 
  #                                    type="response"))
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
                        (I_ts_score-I_tN_score))
    
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
                         data_off=data_like_opp[, c("districtId", "date", "Incoming",
                                                    "Outgoing")])
  
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


optimize_Like_wcheck<-function(w_opt_prev=w_start_adj, 
                        data_like_opp=datalist[[l]]){
  
  Off_data_temp<-Off_cal(w_hat=w_opt_prev, 
                         data_off=data_like_opp[, c("districtId", "date", "Incoming",
                                                    "Outgoing")])
  
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
  
  Sum_i<-t(w_opt_prev%*%t(Off_data_temp[, (grepl("Incoming_", names(Off_data_temp)))]))
  Delta_11<-((Off_data_temp[,(grepl("_1_", names(Off_data_temp))&
                               (grepl("Incom", names(Off_data_temp))))]-
               Off_data_temp[,(grepl(paste0("_", length(w_opt_prev),"_"), names(Off_data_temp))&
                                 (grepl("Incom", names(Off_data_temp))))])/ifelse(Sum_i==0, 
                                                                                  Sum_i+0.000001, 
                                                                                  Sum_i))
  for(i in 2:(length(w_opt_prev)-1)){
    Delta_11<-cbind(Delta_11, (Off_data_temp[,(grepl(paste0("_", i,"_"), names(Off_data_temp))&
                      (grepl("Incom", names(Off_data_temp))))]-
    Off_data_temp[,(grepl(paste0("_", length(w_opt_prev),"_"), names(Off_data_temp))&
                      (grepl("Incom", names(Off_data_temp))))])/ifelse(Sum_i==0, 
                                                                       Sum_i+0.000001, 
                                                                       Sum_i))    
  }
  colnames(Delta_11)<-paste0("lag_", 1:(length(w_start)-1))

  
  return(list("w_sol"=Opt$solution, "Infor"=Dmat, "w_sol_uncon"=Opt$unconstrained.solution,
              "Delta_11"=Delta_11))
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
                                max.k=1000,
                                rep_run=TRUE){
  
  
  set.seed(123)
  data_rdiff=as.data.frame(data_rdiff)
  
  if(rep_run){
    Lagged_Data<-rep(data_rdiff[data_rdiff$date==min(data_rdiff$date),], (length(w_sim)+1))%>%
      mutate(date=date-rep(1:(length(w_sim)+1), each=length(unique(data_rdiff$districtId))))%>%
      rbind(data_rdiff)%>%
      arrange(date, districtId)  
  }else{
    Lagged_Data=as.data.frame(data_rdiff)
  }
  
  
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
  
  
  
  for(datum in unique(data_rdiff$date)){
    IN_poss<-matrix(rep(0:max.k, nrow(lambda_out_array)), 
                    nrow=nrow(lambda_out_array), 
                    byrow=TRUE)
    
    Out_poss<-IN_poss- matrix(rep(as.matrix(Diff_sim[, 
                              names(Diff_sim)==as.character(as.Date(datum))]), length(0:max.k)),
                              byrow=FALSE, ncol=length(0:max.k))
    
    lambda_out_array[, names(lambda_out_array)==as.character(as.Date(datum))]<-
      t(w_sim%*%t(In_sim[,as.character(as.Date((datum-1):(datum-length(w_sim))))]))
    lambda_out_array[lambda_out_array[, names(lambda_out_array)==as.character(as.Date(datum))]==0,
                     names(lambda_out_array)==as.character(as.Date(datum))]<-0.000001
    Prob<-apply(IN_poss, 2, dpois, lambda=lambda_in_array[, 
                 names(lambda_in_array)==as.character(as.Date(datum))])*
      apply(Out_poss, 2, dpois, lambda=lambda_out_array[, 
            names(lambda_out_array)==as.character(as.Date(datum))])
    
    Prob<-Prob/ifelse(rowSums(Prob, na.rm=TRUE)==0, 1, rowSums(Prob, na.rm=TRUE))
    In_sim[, names(In_sim)==as.character(as.Date(datum))]<- apply(Prob, 1, sample, x=0:max.k, size=1, replace=TRUE)
    
    Out_sim[, names(Out_sim)==as.character(as.Date(datum))]<-In_sim[, names(In_sim)==as.character(as.Date(datum))]-
      Diff_sim[, names(Diff_sim)==as.character(as.Date(datum))]
    
  }
  
  
  RDIFF_Data<-suppressMessages((gather(In_sim, key="date",
                                       value = "Incoming", -1)%>%
                                  mutate(date=as.Date(date)))%>%
                                 full_join((gather(Out_sim, key="date",
                                                   value = "Outgoing", -1)%>%
                                              mutate(date=as.Date(date))))%>%
                                 filter(date%in% unique(data_rdiff$date))%>%
                                 mutate(date=as.Date(date))%>%
                                 full_join(data_rdiff))
  
  
  RDIFF_Data_all<-suppressMessages((gather(In_sim, key="date",
                                           value = "Incoming", -1)%>%
                                      mutate(date=as.Date(date)))%>%
                                     full_join((gather(Out_sim, key="date",
                                                       value = "Outgoing", -1)%>%
                                                  mutate(date=as.Date(date))))%>%
                                     mutate(date=as.Date(date))%>%
                                     full_join(data_rdiff))
  
  return(list("data"=RDIFF_Data, "data_all"=RDIFF_Data_all))
}




######### Jack Knife Simulation #######

# ############# Coefficients Median  #############################################
# 
# VCOV_coef<-lapply(modellist[200:400], vcov)
# VCOV_coef_arr<-array(0 , dim=c(length(1:400), nrow(VCOV_coef[[1]]), ncol(VCOV_coef[[1]])))
# 
# for(i in 1:200){
#   VCOV_coef_arr[i,,]<-as.matrix(VCOV_coef[[i]])
# }
# 
# VCOV_coef_med<-apply(VCOV_coef_arr, c(2,3), median)
# 
# 
# coef_fix<-cbind(names(coef(modellist[[200]])[1:10]),
#       round(apply(matrix(unlist(lapply(modellist[360:400], coef)), 
#                    byrow=TRUE, ncol=length(coef(modellist[[200]]))), 
#             2,median)[1:10], 3),
#       round(sqrt(diag(VCOV_coef_med))[1:10], 3))
# 
# colnames(coef_fix)<-c("Covariates", "Coefficients", "Std Dev")
#
#coef_fix

sim_data_jk<-function(w=median,
                      ss=123,
                      deterministic=TRUE,
                      coef=Coef_in,
                      Xmat_in_data=train_data_l_start$data%>%
                        dplyr::select(-c("Incoming", "Outgoing")),
                      model_in=modellist[[1]]){
  
  # set seed
  set.seed(ss)
  
  Xmat_newdata<-rbind(rep(Xmat_in_data[Xmat_in_data$date==min(Xmat_in_data$date),],
                          length(w))%>%
                        mutate(date=date-rep(1:length(w),each=length(unique(Xmat_in_data$districtId)))),
                      Xmat_in_data)
  
  lambda.in<-predict(model_in, newdata=Xmat_newdata, type="response")
  Incoming<-rpois(length(lambda.in), lambda=lambda.in)
  
  data_sim<-as.data.frame(cbind(Xmat_newdata, "Incoming"=Incoming))
  
  all_length<-suppressMessages(as.data.frame(cbind("total_sample"=
                                                     sample(1:length(w),
                                                            size=sum(data_sim$Incoming),
                                                            prob=w, replace=TRUE),
                                                   "date"=rep(data_sim$date,
                                                              times=data_sim$Incoming),
                                                   "districtId"=rep(data_sim$districtId,
                                                                    times=data_sim$Incoming)))%>%
                                 dplyr::group_by(date, districtId, total_sample)%>%
                                 dplyr::summarize(count=n()))%>%
    mutate(date=as.Date(as.numeric(date)), districtId=as.numeric(districtId))
  
  data_tot<-data_sim%>%
    mutate(date=as.Date(as.numeric(date)))
  for(i in 1:length(w)){
    lag<-all_length[all_length$total_sample==i,]%>%
      dplyr::mutate(date=date+i)%>%
      dplyr::select(date, count, districtId)%>%
      mutate(date=as.Date(as.numeric(date)))
    names(lag)<-c("date", paste0("Incom_lag",i), "districtId")
    data_tot<- suppressMessages( data_tot%>%
                                   mutate(date=as.Date(as.numeric(date)),
                                          districtId=as.numeric(districtId))%>%
                                   full_join(lag%>%
                                               mutate(date=as.Date(as.numeric(date)),
                                                      districtId=as.numeric(districtId)), by=c("date", "districtId"))%>%
                                   na.replace(replace=0))
  }
  data_tot<-data_tot%>%
    filter(date>=min(data_sim$date)+length(w)&
             date<=max(data_sim$date))%>%
    arrange(districtId, date)
  data_tot$Outgoing<-rowSums(data_tot[,grepl("Incom_lag", names(data_tot))])
  data_tot$Diff<-data_tot$Incoming-data_tot$Outgoing
  
  
  return(data_tot)
}


######################################################
############## Date ##################################
######################################################

train_data_l<-readRDS("training_data.rds")%>%
  mutate(date=as.Date(date))%>%
  filter(date<=as.Date("2021-12-31")&date>=as.Date("2021-08-01"))%>%
  mutate(Incoming=rpois(n=length(exp(0.3*G_5_7)), lambda=exp(0.3*G_5_7)), 
         Outgoing=rpois(n=length(exp(0.3*G_5_7)), lambda=exp(0.2*G_5_7)))

formula_in=as.formula(Incoming ~ G_4_7_lag_1+G_5_7_lag_1+G_6_7_lag_1+
                        weekday+s(long, lat)+s(as.numeric(date))) 

I_Mstep<-mgcv::gam(formula_in, data=train_data_l, 
                   family="poisson")



######################################################
############## EM Algorithm ##########################
######################################################

runs1=500
max_lag=30
### start value
w_start<-exp(-0.1*1:max_lag)/sum(exp(-0.1*1:max_lag))
### Save Results
omegas_org<-rbind(omegas_org, 
                  matrix(0, ncol=max_lag, nrow=runs1-400+1))
omegas_jk<-rbind(omegas_jk, matrix(0, ncol=max_lag, nrow=runs1-400+1))
omegas_adj<-rbind(omegas_adj, matrix(0, ncol=max_lag, nrow=runs1-400+1))
modellist<-list()
IncomOutgo<-list()
VCOV<-array(0, dim=c(runs1, length(w_start)-1, length(w_start)-1))
VCOV<-abind(VCOV,array(0, dim=c(runs1-400+1, length(w_start)-1, length(w_start)-1)), along=1)
logLike<-c()
c_hats<-c()
opti_list<-list()
re_model<-list()
z<-401
for(z in 401:runs1){
  
  #### Step 1 #### Outer E-Step
  Outer1<-rdiffpois_array_paper(data_rdiff=(train_data_l%>%
                                  dplyr::select(-c("Incoming", "Outgoing"))),
                        IMSTEP=I_Mstep,
                        w_sim=w_start,
                        max.k=1000)
  
  if(z/runs1<0.9){max_w_start=100}else{max_w_start=3}
  #### Step 2 #### Outgoing M-Step 
  for(i in 1:max_w_start){
    w_start_i<-optimize_Like(w_opt_prev=w_start, 
                data_like_opp=Outer1$data_all)
    if(sum((c(w_start_i$w_sol, 1-sum(w_start_i$w_sol))-w_start)^2)<0.00001){
      break
      }else{
      w_start<-c(w_start_i$w_sol, 1-sum(w_start_i$w_sol))
      VCOV[z,,]<-w_start_i$Infor
      }
  }
  print(paste("        OMEGA OG: Done with", i, "out of 100 runs"))
  if(any(w_start<0)){
    w_start<-round(w_start, 15)/sum(round(w_start, 15))
    w_start[w_start<0]<-0
    w_start<-w_start/sum(w_start)
    if(w_start[1]==0){
      w_start<-exp(-0.1*1:max_lag)/sum(exp(-0.1*1:max_lag))
      print("WARNING OMEGA HAD TO BE IMPUTED")
    }
    }
  
  VCOV[z,,]<-w_start_i$Infor
  
  if(z/runs1>0.7){

  #### Step 3 #### Outgoing bias correction
  InnerSim<-sim_data_jk(w=w_start,
              ss=123,
              deterministic=TRUE,
              coef=coef(I_Mstep),
              Xmat_in_data=Outer1$data%>%
                dplyr::select(-c("Incoming", "Outgoing")),
              model_in=I_Mstep)
  #### Step 4 #### Outgoing bias correction
  Inner4<-rdiffpois_array_paper(data_rdiff=(InnerSim%>%
                                dplyr::select(-c("Incoming", "Outgoing"))),
                                IMSTEP=I_Mstep,
                                w_sim=w_start,
                                max.k=1000)
  #### Step 5 #### Outgoing exit rate JK
  w_start_jk<-w_start
  for(i in 1:max_w_start){
    w_start_i<-optimize_Like(w_opt_prev=w_start_jk, 
                             data_like_opp=Inner4$data_all)
    if(sum((c(w_start_i$w_sol, 1-sum(w_start_i$w_sol))-w_start_jk)^2)<0.00001){
      w_start_jk<-c(w_start_i$w_sol, 1-sum(w_start_i$w_sol))
      break
    }else{
      w_start_jk<-c(w_start_i$w_sol, 1-sum(w_start_i$w_sol))
    }
  }
  print(paste("        OMEGA JK: Done with", i, "outof 100 runs"))
  
  #### Step 6 #### Bias correction
  c_hat<-coef(lm(I((w_start-1/max_lag)^2)~I((w_start_jk-1/max_lag)^2)-1))
  w_start_adj<-1/max_lag+ifelse(w_start<1/max_lag, -1, 1)*sqrt(c_hat*(w_start-1/max_lag)^2)
  w_start_adj[w_start_adj<0]<-0
  w_start_adj<-w_start_adj/sum(w_start_adj)
    
  
  c_hats[z]<-c_hat
  if(c_hat<1){
    w_start_jk<-apply(omegas_org[(z-11):(z-1),],2 ,median)/sum(apply(omegas_org[(z-11):(z-1),],2 ,median))
    w_start_adj<-apply(omegas_org[(z-11):(z-1),],2 ,median)/sum(apply(omegas_org[(z-11):(z-1),],2 ,median))  
  }
  omegas_jk[z,]<-w_start_jk
  omegas_adj[z,]<-w_start_adj  
  
  }else{
    w_start_jk<-w_start
    w_start_adj<-w_start
  }
  #### Step 7 ####
  Outer7<-rdiffpois_array_paper(data_rdiff=(train_data_l%>%
                                 dplyr::select(-c("Incoming", "Outgoing"))),
                                IMSTEP=I_Mstep,
                                w_sim=w_start_adj,
                                max.k=1000)
  
  #### Step 8 #### Incoming M-Step
  I_Mstep<-mgcv::gam(formula_in, data=Outer7$data, 
                     family="poisson")
  #### Step End ####
  omegas_org[z,]<-w_start
  modellist[[z]]<-I_Mstep
  
  
  
  logLike[z]<-sum(log(pskellam(q=Outer7$data$Diff, 
                               lambda1=Outer7$data$Incoming, 
                               lambda2=Outer7$data$Outgoing)))
  w_start<-w_start_adj
  Sum_z<-Off_cal(w_hat=w_start, 
          data_off=Outer7$data_all[, c("districtId", "date", "Incoming",
                                     "Outgoing")])
  
  IncomOutgo[[z]]<-suppressMessages(as.data.frame(Outer7$data[,c("districtId", "date", 
                                                                 "Incoming", "Outgoing")])%>%
    full_join(cbind(Sum_z[c("date", "districtId")], "Sum_w"=t(w_start%*%t(Sum_z[, 
                                                 grepl("Incoming_", names(Sum_z))])))))
  
  # re_model[[z]]<-glm(as.formula(paste0("Outgoing~-1+", 
  #                   paste0(names(Sum_z)[grepl("Incoming_", names(Sum_z))], collapse = "+ "))), 
  #                   data=Sum_z, family="poisson")
  # sqrt(diag(vcov(re_model[[z]])))
  print(paste("Done with", z, "outof", runs1, "runs"))
}






################################################################################
############# Results ##########################################################
################################################################################


############# Exit Rates ########## Adjusted ###################################

cumsum(apply(omegas_adj[360:400,], 2, median))

Exit_rates<-ggplot()+
  geom_hline(aes(yintercept=1/(max_lag), col="1/30"))+
  geom_ribbon(aes(x=1:(max_lag-1), ymin=apply(omegas_adj[360:400,], 
        2, median)[-30]/sum(apply(omegas_adj[360:400,], 
                             2, median))-1.96*sqrt(diag(solve(apply(VCOV[360:400,,], 
                        c(2,3), median))))[-30], ymax=apply(omegas_adj[360:400,], 
                         2, median)[-30]/sum(apply(omegas_adj[360:400,], 
                           2, median))+1.96*sqrt(diag(solve(apply(VCOV[360:400,,], 
                c(2,3), median))))[-30], 
        fill="95% Confidence interval"))+
    geom_line(aes(x=1:(max_lag-1), y=(apply(omegas_adj[360:400,], 
                                            2, median)/sum(apply(omegas_adj[360:400,], 
                                                                 2, median)))[-30],
                  col="Exit Rate"))+
  theme_pubr()+
  labs(color=" ")+
  ggtitle("Bias adjusted exit rate")+
          # paste0("Median adjustment scale ", round(median(C_hat[360:400]), 
          #                                          2)))+
  ylab("Value of exit rate")+ 
  scale_colour_manual(values = brewer.pal(n = 5, name = "Blues")[4:5])+
  scale_fill_manual(values = brewer.pal(n = 5, name = "Blues")[2:5])+
  xlab("Lag")+
  scale_x_continuous(breaks = seq(1, max_lag, 1))+
  theme(legend.title = element_blank())

Exit_rates




Exit_rates<-ggplot()+
  geom_hline(aes(yintercept=1/(max_lag), col="1/30"))+
  geom_ribbon(aes(x=1:(max_lag-1), ymin=apply(omegas_adj[360:400,], 
        2, median)[-30]/sum(apply(omegas_adj[360:400,], 
                             2, median))-1.96*sqrt(diag(solve(apply(VCOV[360:400,,], 
                        c(2,3), median))))[-30], ymax=apply(omegas_adj[360:400,], 
                         2, median)[-30]/sum(apply(omegas_adj[360:400,], 
                           2, median))+1.96*sqrt(diag(solve(apply(VCOV[360:400,,], 
                c(2,3), median))))[-30], 
        fill="95% Confidence interval"))+
    geom_line(aes(x=1:(max_lag-1), y=(apply(omegas_adj[360:400,], 
                                            2, median)/sum(apply(omegas_adj[360:400,], 
                                                                 2, median)))[-30],
                  col="Exit Rate"))+
  theme_pubr()+
  labs(color=" ")+
  ggtitle("Bias adjusted exit rate", 
          paste0("Median estimate pull towards 1/30, c=", round(median(sqrt(c_hats), na.rm=TRUE), 2)))+
          # paste0("Median adjustment scale ", round(median(C_hat[360:400]), 
          #                                          2)))+
  ylab(expression(hat(omega)))+ 
  scale_colour_manual(values = brewer.pal(n = 5, name = "Blues")[4:5])+
  scale_fill_manual(values = brewer.pal(n = 5, name = "Blues")[2:5])+
  xlab("Lag")+
  scale_x_continuous(breaks = seq(1, max_lag, 1))+
  theme(legend.title = element_blank())

ggsave(Exit_rates, file="C:/Users/ra98jiq/Downloads/CovidExitRate.pdf", width=9, height=3.5)


############# Estimated Incoming and Outgoing ##################################


InOutIt<-IncomOutgo[[5]][,c("date", "districtId", "Incoming", "Outgoing")]
names(InOutIt)<-c("date", "districtId",paste0(c("Incoming_", "Outgoing_"), 2))


for(i in 6:length(IncomOutgo)){
  InOutIt<-as.data.frame(InOutIt)%>%
    full_join(IncomOutgo[[i]][,c("date", "districtId", "Incoming", "Outgoing")])
  names(InOutIt)[(ncol(InOutIt)-1):ncol(InOutIt)]<-paste0(c("Incoming_", "Outgoing_"), 
                                                          i)
  
}


InOutIt_long<-InOutIt%>%
  gather(key="Run", value="Sim", -date, -districtId)
InOutIt_long$Run_num<-as.numeric(substr(InOutIt_long$Run, 10,12))


InOut_est<-cbind(InOutIt[,c(1,2)], "Incoming_est"=apply(InOutIt[,names(InOutIt)%in%paste0("Incoming_", 200:400)],
                                                        1, median))%>%
  full_join(cbind(InOutIt[,c(1,2)], "Outgoing_est"=apply(InOutIt[,names(InOutIt)%in%paste0("Outgoing_", 200:400)],
                                                         1, median)))%>%
  mutate(date=as.Date(date), districtId=as.numeric(districtId))%>%
  full_join(train_data_l[,c("date", "districtId", "Diff")]%>%
              mutate(date=as.Date(date), districtId=as.numeric(districtId)))%>%
  filter(date%in%unique(as.Date(InOutIt$date)))


TotalCases<-readRDS("DIVI_2022-10-27.rds")


Total_Cases<-TotalCases%>%
  mutate(date=as.Date(date), districtId=gemeindeschluessel)%>%
  filter(date%in%unique(train_data_l$date))%>%
  full_join(InOut_est)

Munchen<-ggplot(Total_Cases[Total_Cases$gemeindeschluessel==9162,])+
  geom_ribbon(aes(x=date, ymax=cumsum(Incoming_est), fill="cum. Incoming"), ymin=0)+
  geom_ribbon(aes(x=date, ymax=cumsum(Incoming_est)-cumsum(Outgoing_est), fill="cum. Incoming- cum. Outgoing"), ymin=0)+
  geom_line(aes(x=date, y=faelle_covid_aktuell, col="Total Covid ICU"), linewidth=0.8)+
  theme(legend.position="none")+
  theme_pubr()+
  scale_colour_manual(values = brewer.pal(n = 8, name = "Blues")[8])+
  scale_fill_manual(values = brewer.pal(n = 8, name = "Blues")[4:7])+
  ggtitle("District SK München")+
  theme(legend.title = element_blank())+
  xlab("Date")+ylab("Estimated cases")


Flensburg<-ggplot(Total_Cases[Total_Cases$gemeindeschluessel==1001,])+
  geom_ribbon(aes(x=date, ymax=cumsum(Incoming_est), fill="cum. Incoming"), ymin=0)+
  geom_ribbon(aes(x=date, ymax=cumsum(Incoming_est)-cumsum(Outgoing_est), fill="cum. Incoming- cum. Outgoing"), ymin=0)+
  geom_line(aes(x=date, y=faelle_covid_aktuell, col="Total Covid ICU"), linewidth=0.8)+
  theme(legend.position="none")+
  theme_pubr()+
  scale_colour_manual(values = brewer.pal(n = 8, name = "Blues")[8])+
  scale_fill_manual(values = brewer.pal(n = 8, name = "Blues")[4:7])+
  ggtitle("District SK Flensburg")+
  theme(legend.title = element_blank())+
  xlab("Date")+ylab("Estimated cases")

Hamburg<-ggplot(Total_Cases[Total_Cases$gemeindeschluessel==2000,])+
  geom_ribbon(aes(x=date, ymax=cumsum(Incoming_est), fill="cum. Incoming"), ymin=0)+
  geom_ribbon(aes(x=date, ymax=cumsum(Incoming_est)-cumsum(Outgoing_est), fill="cum. Incoming- cum. Outgoing"), ymin=0)+
  geom_line(aes(x=date, y=faelle_covid_aktuell, col="Total Covid ICU"), linewidth=0.8)+
  theme(legend.position="none")+
  theme_pubr()+
  scale_colour_manual(values = brewer.pal(n = 8, name = "Blues")[8])+
  scale_fill_manual(values = brewer.pal(n = 8, name = "Blues")[4:7])+
  ggtitle("District SK Hamburg")+
  theme(legend.title = element_blank())+
  xlab("Date")+ylab("Estimated cases")


Berlin<-ggplot(Total_Cases[Total_Cases$gemeindeschluessel==11000,])+
  geom_ribbon(aes(x=date, ymax=cumsum(Incoming_est), fill="cum. Incoming"), ymin=0)+
  geom_ribbon(aes(x=date, ymax=cumsum(Incoming_est)-cumsum(Outgoing_est), fill="cum. Incoming- cum. Outgoing"), ymin=0)+
  geom_line(aes(x=date, y=faelle_covid_aktuell, col="Total Covid ICU"), linewidth=0.8)+
  theme(legend.position="none")+
  theme_pubr()+
  scale_colour_manual(values = brewer.pal(n = 8, name = "Blues")[8])+
  scale_fill_manual(values = brewer.pal(n = 8, name = "Blues")[4:7])+
  ggtitle("District SK Berlin")+
  theme(legend.title = element_blank())+
  xlab("Date")+ylab("Estimated cases")


Munster<-ggplot(Total_Cases[Total_Cases$gemeindeschluessel==5515,])+
  geom_ribbon(aes(x=date, ymax=cumsum(Incoming_est), fill="cum. Incoming"), ymin=0)+
  geom_ribbon(aes(x=date, ymax=cumsum(Incoming_est)-cumsum(Outgoing_est), fill="cum. Incoming- cum. Outgoing"), ymin=0)+
  geom_line(aes(x=date, y=faelle_covid_aktuell, col="Total Covid ICU"), linewidth=0.8)+
  theme(legend.position="none")+
  theme_pubr()+
  scale_colour_manual(values = brewer.pal(n = 8, name = "Blues")[8])+
  scale_fill_manual(values = brewer.pal(n = 8, name = "Blues")[4:7])+
  ggtitle("District SK Münster")+
  theme(legend.title = element_blank())+
  xlab("Date")+ylab("Estimated cases")



plot_grid(Hamburg, Munster, Berlin, Munchen, ncol=2)


############# Spatial Effects  #################################################


smooth_names<-colnames(predict(modellist[[1]], type="lpmatrix"))


coeffcients_smooth<-apply(matrix(unlist(lapply(modellist[200:400], coef)), 
                                 byrow=TRUE, ncol=length(coef(modellist[[200]]))), 
                  2,median)[grepl("long", smooth_names)]
BasisFn_smooth<-predict(modellist[[1]], newdata=train_data_l, type="lpmatrix")[,grepl("long", smooth_names)]
SmoothEst<-unique(cbind(train_data_l[, c("date", "districtId")], 
                 "SmoothEff"=t(coeffcients_smooth%*%t(BasisFn_smooth)))%>%
  full_join(train_data_l)%>%
  dplyr::select(districtId, SmoothEff, geometry))%>%
  drop_na()%>%
  st_as_sf(crs = 4326)

library(ggplot2)
library(sf)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)

# Sicherstellen, dass `shpdatei` das richtige CRS hat
shpdatei <- st_transform(SmoothEst, crs = 4326)
# cities
cities <- data.frame(
  city = c("Hamburg", "Berlin", "Dortmund", "Munich", "Dresden", "Stuttgart"),
  lon = c(10.0, 13.4, 7.47, 11.58, 13.73, 9.18),
  lat = c(53.55, 52.52, 51.51, 48.14, 51.05, 48.78)
)
# Convert to an sf object
cities_sf <- st_as_sf(cities, coords = c("lon", "lat"), crs = 4326)
# Hintergrundkarte mit Natural Earth-Daten (Alternative zu OSM)
germany <- ne_countries(scale = "medium", country="Germany", returnclass = "sf")

mapsmooth<-ggplot() +
  geom_sf(data = germany, fill = "gray90", color = "white") +  # Germany background
  geom_sf(data = shpdatei, aes(fill = SmoothEff), color = "white", size = 0.2) +  # Spatial data with fill
  theme_pubr()+
  scale_fill_distiller(name=expression({{hat(f^1)}}(long, lat))) +  # Handle NA values in color scale
  annotation_scale(location = "bl", width_hint = 0.3) +  # Scale bar
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering()) +
  coord_sf(xlim = st_bbox(germany)[c("xmin", "xmax")], 
           ylim = st_bbox(germany)[c("ymin", "ymax")], 
           expand = FALSE) +  # Crop to Germany
  theme_minimal()+
  xlab("Latitude")+ylab("Longitude")+
  ggtitle("a) Estimated smooth spatial effect")+
  geom_sf(data = cities_sf, color = "black", size = 2) +  # Add points for cities
  geom_text(data = cities_sf, aes(label = city, geometry = geometry), 
            stat = "sf_coordinates", nudge_y = 0.4, size = 4, color = "black")#+
 # scale_fill_gradient2(low="#BFDDF6", mid = "lightgrey", high = "#4A80C0", name = expression({{hat(f^1)}}(long, lat))) 
  
mapsmooth



############# Temporal Effects  #################################################


coeffcients_smooth_t<-apply(matrix(unlist(lapply(modellist[200:400], coef)), 
                                 byrow=TRUE, ncol=length(coef(modellist[[200]]))), 
                          2,median)[grepl("date", smooth_names)]

coeffcients_smooth_t_min<-apply(matrix(unlist(lapply(modellist[200:400], coef)), 
                                   byrow=TRUE, ncol=length(coef(modellist[[200]]))), 
                            2,quantile, 0.025)[grepl("date", smooth_names)]
coeffcients_smooth_t_max<-apply(matrix(unlist(lapply(modellist[200:400], coef)), 
                                       byrow=TRUE, ncol=length(coef(modellist[[200]]))), 
                                2,quantile, 0.975)[grepl("date", smooth_names)]

BasisFn_smooth_t<-predict(modellist[[1]], newdata=train_data_l, type="lpmatrix")[,grepl("date", smooth_names)]
SmoothEst_t<-unique(cbind(train_data_l[, c("date", "districtId")], 
                        "SmoothEff"=t(coeffcients_smooth_t%*%t(BasisFn_smooth_t)),
                        "SmoothEff_2.5"=t(coeffcients_smooth_t_min%*%t(BasisFn_smooth_t)),
                        "SmoothEff_97.5"=t(coeffcients_smooth_t_max%*%t(BasisFn_smooth_t)))%>%
                    full_join(train_data_l)%>%
                    dplyr::select(date, SmoothEff, SmoothEff_2.5, SmoothEff_97.5))%>%
  unique()


Temp_smooth<-ggplot(SmoothEst_t)+
  geom_ribbon(aes(x=as.Date(date), ymin=SmoothEff_2.5, ymax=SmoothEff_97.5, 
                  fill="95% Interquartile range"))+
  geom_line(aes(x=as.Date(date), y=SmoothEff,
                col="Estimated Temporal Effect"))+
  theme_pubr()+
  labs(color=" ")+
  ggtitle("b) Smoothed estimated temporal effect")+
  # paste0("Median adjustment scale ", round(median(C_hat[360:400]), 
  #                                          2)))+
  ylab(expression({{hat(f^2)}}(t)))+ 
  scale_colour_manual(values = brewer.pal(n = 5, name = "Blues")[4:5])+
  scale_fill_manual(values = brewer.pal(n = 5, name = "Blues")[2:5])+
  xlab("Date (2021)")+
  #scale_x_continuous(breaks = seq(min(SmoothEst_t$date), max(SmoothEst_t$date), 30))+
  theme(legend.title = element_blank())





ggsave(plot_grid(mapsmooth, Temp_smooth, rel_heights = c(0.9, 0.1)), file="SmoothEstimate.pdf", width=14, height=5)



############# Coefficients Median  #############################################

VCOV_coef<-lapply(modellist[200:400], vcov)
VCOV_coef_arr<-array(0 , dim=c(length(1:400), nrow(VCOV_coef[[1]]), ncol(VCOV_coef[[1]])))

for(i in 1:200){
  VCOV_coef_arr[i,,]<-as.matrix(VCOV_coef[[i]])
}

VCOV_coef_med<-apply(VCOV_coef_arr, c(2,3), median)


coef_fix<-cbind(names(coef(modellist[[200]])[1:10]),
      round(apply(matrix(unlist(lapply(modellist[360:400], coef)), 
                   byrow=TRUE, ncol=length(coef(modellist[[200]]))), 
            2,median)[1:10], 3),
      round(sqrt(diag(VCOV_coef_med))[1:10], 3))

colnames(coef_fix)<-c("Covariates", "Coefficients", "Std Dev")

coef_fix

################################################################################
############# Variance adjustment ##############################################
################################################################################


############# Fixed Coefficients ###############################################

VCOV_coef<-lapply(modellist[200:400], vcov)
VCOV_coef_arr<-array(0 , dim=c(length(1:400), 
                               nrow(VCOV_coef[[1]]), 
                               ncol(VCOV_coef[[1]])))

for(i in 1:200){
  VCOV_coef_arr[i,,]<-as.matrix(VCOV_coef[[i]])
}

coeffimed<-apply(matrix(unlist(lapply(modellist[360:400], coef)), 
                        byrow=TRUE, ncol=length(coef(modellist[[200]]))), 
                 2,median)

Matrix_Coef<-matrix(unlist(lapply(modellist[360:400], coef)), 
       byrow=TRUE, ncol=length(coef(modellist[[200]])))
Matrix_Coef_med_rep<-matrix(rep(coeffimed, nrow(Matrix_Coef)), 
                            ncol=length(coef(modellist[[200]])),
            byrow=TRUE)

varianceRuns<-360:400







VCOV_sum<-array(0 , dim=c(length(varianceRuns), 
                          nrow(VCOV_coef[[1]]), 
                          ncol(VCOV_coef[[1]])))

for(i in 1:length(varianceRuns)){
  VCOV_sum[i,,]<-(Matrix_Coef[i,]-coeffimed)%*%t((Matrix_Coef[i,]-coeffimed))+
    vcov(modellist[[varianceRuns[i]]])
}

coef_fix_v<-cbind(names(coef(modellist[[200]])[1:10]),
                round(apply(matrix(unlist(lapply(modellist[360:400], coef)), 
                                   byrow=TRUE, ncol=length(coef(modellist[[200]]))), 
                            2,median)[1:10], 3),
                round(sqrt(diag(apply(VCOV_sum, c(2,3), mean)))[1:10], 3))

colnames(coef_fix_v)<-c("Covariates", "Coefficients", "Std Dev")

coef_fix_v

############# Fixed Coefficients ###############################################

varianceRuns<-360:400

coef_med<-apply(matrix(unlist(lapply(modellist[varianceRuns], coef)),
  ncol=length(lapply(modellist[varianceRuns], coef)[[1]]),
  byrow=TRUE), 2, median)


VCOV_Coef_Runs_s<-array(0, dim=c(length(varianceRuns), 
                              ncol=ncol(vcov(modellist[[i]])), 
                              nrow=nrow(vcov(modellist[[i]]))))
VCOV_Coef_Runs_b<-array(0, dim=c(length(varianceRuns), 
                                 ncol=ncol(vcov(modellist[[i]])), 
                                 nrow=nrow(vcov(modellist[[i]]))))

VCOV_Coef_Runs_w_s<-array(0, dim=c(length(varianceRuns), 
                                 ncol=max_lag-1, 
                                 nrow=max_lag-1))
VCOV_Coef_Runs_w_b<-array(0, dim=c(length(varianceRuns), 
                                 ncol=max_lag-1, 
                                 nrow=max_lag-1))
estimated_adj_w<-apply(omegas_adj[varianceRuns,], 2, median)/sum(apply(omegas_adj[varianceRuns,], 2, median))
for(i in varianceRuns){
  VCOV_Coef_Runs_b[(i)-min(varianceRuns)+1,,]<-
    ((unlist(lapply(modellist[i], coef))-coef_med)%*%t(
      unlist(lapply(modellist[i], coef))-coef_med))
  VCOV_Coef_Runs_s[(i)-min(varianceRuns)+1,,]<-
    vcov(modellist[[i]])
  VCOV_Coef_Runs_w_s[(i)-min(varianceRuns)+1,,]<-solve(VCOV[i,,])
  VCOV_Coef_Runs_w_b[(i)-min(varianceRuns)+1,,]<-((omegas_adj[i,]-estimated_adj_w)%*%t((omegas_adj[i,]-estimated_adj_w)))[1:(max_lag-1), 1:(max_lag-1)]
}



stddev_coef<-sqrt(diag(apply(VCOV_Coef_Runs_s, c(2,3), median)+(1+(max(varianceRuns)-min(varianceRuns))^(-1))*apply(VCOV_Coef_Runs_b, c(2,3), sum)/(max(varianceRuns)-min(varianceRuns)-1)))
stddev_w<-sqrt(diag(apply(VCOV_Coef_Runs_w_s, c(2,3), median)+(1+(max(varianceRuns)-min(varianceRuns))^(-1))*apply(VCOV_Coef_Runs_w_b, c(2,3), sum)/(max(varianceRuns)-min(varianceRuns)-1)))
cbind(names(coef(modellist[[1]])), round(stddev_coef, 3))

Exit_rates<-ggplot()+
  geom_hline(aes(yintercept=1/(max_lag), col="1/30"), linetype=2)+
  geom_ribbon(aes(x=1:(max_lag-1), ymin=estimated_adj_w[-max_lag]-1.96*stddev_w, 
                  ymax=estimated_adj_w[-max_lag]+1.96*stddev_w, 
                  fill="95% Confidence interval"))+
  geom_line(aes(x=1:(max_lag-1), y=estimated_adj_w[-max_lag],
                col="Estimated exit rate"))+
  theme_pubr()+
  labs(color=" ")+
  ggtitle("Estimated exit rate", expression(paste(sqrt(hat(c)), "=1.18", sep="")))+
  # paste0("Median adjustment scale ", round(median(C_hat[360:400]), 
  #                                          2)))+
  ylab("Value of exit rate")+ 
  scale_colour_manual(values = c("darkgrey", brewer.pal(n = 5, name = "Blues")[4:5]))+
  scale_fill_manual(values = c(brewer.pal(n = 5, name = "Blues")[2:5]))+
  xlab("Length of stay")+
  scale_x_continuous(breaks = seq(1, max_lag-1, 1))+
  theme(legend.title = element_blank())


Exit_rates

ggsave(Exit_rates, file="ExitRate.pdf", width=7, height=4)



###############################################################################

# Create a data frame with city coordinates

cities <- data.frame(
  city = c("Hamburg", "Berlin", "Dortmund", "Munich", "Dresden", "Stuttgart"),
  lon = c(10.0, 13.4, 7.47, 11.58, 13.73, 9.18),
  lat = c(53.55, 52.52, 51.51, 48.14, 51.05, 48.78)
)


# Convert to an sf object
cities_sf <- st_as_sf(cities, coords = c("lon", "lat"), crs = 4326)
# Hintergrundkarte mit Natural Earth-Daten (Alternative zu OSM)
germany <- ne_countries(scale = "medium", country="Germany", returnclass = "sf")

mapsmooth<-ggplot() +
  geom_sf(data = germany, fill = "gray90", color = "white") +  # Germany background
  geom_sf(data = shpdatei, aes(fill = SmoothEff), color = "white", size = 0.2) +  # Spatial data with fill
  theme_pubr()+
  scale_fill_distiller(name="Spatial Effects") +  # Handle NA values in color scale
  annotation_scale(location = "bl", width_hint = 0.3) +  # Scale bar
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering()) +
  coord_sf(xlim = st_bbox(germany)[c("xmin", "xmax")], 
           ylim = st_bbox(germany)[c("ymin", "ymax")], 
           expand = FALSE) +  # Crop to Germany
  theme_minimal()+
  ggtitle("Estimated Smooth Spatial Effects")+
  geom_sf(data = cities_sf, color = "black", size = 2) +  # Add points for cities
  geom_text(data = cities_sf, aes(label = city, geometry = geometry), 
            stat = "sf_coordinates", nudge_y = 0.4, size = 4, color = "black")


ggsave(mapsmooth, file="SmoothEff.pdf",
       height = 6, width = 5)



TotalCases_min<-as.data.frame(readRDS("DIVI_2022-10-27.rds")%>%
                            filter(as.Date(date)%in%unique(train_data_l$date))%>%
                            na.omit()%>%
                            dplyr::group_by(gemeindeschluessel)%>%
                            dplyr::summarize(min_cap=max((betten_frei+betten_belegt), na.rm = TRUE))%>%
                            mutate(districtId=as.character(gemeindeschluessel))%>%
                            full_join(SmoothEst))%>%
  na.replace(0)%>%
  st_as_sf(crs = 4326)%>%
  arrange(min_cap)





TotalCases<-as.data.frame(readRDS("DIVI_2022-10-27.rds")%>%
  filter(as.Date(date)%in%unique(train_data_l$date))%>%
    na.omit()%>%
  dplyr::group_by(gemeindeschluessel)%>%
  dplyr::summarize(max_cov=max(faelle_covid_aktuell/(betten_frei+betten_belegt+0.00001)*100, na.rm = TRUE))%>%
  mutate(districtId=as.character(gemeindeschluessel))%>%
  full_join(SmoothEst))%>%
  na.replace(0)%>%
  st_as_sf(crs = 4326)

names(train_data_l)
c(TotalCases[TotalCases$gemeindeschluessel=="2000","max_cov"])$max_cov # "SK Hamburg"
# 12.60794
c(TotalCases[TotalCases$gemeindeschluessel=="5913","max_cov"])$max_cov # "Berlin"
# 8.243727
c(TotalCases[TotalCases$gemeindeschluessel=="8111","max_cov"])$max_cov # "SK Dresden"
# 26
c(TotalCases[TotalCases$gemeindeschluessel=="9162","max_cov"])$max_cov # "SK Stuttgart"
# 26.07914
c(TotalCases[TotalCases$gemeindeschluessel=="11000","max_cov"])$max_cov # "SK München"
# 22.12544
c(TotalCases[TotalCases$gemeindeschluessel=="14612","max_cov"])$max_cov # "SK Dortmund"
# 37.93103

sort(unique(train_data_l$district))
unique(train_data_l$districtId[train_data_l$district=="SK Hamburg"|
                          train_data_l$district=="Berlin"|
                          train_data_l$district=="SK Dresden"|
                          train_data_l$district=="SK Stuttgart"|
                          train_data_l$district=="SK München"|
                          train_data_l$district=="SK Dortmund"])


 


og_Dist<-readRDS("DIVI_2022-10-27.rds")

og_Dist%>%
  filter(as.Date(date)%in%unique(train_data_l$date))%>%
  mutate(per_occ=faelle_covid_aktuell/(betten_frei+betten_belegt+0.00000000001)*100)%>%
  arrange(per_occ, decreasing=TRUE)




shpdatei <- st_transform(TotalCases, crs = 4326)
# Compute centroids of each district
district_centroids <- st_centroid(shpdatei)


MapIntro <- ggplot() +
  geom_sf(data = germany, fill = "gray90", color = "white") +  # Germany background
  geom_sf(data = TotalCases, aes(fill = max_cov), color = "white", size = 0.2) +  # Uniform blue fill
  geom_sf(data = district_centroids, color = "white", size = 0.6) +  # District centroids
  theme_pubr() +
  annotation_scale(location = "bl", width_hint = 0.3) +  # Scale bar
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering()) +
  coord_sf(xlim = st_bbox(germany)[c("xmin", "xmax")], 
           ylim = st_bbox(germany)[c("ymin", "ymax")], 
           expand = FALSE) +  # Crop to Germany
  theme_minimal() +
  xlab("")+ylab("")+
  ggtitle("a) Maximum ICU occupancy COVID-19 patients", "As percentage of total capacity, over time period") +
  geom_sf(data = cities_sf, color = "black", size = 2) +  # Add points for cities
  # Text outline (black)
  geom_text(data = cities_sf, aes(label = city, geometry = geometry), 
            stat = "sf_coordinates", nudge_y = 0.4, size = 4, color = "black")+
  scale_fill_gradient(low = "#BFDDF6", high = "#4A80C0", name = "Max Occupied \n (%)")  # Add city labels

MapIntro



###############################################################################



InfectRate<-train_data_l%>%
  group_by(date)%>%
  dplyr::summarize(Infec_3559=mean(G_4_7), Infec_6079=mean(G_5_7), 
            Infec_80=mean(G_6_7))
  

InfectRatePlot<-ggplot(InfectRate)+
  geom_line(aes(x=date, y=Infec_3559, col="Age group 35-59"), size=1.2)+
  geom_line(aes(x=date, y=Infec_6079, col="Age group 60-79"), size=1.2)+
  geom_line(aes(x=date, y=Infec_80, col="Age group 80+"), size=1.2)+
  theme_pubr() +
  scale_colour_manual(values = brewer.pal(n = 5, name = "Blues")[3:5])+
  xlab("Date (2021)")+ ylab("Average Logged Infection Rate")+
  theme(legend.title = element_blank())+
  ggtitle("b) Log 7-day-average infection rate", "Averaged by age group for all districts over time")


ggsave(plot_grid(MapIntro, InfectRatePlot), file="IntroData.pdf",
       height = 5, width = 10)

################################################################################

DIVI_cap<-read.csv("Intensivregister_Bundeslaender_Kapazitaeten.csv")

IncomingGermany_true<-DIVI_cap[,c("datum", "bundesland_id", "bundesland_name", "faelle_covid_erstaufnahmen")]%>%
  na.omit()%>%mutate(date=datum)


IncomingGermany_est<-merge(cbind(InOutIt[, c("date", "districtId")], "EstIn"=apply(InOutIt[, grepl(c("Incoming"), names(InOutIt))][, 356:396], 1, median)),
  cbind(InOutIt[, c("date", "districtId")], "EstOut"=apply(InOutIt[, grepl(c("Outgoing"), names(InOutIt))][, 356:396], 1, median)))%>%
  mutate(bundesland_id=floor(as.numeric(districtId)/1000))%>%
  dplyr::group_by(date, bundesland_id)%>%
  dplyr::summarize(IncomEst=sum(EstIn))

est_vs_true<-merge(IncomingGermany_est, IncomingGermany_true)


est_vs_true_plot<-ggplot(est_vs_true)+
  geom_line(aes(x=date, y=faelle_covid_erstaufnahmen, col="RKI Reported"), size=0.8)+
  geom_line(aes(x=date, y=IncomEst, col="Estimate"), size=0.8)+
  facet_wrap(factor(est_vs_true$bundesland_name))+
  theme_pubr()+
  ylab("Total admitted ICU patients with COVID-19 per county")+
  xlab("Date")+
  scale_colour_manual(values = brewer.pal(n = 5, name = "Blues")[c(3,6)])+
  theme(legend.title = element_blank())+
  ggtitle("RKI reported vs Estimated")
  
ggsave(est_vs_true_plot, file="RKIvsEst.pdf", height=10, width=10)



################################################################################

LogLike<-ggplot()+
  geom_line(aes(x=1:length(logLike), y=logLike[1:length(logLike)], col="Log-Likelihood"))+
  theme_pubr()+
  labs(color=" ")+
  ggtitle("Log-Likelihood over iterations", "ICU Covid-19 Data")+
  ylab("Log-Likelihood")+ 
  scale_colour_manual(values = c(brewer.pal(n = 5, name = "Blues")[4:5]))+
  scale_fill_manual(values = c(brewer.pal(n = 5, name = "Blues")[2:5]))+
  xlab("sEM Iteration")+
  theme(legend.title = element_blank())

ggsave(LogLike, file="C:/Users/ra98jiq/Downloads/LogLikeCov.pdf", height=3, width=5)



