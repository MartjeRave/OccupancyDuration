
######################################################
##### Offset calculation through lagged incoming #####
######################################################
# This function is used to calculate the number of outgoing
# patients based on the estimated exit rate and the simulated 
# incoming patients previously coming in.

Off_cal<-function(data_off=datalist[[j]], w_hat=w_true){
  
  data_tot<-data_off[, !grepl("Incoming_", names(data_off))] #Takes all previously included lags out of the data set
  data_temp<-data_off[, c("date", "districtId", "Incoming")] #Takes only the corrently simulated incoming units
  for(i in 1:length(w_hat)){ # creates a lag, merges it to slimmed data and sums over to a total incoming patients (weighted by the exit rate)
    data_temp<-data_off[, c("date", "districtId", "Incoming")]
    data_temp$date<-data_off$date+i
    names(data_temp)<-c("date", "districtId", paste0("Incoming_", i, "_"))
    data_tot<-data_tot%>%
      dplyr::full_join(data_temp, by=c("date", "districtId"))
  }
  
  Offset_data<-data_tot%>% # drops excess data (created though appending lags which reach beyond the max date)
    drop_na()
  return(Offset_data)
  
  
}

######################################################
############### Simulating data  #####################
######################################################
# This function is used to make simulated data with which 
# we create the simulated data.

# Here we take the true exit rate (assumed or randomly chosen),
# 200 districts, and 200 time points (after substracting the initial x time 
# points with which you calculate the number of outgoing units so we add 
# the length of exit rate), and essentially we set the coefficients in the 
# function, but could also be input factors if so chosen.
# Since we are looking to emulate somewhat resonable data, we just stuck with 
# these coefficients: In= 0.5+ 1* N + 0.2* M. 

sim_data_offset<-function(w=w_true,
                          districtId_sim=1:200, 
                          date_sim=1:(200+ncol(w_true_Grid)),
                          ss=123, with_coef=TRUE){
  # set seed   
  set.seed(ss)
  
  # generate M and N (covariates) 
  # use these with coefficients (set to b0=0.5, b1=1 and b2=0.2) to get incoming units
  
  if(with_coef){
    N<-rgamma(length(districtId_sim)*
                length(date_sim),0.1,0.5)
    M<-rep(rgamma(length(districtId_sim),1,3), each=length(date_sim))
    Incoming<-rpois(length(M), lambda=exp(0.5+1*M+0.2*N))    
    
    data_sim<-as.data.frame(cbind("date"=rep(date_sim,  # emulate data structure of covid data
                                             times=length(districtId_sim)), # emulate data structure of covid data 
                                  "districtId"=rep(districtId_sim, 
                                                   each=length(date_sim)), 
                                  M, N, Incoming))
  }else{
    Incoming<-rpois(length(date_sim)*length(districtId_sim), lambda=10)
    data_sim<-as.data.frame(cbind("date"=rep(date_sim, 
                                             times=length(districtId_sim)), 
                                  "districtId"=rep(districtId_sim, 
                                                   each=length(date_sim)), 
                                  Incoming))
  }
  
 
    ### For each incoming unit we sample a length of stay 
    ## on a given date the number of units sampled to have exited that day 
    ## are aggregated to the total of outgoing units (we call this deterministic) 
  
    # sample length of stay
    all_length<-suppressMessages(as.data.frame(cbind("total_sample"=
                                                       sample(1:length(w), 
                                                              size=sum(data_sim$Incoming), 
                                                              prob=w, replace=TRUE),
                                                     "date"=rep(data_sim$date, 
                                                                times=data_sim$Incoming),
                                                     "districtId"=rep(data_sim$districtId, 
                                                                      times=data_sim$Incoming)))%>%
                                   dplyr::group_by(date, districtId, total_sample)) 
    # count number of units going out at a length of stay, for a given day and given district
    all_length<-sqldf("SELECT COUNT(*) AS count, date, districtId, total_sample FROM all_length GROUP BY total_sample, date, districtId")
    
    # add the number of outgoing units to the date on which they are actually outgoing (given length of stay)
    data_tot<-data_sim 
    for(i in 1:length(w)){
      lag<-all_length[all_length$total_sample==i,]%>%
        dplyr::mutate(date=date+i)%>%
        dplyr::select(date, count, districtId)
      names(lag)<-c("date", paste0("Incom_lag",i), "districtId")
      data_tot<-suppressMessages( data_tot%>%
                                     full_join(lag, by=c("date", "districtId"))%>%
                                     na.replace(replace=0))
    }
    data_tot<-data_tot%>%
      filter(date>=min(data_sim$date)+length(w)&
               date<=max(data_sim$date))%>%
      arrange(districtId, date) # restrict to date 
    data_tot$Outgoing<-rowSums(data_tot[,grepl("lag", names(data_tot))]) # sum all outgoing units
    data_tot$Diff<-data_tot$Incoming-data_tot$Outgoing # Diff function
  
  
  return(data_tot)
}

######################################################
###############  Estimating Exit Rate ################
######################################################
# This function is the set up for the estimation of the exit rate
# Since we are using solve.QP() we use the first derivative of the Poisson 
# model and the second derivative of the Poisson model.
# See appendix in Paper

######## First derivative ######
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



######## Second derivative ######
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

######## Maximizing the likelihood with respect to the exit rate ######
###### Off_cal, Score_func_w and Neg2_der_w are inputs there ##########
## See a more detailed description of the mathematical bg here: https://journals.sagepub.com/doi/full/10.1177/1471082X241235024
optimize_Like<-function(w_opt_prev=w_try4_new, 
                        data_like_opp=datalist[[l]]){
  
  Off_data_temp<-Off_cal(w_hat=w_opt_prev, 
                         data_off=data_like_opp[, 
                                                !grepl("Incom_", names(data_like_opp))])
  
  score_w<-Score_func_w(w_score_fun=w_opt_prev,
                        data_score_fun=Off_data_temp)
  
  FisherInf<-Neg2_der_w(w_derivative=w_opt_prev,
                        data_derivative=Off_data_temp)
  
  Dmat<-FisherInf
  d<-(score_w+ FisherInf%*%w_opt_prev[-length(w_opt_prev)])
  b<-c(rep(0, length(w_opt_prev)-1), rep(-1, length(w_opt_prev)-1))
  
  
  Amat<-cbind(diag(1, length(w_opt_prev)-1), diag(-1, length(w_opt_prev)-1))
  
  Opt<-solve.QP(Dmat=Dmat,
                dvec=d,
                Amat=Amat,
                bvec=b)
  
  return(list("w_sol"=Opt$solution, "Infor"=Dmat, "w_sol_uncon"=Opt$unconstrained.solution))
}


######## Iteratively  optimizing the exit rate until the estimation ###### 
######## has reached convergence #########################################


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
###############  Simulate number of units (in and out) ###
######################################################



rdiffpois_array_paper =function(data_rdiff=preptrain_try, # data
                                IMSTEP=I_Mstep, # incoming gam
                                w_sim=w_start, # exit rate
                                max.k=300){ # probability calculations cut at maximum 300 incoming units   
  
  ### set seed
  set.seed(123) 
  
  ### set up for calculated lag for the first x dates of observation as the lag is not actually observed before the first day (so we impute with first value)
  Lagged_Data<-rep(data_rdiff[data_rdiff$date==min(data_rdiff$date),], (length(w_sim)+1))%>%
    mutate(date=date-rep(1:(length(w_sim)+1), each=length(unique(data_rdiff$districtId))))%>%
    rbind(data_rdiff)%>%
    arrange(date, districtId)
  
  ### the intensity parameter of the incoming units
  Lagged_Data$lambda_in<-predict.gam(IMSTEP, newdata=Lagged_Data, type="response")
  
  #### 2D array for faster calculations date (rows), districts (colums), lambda_in (values)  
  lambda_in_array<-(Lagged_Data%>%
                      dplyr::select(date, districtId, lambda_in)%>%
                      mutate(lambda_in=ifelse(lambda_in==0, 
                                              0.000001, lambda_in))%>%
                      arrange(date, districtId)%>%
                      spread(key="date", value="lambda_in"))
  
  #### set up for simulation an incoming number of units 
  In_sim<-(Lagged_Data%>%
             dplyr::select(date, districtId)%>%
             mutate(in_sim=NA)%>%
             arrange(date, districtId)%>%
             spread(key="date", value="in_sim"))
  
  ### simulate incoming units based on the intensity parameter (for outgoing intesnity parameter)
  In_sim[, 2:(length(w_sim)+2)]<-apply(lambda_in_array[,2:(length(w_sim)+2)], 2, 
                                       rpois, n=nrow(In_sim))
  ### set up outgoing intensity parameter
  lambda_out_array<-In_sim[, -c(2:(length(w_sim)+2))]
  
  ### incoming simulated summed and weighted by exit rate
  lambda_out_array[,2]<-t(w_sim%*%t(In_sim[, (length(w_sim)+2):3]))
  
  #### we essentially make a pivot_wide() (in 3D) to make the calculations a little faster
  Diff_sim<-(data_rdiff%>%
               dplyr::select(date, districtId, Diff)%>%
               arrange(date, districtId)%>%
               spread(key="date", value="Diff"))
  #### set up for simulation an incoming number of units 
  Out_sim<-(data_rdiff%>%
              dplyr::select(date, districtId)%>%
              mutate(out_sim=NA)%>%
              arrange(date, districtId)%>%
              spread(key="date", value="out_sim"))
  datum<-unique(data_rdiff$date)[1]
  
  ### What we have until now: the 2D array for incoming and outgoing intensity
  ### The observed Difference in a 2D array.
  ### Now we calculate the probability of each possible number of incoming unit 
  ### 1:max.k, and outgoing unit 1:max.k-Diff. 
  ### multiply the probabilities, norm them and sample from each vector of probabilities
  ### the number of incoming units and calculate the number of outgoing units
  
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
    In_sim[, names(In_sim)==as.character(datum)]<-apply(Prob, 1, sample, x=0:max.k, size=1, replace=TRUE)
    
    Out_sim[, names(Out_sim)==as.character(datum)]<-In_sim[, names(In_sim)==as.character(datum)]-
      Diff_sim[, names(Diff_sim)==as.character(datum)]
    
  }
  
  ### reformat into the original data format
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




######################################################
############# Bias Adjustment  #######################
######################################################
### Simulate number of incoming units based on estimated coefficients

sim_data_jk<-function(w=median,
                      ss=123, 
                      deterministic=TRUE,
                      coef=Coef_in,
                      covariates=datalist_prep[,c("date", 
                                                  "districtId",
                                                  "M", "N")]){
  ### set seed
  set.seed(ss)
  ### set up for starting data (lag is not observed before the first date so we impute)
  covariates<-covariates%>%
    rbind(rep(covariates[covariates$date==min(covariates$date),], length(w))%>%
            mutate(date=date-rep(1:length(w), each=length(unique(covariates$districtId)))))%>%
    arrange(date, districtId)
  ### calculate incoming intensity parameter
  lambda.in<-exp(t(coef%*%t(cbind(1, covariates%>%dplyr::select(-date, -districtId)))))
  ### use it to simulate incoming patients
  Incoming<-rpois(length(lambda.in), lambda=lambda.in)
  ### append it to og data
  data_sim<-cbind(covariates, "Incoming"=Incoming)
  ### sample length of stay for each incoming unit
  all_length<-suppressMessages(as.data.frame(cbind("total_sample"=
                                                     sample(1:length(w), 
                                                            size=sum(data_sim$Incoming), 
                                                            prob=w, replace=TRUE),
                                                   "date"=rep(data_sim$date, 
                                                              times=data_sim$Incoming),
                                                   "districtId"=rep(data_sim$districtId, 
                                                                    times=data_sim$Incoming)))%>%
                                 dplyr::group_by(date, districtId, total_sample))
  ### sample length of stay for eachi ncoming unit
  all_length<-sqldf("SELECT COUNT(*) AS count, date, districtId, total_sample FROM all_length GROUP BY total_sample, date, districtId")
  ### append to og data (lagged to the sampled length of stay) and sum over sample number of outgoing units per date and district
  data_tot<-data_sim
  for(i in 1:length(w)){
    lag<-all_length[all_length$total_sample==i,]%>%
      dplyr::mutate(date=date+i)%>%
      dplyr::select(date, count, districtId)
    names(lag)<-c("date", paste0("Incom_lag",i), "districtId")
    data_tot<-suppressMessages( data_tot%>%
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



### pretty self explantory, though the function allows for bias adjustment with 
### several functions not just the squared adjustment 
### we correct for the pull towards the probability of 1/L

bias_adj<-function(w_og=omega_org, w_jk=omega_jk, max_l=max_lag, 
                   func=c("none", "sq", "lin", "exp", "abs")){
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
