#####################
#### Packages
#####################
# Datenbearbeitung
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
library(rgeos)
library(skellam)
library(maptools)
library(RColorBrewer)
library(gridExtra)
library(rgeos)
library(rgdal)
library(leaflet)
library(fitdistrplus)
library(lme4)
library(skellam)
library(mgcv)
library(splines)
library(rgdal)
library(sp)
library(dplyr)
library(rgdal)
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

#########FUNCTIONS

smooth_terms<-function(object){
  smoothed_terms<-data.frame(matrix(1:length(names(attr(object@misc$Xvlm.aug,"spar.vlm"))[grepl("s\\(", names(attr(object@misc$Xvlm.aug,"spar.vlm")))]), nrow=1))
  names(smoothed_terms)<-names(attr(object@misc$Xvlm.aug,"spar.vlm"))[grepl("s\\(", names(attr(object@misc$Xvlm.aug,"spar.vlm")))]
  return(smoothed_terms)
}

week_numeric<-function(data_date){
  start_date<-c(4,5,6,0,1,2,3)[as.numeric(wday(data_date))]
  start_of_week<-(as.Date(data_date)-start_date)
  return(start_of_week)
}


######Colourblind friendly Palet
cbPalette <- c("#E69F00","#999999",  "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette



########## Variance estimation ###########################
RubinVar<-function(Varianceruns=Run1$Coeff_var_in[200:500,,], 
                   Coefficientruns=Coeff_in[200:500,-1], runs=500){
  Rub1<-apply(Varianceruns,c(2,3),mean)
  betaminusbetaline<-(Coefficientruns-
                        matrix(rep(apply(Coefficientruns, 2, mean),
                                   each=nrow(Varianceruns)),
                       nrow=nrow(Varianceruns), ncol=ncol(Coefficientruns)))  
  variance_step<-array(NA,dim=c(nrow(betaminusbetaline), ncol(Coefficientruns),
                                ncol(Coefficientruns)))
  for(i in 1:nrow(betaminusbetaline)){
    variance_step[i,,]<-t(as.matrix(betaminusbetaline[i,]))%*%
      as.matrix(betaminusbetaline[i,])
  }
  Rub2<-apply(variance_step,c(2,3),mean)
  Rubinvars<-Rub1+(1+1/nrow(Varianceruns))*(1/(nrow(Varianceruns)-1))*Rub2
  return(Rubinvars)
}


########## E step ###########################
rdiffpois = function( lambda.in=lambda0.in, lambda.out=lambda0.out , D=train$Diff,
                      max_In=2000){
  ## Funktion zur Simulation des E-step. Dabei wird 
  ## I ~ Poisson (lambda.in) und R ~ Poisson(lambda.out)
  ## simuliert, wobei gelten muss I - R = D
  
  k.in = matrix(rep(c(0: max_In), each=length(lambda0.in)), 
                nrow=length(lambda0.in))
  
  D=matrix(rep(D, max_In+1), nrow=length(lambda.in))
  k.out = k.in - D
  index<-(k.in>=0)&(k.out>=0)
  
  
  I<-c()
  R<-c()
  
  for(i in 1:length(lambda0.in)){
    k.in.i = k.in[i,index[i,]]
    k.out.i = k.out[i,index[i,]]
    prob = dpois(k.in.i, lambda= lambda.in[i]) * dpois( k.out.i, lambda = lambda.out[i])
    x=sum(prob)
    if(sum(prob)==0){
      x=1
    }
    prob = prob/x
    j = sample(length(prob), 1, prob=prob)
    
    I[i]<-k.in.i[j]
    R[i]<-k.out.i[j]
  }
  
  return(list( I = I, R = R))
}


########## M step ###########################
Mstep<-function(train1=train_train,
                formula_in1=formula_in, 
                formula_out1=formula_out, 
                data_tot1=data_tot){
  
  I_Mstep<-mgcv::gam(formula_in1, data=train1, family="poisson")
  
  pred_data<-data_tot1%>%
    dplyr::filter(data_tot1$date>=(min(train1$date)-(10))&
                    data_tot1$date<(max(train1$date)+(10)))
    
  sim_start<-cbind(pred_data[,c("districtId", "date")], 
                   "sim"=predict.gam(I_Mstep, newdata = pred_data, 
                                     type="response"))
  
  date_tot<-sim_start
  

  for(i in (1:56)){
    dat_int<-sim_start
    dat_int$date<-sim_start$date+i
    names(dat_int)<-c("districtId", "date", paste0("sim_", i))
    date_tot<-date_tot%>%
      full_join(dat_int, by=c("districtId", "date"))
    
    date_tot[is.na(date_tot[, paste0("sim_", i)]) , paste0("sim_", i)]<-
      date_tot[is.na(date_tot[, paste0("sim_", i)]) , (grep(paste0("sim_", i), 
                                                            names(date_tot))-1)]  
  }
  
  date_tot<-date_tot%>%
    filter(date<=max(train1$date)&
             date>=min(train1$date))
  suv_days_p<-t(suvdays)[2,]
  names(suv_days_p)<-t(suvdays)[1,]
  date_tot$sim_week<-as.matrix(date_tot[, grep("_", 
                               names(date_tot))])%*%as.matrix(1-suv_days_p)
  train_out<-train1%>%
    dplyr::full_join(date_tot[,c("districtId", "date", "sim_week")], 
                     by=c("districtId", "date"))
  
  R_Mstep<-mgcv::gam(formula_out1, offset = log(train_out$sim_week+0.1),
                       data=train1, family="poisson")
    
  
  return(list("I_Mstep"=I_Mstep, "R_Mstep"=R_Mstep))
  
}


########## EM algorithm ###########################
EMsimm<-function(train_train=train, 
                 runs1=500, 
                 data_tot= train_tot,
                 formula_in=as.formula(I_sim ~ weekday + 
                                         s(as.numeric(date), bs = "ps") + 
                                         G_4_7_lag_1 + 
                                         G_5_7_lag_1 + 
                                         G_6_7_lag_1 + 
                                         s(lat, long, bs = "tp")), 
                 formula_out=as.formula(R_sim ~ weekday + 
                                          s(as.numeric(date), bs = "ps") + 
                                          G_4_7_lag_1 + 
                                          G_5_7_lag_1 + 
                                          G_6_7_lag_1 + 
                                          s(lat, long, bs = "tp")),
                 nr_coefficients=48){
  
  lambda0.in<-rep(sample(1:10,1), nrow(train_train))
  lambda0.out<-rep(sample(11:20,1), nrow(train_train))
  
  start_sim<-rdiffpois(lambda0.in,lambda0.out, train_train$Diff)
  
  train_train$I_sim<-start_sim$I
  train_train$R_sim<-start_sim$R
  
  fitlist<-list()
  modellist<-list()
  Coeff_in<-NULL
  Coeff_out<-NULL
  
  nr_coef_in<-nr_coefficients
  nr_coef_out<-nr_coefficients
  
  VCOV_in<-matrix(0, nr_coef_in, nr_coef_in)
  VCOV_out<-matrix(0, nr_coef_out, nr_coef_out)
  
  Coeff_var_in<-array(NA, dim=c(runs1, nrow(VCOV_in),ncol(VCOV_in)))
  Coeff_var_out<-array(NA, dim=c(runs1, nrow(VCOV_out),ncol(VCOV_out)))
  loglike<-c()
  
  # put this before start of loop
  # total = length(1:runs)
  
  # Eine höhere Letalität und lange Beatmungsdauer
  # unterscheiden COVID-19 von schwer verlaufenden
  # Atemwegsinfektionen in Grippewellen
  # Autorinnen und Autoren
  # Kristin Tolksdorf | Dr. Silke Buda | Dr. Ekkehard Schuler |
  #   Prof. Lothar H. Wieler | Prof. Walter Haas
  # Korrespondenz: TolksdorfK@rki.de
  survival_days<-c(1:56)
  survival_prob<-c(0.20, 0.32, 0.4, 0.46, 0.51, 
                   0.55, 0.57, 0.60, 0.63, 0.69,
                   0.71, 0.73, 0.74, 0.77, 0.78,
                   0.79,seq(0.82, 1, by=(1-0.82)/(56-17)))
  
  suvdays<-cbind("days"=survival_days, "prob"=survival_prob)
  
  for(k in 1:runs1){
    train_train$I_sim<- start_sim$I
    train_train$R_sim<- start_sim$R
    
    Msteprun<-Mstep(train1=train_train, 
                    formula_in1=formula_in, 
                    formula_out1=formula_out, 
                    data_tot1=data_tot)
    
    Coeff_var_in[k,,]<-vcov(Msteprun$I_Mstep)
    Coeff_var_out[k,,]<-vcov(Msteprun$R_Mstep)
    lambda.in.new<-predict(Msteprun$I_Mstep, type="response")
    lambda.out.new<-predict(Msteprun$R_Mstep, type="response")
    
    loglike[k]<-sum(log(pskellam(train_train$Diff, lambda.in.new, 
                                 lambda.out.new)))
    
    print(paste("iteration", k, "is done"))
    fitlist[[k]]<-rdiffpois(lambda.in.new,lambda.out.new, train1$Diff)
    start_sim<-as.data.frame(cbind("I"=fitlist[[k]]$I, "R"=fitlist[[k]]$R))
    
    modellist[[k]]<-list("I_mod"=Msteprun$I_Mstep,"R_mod"= Msteprun$R_Mstep)
    if(k!=1){
      Coeff_in<-rbind(Coeff_in, as.data.frame(cbind(k, 
                                                    t(coef(Msteprun$I_Mstep)))))
      Coeff_out<-rbind(Coeff_out, as.data.frame(cbind(k, 
                                                      t(coef(Msteprun$R_Mstep)))))
      # suppressMessages(Coeff_out<-Coeff_out%>%
      #                    full_join(as.data.frame(t(coef(Msteprun$R_Mstep)))))
    }else{
      Coeff_in<-as.data.frame(cbind(k, t(coef(Msteprun$I_Mstep))))
      Coeff_out<-as.data.frame(cbind(k, t(coef(Msteprun$R_Mstep))))
      
    }
    
  }
  
  
  return(list("Coeff_in"=Coeff_in, "Coeff_out"=Coeff_out, 
              "Coeff_var_in"=Coeff_var_in, "Coeff_var_out"=Coeff_var_out,
              # "Varin"=Varin, "Varout"=Varout, 
              "starvaluein"=unique(lambda0.in),
              "starvalueout"=unique(lambda0.out), "Models"=modellist, 
              "Response"=fitlist#, "tav"=tav1, "tot_days"=tot_days2
  ))
}


formula_in<-as.formula(I_sim ~ weekday + 
                         s(as.numeric(date), bs = "ps") + 
                         G_4_7_lag_1 + 
                          G_5_7_lag_1 + 
                         G_6_7_lag_1 + 
                         s(lat, long, bs = "tp"))
formula_out<-as.formula(R_sim~ weekday+ 
                          s(as.numeric(date), bs="ps")+
                          G_4_7_lag_1+
                          G_5_7_lag_1+ 
                          G_6_7_lag_1+
                          s(lat, long, bs="tp"))



# ##########  OBSOLETE:::: EM algorithm Change of offset to try which one makes sensible results ###########################
# 
# Mstep_try<-function(train=train1, offs=FALSE, 
#                 formula_in_mstep=formula_in1, 
#                 formula_out_mstep=formula_out1, 
#                 tav=1, tot_days=30, data_tot1=data_tot){
#   
#   I_Mstep<-mgcv::gam(formula_in_mstep, data=train, family="poisson")
#   
#   pred_data<-data_tot1%>%
#     dplyr::filter(data_tot1$date>=(min(train$date)-(tot_days+7))&
#                     data_tot1$date<(max(train$date)+(tot_days+7)))
#   
#   Incoming_pred<-cbind(pred_data[,c("districtId", "date")], 
#                    "sim"=predict.gam(I_Mstep, newdata = pred_data, type="response"))
#   
#   date_tot<-Incoming_pred
#   
#   for(i in (tav:(tav+tot_days-1))){
#     dat_int<-Incoming_pred
#     dat_int$date<-Incoming_pred$date+i
#     names(dat_int)<-c("districtId", "date", paste0("sim_", i))
#     date_tot<-date_tot%>%
#       full_join(dat_int, by=c("districtId", "date"))
#     
#     date_tot[is.na(date_tot[, paste0("sim_", i)]) , paste0("sim_", i)]<-
#       date_tot[is.na(date_tot[, paste0("sim_", i)]) , (grep(paste0("sim_", i), names(date_tot))-1)]  
#   }
#   
#   date_tot<-date_tot%>%
#     filter(date<=max(train$date)&
#              date>=min(train$date))
#   
#   date_tot$sim_week<-rowSums(date_tot[, grep("_", names(date_tot))])
#   
#   train_out<-train%>%
#     dplyr::full_join(date_tot[,c("districtId", "date", "sim_week")], by=c("districtId", "date"))
#   
#   if(offs){
#     R_Mstep<-mgcv::gam(formula_out_mstep, offset = log(sim_week+0.1),
#                        data=train_out, family="poisson")
#     
#   }else{
#     R_Mstep<-mgcv::gam(formula_out_mstep, 
#                        data=train_out, family="poisson")
#     
#   }
#   
#   return(list("I_Mstep"=I_Mstep, "R_Mstep"=R_Mstep))
#   
# }
# 
# 
# 
# 
# 
# 
# RubinVar<-function(Varianceruns=Run1$Coeff_var_in[200:500,,], Coefficientruns=Coeff_in_[200:500,-1], runs=500){
#   Rub1<-apply(Varianceruns,c(2,3),median)
#   betaminusbetaline<-(Coefficientruns-
#                         matrix(rep(apply(Coefficientruns, 2, median),
#                                    each=nrow(Varianceruns)),
#                                nrow=nrow(Varianceruns), ncol=ncol(Coefficientruns)))  
#   variance_step<-array(NA,dim=c(nrow(betaminusbetaline), ncol(Coefficientruns),ncol(Coefficientruns)))
#   for(i in 1:nrow(betaminusbetaline)){
#     variance_step[i,,]<-t(as.matrix(betaminusbetaline[i,]))%*%
#       as.matrix(betaminusbetaline[i,])
#   }
#   Rub2<-apply(variance_step,c(2,3),median)
#   Rubinvars<-Rub1+(1+1/nrow(Varianceruns))*(1/(nrow(Varianceruns)-1))*Rub2
#   return(Rubinvars)
# }






########## Plotfunktion for output of models ###########################


plotsfunc<-function(Coeff_in_=Run1$Coeff_out[,-1], 
                    k=1, 
                    model_for_basis=Run1$Models[[1]]$R_mod,
                    training_data=train, 
                    label=c("Incoming", "Outgoing"),
                    wdepr=c("yes", "no"), Var_in=Variance_out){
  model_for_space=model_for_basis
  model_for_basis=model_for_basis[[1]]$R_mod
  
  Coeff_in_function=cbind("k"=k:nrow(Coeff_in_), Coeff_in_[k:nrow(Coeff_in_),])
  
  variance1<-cbind(as.data.frame(names(Coeff_in_)), "variance1"=diag(Var_in))
  variance1$conf95<-1.96*sqrt(variance1$variance1)
  medians_1<-apply(Coeff_in_[200:500,], 2, median)
  medians_pervar<-as.data.frame(cbind("medians"=medians_1, 
                                      "lci"=medians_1-variance1$conf95, 
                                      "uci"=medians_1+variance1$conf95))
  
  
  
  ################## Coefficients for linear effects #########################
  gInter<-ggplot(data=Coeff_in_function,aes(x=k, y=`(Intercept)`, col="Intercept"))+
    # geom_rect(aes(xmin=-Inf, xmax=Inf,
    #               ymin=medians_pervar$lci[row.names(medians_pervar)=="(Intercept)"],
    #               ymax=medians_pervar$uci[row.names(medians_pervar)=="(Intercept)"],
    #               fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
    #           col="white")+
    geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)=="(Intercept)"],
                   col="Intercept"))+
    geom_line()+
    #coord_cartesian(ylim=c(-3,-2))+
    #        geom_smooth(se = FALSE)+
    ggtitle(paste(label, "patients"), "Intercept")+theme_pubr()+
    labs(col=" ")+
    ylab("Intercept")+
    scale_color_manual(values=cbPalette)+
    xlab("Iteration")
  # 
  
  gWeek<-ggplot(data=Coeff_in_function)+
    #        geom_smooth(aes(x=k, y=weekdayMon, col="Mon_Effect"))+
    # geom_rect(aes(xmin=-Inf, xmax=Inf, 
    #               ymin=medians_pervar$lci[row.names(medians_pervar)=="weekdayMon"], 
    #               ymax=medians_pervar$uci[row.names(medians_pervar)=="weekdayMon"], 
    #               fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
    #           col="white")+
    # geom_rect(aes(xmin=-Inf, xmax=Inf, 
    #               ymin=medians_pervar$lci[row.names(medians_pervar)=="weekdayTue"], 
    #               ymax=medians_pervar$uci[row.names(medians_pervar)=="weekdayTue"], 
    #               fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
    #           col="white")+
  # #            geom_smooth(aes(x=k, y=weekdayTueday, col="Tueday_Effect"))+
  # geom_rect(aes(xmin=-Inf, xmax=Inf, 
  #               ymin=medians_pervar$lci[row.names(medians_pervar)=="weekdayWed"], 
  #               ymax=medians_pervar$uci[row.names(medians_pervar)=="weekdayWed"], 
  #               fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
  #           col="white")+
  # #            geom_smooth(aes(x=k, y=weekdayWed, col="Wed_Effect"))+
  # geom_rect(aes(xmin=-Inf, xmax=Inf, 
  #               ymin=medians_pervar$lci[row.names(medians_pervar)=="weekdayThu"], 
  #               ymax=medians_pervar$uci[row.names(medians_pervar)=="weekdayThu"], 
  #               fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
  #           col="white")+
  # #            geom_smooth(aes(x=k, y=weekdayThuday, col="Thuday_Effect"))+
  # geom_rect(aes(xmin=-Inf, xmax=Inf, 
  #               ymin=medians_pervar$lci[row.names(medians_pervar)=="weekdaySat"], 
  #               ymax=medians_pervar$uci[row.names(medians_pervar)=="weekdaySat"], 
  #               fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
  #           col="white")+
  # #        geom_smooth(aes(x=k, y=weekdaySat, col="Sat_Effect"))+
  # geom_rect(aes(xmin=-Inf, xmax=Inf, 
  #               ymin=medians_pervar$lci[row.names(medians_pervar)=="weekdaySun"], 
  #               ymax=medians_pervar$uci[row.names(medians_pervar)=="weekdaySun"], 
  #               fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
  #           col="white")+
  #          geom_smooth(aes(x=k, y=weekdaySunday, col="Sunday_Effect"), se = FALSE)+
  ggtitle(" ", "Weekday effect")+
    #coord_cartesian(ylim=c(-0.6,0.4))+
    geom_line(aes(x=k, y=weekdaySat, col="Saturday"))+
    geom_line(aes(x=k, y=weekdaySun, col="Sunday"))+
    geom_line(aes(x=k, y=weekdayMon, col="Monday"))+
    geom_line(aes(x=k, y=weekdayTue, col="Tuesday"))+
    geom_line(aes(x=k, y=weekdayWed, col="Wednesday"))+
    geom_line(aes(x=k, y=weekdayThu, col="Thursday"))+
    geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)=="weekdayMon"], 
                   col="Monday"), alpha=0.5)+
    geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)=="weekdayTue"], 
                   col="Tuesday"), alpha=0.5)+
    geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)=="weekdayWed"], 
                   col="Wednesday"), alpha=0.5)+
    geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)=="weekdayThu"], 
                   col="Thursday"), alpha=0.5)+
    geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)=="weekdaySat"], 
                   col="Saturday"), alpha=0.5)+
    geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)=="weekdaySun"], 
                   col="Sunday"), alpha=0.5)+
    ylab("Weekday effect")+theme_pubr()+
    scale_color_manual(values=cbPalette)+labs(col=" ")+
    xlab("Iteration")
  
  gInfect<-ggplot(data=Coeff_in_function)+
    # geom_rect(aes(xmin=-Inf, xmax=Inf, 
    #               ymin=medians_pervar$lci[row.names(medians_pervar)=="G_4_7_lag_1"], 
    #               ymax=medians_pervar$uci[row.names(medians_pervar)=="G_4_7_lag_1"], 
    #               fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
    #           col="white")+
    # geom_rect(aes(xmin=-Inf, xmax=Inf, 
    #               ymin=medians_pervar$lci[row.names(medians_pervar)=="G_5_7_lag_1"], 
    #               ymax=medians_pervar$uci[row.names(medians_pervar)=="G_5_7_lag_1"], 
    #               fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
    #           col="white")+
    # geom_rect(aes(xmin=-Inf, xmax=Inf, 
  #               ymin=medians_pervar$lci[row.names(medians_pervar)=="G_6_7_lag_1"], 
  #               ymax=medians_pervar$uci[row.names(medians_pervar)=="G_6_7_lag_1"], 
  #               fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
  #           col="white")+
  geom_line(aes(x=k, y=G_4_7_lag_1, col="35-59 yo"))+
    #geom_smooth(aes(x=k, y=G_4_7_lag_1, col="35-59 yo"), se = FALSE)+
    geom_line(aes(x=k, y=G_5_7_lag_1, col="60-79 yo"))+
    # geom_smooth(aes(x=k, y=G_5_7_lag_1, col="60-79 yo"), se = FALSE)+
    geom_line(aes(x=k, y=G_6_7_lag_1, col="80+ yo"))+
    #  geom_smooth(aes(x=k, y=G_6_7_lag_1, col="80+ yo"), se = FALSE)+
    geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)==
                                                       "G_4_7_lag_1"], 
                   col="35-59 yo"))+
    geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)==
                                                       "G_5_7_lag_1"], 
                   col="60-79 yo"))+
    geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)==
                                                       "G_6_7_lag_1"], 
                   col="80+ yo"))+labs(col=" ")+
    
    # coord_cartesian(ylim=c(0,0.6))+
    ggtitle(" ", "Infection effect by age group")+
    ylab("Lagged infection")+theme_pubr()+
    scale_color_manual(values=cbPalette)+
    xlab("Iteration")
  
  if(wdepr=="yes"){
    gdepri<-ggplot(data=Coeff_in_function)+
      geom_rect(aes(xmin=-Inf, xmax=Inf, 
                    ymin=medians_pervar$lci[row.names(medians_pervar)=="GIMD10"], 
                    ymax=medians_pervar$uci[row.names(medians_pervar)=="GIMD10"], 
                    fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
                col="white")+
      geom_line(aes(x=k, y=GIMD10, col="Deprivation Index"))+
      geom_smooth(aes(x=k, y=GIMD10, col="Deprivation Index"), se = FALSE)+
      ggtitle("Estimated effect of deprivation Index",
              paste(label, "patients"))+
      geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)=="GIMD10"], 
                     col="Deprivation Index"))+
      #coord_cartesian(ylim=c(0,0.008))+
      ylab("Deprivation Index")+theme_pubr()+
      scale_color_manual(values=cbPalette)+
      xlab("Iteration")
  }else{
    gdepri<-NULL
  }
  
  
  gTime<-ggplot(data=Coeff_in_function)+
    # geom_rect(aes(xmin=-Inf, xmax=Inf, 
    #               ymin=medians_pervar$lci[row.names(medians_pervar)=="s(as.numeric(date)).1"], 
    #               ymax=medians_pervar$uci[row.names(medians_pervar)=="s(as.numeric(date)).1"], 
    #               fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
    #           col="white")+
    # geom_rect(aes(xmin=-Inf, xmax=Inf, 
    #               ymin=medians_pervar$lci[row.names(medians_pervar)=="s(as.numeric(date)).2"], 
    #               ymax=medians_pervar$uci[row.names(medians_pervar)=="s(as.numeric(date)).2"], 
    #               fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
    #           col="white")+
    # geom_rect(aes(xmin=-Inf, xmax=Inf, 
  #               ymin=medians_pervar$lci[row.names(medians_pervar)=="s(as.numeric(date)).3"], 
  #               ymax=medians_pervar$uci[row.names(medians_pervar)=="s(as.numeric(date)).3"], 
  #               fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
  #           col="white")+
  # geom_rect(aes(xmin=-Inf, xmax=Inf, 
  #               ymin=medians_pervar$lci[row.names(medians_pervar)=="s(as.numeric(date)).4"], 
  #               ymax=medians_pervar$uci[row.names(medians_pervar)=="s(as.numeric(date)).4"], 
  #               fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
  #           col="white")+
  # geom_rect(aes(xmin=-Inf, xmax=Inf, 
  #               ymin=medians_pervar$lci[row.names(medians_pervar)=="s(as.numeric(date)).5"], 
  #               ymax=medians_pervar$uci[row.names(medians_pervar)=="s(as.numeric(date)).5"], 
  #               fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
  #           col="white")+
  # geom_rect(aes(xmin=-Inf, xmax=Inf, 
  #               ymin=medians_pervar$lci[row.names(medians_pervar)=="s(as.numeric(date)).6"], 
  #               ymax=medians_pervar$uci[row.names(medians_pervar)=="s(as.numeric(date)).6"], 
  #               fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
  #           col="white")+
  # geom_rect(aes(xmin=-Inf, xmax=Inf, 
  #               ymin=medians_pervar$lci[row.names(medians_pervar)=="s(as.numeric(date)).7"], 
  #               ymax=medians_pervar$uci[row.names(medians_pervar)=="s(as.numeric(date)).7"], 
  #               fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
  #           col="white")+
  # geom_rect(aes(xmin=-Inf, xmax=Inf, 
  #               ymin=medians_pervar$lci[row.names(medians_pervar)=="s(as.numeric(date)).8"], 
  #               ymax=medians_pervar$uci[row.names(medians_pervar)=="s(as.numeric(date)).8"], 
  #               fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
  #           col="white")+
  # geom_rect(aes(xmin=-Inf, xmax=Inf, 
  #               ymin=medians_pervar$lci[row.names(medians_pervar)=="s(as.numeric(date)).9"], 
  #               ymax=medians_pervar$uci[row.names(medians_pervar)=="s(as.numeric(date)).9"], 
  #               fill="95% Confidence Interval"), alpha=0.5, fill="light grey",
  #           col="white")+
  geom_line(aes(x=k, y=`s(as.numeric(date)).1`, col="Smooth_over_time_1"))+
    #          geom_smooth(aes(x=k, y=`s(as.numeric(date)).1`, col="Smooth_over_time_1"), se = FALSE)+
    geom_line(aes(x=k, y=`s(as.numeric(date)).2`, col="Smooth_over_time_2"))+
    #         geom_smooth(aes(x=k, y=`s(as.numeric(date)).2`, col="Smooth_over_time_2"), se = FALSE)+
    geom_line(aes(x=k, y=`s(as.numeric(date)).3`, col="Smooth_over_time_3"))+
    #        geom_smooth(aes(x=k, y=`s(as.numeric(date)).3`, col="Smooth_over_time_3"), se = FALSE)+
    geom_line(aes(x=k, y=`s(as.numeric(date)).4`, col="Smooth_over_time_4"))+
    #        geom_smooth(aes(x=k, y=`s(as.numeric(date)).4`, col="Smooth_over_time_4"), se = FALSE)+
    geom_line(aes(x=k, y=`s(as.numeric(date)).5`, col="Smooth_over_time_5"))+
    #        geom_smooth(aes(x=k, y=`s(as.numeric(date)).5`, col="Smooth_over_time_5"), se = FALSE)+
    geom_line(aes(x=k, y=`s(as.numeric(date)).6`, col="Smooth_over_time_6"))+
    #        geom_smooth(aes(x=k, y=`s(as.numeric(date)).6`, col="Smooth_over_time_6"), se = FALSE)+
    geom_line(aes(x=k, y=`s(as.numeric(date)).7`, col="Smooth_over_time_7"))+
    #        geom_smooth(aes(x=k, y=`s(as.numeric(date)).7`, col="Smooth_over_time_7"), se = FALSE)+
    geom_line(aes(x=k, y=`s(as.numeric(date)).8`, col="Smooth_over_time_8"))+
    #        geom_smooth(aes(x=k, y=`s(as.numeric(date)).8`, col="Smooth_over_time_8"), se = FALSE)+
    geom_line(aes(x=k, y=`s(as.numeric(date)).9`, col="Smooth_over_time_9"))+
    #        geom_smooth(aes(x=k, y=`s(as.numeric(date)).9`, col="Smooth_over_time_9"), se = FALSE)+
    ylab("Smoothing coefficients")+
    geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)=="s(as.numeric(date)).1"], 
                   col="Smooth_over_time_1"))+
    geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)=="s(as.numeric(date)).2"], 
                   col="Smooth_over_time_2"))+
    geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)=="s(as.numeric(date)).3"], 
                   col="Smooth_over_time_3"))+
    geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)=="s(as.numeric(date)).4"], 
                   col="Smooth_over_time_4"))+
    geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)=="s(as.numeric(date)).5"], 
                   col="Smooth_over_time_5"))+
    geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)=="s(as.numeric(date)).6"], 
                   col="Smooth_over_time_6"))+
    geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)=="s(as.numeric(date)).7"], 
                   col="Smooth_over_time_7"))+
    geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)=="s(as.numeric(date)).8"], 
                   col="Smooth_over_time_8"))+
    geom_hline(aes(yintercept=medians_pervar$medians[row.names(medians_pervar)=="s(as.numeric(date)).9"], 
                   col="Smooth_over_time_9"))+
    #coord_cartesian(ylim=c(-1,1))+
    ggtitle("Smooth time coeefficients", 
            paste(label, "patients"))+theme_pubr()+
    xlab("Iteration")
  
  ################## Smooths over time and space #########################
  ######## Time ########
  
  timeandspace<-Coeff_in_[k:nrow(Coeff_in_),]
  
  Coef_time<-matrix(rep(unlist(ifelse(grepl("date", names(timeandspace)),1,0)), 
                        nrow(timeandspace)), ncol=ncol(timeandspace), 
                    byrow = TRUE)*timeandspace
  
  
  vars_time<-cbind("date"=model_for_basis$model$`as.numeric(date)`, 
                   as.data.frame(t(as.matrix(Coef_time)%*%t(as.matrix(
                     model.matrix(model_for_basis))))))%>%
    gather(key="Runs", value="Smooth", -date)%>%
    mutate(Runs=as.factor(Runs))
  
  avg_time<-cbind("date"=training_data$date , 
                  "Avg"=as.data.frame(t(apply(Coef_time,2, median)%*%t(
                    as.matrix(model.matrix(model_for_basis))))))
  
  time_smoother<-ggplot()+
    geom_line(data=vars_time, aes(x=as.Date(date), y=Smooth, group=Runs), 
              alpha=0.03, size=1.01)+
    geom_line(aes(x=date, y=V1), data=avg_time, col=cbPalette[1], size=1.03)+
    xlab("Date")+ylab("Smoothing over time")+
    ggtitle(paste(label, "patients"), "Smooth over time")+theme_pubr()
  
  ######## Space ########
  load("2. Data/0. Data Extra/district_data.RData")
  load("2. Data/0. Data Extra/bundesland_data.RData")
  
  district_data$bundesland = factor(district_data$bundesland)
  ##############################################################################
  
  # offset<-matrix(0, nrow=nrow(training_data), ncol = (length(model_for_space)-k+1))
  # if(label=="Outgoing"){
  #   for(i in k:(length(model_for_space))){
  #     offset[,(i-k)+1]<-model_for_space[[i]]$R_mod$offset
  #   }
  # }
  

  Coef_lat<-cbind("districtId"=training_data$districtId, 
                  as.data.frame(t(as.matrix(matrix(rep(unlist(ifelse(grepl("lat",
                                            names(timeandspace)),1,0)), 
                                                       nrow(timeandspace)), 
                                                   ncol=ncol(timeandspace), 
                                                   byrow = TRUE)*
                                              timeandspace)%*%t(
                                                model.matrix(model_for_basis)))))%>%
    ####-offset check!
    gather(key="Runs", value="Smooth", -districtId)%>%
    mutate(Runs=as.factor(Runs))%>%
    dplyr::group_by(districtId)%>%
    dplyr::summarize(Median=median(Smooth), Quant025=quantile(Smooth, probs=0.025), 
                     Quant975=quantile(Smooth, probs=0.975))
  
  
  
  district_data = st_set_crs(district_data,value = "EPSG:4326")
  bundesland_data = st_set_crs(bundesland_data,value = "EPSG:4326")
  
  district_data<-merge(district_data, Coef_lat, by.x="lk_id", by.y="districtId", 
                       all.y=TRUE)
  
  
  long_lat1<-ggplot() +
    geom_sf(data = district_data, aes(fill=Median)) +
    theme_pubr() +
    geom_sf(data = bundesland_data, aes(), col = "black", alpha = 0.00001) +
    theme(axis.ticks = element_blank(),
          axis.text =  element_blank(),
          strip.text.y = element_text(size = 17),
          strip.text.x = element_text(size = 17),
          axis.line =  element_blank(),plot.title = element_text(hjust = 0.5)) +
    scale_fill_gradient2(name = "", # option="G"
                         low = cbPalette[6], high = cbPalette[7],na.value="grey",
                         guide = guide_colourbar( barwidth = 8), limits=c(-1.2,
                                                                          1.2))+
    theme(legend.position="bottom")+
    ggtitle(paste(label, "patients"), "Smoothing over space")
  
  long_lat025<-ggplot() +
    geom_sf(data = district_data, aes(fill=Quant025)) +
    theme_pubr() +
    geom_sf(data = bundesland_data, aes(), col = "black", alpha = 0.00001) +
    theme(axis.ticks = element_blank(),
          axis.text =  element_blank(),
          strip.text.y = element_text(size = 17),
          strip.text.x = element_text(size = 17),
          axis.line =  element_blank(),plot.title = element_text(hjust = 0.5)) +
    scale_fill_gradient2(name = "", # option="G"
                         low = cbPalette[6], high = cbPalette[7],na.value="grey",
                         guide = guide_colourbar( barwidth = 8))+     
    theme(legend.position="bottom")+
    ggtitle(paste(label, "patients"), "2.5% Quantile of Smoothing over space")
  
  
  long_lat975<-ggplot() +
    geom_sf(data = district_data, aes(fill=Quant975)) +
    theme_pubr() +
    geom_sf(data = bundesland_data, aes(), col = "black", alpha = 0.00001) +
    theme(axis.ticks = element_blank(),
          axis.text =  element_blank(),
          strip.text.y = element_text(size = 17),
          strip.text.x = element_text(size = 17),
          axis.line =  element_blank(),plot.title = element_text(hjust = 0.5)) +
    scale_fill_gradient2(name = "", # option="G"
                         low = cbPalette[6], high = cbPalette[7],na.value="grey",
                         guide = guide_colourbar( barwidth = 8))+   
    theme(legend.position="bottom")+
    ggtitle(paste(label, "patients"), "97.5% Quantile Smoothing over space")
  
  
  if(wdepr=="yes"){
    plotgrid_in<-plot_grid(gInter, gWeek, gdepri, 
                           gInfect, ncol=1, align="v")
  }else{
    plotgrid_in<-plot_grid(gInter, gWeek, 
                           gInfect, ncol=1, align="v")
  }
  plotgrid_list<-list("Coefficient"=plotgrid_in, "Time"=time_smoother, "Long_Lat_median"=long_lat1,
                      "Long_Lat025"=long_lat025, "Long_Lat975"=long_lat975)
  
  return(plotgrid_list)
}



########## Summary functions of response by date or district for all runs ###########################

Summary_by_run<-function(data_response=Incoming_response, groupId=Incoming_response$date){
  data_response<-cbind(data_response, groupId)
  Summ_data<-data_response%>%
    dplyr::group_by(groupId)%>%
    dplyr::summarize(Run_1=sum(Run_1), Run_2=sum(Run_2), Run_3=sum(Run_3), 
                     Run_4=sum(Run_4), Run_5=sum(Run_5), Run_6=sum(Run_6), 
                     Run_7=sum(Run_7), Run_8=sum(Run_8), Run_9=sum(Run_9), 
                     Run_10=sum(Run_10), Run_11=sum(Run_11), Run_12=sum(Run_12), 
                     Run_13=sum(Run_13), Run_14=sum(Run_14), Run_15=sum(Run_15), 
                     Run_16=sum(Run_16), Run_17=sum(Run_17), Run_18=sum(Run_18), 
                     Run_19=sum(Run_19), Run_20=sum(Run_20), Run_21=sum(Run_21), 
                     Run_22=sum(Run_22), Run_23=sum(Run_23), Run_24=sum(Run_24), 
                     Run_25=sum(Run_25), Run_26=sum(Run_26), Run_27=sum(Run_27), 
                     Run_28=sum(Run_28), Run_29=sum(Run_29), Run_30=sum(Run_30), 
                     Run_31=sum(Run_31), Run_32=sum(Run_32), Run_33=sum(Run_33), 
                     Run_34=sum(Run_34), Run_35=sum(Run_35), Run_36=sum(Run_36), 
                     Run_37=sum(Run_37), Run_38=sum(Run_38), Run_39=sum(Run_39), 
                     Run_40=sum(Run_40), Run_41=sum(Run_41), Run_42=sum(Run_42),
                     Run_43=sum(Run_43), Run_44=sum(Run_44), Run_45=sum(Run_45),
                     Run_46=sum(Run_46), Run_47=sum(Run_47), Run_48=sum(Run_48), 
                     Run_49=sum(Run_49), Run_50=sum(Run_50), Run_51=sum(Run_51), 
                     Run_52=sum(Run_52), Run_53=sum(Run_53), Run_54=sum(Run_54), 
                     Run_55=sum(Run_55), Run_56=sum(Run_56), Run_57=sum(Run_57), 
                     Run_58=sum(Run_58), Run_59=sum(Run_59), Run_60=sum(Run_60), 
                     Run_61=sum(Run_61), Run_62=sum(Run_62), Run_63=sum(Run_63), 
                     Run_64=sum(Run_64), Run_65=sum(Run_65), Run_66=sum(Run_66), 
                     Run_67=sum(Run_67), Run_68=sum(Run_68), Run_69=sum(Run_69), 
                     Run_70=sum(Run_70), Run_71=sum(Run_71), Run_72=sum(Run_72), 
                     Run_73=sum(Run_73), Run_74=sum(Run_74), Run_75=sum(Run_75), 
                     Run_76=sum(Run_76), Run_77=sum(Run_77), Run_78=sum(Run_78), 
                     Run_79=sum(Run_79), Run_80=sum(Run_80), Run_81=sum(Run_81), 
                     Run_82=sum(Run_82), Run_83=sum(Run_83), Run_84=sum(Run_84), 
                     Run_85=sum(Run_85), Run_86=sum(Run_86), Run_87=sum(Run_87), 
                     Run_88=sum(Run_88), Run_89=sum(Run_89), Run_90=sum(Run_90), 
                     Run_91=sum(Run_91), Run_92=sum(Run_92), Run_93=sum(Run_93), 
                     Run_94=sum(Run_94), Run_95=sum(Run_95), Run_96=sum(Run_96), 
                     Run_97=sum(Run_97), Run_98=sum(Run_98), Run_99=sum(Run_99), 
                     Run_100=sum(Run_100), Run_101=sum(Run_101), Run_102=sum(Run_102), 
                     Run_103=sum(Run_103), Run_104=sum(Run_104), Run_105=sum(Run_105), 
                     Run_106=sum(Run_106), Run_107=sum(Run_107), Run_108=sum(Run_108), 
                     Run_109=sum(Run_109), Run_110=sum(Run_110), Run_111=sum(Run_111), 
                     Run_112=sum(Run_112), Run_113=sum(Run_113), Run_114=sum(Run_114), 
                     Run_115=sum(Run_115), Run_116=sum(Run_116), Run_117=sum(Run_117), 
                     Run_118=sum(Run_118), Run_119=sum(Run_119), Run_120=sum(Run_120), 
                     Run_121=sum(Run_121), Run_122=sum(Run_122), Run_123=sum(Run_123), 
                     Run_124=sum(Run_124), Run_125=sum(Run_125), Run_126=sum(Run_126), 
                     Run_127=sum(Run_127), Run_128=sum(Run_128), Run_129=sum(Run_129), 
                     Run_130=sum(Run_130), Run_131=sum(Run_131), Run_132=sum(Run_132), 
                     Run_133=sum(Run_133), Run_134=sum(Run_134), Run_135=sum(Run_135), 
                     Run_136=sum(Run_136), Run_137=sum(Run_137), Run_138=sum(Run_138), 
                     Run_139=sum(Run_139), Run_140=sum(Run_140), Run_141=sum(Run_141), 
                     Run_142=sum(Run_142), Run_143=sum(Run_143), Run_144=sum(Run_144), 
                     Run_145=sum(Run_145), Run_146=sum(Run_146), Run_147=sum(Run_147), 
                     Run_148=sum(Run_148), Run_149=sum(Run_149), Run_150=sum(Run_150), 
                     Run_151=sum(Run_151), Run_152=sum(Run_152), Run_153=sum(Run_153), 
                     Run_154=sum(Run_154), Run_155=sum(Run_155), Run_156=sum(Run_156), 
                     Run_157=sum(Run_157), Run_158=sum(Run_158), Run_159=sum(Run_159), 
                     Run_160=sum(Run_160), Run_161=sum(Run_161), Run_162=sum(Run_162), 
                     Run_163=sum(Run_163), Run_164=sum(Run_164), Run_165=sum(Run_165), 
                     Run_166=sum(Run_166), Run_167=sum(Run_167), Run_168=sum(Run_168), 
                     Run_169=sum(Run_169), Run_170=sum(Run_170), Run_171=sum(Run_171), 
                     Run_172=sum(Run_172), Run_173=sum(Run_173), Run_174=sum(Run_174), 
                     Run_175=sum(Run_175), Run_176=sum(Run_176), Run_177=sum(Run_177), 
                     Run_178=sum(Run_178), Run_179=sum(Run_179), Run_180=sum(Run_180), 
                     Run_181=sum(Run_181), Run_182=sum(Run_182), Run_183=sum(Run_183), 
                     Run_184=sum(Run_184), Run_185=sum(Run_185), Run_186=sum(Run_186), 
                     Run_187=sum(Run_187), Run_188=sum(Run_188), Run_189=sum(Run_189), 
                     Run_190=sum(Run_190), Run_191=sum(Run_191), Run_192=sum(Run_192), 
                     Run_193=sum(Run_193), Run_194=sum(Run_194), Run_195=sum(Run_195), 
                     Run_196=sum(Run_196), Run_197=sum(Run_197), Run_198=sum(Run_198), 
                     Run_199=sum(Run_199), Run_200=sum(Run_200), Run_201=sum(Run_201), 
                     Run_202=sum(Run_202), Run_203=sum(Run_203), Run_204=sum(Run_204), 
                     Run_205=sum(Run_205), Run_206=sum(Run_206), Run_207=sum(Run_207), 
                     Run_208=sum(Run_208), Run_209=sum(Run_209), Run_210=sum(Run_210), 
                     Run_211=sum(Run_211), Run_212=sum(Run_212), Run_213=sum(Run_213), 
                     Run_214=sum(Run_214), Run_215=sum(Run_215), Run_216=sum(Run_216), 
                     Run_217=sum(Run_217), Run_218=sum(Run_218), Run_219=sum(Run_219), 
                     Run_220=sum(Run_220), Run_221=sum(Run_221), Run_222=sum(Run_222), 
                     Run_223=sum(Run_223), Run_224=sum(Run_224), Run_225=sum(Run_225), 
                     Run_226=sum(Run_226), Run_227=sum(Run_227), Run_228=sum(Run_228), 
                     Run_229=sum(Run_229), Run_230=sum(Run_230), Run_231=sum(Run_231), 
                     Run_232=sum(Run_232), Run_233=sum(Run_233), Run_234=sum(Run_234), 
                     Run_235=sum(Run_235), Run_236=sum(Run_236), Run_237=sum(Run_237), 
                     Run_238=sum(Run_238), Run_239=sum(Run_239), Run_240=sum(Run_240), 
                     Run_241=sum(Run_241), Run_242=sum(Run_242), Run_243=sum(Run_243), 
                     Run_244=sum(Run_244), Run_245=sum(Run_245), Run_246=sum(Run_246), 
                     Run_247=sum(Run_247), Run_248=sum(Run_248), Run_249=sum(Run_249), 
                     Run_250=sum(Run_250), Run_251=sum(Run_251), Run_252=sum(Run_252), 
                     Run_253=sum(Run_253), Run_254=sum(Run_254), Run_255=sum(Run_255), 
                     Run_256=sum(Run_256), Run_257=sum(Run_257), Run_258=sum(Run_258), 
                     Run_259=sum(Run_259), Run_260=sum(Run_260), Run_261=sum(Run_261), 
                     Run_262=sum(Run_262), Run_263=sum(Run_263), Run_264=sum(Run_264),
                     Run_265=sum(Run_265), Run_266=sum(Run_266), Run_267=sum(Run_267), 
                     Run_268=sum(Run_268), Run_269=sum(Run_269), Run_270=sum(Run_270), 
                     Run_271=sum(Run_271), Run_272=sum(Run_272), Run_273=sum(Run_273), 
                     Run_274=sum(Run_274), Run_275=sum(Run_275), Run_276=sum(Run_276),
                     Run_277=sum(Run_277), Run_278=sum(Run_278), Run_279=sum(Run_279),
                     Run_280=sum(Run_280), Run_281=sum(Run_281), Run_282=sum(Run_282), 
                     Run_283=sum(Run_283), Run_284=sum(Run_284), Run_285=sum(Run_285), 
                     Run_286=sum(Run_286), Run_287=sum(Run_287), Run_288=sum(Run_288), 
                     Run_289=sum(Run_289), Run_290=sum(Run_290), Run_291=sum(Run_291), 
                     Run_292=sum(Run_292), Run_293=sum(Run_293), Run_294=sum(Run_294), 
                     Run_295=sum(Run_295), Run_296=sum(Run_296), Run_297=sum(Run_297), 
                     Run_298=sum(Run_298), Run_299=sum(Run_299), Run_300=sum(Run_300), 
                     Run_301=sum(Run_301), Run_302=sum(Run_302), Run_303=sum(Run_303), 
                     Run_304=sum(Run_304), Run_305=sum(Run_305), Run_306=sum(Run_306), 
                     Run_307=sum(Run_307), Run_308=sum(Run_308), Run_309=sum(Run_309), 
                     Run_310=sum(Run_310), Run_311=sum(Run_311), Run_312=sum(Run_312), 
                     Run_313=sum(Run_313), Run_314=sum(Run_314), Run_315=sum(Run_315), 
                     Run_316=sum(Run_316), Run_317=sum(Run_317), Run_318=sum(Run_318), 
                     Run_319=sum(Run_319), Run_320=sum(Run_320), Run_321=sum(Run_321), 
                     Run_322=sum(Run_322), Run_323=sum(Run_323), Run_324=sum(Run_324), 
                     Run_325=sum(Run_325), Run_326=sum(Run_326), Run_327=sum(Run_327), 
                     Run_328=sum(Run_328), Run_329=sum(Run_329), Run_330=sum(Run_330), 
                     Run_331=sum(Run_331), Run_332=sum(Run_332), Run_333=sum(Run_333), 
                     Run_334=sum(Run_334), Run_335=sum(Run_335), Run_336=sum(Run_336), 
                     Run_337=sum(Run_337), Run_338=sum(Run_338), Run_339=sum(Run_339), 
                     Run_340=sum(Run_340), Run_341=sum(Run_341), Run_342=sum(Run_342), 
                     Run_343=sum(Run_343), Run_344=sum(Run_344), Run_345=sum(Run_345), 
                     Run_346=sum(Run_346), Run_347=sum(Run_347), Run_348=sum(Run_348), 
                     Run_349=sum(Run_349), Run_350=sum(Run_350), Run_351=sum(Run_351), 
                     Run_352=sum(Run_352), Run_353=sum(Run_353), Run_354=sum(Run_354), 
                     Run_355=sum(Run_355), Run_356=sum(Run_356), Run_357=sum(Run_357), 
                     Run_358=sum(Run_358), Run_359=sum(Run_359), Run_360=sum(Run_360), 
                     Run_361=sum(Run_361), Run_362=sum(Run_362), Run_363=sum(Run_363), 
                     Run_364=sum(Run_364), Run_365=sum(Run_365), Run_366=sum(Run_366), 
                     Run_367=sum(Run_367), Run_368=sum(Run_368), Run_369=sum(Run_369), 
                     Run_370=sum(Run_370), Run_371=sum(Run_371), Run_372=sum(Run_372), 
                     Run_373=sum(Run_373), Run_374=sum(Run_374), Run_375=sum(Run_375), 
                     Run_376=sum(Run_376), Run_377=sum(Run_377), Run_378=sum(Run_378), 
                     Run_379=sum(Run_379), Run_380=sum(Run_380), Run_381=sum(Run_381), 
                     Run_382=sum(Run_382), Run_383=sum(Run_383), Run_384=sum(Run_384), 
                     Run_385=sum(Run_385), Run_386=sum(Run_386), Run_387=sum(Run_387), 
                     Run_388=sum(Run_388), Run_389=sum(Run_389), Run_390=sum(Run_390), 
                     Run_391=sum(Run_391), Run_392=sum(Run_392), Run_393=sum(Run_393), 
                     Run_394=sum(Run_394), Run_395=sum(Run_395), Run_396=sum(Run_396), 
                     Run_397=sum(Run_397), Run_398=sum(Run_398), Run_399=sum(Run_399), 
                     Run_400=sum(Run_400), Run_401=sum(Run_401), Run_402=sum(Run_402), 
                     Run_403=sum(Run_403), Run_404=sum(Run_404), Run_405=sum(Run_405), 
                     Run_406=sum(Run_406), Run_407=sum(Run_407), Run_408=sum(Run_408), 
                     Run_409=sum(Run_409), Run_410=sum(Run_410), Run_411=sum(Run_411), 
                     Run_412=sum(Run_412), Run_413=sum(Run_413), Run_414=sum(Run_414), 
                     Run_415=sum(Run_415), Run_416=sum(Run_416), Run_417=sum(Run_417), 
                     Run_418=sum(Run_418), Run_419=sum(Run_419), Run_420=sum(Run_420), 
                     Run_421=sum(Run_421), Run_422=sum(Run_422), Run_423=sum(Run_423), 
                     Run_424=sum(Run_424), Run_425=sum(Run_425), Run_426=sum(Run_426), 
                     Run_427=sum(Run_427), Run_428=sum(Run_428), Run_429=sum(Run_429), 
                     Run_430=sum(Run_430), Run_431=sum(Run_431), Run_432=sum(Run_432), 
                     Run_433=sum(Run_433), Run_434=sum(Run_434), Run_435=sum(Run_435), 
                     Run_436=sum(Run_436), Run_437=sum(Run_437), Run_438=sum(Run_438), 
                     Run_439=sum(Run_439), Run_440=sum(Run_440), Run_441=sum(Run_441), 
                     Run_442=sum(Run_442), Run_443=sum(Run_443), Run_444=sum(Run_444), 
                     Run_445=sum(Run_445), Run_446=sum(Run_446), Run_447=sum(Run_447), 
                     Run_448=sum(Run_448), Run_449=sum(Run_449), Run_450=sum(Run_450), 
                     Run_451=sum(Run_451), Run_452=sum(Run_452), Run_453=sum(Run_453), 
                     Run_454=sum(Run_454), Run_455=sum(Run_455), Run_456=sum(Run_456), 
                     Run_457=sum(Run_457), Run_458=sum(Run_458), Run_459=sum(Run_459), 
                     Run_460=sum(Run_460), Run_461=sum(Run_461), Run_462=sum(Run_462), 
                     Run_463=sum(Run_463), Run_464=sum(Run_464), Run_465=sum(Run_465), 
                     Run_466=sum(Run_466), Run_467=sum(Run_467), Run_468=sum(Run_468), 
                     Run_469=sum(Run_469), Run_470=sum(Run_470), Run_471=sum(Run_471), 
                     Run_472=sum(Run_472), Run_473=sum(Run_473), Run_474=sum(Run_474), 
                     Run_475=sum(Run_475), Run_476=sum(Run_476), Run_477=sum(Run_477), 
                     Run_478=sum(Run_478), Run_479=sum(Run_479), Run_480=sum(Run_480), 
                     Run_481=sum(Run_481), Run_482=sum(Run_482), Run_483=sum(Run_483), 
                     Run_484=sum(Run_484), Run_485=sum(Run_485), Run_486=sum(Run_486), 
                     Run_487=sum(Run_487), Run_488=sum(Run_488), Run_489=sum(Run_489), 
                     Run_490=sum(Run_490), Run_491=sum(Run_491), Run_492=sum(Run_492), 
                     Run_493=sum(Run_493), Run_494=sum(Run_494), Run_495=sum(Run_495), 
                     Run_496=sum(Run_496), Run_497=sum(Run_497), Run_498=sum(Run_498), 
                     Run_499=sum(Run_499), Run_500=sum(Run_500))
  return(Summ_data)
  
}
