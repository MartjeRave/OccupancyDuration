
######################################################
######### Simulation Model Application- ##############
######### Exploration Bias Adjustment ################
######################################################

function_folder<-"OccupancyDuration"
functionpath<-substr(dirname(rstudioapi::getSourceEditorContext()$path), 
                     1, unlist(gregexpr(function_folder, 
                                        dirname(rstudioapi::getSourceEditorContext()$path)))+(nchar(function_folder)-1))
setwd(functionpath)

######################################################
######### Libraries ##################################
######################################################

source(paste0(functionpath, "/1. Code/0. Functions.R"))



######################################################
######### Homemade Functions #########################
######################################################

source(paste0(functionpath, "/1. Code/0.41 Home_Functions_Sim_Unadj.R"))





######################################################
##### Offset calculation through lagged incoming #####
######################################################


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
  
  modellist[[z]]<-I_Mstep
  
  omegas_org[z,]<-c(w_start)
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

################################################################################
################# Saving Figure 1  Paper #########################################
################################################################################


Exit_rates<-ggplot()+
  geom_segment(aes(x = 1, y = (apply(omegas_org[ varianceRuns,], 
                                     2, median))[1], xend = 1, yend = w_true_Grid[1,1], col="Bias"))+
  geom_segment(aes(x = 2, y = (apply(omegas_org[ varianceRuns,], 
                                     2, median))[2], xend = 2, yend = w_true_Grid[1,2], col="Bias"))+
  geom_segment(aes(x = 3, y = (apply(omegas_org[ varianceRuns,], 
                                     2, median))[3], xend = 3, yend = w_true_Grid[1,3], col="Bias"))+
  geom_segment(aes(x = 4, y = (apply(omegas_org[ varianceRuns,], 
                                     2, median))[4], xend = 4, yend = w_true_Grid[1,4], col="Bias"))+
  geom_segment(aes(x = 5, y = (apply(omegas_org[ varianceRuns,], 
                                     2, median))[5], xend = 5, yend = w_true_Grid[1,5], col="Bias"))+
  geom_segment(aes(x = 6, y = (apply(omegas_org[ varianceRuns,], 
                                     2, median))[6], xend = 6, yend = w_true_Grid[1,6], col="Bias"))+
  geom_segment(aes(x = 7, y = (apply(omegas_org[ varianceRuns,], 
                                     2, median))[7], xend = 7, yend = w_true_Grid[1,7], col="Bias"))+
  geom_segment(aes(x = 8, y = (apply(omegas_org[ varianceRuns,], 
                                     2, median))[8], xend = 8, yend = w_true_Grid[1,8], col="Bias"))+
  geom_segment(aes(x = 9, y = (apply(omegas_org[ varianceRuns,], 
                                     2, median))[9], xend = 9, yend = w_true_Grid[1,9], col="Bias"))+
  geom_segment(aes(x = 10, y = (apply(omegas_org[ varianceRuns,], 
                                      2, median))[10], xend = 10, yend = w_true_Grid[1,10], col="Bias"))+
  geom_segment(aes(x = 11, y = (apply(omegas_org[ varianceRuns,], 
                                      2, median))[11], xend = 11, yend = 0, col="Bias"))+
 geom_line(aes(x=1:11, y=1/12, col="1/12"), linetype=2)+
  geom_line(aes(x=1:(max_lag-1), y=(apply(omegas_org[ varianceRuns,], 
                                           2, median))[-12],
                 col="Est. Exit Rate"))+
  theme_pubr()+
  labs(color=" ")+
  ggtitle("a) Exit rate", " ")+
  geom_point(aes(x=1:(max_lag-1), y=(apply(omegas_org[ varianceRuns,], 
                                           2, median))[-12],
                 col="Est. Exit Rate"), shape=15, size=2.5)+
  geom_line(aes(x=1:11, y=c(w_true_Grid[1,], 0), col="Ground Truth"))+
  geom_point(aes(x=1:11, y=c(w_true_Grid[1,], 0), col="Ground Truth"), size=2.5)+
  ylab(expression(hat(omega)))+ 
  scale_colour_manual(values = c("grey", "darkgrey", brewer.pal(n = 5, name = "Blues")[c(3, 5, 6)]))+
  scale_fill_manual(values = c("grey", "darkgrey",brewer.pal(n = 5, name = "Blues")[c(3, 5, 6)]))+
  xlab("Length of stay")+
  scale_x_continuous(breaks = seq(1, max_lag, 1))+
  theme(legend.title = element_blank())

 
# ## Step 3 ######## Jack Knife Adjustment  ##################################
# 
# omega_org<-apply(omegas_org[ varianceRuns,], 2, median)
# 
# w_true_Grid[1,]
# Jack_Knife_data<-sim_data_jk(w=omega_org,
#                              ss=123, 
#                              deterministic=TRUE,
#                              coef=c(0.5, 1, 0.2),
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


C_hat<-coef(lm(I((c(w_true_Grid[1,],0,0)-1/max_lag)^2)~I((omega_org-1/max_lag)^2)-1))
omega_adj<-1/max_lag+ifelse(omega_org<1/max_lag, -1, 1)*sqrt((C_hat)*(omega_org-1/max_lag)^2)
omega_adj[omega_adj<0]=0
omega_adj<-omega_adj/sum(omega_adj)




Exit_rates_adj<-ggplot()+
  geom_segment(aes(x = 1, y = (apply(omegas_org[ varianceRuns,], 
                                     2, median))[1], xend = 1, yend = w_true_Grid[1,1], col="Bias"))+
  geom_segment(aes(x = 2, y = (apply(omegas_org[ varianceRuns,], 
                                     2, median))[2], xend = 2, yend = w_true_Grid[1,2], col="Bias"))+
  geom_segment(aes(x = 3, y = (apply(omegas_org[ varianceRuns,], 
                                     2, median))[3], xend = 3, yend = w_true_Grid[1,3], col="Bias"))+
  geom_segment(aes(x = 4, y = (apply(omegas_org[ varianceRuns,], 
                                     2, median))[4], xend = 4, yend = w_true_Grid[1,4], col="Bias"))+
  geom_segment(aes(x = 5, y = (apply(omegas_org[ varianceRuns,], 
                                     2, median))[5], xend = 5, yend = w_true_Grid[1,5], col="Bias"))+
  geom_segment(aes(x = 6, y = (apply(omegas_org[ varianceRuns,], 
                                     2, median))[6], xend = 6, yend = w_true_Grid[1,6], col="Bias"))+
  geom_segment(aes(x = 7, y = (apply(omegas_org[ varianceRuns,], 
                                     2, median))[7], xend = 7, yend = w_true_Grid[1,7], col="Bias"))+
  geom_segment(aes(x = 8, y = (apply(omegas_org[ varianceRuns,], 
                                     2, median))[8], xend = 8, yend = w_true_Grid[1,8], col="Bias"))+
  geom_segment(aes(x = 9, y = (apply(omegas_org[ varianceRuns,], 
                                     2, median))[9], xend = 9, yend = w_true_Grid[1,9], col="Bias"))+
  geom_segment(aes(x = 10, y = (apply(omegas_org[ varianceRuns,], 
                                      2, median))[10], xend = 10, yend = w_true_Grid[1,10], col="Bias"))+
  geom_segment(aes(x = 11, y = (apply(omegas_org[ varianceRuns,], 
                                      2, median))[11], xend = 11, yend = 0, col="Bias"))+
  geom_line(aes(x=1:11, y=1/12, col="1/12"), linetype=2)+
  geom_line(aes(x=1:(max_lag-1), y=(apply(omegas_org[ varianceRuns,], 
                                          2, median))[-12],
                col="Est. Exit Rate"))+
  theme_pubr()+
  labs(color=" ")+
  ggtitle("b) Exit rate (ex post adjusted)", expression(sqrt(hat(c))*"="*1.52))+
  geom_point(aes(x=1:(max_lag-1), y=(apply(omegas_org[ varianceRuns,], 
                                           2, median))[-12],
                 col="Est. Exit Rate"), shape=15, size=2.5)+
  geom_line(aes(x=1:11, y=c(w_true_Grid[1,], 0), col="Ground Truth"))+
  geom_point(aes(x=1:11, y=c(w_true_Grid[1,], 0), col="Ground Truth"), size=2.5)+
  ylab(expression(hat(omega)))+ 
  geom_line(aes(x=1:(max_lag-1), y=(omega_adj[-12]),
                col="Bias Adj. Exit Rate"))+
  geom_point(aes(x=1:(max_lag-1), y=(omega_adj[-12]),
                 col="Bias Adj. Exit Rate"), shape=17, size=2.5)+
  scale_colour_manual(values = c("grey", "darkgrey", "black", brewer.pal(n = 5, name = "Blues")[c(3, 5)]))+
  scale_fill_manual(values = c("grey", "darkgrey","black", brewer.pal(n = 5, name = "Blues")[c(3, 5)]))+
  xlab("Length of stay")+
  scale_x_continuous(breaks = seq(1, max_lag, 1))+
  theme(legend.title = element_blank())

Exit_rates_adj



TR_vs_Est<-ggplot() +
  geom_point(aes(
    y = (as.numeric(w_true_Grid[1,]) - 1/12)^2,
    x = (omega_org[1:10] - 1/12)^2,
    col = " "   # just a dummy legend key
  )) +
  geom_line(aes(
    y = C_hat * (omega_org[1:10] - 1/12)^2,
    x = (omega_org[1:10] - 1/12)^2,
    col = "sqrt(c)(omega)"   # another dummy key, will be replaced by label
  )) +
  theme_pubr()+
  labs(color=" ")+
  labs("sqrt(c)(omega)"=expression(sqrt(hat(c)) * (hat(omega)[l] - 1/12)^2))+
  ggtitle("c) Ground truth over estimated exit rate", "Square distance to 1/12")+
  ylab(expression((omega-1/12)^2))+ 
  scale_colour_manual(values = c( brewer.pal(n = 5, name = "Blues")[c(3, 5)]),
                      labels=c(" ", expression(hat(c) * (hat(omega)[l] - 1/12)^2)))+
  scale_fill_manual(values = c(brewer.pal(n = 5, name = "Blues")[c(3, 5)]))+
  xlab(expression((hat(omega)-1/12)^2))+
  theme(legend.title = element_blank())


################################################################################
############# Paper: Figure 2  #################################################
################################################################################

Figure_1<-plot_grid(plot_grid(Exit_rates, Exit_rates_adj, nrow=1),
          plot_grid(NA, TR_vs_Est, NA, nrow=1), nrow=2)
          
ggsave(Figure_1, file="3_Results/Figure_5.pdf", width =10, height=6)



