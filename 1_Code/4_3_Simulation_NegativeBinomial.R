
######################################################
######### Setting path to wd  ########################
######################################################

function_folder <- "OccupancyDuration-main"
functionpath<-substr(dirname(rstudioapi::getSourceEditorContext()$path), 
                     1, unlist(gregexpr(function_folder, 
                                        dirname(rstudioapi::getSourceEditorContext()$path)))+(nchar(function_folder)-1))

setwd(functionpath)


######################################################
######### Libraries ##################################
######################################################

source(paste0(functionpath, "/1_Code/0_Libraries.R"))


######################################################
######### Homemade Functions #########################
######################################################

source(paste0(functionpath, "/1_Code/0_43_Home_Functions_Sim_NegativeBinomial.R"))



############ Application ###############################################################################

######################################################
############### Set up Ground Truth ##################
######################################################


w_true_Grid<-matrix(0, 1, 10)
w_true_Grid[1,]<-exp(0.4*seq(10,1, by=-1))/sum(exp(0.4*seq(10,1, by=-1)))

the<-c(0.5, 1, 5, 10)
name_the<-c("05", "1", "5", "10")


StartT<-Sys.time()

for(ta in 1:length(the)){
  

######################################################
############### Simulating data list #################
######################################################

datalist<-list()
for(j in 311){
  for(i in 1:nrow(w_true_Grid)){
    datalist[[length(datalist)+1]]<-sim_data_offset(w=w_true_Grid[i,],
                                                    districtId_sim=1:200, 
                                                    date_sim=1:(200+ncol(w_true_Grid)),
                                                    misspecified = TRUE,
                                                    theta_mis=the[ta],
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

train_data_l=datalist_real[[1]] 
train_data_l_real=datalist[[1]]

runs1=400
max_lag=12
l<-"sq"

train_data_l<-train_data_l%>%
  arrange(districtId, date)%>%
  mutate(Incoming=sample(10:15,1),
         Outgoing=sample(16:20,1))

formula_in=as.formula(Incoming ~ M+N) 
I_Mstep<-mgcv::gam(formula_in, data=train_data_l, 
                   family="poisson")


w_start<-c()
for(i in (1:(max_lag-1))){
  w_start[i]<-1/i  
}
w_start[1:i]<-w_start[1:i]/sum(w_start[1:i])
w_start[max_lag]=1-sum(w_start)
w_start<-sort(w_start, decreasing=TRUE)

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
  
  w_start<-omega_org<-round(w_start, 10)
  
  if(length(omega)>max_lag){
    print("The omega is longer than the maximum length 
                  of stay you have chosen")}
  
  
  
  ##### Incoming intensity parameter ####
  
  sim_start<-train_data_l_start$data[,c("districtId", "date", 
                                        "Incoming", 
                                        "Outgoing", 
                                        "Diff", "M", "N")]
  
  
  ## Step 3 ######## Jack Knife Adjustment  ##################################
  
  Jack_Knife_data<-sim_data_jk(w=omega_org,
                               ss=123, 
                               deterministic=TRUE,
                               coef=coef(I_Mstep),
                               covariates=sim_start[,
                                                    c("date", "districtId","M", "N")])
  
  train_JK_temp<-Jack_Knife_data[,
                                 c("date", "districtId","M", "N", "Diff")]
  
  ## Step 4 ######## Jack Knife Adjustment Simulation  #######################
  
  start_sim_jk<-rdiffpois_array_paper(data_rdiff=train_JK_temp,
                                      IMSTEP=I_Mstep,
                                      w_sim=omega_org)
  
  ## Step 5 ######## Jack Knife Adjustment Omega  ############################
  w_start_jk<-omega_org
  ##### Seesaw to find best omega for JK ####
  for(i in 1:100){
    new_sol_jk<-optimize_Like(w_opt_prev=w_start_jk,
                              data_like_opp=start_sim_jk$data_all)
    
    if(sum(((round(c(new_sol_jk$w_sol, 
                     1-sum(new_sol_jk$w_sol)), 10)-w_start_jk)^2))<0.0001){
      VCOV_w<-solve(new_sol_jk$Infor)
      break
    }else{
      w_start_jk<-round(c(new_sol_jk$w_sol, 
                          1-sum(new_sol_jk$w_sol)), 10)
    }
  }
  if(i==100){print("We have maxed out the runs to find omega")}
  
  omega_jk<-c(w_start_jk, rep(0, max_lag-length(w_start_jk)))
  
  ## Step 6 ######## Jack Knife Adjustment Omega  ############################
  
  biasadj<-bias_adj(w_og=omega_org,
                    w_jk=omega_jk,
                    max_l=max_lag,
                    func="sq")
  
  
  C_hat[z]<-coef(lm(I((omega_org-1/max_lag)^2)~I((omega_jk-1/max_lag)^2)-1))
  omega_adj<-1/max_lag+ifelse(omega_org<1/max_lag, -1, 1)*sqrt(C_hat[z]*(omega_org-1/max_lag)^2)
  omega_adj[omega_adj<0]=0
  omega_adj<-omega_adj/sum(omega_adj)
  ##########################
  
  w_start<-omega_adj   # AENDERN ZU OMEGA ADJ MIT JK ADJUSTMENT 
  ##########################
  ## Step 7 ######## Simulation Inner E-Step Omega  ############################
  
  Inner_E_Step<-rdiffpois_array_paper(data_rdiff=sim_start%>%
                                        dplyr::select(-c("Incoming", "Outgoing")),
                                      IMSTEP=I_Mstep,
                                      w_sim=omega_adj)
  
  ## Step 8 ######## M-Step for adjusted Simulated Data  #####################
  
  sim_start<-cbind(Inner_E_Step$data[,c("districtId", "date", "Incoming", 
                                        "Outgoing", "Diff", "M", "N")])
  
  
  I_Mstep<-mgcv::gam(formula_in, data=sim_start, 
                     family="poisson")
  
  
  ############################################################################
  
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
  
  omegas_adj[z,]<-omega_adj
  omegas_org[z,]<-omega_org
  omegas_jk[z,]<-omega_jk

  VCOV[z,,]<-new_sol$Infor
  
  if(z!=1){
    Coeff_in<-rbind(Coeff_in, as.data.frame(cbind(z, t(coef(I_Mstep)))))
    
  }else{
    Coeff_in<-as.data.frame(cbind(z, t(coef(I_Mstep))))
    
  }
  print(paste("Inner iteration", z, "out of", runs1, "is done"))
  
}

############ Results Plots ###############################################################################

######################################################
############## Variance for CI #######################
######################################################


varianceRuns<-200:400

coefratemed<-apply(Coeff_in[varianceRuns, 2:4], 2, median)
VCOV_Coef_Runs_b<-array(0, dim=c(length(varianceRuns), 
                                 ncol=ncol(vcov(modellist[[1]])), 
                                 nrow=nrow(vcov(modellist[[1]]))))
VCOV_Coef_Runs_s<-array(0, dim=c(length(varianceRuns), 
                                 ncol=ncol(vcov(modellist[[1]])), 
                                 nrow=nrow(vcov(modellist[[1]]))))

for(i in varianceRuns){
  VCOV_Coef_Runs_b[(i)-min(varianceRuns)+1,,]<-(as.numeric(Coeff_in[i, 2:4]-coefratemed))%*%
    t(as.numeric((Coeff_in[i, 2:4]-coefratemed)))
  VCOV_Coef_Runs_s[(i)-min(varianceRuns)+1,,]<-
    vcov(modellist[[i]])
}

write.csv(Coeff_in, file=paste0("3_Results/Coef_In_Mis_", name_the[ta], ".csv"))
save(modellist, file=paste0("3_Results/modellist_Mis", name_the[ta], ".R"))
write.csv(omegas_adj, file=paste0("3_Results/Exitrate_In_Mis_", name_the[ta], ".csv"))
save(fitlist, file=paste0("3_Results/fitlist_Mis", name_the[ta], ".R"))

######################################################
############## Incoming and Outgoing #################
######################################################

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
  full_join(datalist[[1]][, c("date", "districtId", "Incoming", "Outgoing")])

Max_tot<-max(Inc_Out_Data[, 3:6])
IncomTrue<-ggplot()+
  geom_line(aes(x=1:Max_tot, y=1:Max_tot, col="Perfect fit"), size=1.1)+
  geom_point(aes(y=Inc_Out_Data$Incoming_Est, 
                 x=Inc_Out_Data$Incoming, col="Incoming"), size=1.1, alpha=0.5)+
  theme_pubr()+
  labs(color=" ")+
  scale_colour_manual(values = brewer.pal(n = 5, name = "Blues")[4:6])+ 
  ggtitle("Incoming",  "theta=0.5")+
  ylab("Estimated Incoming")+ 
  xlab("True Incoming") 

OutgoTrue<-ggplot()+
  geom_line(aes(x=1:Max_tot, y=1:Max_tot, col="Perfect fit"), 
            size=1.1)+
  geom_point(aes(y=Inc_Out_Data$Outgoing_Est, 
                 x=Inc_Out_Data$Outgoing, col="Outgoing"), 
             size=1.1, alpha=0.5)+
  theme_pubr()+
  scale_colour_manual(values = brewer.pal(n = 5, name = "Blues")[4:6])+ 
  labs(color=" ")+
  ggtitle("Outgoing",  "theta=0.5")+
  ylab("Estimated Outgoing")+ 
  xlab("True Outgoing") 


Est_True<-plot_grid(IncomTrue, OutgoTrue, nrow=1)


ggsave(Est_True, file=paste0(functionpath, "/3_Results/Figure_", 10+ta, ".pdf"),
       width=8, height=4)


}

EndTime<-Sys.time()

################################################################################
############# Paper: Table 2  ##################################################
################################################################################

Coeff_in_05<-read.csv("3_Results/Coef_In_Mis_05.csv")
Coeff_in_1<-read.csv("3_Results/Coef_In_Mis_1.csv")
Coeff_in_5<-read.csv("3_Results/Coef_In_Mis_5.csv")
Coeff_in_10<-read.csv("3_Results/Coef_In_Mis_10.csv")
Coeff_in_Pois<-read.csv("3_Results/Coef_In_Pois.csv")


# Create the matrix
res <- as.data.frame(rbind(
  "True" = c(0.5, 1, 0.2),
  "Theta= 0.5" = apply(Coeff_in_05[200:400, 3:5], 2, median), 
  "Theta= 1" = apply(Coeff_in_1[200:400, 3:5], 2, median), 
  "Theta= 5" = apply(Coeff_in_5[200:400, 3:5], 2, median), 
  "Theta= 10" = apply(Coeff_in_10[200:400, 3:5], 2, median),
  "Poisson" = apply(Coeff_in_Pois[200:400, 3:5], 2, median))
  )

# Round the values for display
res_round <- round(res, 3)

write.csv(res_round, "3_Results/Table_2.csv")


################################################################################
############# Paper: Figure 5  ##################################################
################################################################################
Exitrate_in_Pois<-read.csv("3_Results/Exitrate_In_Pois")
Exitrate_in_05<-read.csv("3_Results/Exitrate_In_Mis_05.csv")
Exitrate_in_1<-read.csv("3_Results/Exitrate_In_Mis_1.csv")
Exitrate_in_5<-read.csv("3_Results/Exitrate_In_Mis_5.csv")
Exitrate_in_10<-read.csv("3_Results/Exitrate_In_Mis_10.csv")


# Create the matrix (copy your data here)
df <- rbind(
  "True"  = c(w_true_Grid[1, ], 0, 0),
  "Theta = 0.5" = apply(Exitrate_in_05[200:400, -1], 2, median),
  "Theta = 1"  = apply(Exitrate_in_1[200:400, -1], 2, median),
  "Theta = 5"  = apply(Exitrate_in_5[200:400, -1], 2, median),
  "Theta = 10" = apply(Exitrate_in_10[200:400, -1], 2, median),
  "Poisson" = apply(Exitrate_in_Pois[200:400, 3:5], 2, median)
)

# Normalize all but the first row
df[-1, ] <- df[-1, ] / rowSums(df[-1, ])

# Convert to data frame
df_long <- as.data.frame(df) %>%
  mutate(Coef = rownames(.)) %>%
  pivot_longer(-Coef, names_to = "Position", values_to = "Value") %>%
  mutate(Position = as.numeric(gsub("V", "", Position)))


# Define custom blue colours
blue_scale <- c(
  "True" = "#08306B",        # darkest blue
  "Poisson" = "#2171B5",
  "Theta = 0.5" = "#4292C6",
  "Theta = 1"   = "#6BAED6",
  "Theta = 5"   = "#9ECAE1",
  "Theta = 10"  = "#C6DBEF"  # lightest
)

# Create the plot
diff_leng_stay <- ggplot(df_long, aes(x = Position, y = Value, 
                                      group = Coef, 
                                      colour = Coef)) +
  theme_pubr() +
  geom_line(data = df_long %>% filter(Coef == "True"), size = 1.5) +
  geom_line(data = df_long %>% filter(Coef != "True"), size = 0.7) +
  geom_point() +
  scale_x_continuous(breaks = 1:12) +
  scale_colour_manual(values = blue_scale) +
  labs(x = "Lag", y = "Estimated exit rate") +
  ggtitle("Estimated exit rates", subtitle = "For different simulated data")


ggsave(diff_leng_stay, file="3_Results/Figure_5.pdf", width =6, height=3 )











