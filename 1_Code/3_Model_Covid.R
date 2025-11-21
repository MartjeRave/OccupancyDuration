######################################################
######### Model sEM applied to COVID-19 ##############
######################################################

function_folder<-"OccupancyDuration"
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

source(paste0(functionpath, "/1_Code/0_3_Home_Functions_Covid.R"))

######################################################
############## Data ##################################
######################################################

train_data_l<-readRDS("2_Data/training_data.rds")%>%
  mutate(date=as.Date(date))%>%
  filter(date<=as.Date("2021-12-31")&date>=as.Date("2021-08-01"))%>%
  mutate(Incoming=rpois(n=length(exp(0.3*G_5_7)), lambda=exp(0.1*G_5_7)), 
         Outgoing=rpois(n=length(exp(0.3*G_5_7)), lambda=exp(0.1*G_5_7)))

formula_in=as.formula(Incoming ~ G_4_7_lag_1+G_5_7_lag_1+G_6_7_lag_1+
                        weekday+s(long, lat)+s(as.numeric(date))) 

I_Mstep<-mgcv::gam(formula_in, data=train_data_l, 
                   family="poisson")



######################################################
############## EM Algorithm ##########################
######################################################

set.seed(123)

#### Set up
runs1=400
max_lag=30
z<-1



### start value
w_start<-exp(-0.1*1:max_lag)/sum(exp(-0.1*1:max_lag))


### Save Results
omegas_org<-rbind(#omegas_org, 
                  matrix(0, ncol=max_lag, nrow=runs1))
omegas_jk<-rbind(#omegas_jk, 
                 matrix(0, ncol=max_lag, nrow=runs1))
omegas_adj<-rbind(#omegas_adj, 
  matrix(0, ncol=max_lag, nrow=runs1))
modellist<-list()
IncomOutgo<-list()
VCOV<-array(0, dim=c(runs1, length(w_start)-1, length(w_start)-1))
logLike<-c()
c_hats<-c()
opti_list<-list()
re_model<-list()
Global_runs<-vector(length=runs1)




for(z in z:runs1){
  gr<-FALSE
  #### Step 1 #### Outer E-Step
  Outer1<-rdiffpois_array_paper(data_rdiff=(train_data_l%>%
                                  dplyr::select(-c("Incoming", "Outgoing"))),
                        IMSTEP=I_Mstep,
                        w_sim=w_start,
                        max.k=1000)
  
  max_w_start=100
  #### Step 2 #### Outgoing M-Step 
  for(i in 1:max_w_start){
    w_start_i<-tryCatch({optimize_Like(w_opt_prev=w_start, 
                data_like_opp=Outer1$data_all)},
      error=function(e){
        message(paste("Fell into local maximum and have to replace"))
        interim_res<-list()
        interim_res$w_sol<-omegas_org[z-1,-max_lag]
        interim_res$Infor<-as.matrix(VCOV[z-1,,])
        Global_runs[i]<<-TRUE
        return(interim_res)
      })
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
      gr<<-TRUE
      w_start<-exp(-0.2*1:max_lag)/sum(exp(-0.2*1:max_lag))
      print("We have reached a local max")
    }
  }
  if(i==max_w_start){ #### This is a failsafe in case the starting values for incoming and outgoing are poorly sampled
    if(z>20){
      w_start<-omegas_org[z-1,-max_lag]
    }else{
      w_start<-exp(-0.2*1:max_lag)/sum(exp(-0.2*1:max_lag))      
    }
    print("i=max runs QuadProg (1) We have run into a local max")
    gr<<-TRUE
  }
  
  
  # VCOV[z,,]<-w_start_i$Infor
  
  if(z>200){

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
    w_start_i<-tryCatch({
      optimize_Like(w_opt_prev=w_start_jk, 
                             data_like_opp=Inner4$data_all)},
      error=function(e){
        message(paste("Fell into local maximum and have to replace"))
        interim_res<-list()
        interim_res$w_sol<-omegas_jk[z-1,-max_lag]
        gr<<-TRUE
        return(interim_res)
        })
    if(sum((c(w_start_i$w_sol, 1-sum(w_start_i$w_sol))-w_start_jk)^2)<0.00001){
      w_start_jk<-c(w_start_i$w_sol, 1-sum(w_start_i$w_sol))
      break
    }else{
      w_start_jk<-c(w_start_i$w_sol, 1-sum(w_start_i$w_sol))
    }
  }
  if(i==max_w_start){
  w_start_jk<-w_start
  print("i=max runs QuadProg (2) We have run into a local max")
  gr<<-TRUE}
  
  print(paste("        OMEGA JK: Done with", i, "outof 100 runs"))
  
  if(length(w_start)!=max_lag){
    w_start<-c(w_start, 1-sum(w_start))
  }
  if(length(w_start_jk)!=max_lag){
  w_start_jk<-c(w_start_jk, 1-sum(w_start_jk))
  }
  #### Step 6 #### Bias correction
  c_hat<-coef(lm(I((w_start-1/max_lag)^2)~I((w_start_jk[1:length(w_start)]-1/max_lag)^2)-1))
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
  
  Global_runs[i]<-gr
  print(paste("Done with", z, "outof", runs1, "runs"))
}






################################################################################
############# Results ##########################################################
################################################################################


VarianceRuns<-250:350


################################################################################
############# Paper: Figure 6  #################################################
################################################################################


smooth_names<-colnames(predict(modellist[[1]], type="lpmatrix"))


coeffcients_smooth<-apply(matrix(unlist(lapply(modellist[VarianceRuns], coef)), 
                                 byrow=TRUE, ncol=length(coef(modellist[[200]]))), 
                  2,median)[grepl("long", smooth_names)]

BasisFn_smooth<-predict(modellist[[1]], newdata=train_data_l, type="lpmatrix")[,
                 grepl("long", smooth_names)]

SmoothEst<-unique(cbind(train_data_l[, c("date", "districtId")], 
                 "SmoothEff"=t(coeffcients_smooth%*%t(BasisFn_smooth)))%>%
  full_join(train_data_l)%>%
  dplyr::select(districtId, SmoothEff, geometry))%>%
  drop_na()%>%
  st_as_sf(crs = 4326)


# Sicherstellen, dass `shpdatei` das richtige CRS hat
shpdatei<-st_transform(SmoothEst, crs = 4326)
# cities
cities<-data.frame(
  city = c("Hamburg", "Berlin", "Dortmund", "Munich", "Dresden", "Stuttgart"),
  lon = c(10.0, 13.4, 7.47, 11.58, 13.73, 9.18),
  lat = c(53.55, 52.52, 51.51, 48.14, 51.05, 48.78)
)
# Convert to an sf object
cities_sf<-st_as_sf(cities, coords = c("lon", "lat"), crs = 4326)
# Hintergrundkarte mit Natural Earth-Daten (Alternative zu OSM)
germany<-ne_countries(scale = "medium", country="Germany", returnclass = "sf")

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
            stat = "sf_coordinates", nudge_y = 0.4, size = 4, color = "black")

############# Temporal Effects  #################################################


coeffcients_smooth_t<-apply(matrix(unlist(lapply(modellist[VarianceRuns], coef)), 
                                 byrow=TRUE, ncol=length(coef(modellist[[200]]))), 
                          2,median)[grepl("date", smooth_names)]

coeffcients_smooth_t_min<-apply(matrix(unlist(lapply(modellist[VarianceRuns], coef)), 
                                   byrow=TRUE, ncol=length(coef(modellist[[200]]))), 
                            2,quantile, 0.025)[grepl("date", smooth_names)]
coeffcients_smooth_t_max<-apply(matrix(unlist(lapply(modellist[VarianceRuns], coef)), 
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
  # paste0("Median adjustment scale ", round(median(C_hat[VarianceRuns]), 
  #                                          2)))+
  ylab(expression({{hat(f^2)}}(t)))+ 
  scale_colour_manual(values = brewer.pal(n = 5, name = "Blues")[4:5])+
  scale_fill_manual(values = brewer.pal(n = 5, name = "Blues")[2:5])+
  xlab("Date (2021)")+
  theme(legend.title = element_blank())


ggsave(plot_grid(mapsmooth, Temp_smooth, rel_heights = c(0.9, 0.1)), 
       file="3_Results/Figure_6.pdf", width=14, height=5)


################################################################################
############# Paper: Table 3  ##################################################
################################################################################

varianceRuns<-VarianceRuns

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
estimated_adj_w<-apply(omegas_adj[varianceRuns,], 2, median)/
  sum(apply(omegas_adj[varianceRuns,], 2, median))
for(i in varianceRuns){
  VCOV_Coef_Runs_b[(i)-min(varianceRuns)+1,,]<-
    ((unlist(lapply(modellist[i], coef))-coef_med)%*%t(
      unlist(lapply(modellist[i], coef))-coef_med))
  VCOV_Coef_Runs_s[(i)-min(varianceRuns)+1,,]<-
    vcov(modellist[[i]])
  VCOV_Coef_Runs_w_s[(i)-min(varianceRuns)+1,,]<-solve(VCOV[i,,])
  VCOV_Coef_Runs_w_b[(i)-min(varianceRuns)+1,,]<-
    ((omegas_adj[i,]-estimated_adj_w)%*%t((omegas_adj[i,]-estimated_adj_w)))[1:(max_lag-1), 
                                                                             1:(max_lag-1)]
}



stddev_coef<-sqrt(diag(apply(VCOV_Coef_Runs_s, c(2,3), median)+(1+(max(varianceRuns)-min(varianceRuns))^(-1))*apply(VCOV_Coef_Runs_b, c(2,3), sum)/(max(varianceRuns)-min(varianceRuns)-1)))
stddev_w<-sqrt(diag(apply(VCOV_Coef_Runs_w_s, c(2,3), median)+(1+(max(varianceRuns)-min(varianceRuns))^(-1))*apply(VCOV_Coef_Runs_w_b, c(2,3), sum)/(max(varianceRuns)-min(varianceRuns)-1)))
Table_3_prep<-as.data.frame(cbind(names(coef(modellist[[1]])),round(coef_med,3), round(stddev_coef, 3))[c(1,2,3,4,5,9,10,8,6,7),])
names(Table_3_prep)<-c("Covariates", "Estimates", "Std. Dev.")
Table_3<-Table_3_prep
Table_3$Covariates<-c("Intercept", "Infection Rate 35-59", "Infection Rate 60-79",
                      "Infection Rate 80+","Monday","Tuesday","Wednesday",
                      "Thursday", "Saturday", "Sunday")

write.csv(Table_3, file="3_Results/Table_3.csv")
# Here we usually observe some differences in rounding

################################################################################
############# Paper: Figure 7  #################################################
################################################################################


Exit_rates<-ggplot()+
  geom_hline(aes(yintercept=1/(max_lag), col="1/30"), linetype=2)+
  geom_smooth(aes(x=1:(max_lag-1), y=estimated_adj_w[-max_lag],
                  col="Smoothed estimated exit rate"), linewidth=1.1, linetype=2)+
  geom_ribbon(aes(x=1:(max_lag-1), ymin=estimated_adj_w[-max_lag]-1.96*stddev_w, 
                  ymax=estimated_adj_w[-max_lag]+1.96*stddev_w, 
                  fill="95% Confidence interval"), alpha=0.5)+
  geom_line(aes(x=1:(max_lag-1), y=estimated_adj_w[-max_lag],
                col="Estimated exit rate"), linewidth=1.1)+
  theme_pubr()+
  labs(color=" ")+
  ggtitle("Estimated exit rate", expression(paste(sqrt(hat(c)), "=1.18", sep="")))+
  # paste0("Median adjustment scale ", round(median(C_hat[VarianceRuns]), 
  #                                          2)))+
  ylab("Value of exit rate")+ 
  scale_colour_manual(values = c("darkgrey", brewer.pal(n = 5, name = "Blues")[4:5]))+
  scale_fill_manual(values = c(brewer.pal(n = 5, name = "Blues")[2:5]))+
  xlab("Length of stay")+
  scale_x_continuous(breaks = seq(1, max_lag-1, 1))+
  theme(legend.title = element_blank())


ggsave(Exit_rates, file="3_Results/Figure_7.pdf", width=7, height=4)


################################################################################
############# Paper: Figure 8  #################################################
################################################################################
## Berlin and Brandenburg highlights added in post (blue lines around title)

DIVI_cap<-read.csv("2. Data/Intensivregister_Bundeslaender_Kapazitaeten.csv")

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

ggsave(est_vs_true_plot, file="3_Results/Figure_8.pdf", height=10, width=10)


################################################################################
############# Paper: Figure 10  ################################################
################################################################################


LogLike<-ggplot()+
  geom_vline(aes(xintercept=200, col="Beginning Bias Correction"), linetype=2)+
  geom_line(aes(x=1:350, y=logLike[1:350], col="Log-Likelihood"))+
  theme_pubr()+
  labs(color=" ")+
  ggtitle("a) Log-Likelihood over all 350 runs", "ICU Covid-19 Data")+
  ylab("Log-Likelihood")+ 
  scale_colour_manual(values = c(brewer.pal(n = 5, name = "Blues")[4:5]))+
  scale_fill_manual(values = c(brewer.pal(n = 5, name = "Blues")[2:5]))+
  xlab("sEM Iteration")+
  theme(legend.title = element_blank())


LogLike150<-ggplot()+
  # geom_vline(aes(xintercept=200, col="Beginning Bias Correction"), linetype=2)+
  geom_line(aes(x=201:350, y=logLike[201:350], col="Log-Likelihood"))+
  theme_pubr()+
  labs(color=" ")+
  ggtitle("b) Log-Likelihood over the last 150 runs", "ICU Covid-19 Data")+
  ylab("Log-Likelihood")+ 
  scale_colour_manual(values = c(brewer.pal(n = 5, name = "Blues")[4:5]))+
  scale_fill_manual(values = c(brewer.pal(n = 5, name = "Blues")[2:5]))+
  xlab("sEM Iteration")+
  theme(legend.title = element_blank())


ggsave(plot_grid(LogLike, LogLike150, nrow=1), file="3_Results/Figure_10.pdf", height=3, width=9)

