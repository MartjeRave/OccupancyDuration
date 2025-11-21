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

train_data_l<-readRDS("2. Data/training_data.rds")%>%
  mutate(date=as.Date(date))%>%
  filter(date<=as.Date("2021-12-31")&date>=as.Date("2021-01-01"))
train_data_l<-train_data_l%>%
  mutate(Incoming=rpois(n=length(exp(0.3*train_data_l$G_5_7)), lambda=exp(0.1*train_data_l$G_5_7)), 
         Outgoing=rpois(n=length(exp(0.3*train_data_l$G_5_7)), lambda=exp(0.1*train_data_l$G_5_7)))


formula_in=as.formula(Incoming ~ G_4_7_lag_1+G_5_7_lag_1+G_6_7_lag_1+
                        weekday+s(long, lat)+s(as.numeric(date))) 

I_Mstep<-mgcv::gam(formula_in, data=train_data_l, 
                   family="poisson")



######################################################
############## EM Algorithm ##########################
######################################################

runs1=400
max_lag=40
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
#VCOV<-abind(VCOV,array(0, dim=c(runs1-400+1, length(w_start)-1, length(w_start)-1)), along=1)
logLike<-c()
c_hats<-c()
opti_list<-list()
re_model<-list()
Global_runs<-vector(length=runs1)
z<-1
set.seed(123)
for(z in 1:runs1){
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


InOut_est<-cbind(InOutIt[,c(1,2)], "Incoming_est"=apply(InOutIt[,names(InOutIt)%in%paste0("Incoming_", VarianceRuns)],
                                                        1, median))%>%
  full_join(cbind(InOutIt[,c(1,2)], "Outgoing_est"=apply(InOutIt[,names(InOutIt)%in%paste0("Outgoing_", VarianceRuns)],
                                                         1, median)))%>%
  mutate(date=as.Date(date), districtId=as.numeric(districtId))%>%
  full_join(train_data_l[,c("date", "districtId", "Diff")]%>%
              mutate(date=as.Date(date), districtId=as.numeric(districtId)))%>%
  filter(date%in%unique(as.Date(InOutIt$date)))


TotalCases<-readRDS("2. Data/DIVI_2022-10-27.rds")


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


coeffcients_smooth<-apply(matrix(unlist(lapply(modellist[VarianceRuns], coef)), 
                                 byrow=TRUE, ncol=length(coef(modellist[[200]]))), 
                  2,median)[grepl("long", smooth_names)]
BasisFn_smooth<-predict(modellist[[1]], newdata=train_data_l, type="lpmatrix")[,grepl("long", smooth_names)]
SmoothEst<-unique(cbind(train_data_l[, c("date", "districtId")], 
                 "SmoothEff"=t(coeffcients_smooth%*%t(BasisFn_smooth)))%>%
  full_join(train_data_l)%>%
  dplyr::select(districtId, SmoothEff, geometry))%>%
  drop_na()%>%
  st_as_sf(crs = 4326)


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
  #scale_x_continuous(breaks = seq(min(SmoothEst_t$date), max(SmoothEst_t$date), 30))+
  theme(legend.title = element_blank())



ggsave(plot_grid(mapsmooth, Temp_smooth, rel_heights = c(0.9, 0.1)), file="5_Supplementary Material/Supp_Figure_2.pdf", width=14, height=5)




############# Fixed Coefficients ###############################################

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
cbind(names(coef(modellist[[1]])),round(coef_med,3), round(stddev_coef, 3))



Exit_rates<-ggplot()+
  geom_hline(aes(yintercept=1/(max_lag), col="1/40"), linetype=2)+
  geom_ribbon(aes(x=1:(max_lag-1), ymin=estimated_adj_w[-max_lag]-1.96*stddev_w, 
                  ymax=estimated_adj_w[-max_lag]+1.96*stddev_w, 
                  fill="95% Confidence interval"))+
  geom_line(aes(x=1:(max_lag-1), y=estimated_adj_w[-max_lag],
                col="Estimated exit rate"))+
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



ggsave(Exit_rates, file="5_Supplementary Material/Supp_Figure_3.pdf", width=9, height=4)





################################################################################

DIVI_cap<-read.csv("2_Data/Intensivregister_Bundeslaender_Kapazitaeten.csv")

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

ggsave(est_vs_true_plot, file="5_Supplementary Material/Supp_Figure_1.pdf", height=10, width=10)



write.csv(Coeff_in, file="5_Supplementary Material/Coef_In_Lag_40.csv")
save(modellist, file="5_Supplementary Material/modellist_Lag_40.R")
save(fitlist, file="5_Supplementary Material/fitlist_Lag_40.R")
save(VCOV, file="5_Supplementary Material/VCOV_Coef_Runs_w_s_Lag_40.R")
write.csv(omegas_adj, file="5_Supplementary Material/Exitrate_In_Lag_40.csv")
