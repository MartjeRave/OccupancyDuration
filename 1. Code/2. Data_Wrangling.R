
####################
### Set working directory to current folder
####################

function_folder <- "Skellam_share"
functionpath<-substr(dirname(rstudioapi::getSourceEditorContext()$path), 
                     1, unlist(gregexpr(function_folder, 
                                        dirname(rstudioapi::getSourceEditorContext()$path)))+(nchar(function_folder)-1))

setwd(functionpath)

####################
### Response Divi
####################


source("0. Functions/Functions.R")


#### Read in ####

Divi<-readRDS("2. Data/1. Data Download/RKI DIVI in analysis/DIVI_2022-10-27.rds")

Divi$daten_stand<-as.Date(Divi$date)
Divi$betten_belegt<-Divi$betten_belegt-Divi$faelle_covid_aktuell
Divi<-Divi[,!names(Divi)%in%"date"]
Divi<-unique(Divi)


#### NA imputation ####

if(any(sqldf("SELECT gemeindeschluessel, daten_stand, Count(*) AS COUNT_ FROM Divi GROUP BY 1,2 ORDER BY 3 DESC")$COUNT_>1)){
  print("multiple observations in Divi file")}


#####################
#### Covariates RKI
#####################

#### Read in ####


RKI<-read_rds("2. Data/1. Data Download/RKI DIVI in analysis/RKI_2022-10-27.rds")

RKI$age<-as.factor(RKI$age)
RKI$gender<-as.factor(RKI$gender)
lookupage<-cbind(levels(RKI$age),c("G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_N"))

levels(RKI$age)<-c("G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_N")
levels(RKI$gender)<-c("m", "n", "f")
# table(RKI$age)
RKI$cases<-abs(RKI$cases)
new_rows<-sqldf("SELECT landId, districtId, date, age, gender,
                sum(cases) AS cases FROM RKI GROUP BY landId, districtId, date, age, gender")

new_vars<-levels(as.factor(paste0(RKI$gender, "_", RKI$age)))

for(i in 1:length(new_vars)){
  new_rows$new<-ifelse(grepl(new_vars[i], paste0(new_rows$gender, "_", new_rows$age)), 1, 0)*new_rows$cases
  names(new_rows)<-c(names(new_rows)[-length(names(new_rows))],new_vars[i])
}

### Divi only have one value for Berlin ###
new_rows$districtId[new_rows$landId==11]<-11000


###########################
#############DIVI##########
###########################


#### Cleaning of Data  ####

complete_data<-as.data.frame(rep(as.Date(min(Divi$daten_stand):max(Divi$daten_stand),  origin = "1970-01-01"),
                                 times=length(unique(new_rows$districtId))))
names(complete_data)<-"daten_stand"
complete_data$districtId<-sort(rep(unique(new_rows$districtId),
                                   times=length(unique(complete_data$daten_stand))))
names(complete_data)<-c("daten_stand","gemeindeschluessel")

Divi<-as_tibble(Divi)%>%
  full_join(complete_data, by=c("daten_stand","gemeindeschluessel"))%>%
  arrange(gemeindeschluessel,daten_stand)


###setoffnumbers nills so we can go on and imputing
for(i in 1:ncol(Divi)){
  Divi[Divi$daten_stand==min(Divi$daten_stand),i][is.na(Divi[Divi$daten_stand==min(Divi$daten_stand),i])]<-0
}
Divi$bundesland<-floor(Divi$gemeindeschluessel/1000)



for(k in unique(Divi$gemeindeschluessel)){
  if(any(is.na(Divi[Divi$gemeindeschluessel==k,]))){
    for(i in 1:ncol(Divi[Divi$gemeindeschluessel==k,])){
      for(j in 2:nrow(Divi[Divi$gemeindeschluessel==k,])){
        if(is.na(Divi[Divi$gemeindeschluessel==k,i][j,])&((k==Divi[Divi$gemeindeschluessel==k,"gemeindeschluessel"][j-1,]))){
          Divi[Divi$gemeindeschluessel==k,i][j,]<-Divi[Divi$gemeindeschluessel==k,i][j-1,]
        }
      }
    }
  }
}




new_rows_old<-new_rows

new_rows<-sqldf("SELECT landId, districtId, date, SUM(cases) as cases_tot,
                  SUM(f_G_1+m_G_1+n_G_1) AS G_1,
                  SUM(f_G_2+m_G_2+n_G_2) AS G_2,
                  SUM(f_G_3+m_G_3+n_G_3) AS G_3,
                  SUM(f_G_4+m_G_4+n_G_4) AS G_4,
                  SUM(f_G_5+m_G_5+n_G_5) AS G_5,
                  SUM(f_G_6+m_G_6+n_G_6) AS G_6
                FROM new_rows GROUP BY landId, districtId, date ORDER BY districtId, date")

# new_rows_bund<-sqldf("SELECT landId, date, SUM(cases) as cases_tot,
#                   SUM(f_G_1) AS f_G_1, SUM(f_G_2) AS f_G_2, SUM(f_G_3) AS f_G_3,
#                   SUM(f_G_4) AS f_G_4, SUM(f_G_5) AS f_G_5,
#                   SUM(f_G_6) AS f_G_6, SUM(f_G_N) AS f_G_N,
#                   SUM(m_G_1) AS m_G_1, SUM(m_G_2) AS m_G_2, SUM(m_G_3) AS m_G_3,
#                   SUM(m_G_4) AS m_G_4, SUM(m_G_5) AS m_G_5, SUM(m_G_6) AS m_G_6, SUM(m_G_N) AS m_G_N,
#                   SUM(n_G_1) AS n_G_1, SUM(n_G_2) AS n_G_2, SUM(n_G_3) AS n_G_3,
#                   SUM(n_G_4) AS n_G_4, SUM(n_G_5) AS n_G_5, SUM(n_G_6) AS n_G_6,
#                   SUM(n_G_N) AS n_G_N FROM new_rows_old GROUP BY landId,  date ORDER BY landId, date")
# sum(new_rows[,5:ncol(new_rows)])==sum(new_rows$cases_tot)


new_rows$date<- as.Date(new_rows$date)
RKI_old<-RKI
RKI<-as_tibble(new_rows)
# RKI_bund<-as_tibble(new_rows_bund)



#####################
#### Join for Total
#####################

RKI_Divi<-merge(RKI, Divi, by.y=c("daten_stand","gemeindeschluessel"),
                by.x=c("date",  "districtId"), all=T)%>%
  filter(date>=min(Divi$daten_stand)&date>=min(RKI$date))%>%
  arrange(districtId, date)



# some days are missing for some districts so we join is onto a complete date table for each district



RKI_Divi<-RKI_Divi%>%
  full_join(tibble(date=rep(min(RKI_Divi$date)+seq(0:(max(RKI_Divi$date)-min(RKI_Divi$date)))-1,
                            length(levels(as.factor(RKI_Divi$districtId)))),
                   districtId=sort(as.numeric(rep(levels(as.factor(RKI_Divi$districtId)),
                                                  (max(RKI_Divi$date)-min(RKI_Divi$date)+1))))),
            by =c("date","districtId"))%>%
  filter(date>(min(Divi$daten_stand)-1)& date<(max(Divi$daten_stand)+1))%>%
  dplyr::select(-c("landId"))

#####################
#### Imputation
#####################

# Missing RKI data is missing observations whereas in Divi data it just seems like that there is some data not recorded where there
# should be data, just because we think that the hospital data is a lot less variable that the RKI data


RKI_cols<-names(new_rows)[4:length(names(new_rows))]
RKI_Divi[is.na(RKI_Divi$cases_tot), RKI_cols]<-0

for(i in 1:ncol(RKI_Divi)){
  if(any(is.na(RKI_Divi[,i]))){print(paste(i, "HAS NA's"))}
}


# sqldf("SELECT COUNT(*) AS COUNT_, date,COUNT(DISTINCT districtId) AS COUNT_GEM  FROM
#      RKI_Divi GROUP BY 2 HAVING COUNT_>401 OR  COUNT_GEM<>COUNT_")



#####################
#### Lagged Variables
#####################


lag_set_up<-RKI_Divi[,c("districtId", "date", "cases_tot", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6")]


lag_set_up_1<-lag_set_up
lag_set_up_1$date<-lag_set_up_1$date+1
names(lag_set_up_1)[grepl("G_", names(lag_set_up_1))]<-paste0(names(lag_set_up_1)[grepl("G_", names(lag_set_up_1))], "_d_1")
lag_set_up_1<-lag_set_up_1[,!grepl("cases_tot", names(lag_set_up_1))]

lag_set_up_2<-lag_set_up
lag_set_up_2$date<-lag_set_up_2$date+2
names(lag_set_up_2)[grepl("G_", names(lag_set_up_2))]<-paste0(names(lag_set_up_2)[grepl("G_", names(lag_set_up_2))], "_d_2")
lag_set_up_2<-lag_set_up_2[,!grepl("cases_tot", names(lag_set_up_2))]

lag_set_up_3<-lag_set_up
lag_set_up_3$date<-lag_set_up_3$date+3
names(lag_set_up_3)[grepl("G_", names(lag_set_up_3))]<-paste0(names(lag_set_up_3)[grepl("G_", names(lag_set_up_3))], "_d_3")
lag_set_up_3<-lag_set_up_3[,!grepl("cases_tot", names(lag_set_up_3))]

lag_set_up_4<-lag_set_up
lag_set_up_4$date<-lag_set_up_4$date+4
names(lag_set_up_4)[grepl("G_", names(lag_set_up_4))]<-paste0(names(lag_set_up_4)[grepl("G_", names(lag_set_up_4))], "_d_4")
lag_set_up_4<-lag_set_up_4[,!grepl("cases_tot", names(lag_set_up_4))]

lag_set_up_5<-lag_set_up
lag_set_up_5$date<-lag_set_up_5$date+5
names(lag_set_up_5)[grepl("G_", names(lag_set_up_5))]<-paste0(names(lag_set_up_5)[grepl("G_", names(lag_set_up_5))], "_d_5")
lag_set_up_5<-lag_set_up_5[,!grepl("cases_tot", names(lag_set_up_5))]

lag_set_up_6<-lag_set_up
lag_set_up_6$date<-lag_set_up_6$date+6
names(lag_set_up_6)[grepl("G_", names(lag_set_up_6))]<-paste0(names(lag_set_up_6)[grepl("G_", names(lag_set_up_6))], "_d_6")
lag_set_up_6<-lag_set_up_6[,!grepl("cases_tot", names(lag_set_up_6))]

lag_set_up_7<-lag_set_up
lag_set_up_7$date<-lag_set_up_7$date+7
names(lag_set_up_7)[grepl("G_", names(lag_set_up_7))]<-paste0(names(lag_set_up_7)[grepl("G_", names(lag_set_up_7))], "_d_7")
lag_set_up_7<-lag_set_up_7[,!grepl("cases_tot", names(lag_set_up_7))]

# library("tidyverse")
# install.packages("cli", dependencies = TRUE)
# detach("package:cli", unload=TRUE)

# cbind(as.data.frame(as.Date(rep(min(lag_set_up$date):max(lag_set_up$date), each=length(unique(as.factor(lag_set_up$districtId)))))),
#       as.data.frame(unique(lag_set_up$districtId)))

all_data<-as.tbl(cbind(as.data.frame(rep(unique(lag_set_up$districtId), each=length(as.Date(min(lag_set_up$date)+seq(0, min(lag_set_up$date):max(min(lag_set_up$date))))))),
                       date=rep(as.Date(min(lag_set_up$date)+seq(0, min(lag_set_up$date):max(min(lag_set_up$date)))), times=length(unique(lag_set_up$districtId)))))

names(all_data)<-c("districtId", "date")

unique(all_data$date)

lag_set_up_com<-lag_set_up%>%
  full_join(all_data)%>%
  filter(date>as.Date("2020-05-01")&date<Sys.Date())%>%
  arrange(districtId, date)%>%
  full_join(lag_set_up_1)%>%
  filter(date>as.Date("2020-05-01")&date<Sys.Date())%>%
  full_join(lag_set_up_2)%>%
  filter(date>as.Date("2020-05-01")&date<Sys.Date())%>%
  full_join(lag_set_up_3)%>%
  filter(date>as.Date("2020-05-01")&date<Sys.Date())%>%
  full_join(lag_set_up_4)%>%
  filter(date>as.Date("2020-05-01")&date<Sys.Date())%>%
  full_join(lag_set_up_5)%>%
  filter(date>as.Date("2020-05-01")&date<Sys.Date())%>%
  full_join(lag_set_up_6)%>%
  filter(date>as.Date("2020-05-01")&date<Sys.Date())%>%
  full_join(lag_set_up_7)%>%
  filter(date>as.Date("2020-05-01")&date<Sys.Date())%>%
  filter(date>as.Date("2020-05-01"))


#####################
#### Imputation
#####################

new_rows_com<-lag_set_up_com

for(i in which(grepl("G_", names(new_rows_com))|grepl("cases", names(new_rows_com)))){
  new_rows_com[is.na(new_rows_com[,i]),i]<-0
}


#####################
#### Summed and Inzidences
#####################


districts<-preprocess.districts()

districts$pop1<-districts$pop.m.0.4+districts$pop.w.0.4
districts$pop2<-districts$pop.m.5.14+districts$pop.w.5.14
districts$pop3<-districts$pop.m.15.34+districts$pop.w.15.34
districts$pop4<-districts$pop.m.35.59+districts$pop.w.35.59
districts$pop5<-districts$pop.m.60.79+districts$pop.w.60.79
districts$pop6<-districts$pop.m.80+districts$pop.w.80


### We add up all incidences from 6 days prior to today
new_rows_sum<-sqldf("SELECT districtId, date,
      G_1+G_1_d_1+G_1_d_2+G_1_d_3+G_1_d_4+G_1_d_5+G_1_d_6 AS G_1_7,
      G_2+G_2_d_1+G_2_d_2+G_2_d_3+G_2_d_4+G_2_d_5+G_2_d_6 AS G_2_7,
      G_3+G_3_d_1+G_3_d_2+G_3_d_3+G_3_d_4+G_3_d_5+G_3_d_6 AS G_3_7,
      G_4+G_4_d_1+G_4_d_2+G_4_d_3+G_4_d_4+G_4_d_5+G_4_d_6 AS G_4_7,
      G_5+G_5_d_1+G_5_d_2+G_5_d_3+G_5_d_4+G_5_d_5+G_5_d_6 AS G_5_7,
      G_6+G_6_d_1+G_6_d_2+G_6_d_3+G_6_d_4+G_6_d_5+G_6_d_6 AS G_6_7
      FROM new_rows_com")%>%
  full_join(unique(districts[, 
                                   c("districtId", "pop1", "pop2", "pop3", "pop4", "pop5",  "pop6")]))%>%
  #We take weekly averages and the log of it
  mutate(G_1_7=log(G_1_7/(pop1*7)*100000+0.1),
         G_2_7=log(G_2_7/(pop2*7)*100000+0.1),
         G_3_7=log(G_3_7/(pop3*7)*100000+0.1),
         G_4_7=log(G_4_7/(pop4*7)*100000+0.1),
         G_5_7=log(G_5_7/(pop5*7)*100000+0.1),
         G_6_7=log(G_6_7/(pop6*7)*100000+0.1))%>%
  #the variables we end up with
  dplyr::select(districtId, date, G_1_7, G_2_7, G_3_7, G_4_7, G_5_7,G_6_7)

# 7 day lag to be added
new_rows_lag<-new_rows_sum
new_rows_lag$date<-new_rows_lag$date+1
names(new_rows_lag)[3:length(names(new_rows_lag))]<-paste0(names(new_rows_lag)[3:length(names(new_rows_lag))], "_lag_1")

Inz<-new_rows_sum%>%
  dplyr::select(districtId, date, G_4_7,G_5_7,G_6_7)%>%
  full_join(new_rows_lag)%>%
  dplyr::select(districtId, date, G_4_7,G_5_7,G_6_7, G_4_7_lag_1, G_5_7_lag_1,G_6_7_lag_1)%>%
  mutate(date=as.Date(date))%>%
  dplyr::filter(date > as.Date("2020-12-31") & date < as.Date("2022-01-01"))


#####################
#### Adding back DIVI
#####################

Divi_diff<-RKI_Divi%>%
  dplyr::select(districtId, date, faelle_covid_aktuell)%>%
  mutate(date=date+1, faelle_covid_aktuell_1=faelle_covid_aktuell)%>%
  dplyr::select(districtId, date, faelle_covid_aktuell_1)


RKI_Divi_lag<-RKI_Divi%>%
  dplyr::select(districtId, date, betten_frei, betten_belegt, faelle_covid_aktuell)%>%
  full_join(Inz)%>%
  dplyr::filter(date > as.Date("2020-12-31") & date < as.Date("2022-01-01"))%>%
  full_join(Divi_diff)%>%
  mutate(Diff=faelle_covid_aktuell-faelle_covid_aktuell_1,
         bed_tot=faelle_covid_aktuell+betten_frei+betten_belegt)%>%
  dplyr::filter(date > as.Date("2020-12-31") & date < as.Date("2022-01-01"))%>%
  dplyr::select(districtId, date, Diff, bed_tot,
                G_4_7, G_5_7, G_6_7, G_4_7_lag_1, G_5_7_lag_1, G_6_7_lag_1)%>%
  arrange(date,districtId)


RKI_Divi_lag<-merge(RKI_Divi_lag, unique(RKI_old[,c("district", "districtId")]), by.x="districtId",
                    by.y="districtId", all.x=TRUE)%>%
  arrange(districtId, date)

############################################################################################
##### adding the structure of the districts to the data set. Data was received by Dr Cornelius Fritz, 2021.

load("C:/Users/ra98jiq/Documents/Papers/Skellam_share/2. Data/0. Data Extra/district_data.RData")



RKI_Divi_lag<-merge(RKI_Divi_lag, district_data[,c("lk_id", "long", "lat")], by.x = c("districtId"), by.y=c("lk_id"),
                    all=TRUE)

RKI_Divi_lag$district[RKI_Divi_lag$districtId=="11000"]<-"SK Berlin"

# save.image(file="C:/Users/ra98jiq/Documents/Papers/Skellam/Skellam/dec_01_12.Rdata")


depri<-read_excel("C:/Users/ra98jiq/LRZ Sync+Share/Corona/Data/Demographics/deprivation.xlsx")
#### ABGLEICHUNG VON LANDKREISEN!!!!! ########
depri$AGS[depri$AGS=="13002000"|
            depri$AGS=="13052000"|
            depri$AGS=="13055000"]<-"13071000" #"LK Mecklenburgische Seenplatte"
depri$AGS[depri$AGS=="13053000"|
            depri$AGS=="13051000"|
            depri$AGS=="13056000"]<-"13072000"  #"LK Rostock" 
depri$AGS[depri$AGS=="13061000"|
            depri$AGS=="13005000"|
            depri$AGS=="13057000"] <-"13073000" #"LK Vorpommern-R\xfcgen"
depri$AGS[depri$AGS=="13058000"|
            depri$AGS=="13006000"]<-"13074000" #"LK Nordwestmecklenburg"
depri$AGS[depri$AGS=="13062000"|
            depri$AGS=="13059000"|
            depri$AGS=="13001000"]<-"13075000" #"LK Vorpommern-Greifswald"
depri$AGS[depri$AGS=="13060000"|
            depri$AGS=="13054000"]<-"13076000" #"LK Ludwigslust-Parchim"


RKI_Divi_lag<-depri%>%
  mutate(AGS=as.numeric(AGS)/1000)%>%
  mutate(AGS=ifelse(floor(AGS/1000)==11, 11000, AGS))%>%
  dplyr::group_by(AGS)%>%
  dplyr::summarize(GIMD10=mean(GIMD10))%>%
  dplyr::mutate(districtId= as.character(AGS))%>%
  dplyr::select(districtId, GIMD10)%>%
  dplyr::full_join(RKI_Divi_lag)%>%
  dplyr::mutate(district=ifelse(districtId=="11000","Berlin",district))%>%
  na.omit()

RKI_Divi_lag<-RKI_Divi_lag%>%
  na.omit()


RKI_Divi_lag$district<-gsub( "\xdf", "ß", RKI_Divi_lag$district)
RKI_Divi_lag$district<-gsub("\xf6" ,"ö", RKI_Divi_lag$district)
RKI_Divi_lag$district<-gsub( "\xe4", "ä", RKI_Divi_lag$district)
RKI_Divi_lag$district<-gsub( "\xfc", "ü", RKI_Divi_lag$district)


RKI_Divi_lag<- RKI_Divi_lag[RKI_Divi_lag$districtId!=16056,]

RKI_Divi_lag$weekday<-as.character(wday(RKI_Divi_lag$date, label=TRUE))

saveRDS(RKI_Divi_lag, "C:/Users/ra98jiq/Documents/Papers/Skellam_share/2. Data/2. Training Data/training_data.rds")

# 
# write.csv(RKI_Divi_lag, file="C:/Users/ra98jiq/Documents/Papers/Skellam/Skellam/Data/training_data.csv", fileEncoding = "UTF-8")
