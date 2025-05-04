#####################################################
######### Comparing Bias Correction ##################
######################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
#library(taskscheduleR)
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
library(ggplot2)
library(sf)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)


train_data_l<-readRDS("2. Data/training_data.rds")%>%
  mutate(date=as.Date(date))%>%
  filter(date<=as.Date("2021-12-31")&date>=as.Date("2021-08-01"))%>%
  mutate(Incoming=rpois(n=length(exp(0.3*G_5_7)), lambda=exp(0.1*train_data_l$G_5_7)), 
         Outgoing=rpois(n=length(exp(0.3*G_5_7)), lambda=exp(0.1*train_data_l$G_5_7)))



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


ggsave(mapsmooth, file="3. Model/SmoothEff.pdf",
       height = 6, width = 5)



TotalCases_min<-as.data.frame(readRDS("2. Data/DIVI_2022-10-27.rds")%>%
                                filter(as.Date(date)%in%unique(train_data_l$date))%>%
                                na.omit()%>%
                                dplyr::group_by(gemeindeschluessel)%>%
                                dplyr::summarize(min_cap=max((betten_frei+betten_belegt), na.rm = TRUE))%>%
                                mutate(districtId=as.character(gemeindeschluessel))%>%
                                full_join(SmoothEst))%>%
  na.replace(0)%>%
  st_as_sf(crs = 4326)%>%
  arrange(min_cap)



TotalCases<-as.data.frame(readRDS("2. Data/DIVI_2022-10-27.rds")%>%
                            filter(as.Date(date)%in%unique(train_data_l$date))%>%
                            na.omit()%>%
                            dplyr::group_by(gemeindeschluessel)%>%
                            dplyr::summarize(max_cov=max(faelle_covid_aktuell/(betten_frei+betten_belegt+0.00001)*100, na.rm = TRUE))%>%
                            mutate(districtId=as.character(gemeindeschluessel))%>%
                            full_join(SmoothEst))%>%
  na.replace(0)%>%
  st_as_sf(crs = 4326)

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

og_Dist<-readRDS("2. Data/DIVI_2022-10-27.rds")

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


ggsave(plot_grid(MapIntro, InfectRatePlot), file="2. Data/IntroData.pdf",
       height = 5, width = 10)


