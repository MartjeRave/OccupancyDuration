####################
### Set working directory to current folder
####################

function_folder <- "Skellam_share"
functionpath<-substr(dirname(rstudioapi::getSourceEditorContext()$path), 
                     1, unlist(gregexpr(function_folder, 
                                        dirname(rstudioapi::getSourceEditorContext()$path)))+(nchar(function_folder)-1))


source(paste0(functionpath, "/0. Functions/Functions.R"))

################################################################################
############## Downloaded on the 27.10.2022 ####################################
################################################################################


path.data <- paste0(functionpath, "/2. Data/1. Data Download/RKI DIVI in analysis")  

#### DIVI ######################################################################


setwd(path.data)

Divi<-read_csv("https://diviexchange.blob.core.windows.net/%24web/zeitreihe-tagesdaten.csv")
saveRDS(Divi, file=paste0("DIVI_", Sys.Date(), ".rds"))

# Divi<-readRDS(paste0("DIVI_", Sys.Date(), ".rds"))

#### RKI #######################################################################


data.RKI <- read_csv("https://www.arcgis.com/sharing/rest/content/items/f10774f1c63e40168479a1feb6c7ca74/data")

data.RKI<-data.RKI[,-1]

data.RKI$districtId<-as.numeric(data.RKI$IdLandkreis)
data.RKI <- data.RKI[which(data.RKI$districtId != 3152), ] ### 3152 is not a district anymore

saveRDS(data.RKI, paste0("RKI_", Sys.Date(), ".rds"))


################################################################################
################## Old code ####################################################
################################################################################
# !!!!!!!!!!!!!!!!!!!!!!!!! This code does not work anymore, but this is how the 
# original data was downloaded. We include it here for your information. 
# Please see the updated download code below.!!!!!!!!!!!!!!!!!!!!!!!


# # set the path where the pdf reports should be stored
# path.data <- "Your path"
# 
# setwd(path.data)
# 
# days <- as.character(seq.Date(from = as.Date("2022-09-08"), to = as.Date(Sys.Date()), by = "day"))
# 
# for(i in days){
#   download_url <- paste0("https://www.divi.de/joomlatools-files/docman-files/divi-intensivregister-tagesreports-csv/DIVI-Intensivregister_",
#                          i, "_09-15.csv") #versuche Uhrzeit 9:15
#   dat <- tryCatch({
#     suppressWarnings(read.csv(download_url))
#   },
#   error = function(e){ #versuche andere Uhrzeit: 12:15 Uhr
#     download_url2 <- paste0("https://www.divi.de/joomlatools-files/docman-files/divi-intensivregister-tagesreports-csv/DIVI-Intensivregister_",
#                             i, "_12-15.csv")
#     dat <- tryCatch({
#       suppressWarnings(read.csv(download_url2))
#     },
#     error = function(e){ #versuche andere Uhrzeit: 14:15 Uhr
#       download_url3 <- paste0("https://www.divi.de/joomlatools-files/docman-files/divi-intensivregister-tagesreports-csv/DIVI-Intensivregister_",
#                               i, "_14-15.csv")
#       dat <- tryCatch({
#         suppressWarnings(read.csv(download_url3))
#       },
#       error = function(e){
#         print(paste0("keine Daten gefunden, Tag: ", i))
#         return(NULL)
#       })
#     })
#   })
#   
#   if("faelle_covid_aktuell" %in% colnames(dat)){
#     dat$daten_stand <- as.Date(i)
#     saveRDS(dat, file = paste0("data_divi_", i, ".rds"))
#   } 
# }
# 
# 
# # This function reads the raw RKI data, gives suitable column names
# # and saves the table as tibble
# # Caution: The style of the RKI datasets has changed over time and probably
# # will change again, i.e. new columns are added or columns are rearranged
# 
# # Input:
# # - all: if TRUE, all raw data available are being read and saved,
# # if FALSE (default), only datasets which have not been read yet
# # will be read. Setting "all" to TRUE should be avoided.
# 
# # Output: none
# 
# read.RKI <- function(all = FALSE){
#   
#   # get all filenames and the respective dates from the original RKI datasets
#   files.RKI <- list.files(path= paste0(path.LRZ, "Data/RKI"))
#   dates.RKI <- as.POSIXct(sub("\\..*" , "", sub("^[^2]*", "" , files.RKI)), 
#                           hour = 0, tz = "GMT")
#   
#   # get all filenames and the respective dates from the data that already has been preprocessed
#   files.read <- list.files(path= paste0(path.LRZ, "Data/Formatted"))
#   dates.read <- as.POSIXct(sub("\\..*" , "", sub("^[^2]*", "" , files.read)), 
#                            hour = 0, tz = "GMT")
#   
#   if (!all){
#     # only preprocess if data has not been preprocessed yet
#     files.RKI <- files.RKI[which(!is.element(dates.RKI, dates.read))]
#   }
#   
#   districts <- preprocess.districts()
#   
#   for (file in files.RKI) {
#     # read RKI dataset and get reporting data
#     data.RKI <- read_csv(file = paste0(path.LRZ, "Data/RKI/", file))
#     reporting.date <- as.POSIXct(sub("\\..*" , "", sub("^[^2]*", "" , file)), 
#                                  hour = 0, tz = "GMT")
#     
#     # number of columns in the RKI dataset
#     cols <- as.character(ncol(data.RKI))
#     
#     # label columns names depending on cols
#     names(data.RKI) <- switch(cols,
#                               "18" = c("id", "landId","land","district","age","gender",
#                                        "cases","deaths","date","districtId",
#                                        "updated","newcase","newdeath","date.desease",
#                                        "cured","newcure", "is.beginning.desease", "age2"),
#                               "17" =  c("landId","land","district","age","gender",
#                                         "cases","deaths","id","date","districtId",
#                                         "updated","newcase","newdeath","date.desease",
#                                         "cured","newcure", "is.beginning.desease"),
#                               "14" = c("landId","land","district","age","gender",
#                                        "cases","deaths","id","date","districtId",
#                                        "updated","newcase","newdeath","date.desease"))
#     
#     data <- data.RKI %>% mutate(gender = as.character(gender))
#     
#     # change data format of districtId from character to numeric
#     if (is.character(data$districtId)){
#       data$districtId <- as.numeric(data$districtId)
#     }
#     
#     # 3152 is not a districts any more
#     if (length(which(data$districtId == 3152)) > 0) {
#       data <- data[-which(data$districtId == 3152), ]
#     }
#     
#     # delete entries with a non valid district Id
#     ind.delete <- which(!data$districtId %in% districts$districtId)
#     if (length(ind.delete > 0)){
#       data <- data[-ind.delete, ]
#     }
#     
#     # if date is not in desired format
#     if (!is.element("POSIXct", class(data$date))){
#       data <- data %>% mutate(date = as.POSIXct(x = date, tryFormats = c("%Y-%m-%d", "%m/%d/%Y", "%Y/%m/%d"),
#                                                 tz = "UTC")) %>% 
#         mutate(date = as_datetime(x = date))
#     }
#     
#     # rename age variable of necessary
#     if (is.element("age_group", names(data))){
#       data <- rename(data, age = age_group) %>% mutate(age = as.character(age))
#     }
#     
#     # save prepocessed dataset
#     saveRDS(data, file = paste0(path.LRZ, "Data/Formatted/cases_GermanTemporal_", 
#                                 as.character(reporting.date), ".rds"))
#   }
# }
# 
# # This function formats the RKI data read in by read.RKI() to have the data in a
# # common style
# 
# # Input: none
# 
# # Output: none
# 
# format.RKI <- function(){
#   
#   files.preprocessed <- list.files(path= paste0(path.LRZ, "Data/Formatted"))
#   dates.preprocessed <- as.POSIXct(sub("\\..*" , "", sub("^[^2]*", "" , files.preprocessed)), 
#                                    hour = 0, tz = "GMT")
#   
#   for (file in files.preprocessed) {
#     data <- as_tibble(read_rds(paste0(path.LRZ, "Data/Formatted/", file))) %>% 
#       mutate(gender = as.character(gender))
#     
#     # change data format of districtId from character to numeric
#     if (is.character(data$districtId)){
#       data$districtId <- as.numeric(data$districtId)
#     }
#     
#     # 3152 is not a districts any more
#     if (length(which(data$districtId == 3152)) > 0) {
#       data.RKI <- data.RKI[-which(data.RKI$districtId == 3152), ]
#     }
#     
#     # if date is not in desired format
#     if (!is.element("POSIXct", class(data$date))){
#       data <- data %>% mutate(date = as.POSIXct(x = date, tryFormats = c("%Y-%m-%d", "%m/%d/%Y", "%Y/%m/%d"),
#                                                 tz = "UTC")) %>% 
#         mutate(date = as_datetime(x = date))
#     }
#     
#     # rename age variable of necessary
#     if (is.element("age_group", names(data))){
#       data <- rename(data, age = age_group) %>% mutate(age = as.character(age))
#     }
#     
#     saveRDS(data, file = paste0(path.LRZ, "Data/Formatted/", file))
#   }
# }
# 
# # This function creates a data frame that with one row for each district. 
# # The columns contain 
# # - the district names and Ids, 
# # - the gender/age group specific population sizes
# # - the coordinates of the centroids of the districts
# # - the population density of the districts (not used in the analyses)
# 
# # Input: none
# 
# # Output: the data frame described above
# 
# preprocess.districts <- function(){
#   
#   # read population and coordinates of districts
#   coordinates <- read_excel(paste(path.LRZ, "Data/Demographics/coordinates.xlsx", sep = ""))
#   population <- read_excel(paste(path.LRZ, "Data/Demographics/population.xlsx", sep = ""))
#   pop.density <- read_excel(paste(path.LRZ, "Data/Demographics/population_total.xlsx", sep = ""))
#   #deprivation <- read_excel(paste(path.LRZ, "Data/Demographics/deprivation.xlsx", sep = ""))
#   
#   
#   districts <- tibble(districtId = as.numeric(population$districtId[seq(1, nrow(population), 2)]),
#                       pop = round(population$gesamt[seq(1, nrow(population), 2)]), 
#                       pop.m = round(population$gesamt[seq(2, nrow(population), 2)]), 
#                       pop.f = pop - pop.m,
#                       pop.m.0.4 = round(rowSums(population[seq(2, nrow(population), 2), 5:9])),
#                       pop.w.0.4 = round(rowSums(population[seq(1, nrow(population), 2), 5:9])) -
#                         round(rowSums(population[seq(2, nrow(population), 2), 5:9])),
#                       pop.m.5.14 = round(rowSums(population[seq(2, nrow(population), 2), 10:19])),
#                       pop.w.5.14 = round(rowSums(population[seq(1, nrow(population), 2), 10:19])) - 
#                         round(rowSums(population[seq(2, nrow(population), 2), 10:19])),
#                       pop.m.15.34 = round(rowSums(population[seq(2, nrow(population), 2), 20:39])),
#                       pop.w.15.34 = round(rowSums(population[seq(1, nrow(population), 2), 20:39])) - 
#                         round(rowSums(population[seq(2, nrow(population), 2), 20:39])),
#                       pop.m.35.59 = round(rowSums(population[seq(2, nrow(population), 2), 40:64])),
#                       pop.w.35.59 = round(rowSums(population[seq(1, nrow(population), 2), 40:64])) - 
#                         round(rowSums(population[seq(2, nrow(population), 2), 40:64])),
#                       pop.m.60.79 = round(rowSums(population[seq(2, nrow(population), 2), 65:84])),
#                       pop.w.60.79 = round(rowSums(population[seq(1, nrow(population), 2), 65:84])) - 
#                         round(rowSums(population[seq(2, nrow(population), 2), 65:84])),
#                       pop.m.80 = round(rowSums(population[seq(2, nrow(population), 2), 85:95])),
#                       pop.w.80 = round(rowSums(population[seq(1, nrow(population), 2), 85:95])) - 
#                         round(rowSums(population[seq(2, nrow(population), 2), 85:95]))) 
#   
#   # add coordinate information
#   districts <- districts[order(districts$districtId), ] %>% 
#     mutate(name = coordinates$name, 
#            lon = as.numeric(coordinates$longitude), 
#            lat = as.numeric(coordinates$latitude), 
#            density = pop.density$perkm2)
#   
#   # add deprivation information
#   #districts <-  cbind(districts, deprivation[,3:13])
#   
#   return(districts)
# }
# 
# library(readr)
# library(readxl)
# library(dplyr)
# library(lubridate)
# 
# path.LRZ <<- "C:/Users/ra98jiq/LRZ Sync+Share/Corona/"  
# 
# data.RKI <- read.csv("https://www.arcgis.com/sharing/rest/content/items/f10774f1c63e40168479a1feb6c7ca74/data", sep = ",",
#                      fileEncoding = "UTF-8")
# 
# write_excel_csv(data.RKI, file = paste0(path.LRZ, "Data/RKI/RKI_COVID19_", Sys.Date(), ".csv"), delim = ",")
# 
# source("C:/Users/ru58paj/Documents/codag-lmu/R/Code/Functions/Preprocessing.R")
# 
# read.RKI()
