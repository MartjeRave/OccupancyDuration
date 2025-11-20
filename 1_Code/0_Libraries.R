#####################
#### Packages
#####################

#install.packages("rgdal")
#library(maptools)
#library(rgeos)
#library(taskscheduleR)

library(abind)
library(cluster)
library(colorspace)
library(corrplot)
library(cowplot)
library(data.table)
library(dplyr)
library(fitdistrplus)
library(forecast)
library(gam)
library(ggdendro)
library(ggplot2)
library(ggpubr)
library(ggspatial)
library(glmmLasso)
library(glmnet)
library(gridExtra)
library(gtable)
library(gtools)
library(Hmisc)
library(ISLR)
library(kableExtra)
library(lazyWeave)
library(leaflet)
library(lme4)
library(lubridate)
library(mefa)
library(mgcv)
library(miceadds)
library(MLmetrics)
library(plotmo)
library(plyr)
library(psych)
library(quadprog)
library(RColorBrewer)
library(RCurl)
library(readr)
library(readxl)
library(Rmisc)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)
library(sf)
library(shapefiles)
library(skellam)
library(sp)
library(splines)
library(sqldf)
library(terra)
library(tidyverse)
library(VGAM)
library(zoo)



#####################
#### Custom Functions:
##################### 
# written by Dr. Marc Schneble
# This is to grab the population size of each district
# (Definitely not the most efficient to go about it but it's legacy)

preprocess.districts <- function(path.LRZ){
  
  # read population and coordinates of districts
  coordinates <- read_excel(paste(path.LRZ, "/2. Data/coordinates.xlsx", sep = ""))
  population <- read_excel(paste(path.LRZ, "/2. Data/population.xlsx", sep = ""))
  pop.density <- read_excel(paste(path.LRZ, "/2. Data/population_total.xlsx", sep = ""))
  #deprivation <- read_excel(paste(path.LRZ, "Data/Demographics/deprivation.xlsx", sep = ""))
  
  
  districts <- tibble(districtId = as.numeric(population$districtId[seq(1, nrow(population), 2)]),
                      pop = round(population$gesamt[seq(1, nrow(population), 2)]),
                      pop.m = round(population$gesamt[seq(2, nrow(population), 2)]),
                      pop.f = pop - pop.m,
                      pop.m.0.4 = round(rowSums(population[seq(2, nrow(population), 2), 5:9])),
                      pop.w.0.4 = round(rowSums(population[seq(1, nrow(population), 2), 5:9])) -
                        round(rowSums(population[seq(2, nrow(population), 2), 5:9])),
                      pop.m.5.14 = round(rowSums(population[seq(2, nrow(population), 2), 10:19])),
                      pop.w.5.14 = round(rowSums(population[seq(1, nrow(population), 2), 10:19])) -
                        round(rowSums(population[seq(2, nrow(population), 2), 10:19])),
                      pop.m.15.34 = round(rowSums(population[seq(2, nrow(population), 2), 20:39])),
                      pop.w.15.34 = round(rowSums(population[seq(1, nrow(population), 2), 20:39])) -
                        round(rowSums(population[seq(2, nrow(population), 2), 20:39])),
                      pop.m.35.59 = round(rowSums(population[seq(2, nrow(population), 2), 40:64])),
                      pop.w.35.59 = round(rowSums(population[seq(1, nrow(population), 2), 40:64])) -
                        round(rowSums(population[seq(2, nrow(population), 2), 40:64])),
                      pop.m.60.79 = round(rowSums(population[seq(2, nrow(population), 2), 65:84])),
                      pop.w.60.79 = round(rowSums(population[seq(1, nrow(population), 2), 65:84])) -
                        round(rowSums(population[seq(2, nrow(population), 2), 65:84])),
                      pop.m.80 = round(rowSums(population[seq(2, nrow(population), 2), 85:95])),
                      pop.w.80 = round(rowSums(population[seq(1, nrow(population), 2), 85:95])) -
                        round(rowSums(population[seq(2, nrow(population), 2), 85:95])))
  
  # add coordinate information
  districts <- districts[order(districts$districtId), ] %>%
    mutate(name = coordinates$name,
           lon = as.numeric(coordinates$longitude),
           lat = as.numeric(coordinates$latitude),
           density = pop.density$perkm2)
  
  # add deprivation information
  #districts <-  cbind(districts, deprivation[,3:13])
  
  return(districts)
}


################################################################################
######### R-session of author ##################################################
################################################################################
# 
# 
# R version 4.5.1 (2025-06-13 ucrt) -- "Great Square Root"
# Copyright (C) 2025 The R Foundation for Statistical Computing
# Platform: x86_64-w64-mingw32/x64
# 
# R is free software and comes with ABSOLUTELY NO WARRANTY.
# You are welcome to redistribute it under certain conditions.
# Type 'license()' or 'licence()' for distribution details.
# 
# R is a collaborative project with many contributors.
# Type 'contributors()' for more information and
# 'citation()' on how to cite R or R packages in publications.
# 
# Type 'demo()' for some demos, 'help()' for on-line help, or
# 'help.start()' for an HTML browser interface to help.
# Type 'q()' to quit R.
# 
# > library(abind)
# > library(cluster)
# > library(colorspace)
# > library(corrplot)
# corrplot 0.95 loaded
# > library(cowplot)
# > library(data.table)
# data.table 1.17.6 using 4 threads (see ?getDTthreads).  Latest news: r-datatable.com
# > library(dplyr)
# 
# Attaching package: ‘dplyr’
# 
# The following objects are masked from ‘package:data.table’:
#   
#   between, first, last
# 
# The following objects are masked from ‘package:stats’:
#   
#   filter, lag
# 
# The following objects are masked from ‘package:base’:
#   
#   intersect, setdiff, setequal, union
# 
# > library(fitdistrplus)
# Loading required package: MASS
# 
# Attaching package: ‘MASS’
# 
# The following object is masked from ‘package:dplyr’:
#   
#   select
# 
# Loading required package: survival
# > library(forecast)
# Registered S3 method overwritten by 'quantmod':
#   method            from
# as.zoo.data.frame zoo 
# > library(gam)
# Loading required package: splines
# Loading required package: foreach
# Loaded gam 1.22-5
# 
# > library(ggdendro)
# > library(ggplot2)
# > library(ggpubr)
# 
# Attaching package: ‘ggpubr’
# 
# The following object is masked from ‘package:forecast’:
#   
#   gghistogram
# 
# The following object is masked from ‘package:cowplot’:
#   
#   get_legend
# 
# > library(ggspatial)
# > library(glmmLasso)
# > library(glmnet)
# Loading required package: Matrix
# Loaded glmnet 4.1-9
# > library(gridExtra)
# 
# Attaching package: ‘gridExtra’
# 
# The following object is masked from ‘package:dplyr’:
#   
#   combine
# 
# > library(gtable)
# > library(gtools)
# 
# Attaching package: ‘gtools’
# 
# The following object is masked from ‘package:glmnet’:
#   
#   na.replace
# 
# > library(Hmisc)
# 
# Attaching package: ‘Hmisc’
# 
# The following object is masked from ‘package:ggdendro’:
#   
#   label
# 
# The following objects are masked from ‘package:dplyr’:
#   
#   src, summarize
# 
# The following objects are masked from ‘package:base’:
#   
#   format.pval, units
# 
# > library(ISLR)
# > library(kableExtra)
# 
# Attaching package: ‘kableExtra’
# 
# The following object is masked from ‘package:dplyr’:
#   
#   group_rows
# 
# > library(lazyWeave)
# > library(leaflet)
# > library(lme4)
# > library(lubridate)
# 
# Attaching package: ‘lubridate’
# 
# The following objects are masked from ‘package:data.table’:
#   
#   hour, isoweek, mday, minute, month, quarter, second, wday, week, yday, year
# 
# The following object is masked from ‘package:cowplot’:
#   
#   stamp
# 
# The following objects are masked from ‘package:base’:
#   
#   date, intersect, setdiff, union
# 
# > library(mefa)
# mefa 3.2-9 	 2024-05-19
# 
# Attaching package: ‘mefa’
# 
# The following objects are masked from ‘package:Hmisc’:
#   
#   label, label<-
#   
#   The following object is masked from ‘package:ggdendro’:
#   
#   label
# 
# The following object is masked from ‘package:data.table’:
#   
#   melt
# 
# > library(mgcv)
# Loading required package: nlme
# 
# Attaching package: ‘nlme’
# 
# The following object is masked from ‘package:lme4’:
#   
#   lmList
# 
# The following object is masked from ‘package:forecast’:
#   
#   getResponse
# 
# The following object is masked from ‘package:dplyr’:
#   
#   collapse
# 
# This is mgcv 1.9-3. For overview type 'help("mgcv-package")'.
# 
# Attaching package: ‘mgcv’
# 
# The following object is masked from ‘package:gtools’:
#   
#   scat
# 
# The following objects are masked from ‘package:gam’:
#   
#   gam, gam.control, gam.fit, s
# 
# > library(miceadds)
# Loading required package: mice
# 
# Attaching package: ‘mice’
# 
# The following object is masked from ‘package:stats’:
#   
#   filter
# 
# The following objects are masked from ‘package:base’:
#   
#   cbind, rbind
# 
# * miceadds 3.17-44 (2024-01-08 19:08:24)
# > library(MLmetrics)
# 
# Attaching package: ‘MLmetrics’
# 
# The following object is masked from ‘package:base’:
#   
#   Recall
# 
# > library(plotmo)
# Loading required package: Formula
# Loading required package: plotrix
# > library(plyr)
# ------------------------------------------------------------------------------------------------------------
#   You have loaded plyr after dplyr - this is likely to cause problems.
# If you need functions from both plyr and dplyr, please load plyr first, then dplyr:
#   library(plyr); library(dplyr)
# ------------------------------------------------------------------------------------------------------------
#   
#   Attaching package: ‘plyr’
# 
# The following objects are masked from ‘package:Hmisc’:
#   
#   is.discrete, summarize
# 
# The following object is masked from ‘package:ggpubr’:
#   
#   mutate
# 
# The following objects are masked from ‘package:dplyr’:
#   
#   arrange, count, desc, failwith, id, mutate, rename, summarise, summarize
# 
# > library(psych)
# 
# Attaching package: ‘psych’
# 
# The following object is masked from ‘package:plotrix’:
#   
#   rescale
# 
# The following object is masked from ‘package:MLmetrics’:
#   
#   AUC
# 
# The following object is masked from ‘package:Hmisc’:
#   
#   describe
# 
# The following object is masked from ‘package:gtools’:
#   
#   logit
# 
# The following objects are masked from ‘package:ggplot2’:
#   
#   %+%, alpha
# 
# > library(quadprog)
# > library(RColorBrewer)
# > library(RCurl)
# 
# Attaching package: ‘RCurl’
# 
# The following object is masked from ‘package:mice’:
#   
#   complete
# 
# > library(readr)
# > library(readxl)
# > library(Rmisc)
# Loading required package: lattice
# > library(rnaturalearth)
# > library(rnaturalearthdata)
# 
# Attaching package: ‘rnaturalearthdata’
# 
# The following object is masked from ‘package:rnaturalearth’:
#   
#   countries110
# 
# > library(scales)
# 
# Attaching package: ‘scales’
# 
# The following object is masked from ‘package:readr’:
#   
#   col_factor
# 
# The following objects are masked from ‘package:psych’:
#   
#   alpha, rescale
# 
# The following object is masked from ‘package:plotrix’:
#   
#   rescale
# 
# > library(sf)
# Linking to GEOS 3.13.1, GDAL 3.11.0, PROJ 9.6.0; sf_use_s2() is TRUE
# > library(shapefiles)
# Loading required package: foreign
# 
# Attaching package: ‘shapefiles’
# 
# The following objects are masked from ‘package:foreign’:
#   
#   read.dbf, write.dbf
# 
# > library(skellam)
# > library(sp)
# > library(splines)
# > library(sqldf)
# Loading required package: gsubfn
# Loading required package: proto
# Loading required package: RSQLite
# > library(terra)
# terra 1.8.54
# 
# Attaching package: ‘terra’
# 
# The following object is masked from ‘package:scales’:
#   
#   rescale
# 
# The following objects are masked from ‘package:psych’:
#   
#   describe, distance, rescale
# 
# The following object is masked from ‘package:plotrix’:
#   
#   rescale
# 
# The following objects are masked from ‘package:Hmisc’:
#   
#   describe, mask, zoom
# 
# The following object is masked from ‘package:ggpubr’:
#   
#   rotate
# 
# The following object is masked from ‘package:MASS’:
#   
#   area
# 
# The following object is masked from ‘package:data.table’:
#   
#   shift
# 
# The following object is masked from ‘package:colorspace’:
#   
#   RGB
# 
# > library(tidyverse)
# ── Attaching core tidyverse packages ────────────────────────────────────────────────────── tidyverse 2.0.0 ──
# ✔ forcats 1.0.0     ✔ tibble  3.3.0
# ✔ purrr   1.0.4     ✔ tidyr   1.3.1
# ✔ stringr 1.5.1     
# ── Conflicts ──────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
# ✖ psych::%+%()             masks ggplot2::%+%()
# ✖ purrr::accumulate()      masks foreach::accumulate()
# ✖ scales::alpha()          masks psych::alpha(), ggplot2::alpha()
# ✖ plyr::arrange()          masks dplyr::arrange()
# ✖ dplyr::between()         masks data.table::between()
# ✖ scales::col_factor()     masks readr::col_factor()
# ✖ nlme::collapse()         masks dplyr::collapse()
# ✖ gridExtra::combine()     masks dplyr::combine()
# ✖ purrr::compact()         masks plyr::compact()
# ✖ tidyr::complete()        masks RCurl::complete(), mice::complete()
# ✖ plyr::count()            masks dplyr::count()
# ✖ plyr::desc()             masks dplyr::desc()
# ✖ purrr::discard()         masks scales::discard()
# ✖ tidyr::expand()          masks Matrix::expand()
# ✖ tidyr::extract()         masks terra::extract()
# ✖ plyr::failwith()         masks dplyr::failwith()
# ✖ mice::filter()           masks dplyr::filter(), stats::filter()
# ✖ dplyr::first()           masks data.table::first()
# ✖ kableExtra::group_rows() masks dplyr::group_rows()
# ✖ lubridate::hour()        masks data.table::hour()
# ✖ plyr::id()               masks dplyr::id()
# ✖ lubridate::isoweek()     masks data.table::isoweek()
# ✖ dplyr::lag()             masks stats::lag()
# ✖ dplyr::last()            masks data.table::last()
# ✖ lubridate::mday()        masks data.table::mday()
# ✖ lubridate::minute()      masks data.table::minute()
# ✖ lubridate::month()       masks data.table::month()
# ✖ plyr::mutate()           masks ggpubr::mutate(), dplyr::mutate()
# ✖ tidyr::pack()            masks Matrix::pack()
# ✖ lubridate::quarter()     masks data.table::quarter()
# ✖ plyr::rename()           masks dplyr::rename()
# ✖ lubridate::second()      masks data.table::second()
# ✖ MASS::select()           masks dplyr::select()
# ✖ Hmisc::src()             masks dplyr::src()
# ✖ lubridate::stamp()       masks cowplot::stamp()
# ✖ plyr::summarise()        masks dplyr::summarise()
# ✖ plyr::summarize()        masks Hmisc::summarize(), dplyr::summarize()
# ✖ purrr::transpose()       masks data.table::transpose()
# ✖ tidyr::unpack()          masks Matrix::unpack()
# ✖ lubridate::wday()        masks data.table::wday()
# ✖ lubridate::week()        masks data.table::week()
# ✖ purrr::when()            masks foreach::when()
# ✖ lubridate::yday()        masks data.table::yday()
# ✖ lubridate::year()        masks data.table::year()
# ℹ Use the conflicted package to force all conflicts to become errors
# > library(VGAM)
# Loading required package: stats4
# 
# Attaching package: ‘VGAM’
# 
# The following objects are masked from ‘package:skellam’:
#   
#   dskellam, rskellam
# 
# The following objects are masked from ‘package:psych’:
#   
#   fisherz, logistic, logit
# 
# The following object is masked from ‘package:miceadds’:
#   
#   round2
# 
# The following object is masked from ‘package:mgcv’:
#   
#   s
# 
# The following object is masked from ‘package:gtools’:
#   
#   logit
# 
# The following objects are masked from ‘package:glmmLasso’:
#   
#   acat, cumulative
# 
# The following object is masked from ‘package:gam’:
#   
#   s
# 
# > library(zoo)
# 
# Attaching package: ‘zoo’
# 
# The following object is masked from ‘package:terra’:
#   
#   time<-
#   
#   The following objects are masked from ‘package:data.table’:
#   
#   yearmon, yearqtr
# 
# The following objects are masked from ‘package:base’:
#   
#   as.Date, as.Date.numeric
