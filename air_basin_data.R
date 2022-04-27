##################################################################
##########      Load library and data
##################################################################

rm(list=ls())
library(foreign)
library(data.table)
library(ggplot2)
library(stringr)
library(tidyselect)
library(rlist)
library(geoR)
library(usmap)
library(raster)
library(sp)
library(maps)
library(maptools)
library(dplyr)
library(spNNGP)
library(geosphere)
library(rgdal)
library(sf)
library(RgoogleMaps)
library(ggmap)
library(rgeos)
library(lubridate)
library(geodist)

# load GEOS,GDAL,PROJ,UDUNITS


#####################################################################
##########      set directory
#####################################################################

# path to home directory
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

home_dir <- paste(code.path.splitted[1: (length(code.path.splitted)-2)], collapse = "/")
data_dir<-paste0(home_dir,'/data/')
res_dir<-paste0(home_dir,'/results/')


#####################################################################
##########      process SCAB shapefile
#####################################################################

## read SCAB boundary file
setwd(data_dir)

poly.layer <- paste('SCAB')
poly.path<-'SCAB'

SCAB_shp <- readOGR(dsn = poly.path,verbose = FALSE,
                  layer = as.character(poly.layer))

#plot(SCAB_shp[1,])

## smooth edges
tmp.shp <- gSimplify(SCAB_shp,tol=0.01)

SCAB_shp_clean <- tmp.shp[1,]

geo_proj = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
SCAB_shp_clean = spTransform(SCAB_shp_clean,geo_proj)


#plot(SCAB_shp_clean)

## make googlemap

scab_stamen <- get_stamenmap(bbox = c(left = bbox(SCAB_shp_clean)[1,1] - 0.1,
                                             
                                             bottom =  bbox(SCAB_shp_clean)[2,1] - 0.1,
                                             
                                             right =  bbox(SCAB_shp_clean)[1,2] + 0.1,
                                             
                                             top =  bbox(SCAB_shp_clean)[2,2] + 0.1), maptype = 'terrain',
                           zoom=8)


#ggmap(scab_stamen)

scab_map <- ggmap(scab_stamen) + geom_polygon(data=fortify(SCAB_shp_clean),
                                                  aes(long, lat,group=piece),colour='blue', fill=NA)+
  ggtitle('South carolina air basin map')+
  theme(plot.title = element_text(hjust = 0.5))

##################################################################
##########    prepare 2-week wednesday
##################################################################
mesa_interval <- seq(as.Date("1998-12-09"), by = "14 days", length.out = 600)
mesa_interval2 <- seq(as.Date("1998-12-16"), by = "14 days", length.out = 600)


mesa_wednesday <- mesa_interval + 7

tmp.t <-  c(as.Date("2008-01-01"), as.Date("2010-12-31"))
long_range.t <-  c(as.Date("2000-01-01"), as.Date("2021-12-31"))
early_range.t <-  c(as.Date("2000-01-01"), as.Date("2010-12-31"))

##################################################################
##########    prepare mesa-agency no2 data
##################################################################


setwd(paste0(data_dir,"/dr0328/"))

### mesa agency 
mesa_ag_cov <- read.delim('dr0328_mesa_agency_covars_NO2.txt',sep=',')
mesa_ag_no2 <- read.delim('dr0328_mesa_agency_NO2.txt',sep=',')

# prepare locations within SCAB
mesa_ag_loc <- mesa_ag_cov[,c('native_id','longitude','latitude')]
mesa_ag_loc <- mesa_ag_loc[complete.cases(mesa_ag_loc),]
coordinates(mesa_ag_loc) <- ~ longitude+latitude
proj4string(mesa_ag_loc) <- proj4string(SCAB_shp_clean)

mesa_ag_loc_inscab <- over(mesa_ag_loc,SCAB_shp_clean)
mesa_ag_loc$scab <- mesa_ag_loc_inscab

mesa_ag_loc_scab <- data.frame(mesa_ag_loc)
mesa_ag_loc_scab <- mesa_ag_loc_scab[complete.cases(mesa_ag_loc_scab),]
mesa_ag_loc_scab <- mesa_ag_loc_scab[,c(1:3)]


mesa_ag_loc_scab_sp <- mesa_ag_loc_scab
coordinates(mesa_ag_loc_scab_sp) <- ~ longitude+latitude
proj4string(mesa_ag_loc_scab_sp) <- proj4string(SCAB_shp_clean)

mesa_ag_dist <- geodist(mesa_ag_loc_scab[,c(3:2)])

#plot(SCAB_shp_clean)
#plot(mesa_ag_loc_scab_sp,add=T)




# prepare N02 within SCAB
mesa_ag_no2_scab <- mesa_ag_no2[mesa_ag_no2$native_id%in%mesa_ag_loc_scab$native_id, ]
mesa_ag_no2_scab <- merge(mesa_ag_no2_scab,mesa_ag_loc_scab,by='native_id')

mesa_ag_no2_scab$start_wednesday <- mesa_interval[ findInterval(as.Date(mesa_ag_no2_scab$day),
                                                                 mesa_interval) ]

mesa_ag_no2_scab$intended_wednesday <- mesa_ag_no2_scab$start_wednesday +7
  
# keep only 2-week periods with >= 12 observations
mesa_ag_no2_scab <- as.data.table(mesa_ag_no2_scab)
mesa_ag_no2_scab_2wk <- mesa_ag_no2_scab[, .(count=.N, no2_2wk = mean(NO2_conc,na.rm=T)),
                                          by=.(native_id, intended_wednesday)][ count >11 ]

mesa_ag_no2_scab_2wk <- merge(mesa_ag_no2_scab_2wk,mesa_ag_loc_scab,by='native_id')



# visualize 
mesa_ag_tmp <- mesa_ag_no2_scab_2wk[ intended_wednesday >= long_range.t[1] &
                                       intended_wednesday <= long_range.t[2]]

mesa_ag_loc_tmp <- mesa_ag_tmp[, head(.SD, 1), by=.(native_id)]


mesa_ag_monitor_loc <- ggmap(scab_stamen) + geom_polygon(data=fortify(SCAB_shp_clean),
                                                  aes(long, lat,group=piece),colour='blue', fill=NA)+
  geom_point(data = mesa_ag_loc_tmp,
             mapping = aes(x = longitude, y = latitude), color='red',
             shape = 20, size = 3, stroke = 0.7, alpha = 0.7) +
  ggtitle('SCAB monitors locations, MESA agency, 2008-2010')+
  theme(plot.title = element_text(hjust = 0.5))


mesa_ag_monitor_tmp_seg <- ggplot(data = mesa_ag_tmp) +
  geom_point(mapping = aes(x = intended_wednesday, y = as.character(native_id)), size = 2) +
  theme_bw(base_size = 14)+
  xlab('time')+ylab('monitor_id')+ ggtitle('MESA agency monitors in SCAB')+
  theme(plot.title = element_text(hjust = 0.5))


### store coordinates for mesa agency, 2000-2010
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
mesa_ag_loc_scab_sf <- st_as_sf(x = mesa_ag_loc_scab,                         
               coords = c("longitude", "latitude"),
               crs = projcrs)

setwd(paste0(data_dir,"/site_loc/"))

if(!dir.exists(paths = paste0('mesa_ag'))){
  dir.create(path = paste0('mesa_ag'))
}

st_write(mesa_ag_loc_scab_sf, "mesa_ag/mesa_ag.shp")


### 60376002 and 60376012 might be the same site
### mesa_L001 co-locate
##################################################################
##########    prepare mesa no2 data
##################################################################

### mesa
mesa_cov <- read.delim('dr0328_mesa_covars_NO2.txt',sep=',')
mesa_no2 <- read.delim('dr0328_mesa_NO2.txt',sep=',')


### check whether mesa data do have intended wednesday
if(FALSE){
mesa_no2_comp <- mesa_no2_scab[complete.cases(mesa_no2_scab),]
mesa_no2_comp$intended_wednesday <- as.Date(mesa_no2_comp$intended_wednesday)
dim(mesa_no2_comp)
dim(mesa_no2_comp[mesa_no2_comp$intended_wednesday %in% mesa_wednesday,] )

mesa_no2_comp$mark <- mesa_no2_comp$intended_wednesday %in% mesa_wednesday
}
#length(unique(mesa_no2$location_id) %in% mesa_cov$location_id)
# locations
mesa_loc <- mesa_cov[,c('native_id','longitude','latitude')]
mesa_loc <- mesa_loc[complete.cases(mesa_loc),]
coordinates(mesa_loc) <- ~ longitude+latitude
proj4string(mesa_loc) <- proj4string(SCAB_shp_clean)

mesa_loc_inscab <- over(mesa_loc,SCAB_shp_clean)
mesa_loc$scab <- mesa_loc_inscab

mesa_loc_scab <- data.frame(mesa_loc)
mesa_loc_scab <- mesa_loc_scab[complete.cases(mesa_loc_scab),]

mesa_no2$native_id <- mesa_no2$site_id
mesa_no2_scab <- mesa_no2[mesa_no2$native_id%in%mesa_loc_scab$native_id, ]
mesa_no2_scab <- merge(mesa_no2_scab,mesa_loc_scab,by='native_id')


# add intended wednesday for missing 
mesa_no2_scab$sample_start_date <- as.Date(mesa_no2_scab$sample_start_date)
mesa_no2_scab$sample_stop_date <- as.Date(mesa_no2_scab$sample_stop_date)

mesa_no2_scab$mid_day <-  mesa_no2_scab$sample_start_date + 
  floor((mesa_no2_scab$sample_stop_date-mesa_no2_scab$sample_start_date)/2)

mesa_no2_scab$candidate_wed <- round_date(mesa_no2_scab$mid_day, "week",3)

mesa_no2_scab <- as.data.table(mesa_no2_scab)

mesa_no2_scab <- mesa_no2_scab[is.na(intended_wednesday)|intended_wednesday == "", 
                               ind := 'Missing']

mesa_no2_scab$intended_wednesday <- as.Date(mesa_no2_scab$intended_wednesday)

mesa_no2_scab <- mesa_no2_scab[ind == 'Missing', 
                               intended_wednesday := candidate_wed]


# classify by intended wednesdays

mesa_no2_scab$wed_week <- (as.numeric(mesa_no2_scab$intended_wednesday -as.Date("1998-12-09"))%% 14)/7

mesa_no2_scab <- mesa_no2_scab[wed_week == 1, 
                               wed_type := 'L']
mesa_no2_scab <- mesa_no2_scab[wed_week == 0, 
                               wed_type := 'LR']

## check 
site_wed <- unique(as.data.frame(mesa_no2_scab)[c("native_id", "wed_type")])


# visualize 
mesa_tmp <- mesa_no2_scab[ intended_wednesday >= early_range.t[1] &
                             intended_wednesday <= early_range.t[2]]

mesa_tmp <- mesa_tmp [ ,site_LR := ifelse(mean(wed_week)>0.5,'L','LR'), by=native_id]


mesa_loc_tmp <- mesa_tmp[, head(.SD, 1), by=.(native_id)]
#mesa_loc_tmp <- mesa_loc_tmp[mesa_loc_tmp$native_id%in% c('mesa_L001','mesa_L002','mesa_LC001','mesa_LC002'),]

mesa_monitor_loc <-  ggmap(scab_stamen) + geom_polygon(data=fortify(SCAB_shp_clean),
                                                       aes(long, lat,group=piece),colour='blue', fill=NA)+
  geom_point(data = mesa_loc_tmp,
             mapping = aes(x = longitude, y = latitude,color=site_LR), 
             shape = 20, size = 2, stroke = 0.7, alpha = 0.7) +
  ggtitle('SCAB monitors locations, MESA, 2000-2018')+
  theme(plot.title = element_text(hjust = 0.5))


mesa_monitor_tmp_seg <-  ggplot(data = mesa_tmp) +
  geom_point(mapping = aes(x = intended_wednesday, y = as.character(native_id),color=site_LR), size = 1) +
  theme_bw(base_size = 14)+
  xlab('time')+ylab('monitor_id')+ ggtitle('MESA monitors in SCAB')+
  theme(plot.title = element_text(hjust = 0.5))


## early range, L and LR separate well

### store coordinates for mesa agency, 2000-2010
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
mesa_loc_scab_sf <- st_as_sf(x = mesa_loc_tmp,                         
                                coords = c("longitude", "latitude"),
                                crs = projcrs)

setwd(paste0(data_dir,"/site_loc/"))

if(!dir.exists(paths = paste0('mesa'))){
  dir.create(path = paste0('mesa'))
}

st_write(mesa_loc_scab_sf, "mesa/mesa.shp")


##################################################################
##########    prepare spiromics agency no2 data
##################################################################

### spiromics  
spiromics_ag_no2 <- read.delim('dr0328_spiromics_agency_NO2.txt',sep=',')

#plot(SCAB_shp_clean)
#plot(mesa_ag_loc_scab_sp,add=T)

# prepare N02 within SCAB
spiromics_ag_no2_scab <- spiromics_ag_no2[spiromics_ag_no2$native_id%in%mesa_ag_loc_scab$native_id, ]
spiromics_ag_no2_scab <- merge(spiromics_ag_no2_scab,mesa_ag_loc_scab,by='native_id')




spiromics_ag_no2_scab$start_wednesday <- mesa_interval[ findInterval(as.Date(spiromics_ag_no2_scab$day),
                                                                mesa_interval) ]

spiromics_ag_no2_scab$intended_wednesday <- spiromics_ag_no2_scab$start_wednesday +7

# keep only 2-week periods with >= 12 observations
spiromics_ag_no2_scab <- as.data.table(spiromics_ag_no2_scab)
spiromics_ag_no2_scab_2wk <- spiromics_ag_no2_scab[, .(count=.N, no2_2wk = mean(NO2_conc,na.rm=T)),
                                         by=.(native_id, intended_wednesday)][ count >11 ]

spiromics_ag_no2_scab_2wk <- merge(spiromics_ag_no2_scab_2wk,mesa_ag_loc_scab,by='native_id')



# visualize 
spiromics_ag_long_range <- spiromics_ag_no2_scab_2wk[ intended_wednesday >= long_range.t[1] &
                                       intended_wednesday <= long_range.t[2]]

spiromics_ag_loc_long_range <- spiromics_ag_long_range[, head(.SD, 1), by=.(native_id)]


spiromics_ag_monitor_loc_long_range <- ggmap(scab_stamen) + geom_polygon(data=fortify(SCAB_shp_clean),
                                                         aes(long, lat,group=piece),colour='blue', fill=NA)+
  geom_point(data = spiromics_ag_loc_long_range,
             mapping = aes(x = longitude, y = latitude), color='red',
             shape = 20, size = 3, stroke = 0.7, alpha = 0.7) +
  ggtitle('SCAB monitors locations, spiromics agency, 2008-2010')+
  theme(plot.title = element_text(hjust = 0.5))


spiromics_ag_monitor_seg_long_range <- ggplot(data = spiromics_ag_long_range) +
  geom_point(mapping = aes(x = intended_wednesday, y = as.character(native_id)), size = 2) +
  theme_bw(base_size = 14)+
  xlab('time')+ylab('monitor_id')+ ggtitle('Spiromics agency monitors in SCAB')+
  theme(plot.title = element_text(hjust = 0.5))



### compare two agency 
mesa_comp <- mesa_ag_no2_scab[,c('native_id','NO2_conc','day')]
spiromics_comp <- spiromics_ag_no2_scab[,c('native_id','NO2_conc','day')]
all_equal(mesa_comp, spiromics_comp)


##################################################################
##########    prepare spiromics no2 data
##################################################################

spiromics_no2 <- read.delim('dr0328_spiromics_ogawa_monitoring.txt',sep=',')
spiromics_no2 <- spiromics_no2[spiromics_no2$poll=='NO2',]

#length(unique(mesa_no2$location_id) %in% mesa_cov$location_id)
# locations

spiromics_no2_scab <- spiromics_no2[spiromics_no2$native_id%in%mesa_loc_scab$native_id, ]
#spiromics_no2_scab <- merge(spiromics_no2_scab,mesa_loc_scab,by='native_id')


# add intended wednesday for missing 
spiromics_no2_scab$sample_start_date <- as.Date(spiromics_no2_scab$sample_start_date)
spiromics_no2_scab$sample_stop_date <- as.Date(spiromics_no2_scab$sample_stop_date)

spiromics_no2_scab$mid_day <-  spiromics_no2_scab$sample_start_date + 
  floor((spiromics_no2_scab$sample_stop_date-spiromics_no2_scab$sample_start_date)/2)

spiromics_no2_scab$intended_wednesday <- round_date(spiromics_no2_scab$mid_day, "week",3)

spiromics_no2_scab <- as.data.table(spiromics_no2_scab)

spiromics_no2_scab$intended_wednesday <- as.Date(spiromics_no2_scab$intended_wednesday)



# classify by intended wednesdays

spiromics_no2_scab$wed_week <- (as.numeric(spiromics_no2_scab$intended_wednesday -as.Date("1998-12-09"))%% 14)/7

spiromics_no2_scab <- spiromics_no2_scab[wed_week == 1, 
                               wed_type := 'L']
spiromics_no2_scab <- spiromics_no2_scab[wed_week == 0, 
                               wed_type := 'LR']

## check 
site_wed <- unique(as.data.frame(spiromics_no2_scab)[c("native_id", "wed_type")])


# visualize 
spiromics_tmp <- spiromics_no2_scab[ intended_wednesday >= early_range.t[1] &
                             intended_wednesday <= early_range.t[2]]

spiromics_tmp <- spiromics_tmp [ ,site_LR := ifelse(mean(wed_week)>0.5,'L','LR'), by=native_id]


spiromics_loc_tmp <- spiromics_tmp[, head(.SD, 1), by=.(native_id)]



spiromics_monitor_loc <- ggmap(scab_stamen) + geom_polygon(data=fortify(SCAB_shp_clean),
                                                                 aes(long, lat,group=piece),colour='blue', fill=NA)+
  geom_point(data = spiromics_loc_tmp,
             mapping = aes(x = longitude, y = latitude,color=site_LR), 
             shape = 20, size = 2, stroke = 0.7, alpha = 0.7) +
  ggtitle('SCAB monitors locations, MESA, 2000-2018')+
  theme(plot.title = element_text(hjust = 0.5))

spiromics_monitor_tmp_seg <- ggplot(data = spiromics_tmp) +
  geom_point(mapping = aes(x = intended_wednesday, y = as.character(native_id),color=site_LR), size = 1) +
  theme_bw(base_size = 14)+
  xlab('time')+ylab('monitor_id')+ ggtitle('spiromics monitors in SCAB, 2000-2018')+
  theme(plot.title = element_text(hjust = 0.5))

## early range, L and LR separate well



### compare mesa and spiromics short term

spiromics_tab <- table(spiromics_tmp$native_id)
mesa_tab <- table(mesa_tmp$native_id)

spiromics_tab-mesa_tab
