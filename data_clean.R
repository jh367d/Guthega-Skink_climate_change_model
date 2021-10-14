


library (galah)
library (stringr)
library (rgdal)
library(sf)
library(maps)
library(mapdata)
library(dplyr)
library(tidyr)
library(ggplot2)



#data from ALA for Eucalypts

library(readr)
records.2021.10.07._Liopholis_guthega_update <- read.csv("~/Documents/R/git/records-2021-10-07-_Liopholis_guthega_update.csv")


#plot data on a map to get initial view
library(ggplot2)
wm <- borders("world", colour="slategray1", fill="slategray1")
ggplot()+ coord_fixed()+ wm +
  geom_point(data = records.2021.07.21, aes(x = decimalLongitude, y = decimalLatitude),
             colour = "green", size = 0.5)+
  theme_bw()

#now clean with coordinate cleaner 
library(countrycode)
library(CoordinateCleaner)
library(rgbif)
library(sp)
library(dplyr)

#columns of interest
records_filter_1 <- records.2021.10.07._Liopholis_guthega_update %>%
  dplyr::select(species, lon, lat, countryCode, individualCount,
                family, taxonRank, coordinateUncertaintyInMeters, year,
                basisOfRecord, institutionCode, datasetName)


#filter records to only Australia
records_filter_au <- records_filter_1 %>%
  filter(countryCode == "AU")

#plot data on a map to get  view
library(ggplot2)
wm <- borders("world", colour="slategray1", fill="slategray1")
ggplot()+ coord_fixed()+ wm +
  geom_point(data = records_filter_au , aes(x = lon, y = lat),
             colour = "green", size = 0.5)+
  theme_bw()




# get rid of records without coordinates as ALA did not get rid of all records without coordinates
records_filter_au_1  <- records_filter_au %>%
  filter(!is.na(lon))%>%
  filter(!is.na(lat))





#automatic cleaning algorithm of CoordinateCleaner
library(rnaturalearth)
#conversion of country code from ISO2c to ISO3c
records_filter_au_1$countryCode <-  countrycode(records_filter_au_1$countryCode, origin =  'iso2c', destination = 'iso3c')


#flag problems (won't include 'country' test as it flags all records)
records_filter_au_1 <- data.frame(records_filter_au_1)
flags <- clean_coordinates(x = records_filter_au_1 ,
                           lon = "lon",
                           lat = "lat",
                           countries = "countryCode",
                           species = "species",
                           tests = c("capitals", "centroids", "equal","gbif", "institutions",
                                     "zeros" )) 
summary(flags)
plot(flags, lon = "lon", lat = "lat")

#now remove problematic records
dat_cl_aus <- records_filter_au_1[flags$.summary,]

#The flagged records
dat_fl_aus <- records_filter_au_1[!flags$.summary,]

#plot flagged records

wm <- borders("world", colour="slategray1", fill="slategray1")
ggplot()+ coord_fixed()+ wm +
  geom_point(data = dat_fl_aus , aes(x = lon, y = lat),
             colour = "green", size = 0.5)+
  theme_bw()



#Remove records with low coordinate precision 
#remove all records with a precision below 100 km
hist(dat_cl_aus$coordinateUncertaintyInMeters / 1000, breaks = 20)

dat_cl_aus <- dat_cl_aus %>%
  filter(coordinateUncertaintyInMeters / 1000 <= 100 | is.na(coordinateUncertaintyInMeters))

#Remove unsuitable data sources, 
table(dat_cl_aus$basisOfRecord)



dat_cl_aus <- filter(dat_cl_aus, basisOfRecord == "HUMAN_OBSERVATION" |
                       basisOfRecord == "OBSERVATION" |
                       basisOfRecord == "PRESERVED_SPECIMEN")

#remove records with suspicious individual counts
table(dat_cl_aus$individualCount)



#remove records of absence i.e (individual count = 0) 
dat_cl_aus <- dat_cl_aus%>%
  filter(individualCount > 0 | is.na(individualCount))%>%
  filter(individualCount < 99 | is.na(individualCount)) 

#plot flagged records

wm <- borders("world", colour="slategray1", fill="slategray1")
ggplot()+ coord_fixed()+ wm +
  geom_point(data = dat_cl_aus, aes(x = lon, y = lat),
             colour = "green", size = 0.5)+
  theme_bw()




#export file
write.csv(dat_cl_aus,"/Users/jackhanigan/Documents/R/ecul_data2.csv", row.names = FALSE)

