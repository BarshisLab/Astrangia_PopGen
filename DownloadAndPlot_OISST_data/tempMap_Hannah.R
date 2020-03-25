library(vegan)
library(maptools)
library(dplyr)
library(ggplot2)
library(raster)
library(xts)

setwd('/Users/hannahaichelman/Documents/ODU_MS/Map_Satellite_Temp/FromJP_OISST_data')
#install.packages('gpclib', type='source')
if (!rgeosStatus()) gpclibPermit()
gshhs.f.b <- 'gshhs_f.b'
sf1 <- getRgshhsMap(gshhs.f.b, xlim = c(-85, -65), ylim = c(30, 45)) %>%
  fortify()
sites <- read.csv('GPSCoordinates.csv')

load('oisst_VARI.Rdata')
#~~~~~~~~~
# All you need to do here is just specify the data layer you want to project
# To do that, just change the column id in the var=all.avg$mean line
# Check colnames(all.avg) to see your options
#~~~~~~~~~
all.avg2 <- data.frame(lon=as.numeric(gsub('X\\.?|_\\.?\\d+\\.?\\d+','',rownames(all.avg)))*-1, 
                       lat=as.numeric(gsub('X\\.?\\d+\\.?\\d+_\\.?','',rownames(all.avg))), 
                       var=all.avg$weekRange) #Change the column id here

ggplot.data <- all.avg2

ggplot() + 
  #geom_raster(data = ggplot.data, aes(x=lon, y=lat, fill=var))+
  geom_tile(data = ggplot.data, aes(x=lon, y=lat, fill=var, width=0.3, height=0.3))+
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = 'grey70', color='black', lwd = 0.1) +
  geom_point(data=sites, aes(x=site_long, y=site_lat, shape = Site), size=2, pch = c(21,24), fill = "white",color = "black") +
  theme_bw()+
  ggtitle("Weekly Range (°C)")+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_fixed(ratio=1.1, xlim = c(-80,-65), ylim = c(30,45), expand = 0)+
  #coord_fixed(ratio=1.1, xlim = c(146.1, 147.6), ylim = c(-19.5, -18.2), expand = 0)+
  scale_fill_gradient2(low = 'dodgerblue', high = 'red', mid = 'white', midpoint = mean(ggplot.data$var))

