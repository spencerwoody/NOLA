# Created by Spencer Woody on 02 Jun 2017

library(ggplot2)
library(ggmap)

nola <- read.csv("nola.csv", header = T)
names(nola)[2] <- "lon"


NOLAmap_coords <- geocode("NOLA")

NOLAmap_coords <- c(median(nola$lon), median(nola$lat))

minlat <- min(nola$lat)
minlon <- min(nola$lon)
maxlat <- max(nola$lat)
maxlon <- max(nola$lon)

NOLAmap <- get_map(location = c(minlon, minlat, maxlon, maxlat), zoom = 13, col = "bw")

mymap <- ggmap(NOLAmap) + 
geom_point(
	data = nola, 
	aes(x = lon, y = lat, col = factor(street)), 
	alpha = 0.5
	) +
scale_colour_brewer(
	name = "", 
	labels = c("St Charles Ave", "S Carrollton Ave", "St Claude Ave"), 
	palette = "Set1"
	) +
theme(legend.position="bottom")
mymap



s