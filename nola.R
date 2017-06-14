# Created by Spencer Woody on 04 Jun 2017

library(ggplot2)
library(ggmap)
library(geosphere)
library(RColorBrewer) # display.brewer.all

# Read in data, change name of longitude column
nola <- read.csv("nola.csv", header = T)
names(nola)[2] <- "lon"

### --------------------------------------------------------------------------
### Map the stores
### --------------------------------------------------------------------------

# Set bounds of map
minlat <- min(nola$lat)
minlon <- min(nola$lon)
maxlat <- max(nola$lat)
maxlon <- max(nola$lon)

# Create background for map image
NOLAmap <- get_map(
	location = c(minlon, minlat, maxlon, maxlat), 
	zoom = 13, 
	col = "bw"
	)

# Plot locations of stores, color them by street name
storemap <- ggmap(NOLAmap) + 
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

# Save this map to a PDF
pdf("img/storemap.pdf")
storemap
dev.off()

### --------------------------------------------------------------------------
### Create distance matrices
### --------------------------------------------------------------------------

coords1 <- cbind(nola$lon[nola$street == 1], nola$lat[nola$street == 1])
coords2 <- cbind(nola$lon[nola$street == 2], nola$lat[nola$street == 2])
coords3 <- cbind(nola$lon[nola$street == 3], nola$lat[nola$street == 3])

distmat1 <- distm(coords1) / 1609.34 # haversine distance, converted to miles
distmat2 <- distm(coords2) / 1609.34
distmat3 <- distm(coords3) / 1609.34

distlist <- list(distmat1, distmat2, distmat3)

M52 <- function (dist, params) {
	# ----------------------------------------------------------------------
	# Calculate the Matern(5/2) covariance matrix from matrix of distances
	# ----------------------------------------------------------------------
	# INPUTS:
	# dist - distance matrix
	# params - vector of three hyperparamers
	#            [1] rho ....... relative distance
	#            [2] tau1.sq ... amplitude
	#            [3] tau2.sq ... nugget
	# ----------------------------------------------------------------------
	# OUTPUT (list): 
  	# M52 - Matern(5/2) covariance matrix
	# ----------------------------------------------------------------------
	
	rho     <- params[1]
	tau1.sq <- params[2]
	tau2.sq <- params[3]
	
	M52 <- { tau1.sq * 
		     (   1 + sqrt(5) * dist / rho + 5 / 3 * (dist / rho) ^ 2   ) *
			   exp(- sqrt(5) * dist / rho) } + tau2.sq * (dist == dist)
		   
	return(M52)
}



