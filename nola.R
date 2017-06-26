# Created by Spencer Woody on 04 Jun 2017

library(ggplot2)
library(ggmap)
library(geosphere)
library(Matrix)
library(mvtnorm)
library(BayesLogit)
library(RColorBrewer) # display.brewer.all()

source("nola_funs.R")

### --------------------------------------------------------------------------
### Prepare the data
### --------------------------------------------------------------------------

# Read in data, change name of longitude column
nola <- read.csv("nola.csv", header = T)
names(nola)[2] <- "lon"

# Reorder based on street 
nola <- nola[order(nola$street), ]
rownames(nola) <- 1:nrow(nola)

# Remove duplicate stores
nola <- nola[-which(duplicated(cbind(nola$lon, nola$lat, 
                                     nola$code, nola$street))), ]

# Number of stores on each street
n1 <- sum(nola$street == 1)
n2 <- sum(nola$street == 2)
n3 <- sum(nola$street == 3)

n.list <- list(n1 = n1, n2 = n2, n3 = n3)

# Total number of stores
n <- nrow(nola)

# Number of streets
m <- 3

# Entire design matrix
X <- as.matrix(nola[, -c(1:3, 12:15)])

# Number of fixed-effect covariates
p <- ncol(X)

# Separate design matrices for stores on each street
X.i <- list()
X.i[[1]] <- X[nola$street == 1, ]
X.i[[2]] <- X[nola$street == 2, ]
X.i[[3]] <- X[nola$street == 3, ]

# Response vectors
y1 <- nola$y1
y2 <- nola$y2
y3 <- nola$y3

# Response vector y3 separated for each store
y3.1 <- y3[nola$street == 1]
y3.2 <- y3[nola$street == 2]
y3.3 <- y3[nola$street == 3]

# For using PG data augmentation
eta3.3 <- y3.3 - 1 / 2
eta3.2 <- y3.2 - 1 / 2
eta3.1 <- y3.1 - 1 / 2

eta3 <- list(eta3.1, eta3.2, eta3.3)

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
    col = "bw")

# Plot locations of stores, color them by street name
storemap <- ggmap(NOLAmap) + 
	geom_point(
		data = nola, 
        aes(x = lon, y = lat, col = factor(street)), 
        alpha = 0.5) +
    scale_colour_brewer(
        name = "", 
        labels = c("St Charles Ave", "S Carrollton Ave", "St Claude Ave"), 
        palette = "Set1") +
	ggtitle("Location of all stores in dataset") +
    theme(
		axis.text = element_text(family="Avenir Next", size=11),
		axis.title = element_text(family="Avenir Next", size=11),
		plot.title = element_text(family="Avenir Next", size=16, hjust = 0.5),
		  legend.position="bottom")

# Save this map to a PDF
pdf("img/storemap.pdf", fonts = "Avenir Next")
storemap
dev.off()

### --------------------------------------------------------------------------
### Create distance matrices
### --------------------------------------------------------------------------

# Coordinates for each steets
coords1 <- cbind(nola$lon[nola$street == 1], nola$lat[nola$street == 1])
coords2 <- cbind(nola$lon[nola$street == 2], nola$lat[nola$street == 2])
coords3 <- cbind(nola$lon[nola$street == 3], nola$lat[nola$street == 3])

# Haversine distance, converted to miles
# https://en.wikipedia.org/wiki/Haversine_formula
distmat1 <- distm(coords1)  / 1609.34 
distmat2 <- distm(coords2)  / 1609.34
distmat3 <- distm(coords3)  / 1609.34

# Create a list
distlist <- list(distmat1, distmat2, distmat3)


### --------------------------------------------------------------------------
### Gibbs sampler
### --------------------------------------------------------------------------

# Number of iterations and burn-in
niter <- 200
nburn <- 100 + 1

# Create empty vectors / matrices for iteration of the Gibbs
b   <- rep(NA, niter + nburn)
tau <- rep(NA, niter + nburn)

beta <- matrix(nrow = niter + nburn, ncol = p)

omega <- matrix(nrow = niter + nburn, ncol = n)

f <- matrix(nrow = niter + nburn, ncol = n)

### Initialize stuff

# psi
psi.n <- mapply(FUN = rnorm, n = n.list)

# Latent variable
omega.n <- mapply(
	FUN = rpg,
	num = n.list,
	h = mapply(FUN = rep, rep(1, m), n.list),
	z = mapply(FUN = rep, rep(0, m), n.list))

# Make diagonal matrix
Omega.n <- mapply(
	FUN = Diagonal,
	n = n.list,
	x = omega.n)

omega[1, ] <- unlist(omega.n)

# beta
beta.n <- rnorm(p)

beta[1, ] <- beta.n

# f
f[1, ] <- unlist(psi.n) - X %*% beta.n

# b
b.n <- 0.01

b[1] <- b.n

# tau
tau.n <- 2

tau[1] <- tau.n

# Initialize covariance matrix
K.i <- lapply(distlist, FUN = M52, params = c(b.n, tau.n , 0))
K.i.inv <- lapply(K.i, solve)

# Standard deviation of log of proposal distribution for MH step
sd1 <- 0.5 # b
sd2 <- 0.5 # tau

for (iter in 2:(niter + nburn)) {
	
	# Reassign "current" hyperparameters
	b.c <- b[iter - 1]
	tau.c <- tau[iter - 1]
	
    # Update psi -----------------------------------------------------------
    psi.n.covmat <- lapply(mapply(FUN = '+', Omega.n, K.i.inv), solve)
	
	psi.n <- mapply(
		FUN = rmvnorm,
		n = rep(1, m),
		mean = mapply(FUN = "%*%", psi.n.covmat, eta3),
		sigma = lapply(psi.n.covmat, as.matrix))
	
	psi.n.unlist <- unlist(psi.n)
	
	f.n <- psi.n.unlist - X %*% beta.n
	
    # Update omegas --------------------------------------------------------
	
	omega.n <- mapply(
		FUN = rpg,
		num = n.list,
		h = mapply(FUN = rep, rep(1, m), n.list),
		z = psi.n)
	
	# Make diagonal matrix
	Omega.n <- mapply(
		FUN = Diagonal,
		n = n.list,
		x = omega.n)
	
    # Update beta ----------------------------------------------------------
	
	K.i.inv <- lapply(K.i, solve)
    K.inv <- bdiag(K.i.inv)
	
    beta.covmat <- solve(as.matrix(crossprod(X, K.inv) %*% X))
	beta.mean <- beta.covmat %*% crossprod(X, K.inv) %*% psi.n.unlist
	
    beta.n <- as.vector(rmvnorm(1, beta.mean, beta.covmat))
	
	Xbeta.n <- mapply(FUN = "%*%", X.i, rep(list(beta.n), m))
	
	#### Update hyperparameters --------------------------------------------
	# b
	
	# Draw from candidate distribution
	b.prime <- rlnorm(1, meanlog = log(b.c), sdlog = sd1)
	
	sumll.b.prime <- mapply(
		FUN = dmvnorm,
		x = psi.n,
		mean = Xbeta.n,
		sigma = lapply(distlist, FUN = M52, params = c(b.prime, tau.c, 0)),
		log = rep(TRUE, 3))
	
	sumll.b.c <- mapply( 
		FUN = dmvnorm,
		x = psi.n,
		mean = Xbeta.n,
		sigma = lapply(distlist, FUN = M52, params = c(b.c, tau.c, 0)),
		log = rep(TRUE, 3))
	
	d.b.c <- dlnorm(b.c, meanlog = log(b.prime), sdlog = sd1, log = TRUE)
	
	d.b.prime <- dlnorm(b.prime, meanlog = log(b.c), sdlog = sd1, log = TRUE)
	
	# Acceptance prob 1
	A1 <- exp( sum(unlist(sumll.b.prime)) + d.b.c - 
			   sum(unlist(sumll.b.c))     - d.b.prime  )
	
	U1 <- runif(1)
	
	if (U1 < A1) {
		b.n <- b.prime
	} else {
		b.n <- b.c
	}
	
	K.i <- lapply(distlist, FUN = M52, params = c(b.n, tau.c ,0))
	
	# tau
	tau.prime <- rlnorm(1, meanlog = log(tau.c), sdlog = sd2)
	
	sumll.tau.prime <- mapply(
		FUN = dmvnorm,
		x = psi.n,
		mean = Xbeta.n,
		sigma = lapply(distlist, FUN = M52, params = c(b.n, tau.prime, 0)),
		log = rep(TRUE, 3))
	
	sumll.tau.c <- mapply( 
		FUN = dmvnorm,
		x = psi.n,
		mean = Xbeta.n,
		sigma = lapply(distlist, FUN = M52, params = c(b.n, tau.c, 0)),
		log = rep(TRUE, 3))
	
	d.tau.c <- dlnorm(
		tau.c, 
		meanlog = log(tau.prime), 
		sdlog = sd1, 
		log = TRUE)
	
	d.tau.prime <- dlnorm(
		tau.prime, 
		meanlog = log(tau.c), 
		sdlog = sd1, 
		log = TRUE)
	
	# Acceptance prob 2
	A2 <- exp( sum(unlist(sumll.tau.prime)) + d.tau.c - 
			   sum(unlist(sumll.tau.c))     - d.tau.prime  )
	
	U2 <- runif(1)
	
	if (U2 < A2) {
		tau.n <- tau.prime
	} else {
		tau.n <- tau.c
	}
 	
	K.i <- lapply(distlist, FUN = M52, params = c(b.n, tau.n ,0))
	
	# Store new updates  ---------------------------------------------------
	beta[iter, ] <- beta.n
	
	omega[iter, ] <- unlist(omega.n)
	
	f[iter, ] <- f.n
	
	b[iter] <- b.n
	tau[iter] <- tau.n

	# Iteration update  ---------------------------------------------------
	if ((iter-1) %% 10 == 0) {
		print(
			sprintf(
				"Iteration %i out of %i...",
				(iter - 1), 
				(niter + nburn - 1)))
	}
}


# Remove burn-in

b <- b[-(1:nburn)]
tau <- tau[-(1:nburn)]

beta <- beta[-(1:nburn), ]

omega <- omega[-(1:nburn), ]

f <- f[-(1:nburn), ]


plot(b, type = "l")

plot(tau, type = "l")



