# Created by Spencer Woody on 25 Jun 2017

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
            exp(- sqrt(5) * dist / rho) } + tau2.sq * (dist == 0)
    
    return(M52)
}

