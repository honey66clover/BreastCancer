mylasso <- function(x,y,lambda){
    x <- as.matrix(x) # convert dataframe
    y <- as.vector(y)
    p <- ncol(x) # number of columns; number of predictors
    tol <- 1e-5 #tolerance; accuracy
    niter <- 100 #number of iteration
    iiter <- 1 #counter of iteration
    distance <- 1
	beta.ols <- coefficients(lm(y~x))  #solution of OLS, including intercept
    beta.old <- beta.ols
    beta.new <- beta.old #initiation
    while(iiter <= niter & distance > tol) {
        for(j in 1:(p+1)){
			cj <- beta.ols[j]
			if (cj>lambda) beta.new[j] <- cj-lambda
			else if (cj<(-lambda)) beta.new[j] <- cj+lambda
			else beta.new[j] <- 0 #soft threshold function
        }
        distance <- norm(as.matrix(beta.old-beta.new),"F") #L2 norm of the residue
        iiter <- iiter + 1
        beta.old <- beta.new
    }
    return(beta.new)
}

cj <- sum(x[,j]*(y-x%*%beta.new))
			if (cj>lambda) beta.new[j] <- cj-lambda
			else if (cj<(-lambda)) beta.new[j] <- cj+lambda
			else beta.new[j] <- 0 #soft threshold function