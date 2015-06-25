mylasso <- function(x,y,lambda){
    x <- as.matrix(x) # convert dataframe
    y <- as.vector(y)
    p <- ncol(x) # number of columns; number of predictors
	n <- nrow(x) #number of rows; number of observations
	
    tol <- 1e-5 #tolerance; accuracy
    niter <- 1000 #number of iteration
    iiter <- 1 #counter of iteration
    distance <- 1
	
	beta <- coefficients(lm(y~x)) #betas are initialized as solution from OLS rep(0,9)#
	cat("OLS coefficients are:", beta.old,'\n')
    beta <- as.vector(beta) #covert to be a column
	const <- rep(1,n) #for intersect
	x<-cbind(const,x)
    while(iiter <= niter & distance > tol) { 
        for(j in 1:(p+1)){
			cj <- t(x[,j])%*%(y-x[,-j]%*%beta[-j])
			if (cj>lambda) beta[j] <- cj-lambda
			else if (cj<(-lambda)) beta[j] <- cj+lambda
			else beta[j] <- 0 #soft threshold function
			cat("cj",cj,"lambda",lambda,"beta",beta,"\n")
        }
        distance <- norm(as.matrix(y-x%*%beta),"F") #L2 norm of the residue
		cat("distance",distance,"\n")
        iiter <- iiter + 1
    }
    return(beta)
}