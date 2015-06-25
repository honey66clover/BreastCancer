library(EbayesThresh) # to introduce soft threshold function, otherwise use three if statements 
mylasso <- function(x,y,lambda){
    x <- as.matrix(x) # convert dataframe
    y <- as.vector(y)
    p <- ncol(x) # number of columns; number of predictors
    tol <- 1e-5 #tolerance; accuracy
    niter <- 100 #number of iteration
    iiter <- 1 #counter of iteration
    dist <- 1
	beta.old <- solve(t(x)%*%x)%*%t(x)%*%y  #solution of OLS
    beta.new <- beta.old #initiation
    while(iiter <= niter & dist > tol) {
        for(j in 1:p){
			beta.new[j] <- threshld(t(x[,j])%*%(y-x[,-j]%*%beta.new[-j]),lambda, hard = FALSE)
        }
        dist <- norm(as.matrix(beta.old-beta.new),"F")		#L2 norm of the residue
        iiter <- iiter + 1
        beta.old <- beta.new
    }
    return(beta.new)
}