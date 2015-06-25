# Example of a function
myfunction <- function(x){
    y <- x^2 + 1
    return(y)
}

###########################################################################

# write your own ridge
myridge <- function(x,y,lambda){
    x <- as.matrix(x) # convert dataframe
    y <- as.vector(y)
    p <- ncol(x) # number of columns; number of predictors
    beta <- solve(t(x)%*%x + lambda * diag(p)) %*% t(x) %*% y
    # t(): transpose
    # diag(): diagonal matrix
    # matrix multiplication %*%
    # solve(): inverse
    return(beta)
}
x <- matrix(rnorm(10),5,2) # rnorm(): generate norm random variables
y <- rnorm(5)
x <- scale(x)
y <- y - mean(y)
beta.hat <- myridge(x,y,0); beta.hat
beta.hat <- myridge(x,y,1); beta.hat
beta.hat <- myridge(x,y,100); beta.hat

###########################################################################

# Coordinate Descent for LASSO in R: skeleton
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

###########################################################################

# remove all data
rm(list=ls(all=TRUE))
# load data
x <- read.table("C:\\Users\\honey66clover\\Desktop\\Weijuan\\Courses\\Statistics\\588 Data Mining\\Dan Yang\\prostate.data")
# normalize predictor
x <- cbind(scale(x[,1:8],TRUE,TRUE),x[,9:10])
#x <- as.matrix(x)
mylasso(x=as.matrix(x[x[,10],1:8]),y=x[x[,10],9],0)

#OLS solutions:
# (Intercept)  x[, 1:8]lcavol x[, 1:8]lweight     x[, 1:8]age    x[, 1:8]lbph 
#     2.47838688      0.66514667      0.26648026     -0.15819522      0.14031117 
#    x[, 1:8]svi     x[, 1:8]lcp x[, 1:8]gleason   x[, 1:8]pgg45 
#     0.31532888     -0.14828568      0.03554917      0.12571982 

###########################################################################

# LARS path

library(lars)

lars.prostate <- lars(x=as.matrix(x[x[,10],1:8]),y=x[x[,10],9],type="lar")
plot(lars.prostate)

lars.prostate <- lars(x=as.matrix(x[x[,10],1:8]),y=x[x[,10],9],type="forward.stagewise")
#> lars.prostate <- lars(x=as.matrix(x[x[,10],1:8]),y=x[x[,10],9],type="lasso")
#> coef(lars.prostate)
#         lcavol    lweight         age        lbph        svi        lcp
# [1,] 0.0000000 0.00000000  0.00000000 0.000000000 0.00000000  0.0000000
# [2,] 0.4059190 0.00000000  0.00000000 0.000000000 0.00000000  0.0000000
# [3,] 0.4756800 0.06611066  0.00000000 0.000000000 0.00000000  0.0000000
# [4,] 0.5321290 0.16878914  0.00000000 0.000000000 0.09162986  0.0000000
# [5,] 0.5332414 0.16987293  0.00000000 0.003547863 0.09488749  0.0000000
# [6,] 0.5498506 0.22078661  0.00000000 0.142241597 0.19744053  0.0000000
# [7,] 0.5562480 0.23122248 -0.03184273 0.162425741 0.20525620  0.0000000
# [8,] 0.6633830 0.26157445 -0.13272669 0.204414190 0.29602406 -0.2560034
# [9,] 0.6795281 0.26305307 -0.14146483 0.210146557 0.30520060 -0.2884928
#          gleason      pgg45
# [1,]  0.00000000 0.00000000
# [2,]  0.00000000 0.00000000
# [3,]  0.00000000 0.00000000
# [4,]  0.00000000 0.00000000
# [5,]  0.00000000 0.00000000
# [6,]  0.00000000 0.08670586
# [7,]  0.00000000 0.10265343
# [8,]  0.00000000 0.23602025
# [9,] -0.02130504 0.26695576












