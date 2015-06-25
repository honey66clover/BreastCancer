mylasso <- function(x,y,lambda){
    x <- as.matrix(x) # convert dataframe
    y <- as.vector(y)
    p <- ncol(x) # number of columns; number of predictors
	n <- nrow(x) #number of rows; number of observations
	
    tol <- 1e-5 #tolerance; accuracy
    niter <- 1000 #number of iteration
    iiter <- 1 #counter of iteration
    distance <- 1
	
	betaols <- coefficients(lm(y~x))  #solution of OLS
	beta <- betaols #betas are initialized as solution from OLS
	const <- rep(1,n) 
	x<-cbind(const,x) #add in one column for intersect
	x<- scale(x,TRUE,TRUE)
    while(iiter <= niter & distance > tol) { 
        for(j in 1:(p+1)){
			cj <- betaols[j]#t(x[,j])%*%(y-x[,-j]%*%beta[-j])
			#cat(t(x[,j])%*%y,t(x[,j])%*%x[,-j]%*%beta[-j],"\n")
			if (cj>lambda) beta[j] <- cj-lambda
			else if (cj<(-lambda)) beta[j] <- cj+lambda
			else beta[j] <- 0 #soft threshold function
			#cat("cj",cj,"lambda",lambda,"beta",beta,"\n")
        }
        distance <- norm(as.matrix(y-x%*%beta),"F") #L2 norm of the residue
		#cat("distance",distance,"\n")
        iiter <- iiter + 1
    }
    return(beta)
}

#plot
betas<-vector()
for (lambda in seq(0,1,0.1)){
	x <- read.table("C:\\Users\\honey66clover\\Desktop\\Weijuan\\Courses\\Statistics\\588 Data Mining\\Dan Yang\\prostate.data")
	# normalize predictor
	x <- cbind(scale(x[,1:8],TRUE,TRUE),x[,9:10])
	beta<-as.vector(mylasso(x[,1:8],x[,9],lambda))
	betas<-rbind(betas,beta)
}
lambda<-seq(0,1,0.01)
plot(lambda,betas[,1])
for (i in 2:9){
	lines(lambda,betas[,i])
}

library(lars)
lars.prostate <- lars(x=as.matrix(x[x[,10],1:8]),y=x[x[,10],9],type="lasso")
coef(lars.prostate,s=0,mode="lambda")	
lcavol     lweight         age        lbph         svi         lcp 
 0.67952814  0.26305307 -0.14146483  0.21014656  0.30520060 -0.28849277 
    gleason       pgg45 
-0.02130504  0.26695576 
beta.old 0.5998959 0.1858765 0.2808051 0.1107597 0.4003204 -0.5932076 -0.6132502 0.9168642 

	