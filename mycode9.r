setwd("C:/Users/dyang/Dropbox/588/myself/notes/lecture8")
setwd("D:/Dropbox/588/myself/notes/lecture9")
###########################################################################
# Three Gaussian example
#FDA
rm(list=ls())
p <- 128 #ncol	
n <- 300 #nrow
# mean for three classes p*1
theta1 <- c(3,0,rep(0,p-2)) #(3,0,0,0,0,0,0,0...)
theta2 <- c(0,-3,rep(0,p-2)) #(0,-3,0,0,0,0,0...)
theta3 <- c(-1,0,rep(0,p-2)) #(-1,0,0,0,0,0,0....)
# the mean matrix, n*p
mu <- rbind(matrix(rep(theta1,n/3),nrow=n/3,byrow=TRUE),
            matrix(rep(theta2,n/3),nrow=n/3,byrow=TRUE),
            matrix(rep(theta3,n/3),nrow=n/3,byrow=TRUE)) 
# observed data, with mena=mu and cov=diag(.5,2,...,2), n*p
x <- mu + matrix(rnorm(n*p),n) %*% diag(c(.5,2,rep(2,p-2)))

# oracle projection; the first two coordinates
pdf(file="LDA_oracle.pdf",height=8,width=8)
plot(x[,1],x[,2],col=rep(2:4,rep(n/3,3)),pch=rep(2:4,rep(n/3,3)))
points(mu[,1],mu[,2],pch=19,cex=2) #add points to plot
dev.off()

# PCA by SVD
x.svd <- svd(x)
# PCA projection
pdf(file="LDA_PCA.pdf",height=8,width=8)
#scores
plot(x.svd$u[,1]*x.svd$d[1],x.svd$u[,2]*x.svd$d[2],col=rep(2:4,rep(n/3,3)),pch=rep(2:4,rep(n/3,3)))
dev.off()
# PCA loadings
pdf(file="LDA_PCA_loading.pdf",height=8,width=8)
par(mfrow=c(2,2))
plot(x.svd$v[,1])
plot(x.svd$v[,2])
plot(x.svd$v[,3])
plot(x.svd$v[,4])
dev.off()
# PCA spectral, sigular values
pdf(file="LDA_PCA_spectral.pdf",height=8,width=8)
plot(x.svd$d)
dev.off()

#--------------------------------------------------------------------------
# compute mean estimate: overall and group means
mu.hat <- apply(x,2,mean)
mu.1.hat <- apply(x[1:(n/3),],2,mean)
mu.2.hat <- apply(x[(n/3+1):(2*n/3),],2,mean)
mu.3.hat <- apply(x[(2*n/3+1):n,],2,mean)

# between class covariance
S.b <- ((n/3)*(mu.1.hat-mu.hat)%*%t(mu.1.hat-mu.hat)+
       (n/3)*(mu.2.hat-mu.hat)%*%t(mu.2.hat-mu.hat)+
       (n/3)*(mu.3.hat-mu.hat)%*%t(mu.3.hat-mu.hat))/(n-1)
# within class covariance
S.w <- (t(x[1:(n/3),] - rep(1,n/3)%*% t(mu.1.hat)) %*% (x[1:(n/3),] - rep(1,n/3)%*% t(mu.1.hat)) +
       t(x[(n/3+1):(2*n/3),] - rep(1,n/3)%*% t(mu.2.hat)) %*% (x[(n/3+1):(2*n/3),] - rep(1,n/3)%*% t(mu.2.hat)) +
       t(x[(2*n/3+1):n,] - rep(1,n/3)%*% t(mu.3.hat)) %*% (x[(2*n/3+1):n,] - rep(1,n/3)%*% t(mu.3.hat)))/(n-3)
# total variance
S.t <- t(x - rep(1,n)%*% t(mu.hat)) %*% (x - rep(1,n)%*% t(mu.hat))
# relation
# S.t - S.b*(n-1) - S.w*(n-K) = 0
max(abs(S.t - S.b*(n-1) - S.w*(n-3)))

# define relative matrix
S <- solve(S.w) %*% S.b
# eigen decomp. of S
S.eig <- svd(S)

# FDA projection
pdf(file="LDA_LDA.pdf",height=8,width=8)
plot(x %*% S.eig$v[,1], x%*% S.eig$v[,2], col=rep(2:4,rep(n/3,3)),pch=rep(2:4,rep(n/3,3)))
dev.off()

# FDA loadings
pdf(file="LDA_LDA_loading.pdf",height=8,width=8)
par(mfrow=c(2,2))
plot(S.eig$v[,1])
plot(S.eig$v[,2])
plot(S.eig$v[,3])
plot(S.eig$v[,4])
dev.off()

# FDA spectral
pdf(file="LDA_LDA_spectral.pdf",height=8,width=8)
plot(S.eig$d)
dev.off()

# three projections together
pdf(file="LDA.pdf",height=3,width=8)
par(mfrow=c(1,3))
plot(x[,1],x[,2],col=rep(2:4,rep(n/3,3)),pch=rep(2:4,rep(n/3,3)),main="oracle projection")
points(mu[,1],mu[,2],pch=19,cex=2)
plot(x.svd$u[,1]*x.svd$d[1],x.svd$u[,2]*x.svd$d[2],col=rep(2:4,rep(n/3,3)),pch=rep(2:4,rep(n/3,3)),main="PCA projection")
plot(x %*% S.eig$v[,1], x%*% S.eig$v[,2], col=rep(2:4,rep(n/3,3)),pch=rep(2:4,rep(n/3,3)),main = "LDA projection")
dev.off()

###########################################################################
# Three Gaussian example
# LDA and DR
rm(list=ls())
p <- 2
n <- 450
# mean for three classes
theta1 <- c(0,0)
theta2 <- c(-3,2)
theta3 <- c(-1,-3)
# the mean matrix
mu <- rbind(matrix(rep(theta1,n/3),nrow=n/3,byrow=TRUE),
            matrix(rep(theta2,n/3),nrow=n/3,byrow=TRUE),
            matrix(rep(theta3,n/3),nrow=n/3,byrow=TRUE))
# observed data
Sigma <- matrix(c(4, -1,-1,.6),2)
Sigma.eig <- svd(Sigma)
Sigma.half <- Sigma.eig$v %*% diag(sqrt(Sigma.eig$d)) %*% t(Sigma.eig$v)
#Sigma.half = diag(rep(1,2))
set.seed(1)
x <- mu + matrix(rnorm(n*p),n) %*% Sigma.half
x.test <- mu + matrix(rnorm(n*p),n) %*% Sigma.half
g <- c(rep(1,n/3),rep(2,n/3),rep(3,n/3))

#--------------------------------------------------------------------------
# LDA with two dim predictors
pi1.hat <- mean(g==1)
pi2.hat <- mean(g==2)
pi3.hat <- mean(g==3)

mu.1.hat <- apply(x[1:(n/3),],2,mean)
mu.2.hat <- apply(x[(n/3+1):(2*n/3),],2,mean)
mu.3.hat <- apply(x[(2*n/3+1):n,],2,mean)

S.w <- (t(x[1:(n/3),] - rep(1,n/3)%*% t(mu.1.hat)) %*% (x[1:(n/3),] - rep(1,n/3)%*% t(mu.1.hat)) +
       t(x[(n/3+1):(2*n/3),] - rep(1,n/3)%*% t(mu.2.hat)) %*% (x[(n/3+1):(2*n/3),] - rep(1,n/3)%*% t(mu.2.hat)) +
       t(x[(2*n/3+1):n,] - rep(1,n/3)%*% t(mu.3.hat)) %*% (x[(2*n/3+1):n,] - rep(1,n/3)%*% t(mu.3.hat)))/(n-3)


slope1 <- solve(S.w) %*% (mu.2.hat-mu.1.hat)
intercept1 <- log(pi2.hat/pi1.hat) - t(mu.2.hat+mu.1.hat) %*% solve(S.w) %*% (mu.2.hat-mu.1.hat)/2
slope1
intercept1

slope2 <- solve(S.w) %*% (mu.3.hat-mu.1.hat)
intercept2 <- log(pi3.hat/pi1.hat) - t(mu.3.hat+mu.1.hat) %*% solve(S.w) %*% (mu.3.hat-mu.1.hat)/2
slope2
intercept2

slope3 <- solve(S.w) %*% (mu.3.hat-mu.2.hat)
intercept3 <- log(pi3.hat/pi2.hat) - t(mu.3.hat+mu.2.hat) %*% solve(S.w) %*% (mu.3.hat-mu.2.hat)/2
slope3
intercept3

#windows()
pdf(file="Identity_LDA.pdf")
plot(x[,1],x[,2],col= c(2:4)[g])
abline(b=-slope1[1]/slope1[2],a=-intercept1/slope1[2],col=5)
abline(b=-slope2[1]/slope2[2],a=-intercept2/slope2[2],col=5)
abline(b=-slope3[1]/slope3[2],a=-intercept3/slope3[2],col=5)
dev.off()

delta.hat <-
cbind(log(pi1.hat)-diag((x.test-rep(1,n)%*% t(mu.1.hat))%*% solve(S.w)%*%t(x.test-rep(1,n)%*% t(mu.1.hat))/2),
      log(pi2.hat)-diag((x.test-rep(1,n)%*% t(mu.2.hat))%*% solve(S.w)%*%t(x.test-rep(1,n)%*% t(mu.2.hat))/2),
      log(pi3.hat)-diag((x.test-rep(1,n)%*% t(mu.3.hat))%*% solve(S.w)%*%t(x.test-rep(1,n)%*% t(mu.3.hat))/2))
g.hat <- apply(delta.hat,1,which.max)
mean(g.hat != g)
# 0.05777778

#--------------------------------------------------------------------------
# LDA with one PC
pi1.hat <- mean(g==1)
pi2.hat <- mean(g==2)
pi3.hat <- mean(g==3)

mu.hat <- apply(x,2,mean)
x.svd <- svd(x-rep(1,n)%*%t(mu.hat))
z.pca <- x %*% x.svd$v[,1]

#windows()
pdf(file="Identity_PCA_direction.pdf")
plot(x[,1],x[,2],col= c(2:4)[g])
abline(b = x.svd$v[2,1]/x.svd$v[1,1], a = mu.hat[2] - mu.hat[1]*x.svd$v[2,1]/x.svd$v[1,1],col=5,lwd=3,lty=2)
abline(b = x.svd$v[2,2]/x.svd$v[1,2], a = mu.hat[2] - mu.hat[1]*x.svd$v[2,2]/x.svd$v[1,2],col=6,lwd=2,lty=3)
dev.off()

mu.1.hat <- apply(z.pca[1:(n/3),,drop=FALSE],2,mean)
mu.2.hat <- apply(z.pca[(n/3+1):(2*n/3),,drop=FALSE],2,mean)
mu.3.hat <- apply(z.pca[(2*n/3+1):n,,drop=FALSE],2,mean)

S.w <- (t(z.pca[1:(n/3),,drop=FALSE] - rep(1,n/3)%*% t(mu.1.hat)) %*% (z.pca[1:(n/3),,drop=FALSE] - rep(1,n/3)%*% t(mu.1.hat)) +
       t(z.pca[(n/3+1):(2*n/3),,drop=FALSE] - rep(1,n/3)%*% t(mu.2.hat)) %*% (z.pca[(n/3+1):(2*n/3),,drop=FALSE] - rep(1,n/3)%*% t(mu.2.hat)) +
       t(z.pca[(2*n/3+1):n,,drop=FALSE] - rep(1,n/3)%*% t(mu.3.hat)) %*% (z.pca[(2*n/3+1):n,,drop=FALSE] - rep(1,n/3)%*% t(mu.3.hat)))/(n-3)


slope1 <- solve(S.w) %*% (mu.2.hat-mu.1.hat)
intercept1 <- log(pi2.hat/pi1.hat) - t(mu.2.hat+mu.1.hat) %*% solve(S.w) %*% (mu.2.hat-mu.1.hat)/2
slope1
intercept1

slope2 <- solve(S.w) %*% (mu.3.hat-mu.1.hat)
intercept2 <- log(pi3.hat/pi1.hat) - t(mu.3.hat+mu.1.hat) %*% solve(S.w) %*% (mu.3.hat-mu.1.hat)/2
slope2
intercept2

slope3 <- solve(S.w) %*% (mu.3.hat-mu.2.hat)
intercept3 <- log(pi3.hat/pi2.hat) - t(mu.3.hat+mu.2.hat) %*% solve(S.w) %*% (mu.3.hat-mu.2.hat)/2
slope3
intercept3

#windows()
pdf(file="Identity_LDA_on_PC.pdf")
plot(z.pca,g,col= c(2:4)[g])
abline(v=-intercept1/slope1[1],col=5)
abline(v=-intercept2/slope2[1],col=5)
dev.off()

z.pca.test <- x.test %*% x.svd$v[,1]
delta.hat <-
cbind(log(pi1.hat)-diag((z.pca.test-rep(1,n)%*% t(mu.1.hat))%*% solve(S.w)%*%t(z.pca.test-rep(1,n)%*% t(mu.1.hat))/2),
      log(pi2.hat)-diag((z.pca.test-rep(1,n)%*% t(mu.2.hat))%*% solve(S.w)%*%t(z.pca.test-rep(1,n)%*% t(mu.2.hat))/2),
      log(pi3.hat)-diag((z.pca.test-rep(1,n)%*% t(mu.3.hat))%*% solve(S.w)%*%t(z.pca.test-rep(1,n)%*% t(mu.3.hat))/2))
g.hat <- apply(delta.hat,1,which.max)
mean(g.hat != g)
# [1] 0.1288889

#--------------------------------------------------------------------------
# LDA with one FD
pi1.hat <- mean(g==1)
pi2.hat <- mean(g==2)
pi3.hat <- mean(g==3)

mu.hat <- apply(x,2,mean)
mu.1.hat <- apply(x[1:(n/3),],2,mean)
mu.2.hat <- apply(x[(n/3+1):(2*n/3),],2,mean)
mu.3.hat <- apply(x[(2*n/3+1):n,],2,mean)

# between class covariance
S.b <- ((n/3)*(mu.1.hat-mu.hat)%*%t(mu.1.hat-mu.hat)+
       (n/3)*(mu.2.hat-mu.hat)%*%t(mu.2.hat-mu.hat)+
       (n/3)*(mu.3.hat-mu.hat)%*%t(mu.3.hat-mu.hat))/(n-1)
# within class covariance
S.w <- (t(x[1:(n/3),] - rep(1,n/3)%*% t(mu.1.hat)) %*% (x[1:(n/3),] - rep(1,n/3)%*% t(mu.1.hat)) +
       t(x[(n/3+1):(2*n/3),] - rep(1,n/3)%*% t(mu.2.hat)) %*% (x[(n/3+1):(2*n/3),] - rep(1,n/3)%*% t(mu.2.hat)) +
       t(x[(2*n/3+1):n,] - rep(1,n/3)%*% t(mu.3.hat)) %*% (x[(2*n/3+1):n,] - rep(1,n/3)%*% t(mu.3.hat)))/(n-3)

# define relative matrix
S <- solve(S.w) %*% S.b
# eigen decomp. of S
S.eig <- svd(S)

#windows()
pdf(file="Identity_FDA_direction.pdf")
plot(x[,1],x[,2],col= c(2:4)[g])
abline(b = S.eig$v[2,1]/S.eig$v[1,1], a = mu.hat[2] - mu.hat[1]*S.eig$v[2,1]/S.eig$v[1,1],col=5,lwd=3,lty=2)
abline(b = S.eig$v[2,2]/S.eig$v[1,2], a = mu.hat[2] - mu.hat[1]*S.eig$v[2,2]/S.eig$v[1,2],col=6,lwd=2,lty=3)
dev.off()

z.lda <- x %*% S.eig$v[,1]


mu.1.hat <- apply(z.lda[1:(n/3),,drop=FALSE],2,mean)
mu.2.hat <- apply(z.lda[(n/3+1):(2*n/3),,drop=FALSE],2,mean)
mu.3.hat <- apply(z.lda[(2*n/3+1):n,,drop=FALSE],2,mean)

S.w <- (t(z.lda[1:(n/3),,drop=FALSE] - rep(1,n/3)%*% t(mu.1.hat)) %*% (z.lda[1:(n/3),,drop=FALSE] - rep(1,n/3)%*% t(mu.1.hat)) +
       t(z.lda[(n/3+1):(2*n/3),,drop=FALSE] - rep(1,n/3)%*% t(mu.2.hat)) %*% (z.lda[(n/3+1):(2*n/3),,drop=FALSE] - rep(1,n/3)%*% t(mu.2.hat)) +
       t(z.lda[(2*n/3+1):n,,drop=FALSE] - rep(1,n/3)%*% t(mu.3.hat)) %*% (z.lda[(2*n/3+1):n,,drop=FALSE] - rep(1,n/3)%*% t(mu.3.hat)))/(n-3)


slope1 <- solve(S.w) %*% (mu.2.hat-mu.1.hat)
intercept1 <- log(pi2.hat/pi1.hat) - t(mu.2.hat+mu.1.hat) %*% solve(S.w) %*% (mu.2.hat-mu.1.hat)/2
slope1
intercept1

slope2 <- solve(S.w) %*% (mu.3.hat-mu.1.hat)
intercept2 <- log(pi3.hat/pi1.hat) - t(mu.3.hat+mu.1.hat) %*% solve(S.w) %*% (mu.3.hat-mu.1.hat)/2
slope2
intercept2

slope3 <- solve(S.w) %*% (mu.3.hat-mu.2.hat)
intercept3 <- log(pi3.hat/pi2.hat) - t(mu.3.hat+mu.2.hat) %*% solve(S.w) %*% (mu.3.hat-mu.2.hat)/2
slope3
intercept3

#windows()
pdf(file="Identity_LDA_on_FDA.pdf")
plot(z.lda,g,col= c(2:4)[g])
abline(v=-intercept1/slope1[1],col=5)
abline(v=-intercept2/slope2[1],col=5)
dev.off()

z.lda.test <- x.test %*% S.eig$v[,1]
delta.hat <-
cbind(log(pi1.hat)-diag((z.lda.test-rep(1,n)%*% t(mu.1.hat))%*% solve(S.w)%*%t(z.lda.test-rep(1,n)%*% t(mu.1.hat))/2),
      log(pi2.hat)-diag((z.lda.test-rep(1,n)%*% t(mu.2.hat))%*% solve(S.w)%*%t(z.lda.test-rep(1,n)%*% t(mu.2.hat))/2),
      log(pi3.hat)-diag((z.lda.test-rep(1,n)%*% t(mu.3.hat))%*% solve(S.w)%*%t(z.lda.test-rep(1,n)%*% t(mu.3.hat))/2))
g.hat <- apply(delta.hat,1,which.max)
mean(g.hat != g)
# [1] 0.1222222


###########################################################################
# South African Heart Disease
# logistic regression

rm(list=ls())
x <- read.table("C:/Users/honey66clover/Desktop/Weijuan/Courses/Statistics/588 Data Mining/Dan Yang/SAheart.data", 	sep=",",head=T,row.names=1)
x[,"famhist"] <- as.numeric(x[,"famhist"])
x <- as.matrix(x)
y <- x[,10]
x <- scale(x[,-c(4,6,10)])
library(glmnet)
ans.logistic <- glmnet(x,y, family=c("binomial"),lambda=0)
coef(ans.logistic)
predict(ans.logistic,type = "coefficients") #same output as previous line
a <- predict(ans.logistic,x,type = "link") 
p <- exp(a)/(exp(a)+1) #probability to be in class 1
predict(ans.logistic,x,type = "response") #same output as previous line
classfication <- table(p>.5,y)
1-sum(diag(classfication))/sum(classfication)
#[1] 0.2705628

ans.logistic.regularized <- glmnet(x, y, family=c("binomial"))
plot(ans.logistic.regularized)

library(glmpath)
ans.logistic.regularized2 <- glmpath(x, y, family=binomial)
plot(ans.logistic.regularized2)








