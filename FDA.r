###################FDA###############################
# compute mean estimate: overall and group means
mu.hat <- apply(x,2,mean)
mu.1.hat <- apply(x[1:154,],2,mean)
mu.2.hat <- apply(x[155:300,],2,mean)
# between class covariance
S.b <- (n1*(mu.1.hat-mu.hat)%*%t(mu.1.hat-mu.hat)+
       n2*(mu.2.hat-mu.hat)%*%t(mu.2.hat-mu.hat))/(n-1)
# within class covariance
S.w <- (t(x[1:n1,] - rep(1,n1)%*% t(mu.1.hat)) %*% (x[1:n1,] - rep(1,n1)%*% t(mu.1.hat)) +
       t(x[(n1+1):n,] - rep(1,n2)%*% t(mu.2.hat)) %*% (x[(n1+1):n,] - rep(1,n2)%*% t(mu.2.hat)))/(n-2)
# total variance
S.t <- t(x - rep(1,n)%*% t(mu.hat)) %*% (x - rep(1,n)%*% t(mu.hat))
#max(abs(S.t - S.b*(n-1) - S.w*(n-3)))=1.153921e-11  for test
# define relative matrix
library(MASS)
S <- ginv(S.w) %*% S.b
# svd decomp. of S
S.svd <- svd(S)
# eigen decomp. of S
S.eig <- eigen(S)
S.eig.values <- Re(S.eig$values)
S.eig.vectors <- Re(S.eig$vectors)
plot(S.eig.values)
#scores
x.fda<-x%*%S.svd$v
plot(x.fda[,1],x.fda[,2],col = c(g)) #scatterplot of the first score vs the second score
#plot(x.svd$u[,1]*x.svd$d[1],x.svd$u[,2]*x.svd$d[2],col = c(rep(1,n1/2),rep(2,n2/2),rep(3,n3/2)),pch=c(rep(1,n1/2),rep(2,n2/2),rep(3,n3/2))) #scatterplot of the first score vs the second score

x.test.fda<-x.test%*%S.svd$v
#x.test.mean <- apply(x.test,2,mean)
#x.test.demean <- x.test-rep(1,nrow(x.test))%*%t(x.test.mean)
#x.test.svd <- svd(x.test.demean)
plot(x.test.fda[,1],x.test.fda[,2],col = c(g.test)) #scatterplot of the first score vs the second score

########Regression on FDA scores################
###linear
library(MASS)
x <- cbind(x.fda[,1],x.fda[,2]) #two predictors
x <- cbind(rep(1,nrow(x)),x) # standardize predictors  and include intercept
B <- ginv(t(x)%*%x)%*%t(x)%*%y
plot(x.fda[,1],x.fda[,2],col= c(1,2)[g])
abline(a=(B[1,1]-B[1,2])/(B[3,2]-B[3,1]),b=(B[2,1]-B[2,2])/(B[3,2]-B[3,1]),col=4) #a:intercept,b:slope -0.001434887,  3.417291 boundary between 4 and 9
yhat<-x%*%B
g.hat<-which.max(yhat[1,])
for (i in 2:nrow(yhat)){
	g.hat<-append(g.hat,which.max(yhat[i,]))
}
t=table(g.hat,g)
trainerr=1-(t[1,1]+t[2,2])/n #error rate 0.1366667
x.test <- cbind(x.test.fda[,1],x.test.fda[,2]) #two predictors
x.test <- cbind(rep(1,nrow(x.test)),x.test) # standardize predictors  and include intercept
plot(x.test.fda[,1],x.test.fda[,2],col= c(1,2,3)[g.test])
abline(a=(B[1,1]-B[1,2])/(B[3,2]-B[3,1]),b=(B[2,1]-B[2,2])/(B[3,2]-B[3,1]),col=4) #boundary between 4 and 9
yhat.test<-x.test%*%B
g.hat.test<-which.max(yhat.test[1,])
for (i in 2:nrow(yhat.test)){
	g.hat.test<-append(g.hat.test,which.max(yhat.test[i,]))
}
t=table(g.hat.test,g.test)
testerr=1-(t[1,1]+t[2,2])/n #0.15

###LDA
x <- cbind(x.fda[,1],x.fda[,2])
pi1.hat <- mean(g==1) 
pi2.hat <- mean(g==2) 
mu.1.hat <- apply(x[g==1,],2,mean) 
mu.2.hat <- apply(x[g==2,],2,mean) 
S.w <- (t(x[g==1,] - rep(1,n1)%*% t(mu.1.hat)) %*% (x[g==1,] - rep(1,n1)%*% t(mu.1.hat)) +
        t(x[g==2,] - rep(1,n2)%*% t(mu.2.hat)) %*% (x[g==2,] - rep(1,n2)%*% t(mu.2.hat)))/(n-2)
delta1<-x%*%ginv(S.w)%*%mu.1.hat + as.vector(log(pi1.hat)-t(mu.1.hat)%*%ginv(S.w)%*%(mu.1.hat)/2)
delta2<-x%*%ginv(S.w)%*%mu.2.hat + as.vector(log(pi2.hat)-t(mu.2.hat)%*%ginv(S.w)%*%(mu.2.hat)/2)
yhat<-cbind(delta1,delta2)	#delta matrix is like yhat in linear regression
g.hat<-which.max(yhat[1,])
for (i in 2:nrow(yhat)){
	g.hat<-append(g.hat,which.max(yhat[i,]))
}
t=table(g.hat,g)
trainerr=1-(t[1,1]+t[2,2])/n #train error rate  0.1333333
slope1 <- ginv(S.w) %*% (mu.2.hat-mu.1.hat) #beta1
intercept1 <- log(pi2.hat/pi1.hat) - t(mu.2.hat+mu.1.hat) %*% ginv(S.w) %*% (mu.2.hat-mu.1.hat)/2 #a1
plot(x[,1],x[,2],col= c(1,2)[g])
abline(b=-slope1[1]/slope1[2],a=-intercept1/slope1[2],col=4) #boundary between 4 and 9

x.test <- cbind(x.test.fda[,1],x.test.fda[,2]) #two predictors
delta1.test<-x.test%*%ginv(S.w)%*%mu.1.hat + as.vector(log(pi1.hat)-t(mu.1.hat)%*%ginv(S.w)%*%(mu.1.hat)/2)
delta2.test<-x.test%*%ginv(S.w)%*%mu.2.hat + as.vector(log(pi2.hat)-t(mu.2.hat)%*%ginv(S.w)%*%(mu.2.hat)/2)
yhat.test<-cbind(delta1.test,delta2.test)	#delta matrix is like yhat in linear regression
g.hat.test<-which.max(yhat.test[1,])
for (i in 2:nrow(yhat.test)){
	g.hat.test<-append(g.hat.test,which.max(yhat.test[i,]))
}
t=table(g.hat.test,g.test)
testerr=1-(t[1,1]+t[2,2])/n #test error rate 0.2
plot(x.test[,1],x.test[,2],col= c(1,2)[g.test])
abline(b=-slope1[1]/slope1[2],a=-intercept1/slope1[2],col=4) 

###Logistic regression
library(glmnet)
x <- cbind(x.fda[,1],x.fda[,2]) #two predictors
ans.logistic <- glmnet(x,g, family=c("multinomial"),lambda=0)
beta1<-coef(ans.logistic)$'1'
beta2<-coef(ans.logistic)$'2'
predict(ans.logistic,type = "coefficients") #same output as previous line
a <- predict(ans.logistic,x,type = "link") 
p <- exp(a)/(exp(a)+1) #probability 
plot(x.fda[,1],x.fda[,2],col= c(1,2)[g])
abline(a=((beta1[1]-beta2[1])/(beta2[3]-beta1[3])), b=((beta1[2]-beta2[2])/(beta2[3]-beta1[3])),col=4)#boundary between 4 and 9
yhat<-cbind(p[,1,1],p[,2,1])	#delta matrix is like yhat in linear regression
g.hat<-which.max(yhat[1,])
for (i in 2:nrow(yhat)){
	g.hat<-append(g.hat,which.max(yhat[i,]))
}
t=table(g.hat,g)
trainerr=1-(t[1,1]+t[2,2])/n #train error rate 0.06666667

x.test <- cbind(x.test.fda[,1],x.test.fda[,2])
a.test <- predict(ans.logistic,x.test,type = "link") 
p.test <- exp(a.test)/(exp(a.test)+1)
plot(x.test.fda[,1],x.test.fda[,2],col= c(1,2)[g.test])
abline(a=((beta1[1]-beta2[1])/(beta2[3]-beta1[3])), b=((beta1[2]-beta2[2])/(beta2[3]-beta1[3])),col=4)#boundary between 4 and 9
yhat.test<-cbind(p.test[,1,1],p.test[,2,1])	#delta matrix is like yhat in linear regression
g.hat.test<-which.max(yhat.test[1,])
for (i in 2:nrow(yhat.test)){
	g.hat.test<-append(g.hat.test,which.max(yhat.test[i,]))
}
t=table(g.hat.test,g.test)
testerr=1-(t[1,1]+t[2,2])/n #0.1766667
