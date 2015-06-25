rm(list=ls(all=TRUE))
x <- read.table("C:\\Users\\honey66clover\\Desktop\\Weijuan\\Courses\\Statistics\\588 Data Mining\\Dan Yang\\project\\breast cancer\\test.csv",sep=",")
x.test <- read.table("C:\\Users\\honey66clover\\Desktop\\Weijuan\\Courses\\Statistics\\588 Data Mining\\Dan Yang\\project\\breast cancer\\train.csv",sep=",")
x <- as.matrix(x)
x.test <- as.matrix(x.test)
n<-nrow(x) #total train sample
n1<-203 #class 1
n2<-66 #class2
g<-x[,1]
x<-x[,2:31]
y1<-cbind(rep(1,203),rep(0,203)) #Bs
y2<-cbind(rep(0,66),rep(1,66)) #Ms
y<-rbind(y1,y2)
g.test<-x.test[,1]
x.test<-x.test[,2:31]
y1.test<-cbind(rep(1,154),rep(0,154)) #Bs
y2.test<-cbind(rep(0,146),rep(1,146)) #Ms
y.test<-rbind(y1.test,y2.test)
##########PCA################33
x.mean <- apply(x,2,mean)
x.demean <- x-rep(1,nrow(x))%*%t(x.mean)
x.svd <- svd(x.demean)
# singular values
plot(x.svd$d)
#scores
x.pca=x%*%x.svd$v
plot(x.pca[,1],x.pca[,2],col = c(g)) #scatterplot of the first score vs the second score

#x.test.mean <- apply(x.test,2,mean)
#x.test.demean <- x.test-rep(1,nrow(x.test))%*%t(x.test.mean)
#x.test.svd <- svd(x.test.demean)
x.test.pca=x.test%*%x.svd$v
plot(x.test.pca[,1],x.test.pca[,2],col = c(g.test)) #scatterplot of the first score vs the second score

########Regression on PCA scores################
###linear
library(MASS)
x <- cbind(x.pca[,1],x.pca[,2]) #two predictors
x <- cbind(rep(1,nrow(x)),x) # standardize predictors  and include intercept
B <- ginv(t(x)%*%x)%*%t(x)%*%y
plot(x.pca[,1],x.pca[,2],col= c(1,2)[g])
abline(a=(B[1,1]-B[1,2])/(B[3,2]-B[3,1]),b=(B[2,1]-B[2,2])/(B[3,2]-B[3,1]),col=4) #a:intercept,b:slope -0.001434887,  3.417291 boundary between 4 and 9
yhat<-x%*%B
g.hat<-which.max(yhat[1,])
for (i in 2:nrow(yhat)){
	g.hat<-append(g.hat,which.max(yhat[i,]))
}
t=table(g.hat,g)
trainerr=1-(t[1,1]+t[2,2])/n #error rate 0.08178439
x.test <- cbind(x.test.pca[,1],x.test.pca[,2]) #two predictors
x.test <- cbind(rep(1,nrow(x.test)),x.test) # standardize predictors  and include intercept
plot(x.test.pca[,1],x.test.pca[,2],col= c(1,2,3)[g.test])
abline(a=(B[1,1]-B[1,2])/(B[3,2]-B[3,1]),b=(B[2,1]-B[2,2])/(B[3,2]-B[3,1]),col=4) #boundary between 4 and 9
yhat.test<-x.test%*%B
g.hat.test<-which.max(yhat.test[1,])
for (i in 2:nrow(yhat.test)){
	g.hat.test<-append(g.hat.test,which.max(yhat.test[i,]))
}
t=table(g.hat.test,g.test)
testerr=1-(t[1,1]+t[2,2])/n #0.1524164

###LDA
x <- cbind(x.pca[,1],x.pca[,2])
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

x.test <- cbind(x.test.pca[,1],x.test.pca[,2]) #two predictors
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
x <- cbind(x.pca[,1],x.pca[,2]) #two predictors
ans.logistic <- glmnet(x,g, family=c("multinomial"),lambda=0)
beta1<-coef(ans.logistic)$'1'
beta2<-coef(ans.logistic)$'2'
predict(ans.logistic,type = "coefficients") #same output as previous line
a <- predict(ans.logistic,x,type = "link") 
p <- exp(a)/(exp(a)+1) #probability 
plot(x.pca[,1],x.pca[,2],col= c(1,2)[g])
abline(a=((beta1[1]-beta2[1])/(beta2[3]-beta1[3])), b=((beta1[2]-beta2[2])/(beta2[3]-beta1[3])),col=4)#boundary between 4 and 9
yhat<-cbind(p[,1,1],p[,2,1])	#delta matrix is like yhat in linear regression
g.hat<-which.max(yhat[1,])
for (i in 2:nrow(yhat)){
	g.hat<-append(g.hat,which.max(yhat[i,]))
}
t=table(g.hat,g)
trainerr=1-(t[1,1]+t[2,2])/n #train error rate 0.06666667

x.test <- cbind(x.test.pca[,1],x.test.pca[,2])
a.test <- predict(ans.logistic,x.test,type = "link") 
p.test <- exp(a.test)/(exp(a.test)+1)
plot(x.test.pca[,1],x.test.pca[,2],col= c(1,2)[g.test])
abline(a=((beta1[1]-beta2[1])/(beta2[3]-beta1[3])), b=((beta1[2]-beta2[2])/(beta2[3]-beta1[3])),col=4)#boundary between 4 and 9
yhat.test<-cbind(p.test[,1,1],p.test[,2,1])	#delta matrix is like yhat in linear regression
g.hat.test<-which.max(yhat.test[1,])
for (i in 2:nrow(yhat.test)){
	g.hat.test<-append(g.hat.test,which.max(yhat.test[i,]))
}
t=table(g.hat.test,g.test)
testerr=1-(t[1,1]+t[2,2])/n #0.1766667


