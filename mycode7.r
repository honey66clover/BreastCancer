setwd("C:/Users/honey66clover/Desktop/Weijuan/Courses/Statistics/588 Data Mining/Dan Yang")
rm(list=ls(all=TRUE))
# read data
# separation ","
x <- read.table("C:\\Users\\honey66clover\\Desktop\\Weijuan\\Courses\\Statistics\\588 Data Mining\\Dan Yang\\train3.data",sep=",")
x <- as.matrix(x)
dim(x)
# [1] 658 256, 658 images, each has 256 pixels


# plot the first image of 3
windows()
image(matrix(x[1,],16))  #16 pixels by 16 pixels image
# gray scale
image(matrix(x[1,],16),col=gray((32:0)/32))
# upside down?
# convert to a 3D array
y <- array(x,c(nrow(x),16,16))
# correct the upside down
y <- y[,,16:1]
windows()
image(y[1,,],col=gray((32:0)/32))
# looks all right
# change the margin!!!
windows()
par(mfrow = c(4,4))
for(i in 1:16){
  image(y[i,,],col=gray((32:0)/32))
}
# convert back to 2D matrix
x <- matrix(y,ncol=256)


# compute the mean image
x.mean <- apply(x,2,mean)
# subtract the mean image !!!!!!!!!! see digit 4 below
x.demean <- x-rep(1,nrow(x))%*%t(x.mean)
# SVD of the centered matrix
x.svd <- svd(x.demean)
# plot mean, the first and second eigen-image
windows(height=3,width=8)
par(mfrow = c(1,3))
image(matrix(x.mean,16),col=gray((32:0)/32))
image(matrix(-x.svd$v[,1],16),col=gray((32:0)/32))
image(matrix(x.svd$v[,2],16),col=gray((32:0)/32))

###########################################################################
# analyze 4 and 9 together
# load 4
rm(list=ls(all=TRUE))
x <- read.table("C:\\Users\\honey66clover\\Desktop\\Weijuan\\Courses\\Statistics\\588 Data Mining\\Dan Yang\\train4.data",sep=",")
x <- as.matrix(x)
y <- array(x,c(nrow(x),16,16))
y <- y[,,16:1]
x1 <- matrix(y,ncol=256)
windows()
par(mfrow = c(4,4),mar=rep(1,4)) #change the margins
for(i in 1:16){
  image(y[i,,],col=gray((32:0)/32))
}
# load 9
x <- read.table("C:\\Users\\honey66clover\\Desktop\\Weijuan\\Courses\\Statistics\\588 Data Mining\\Dan Yang\\train9.data",sep=",")
x <- as.matrix(x)
y <- array(x,c(nrow(x),16,16))
y <- y[,,16:1]
x2 <- matrix(y,ncol=256)
windows()
par(mfrow = c(4,4),mar=rep(1,4))
for(i in 1:16){
  image(y[i,,],col=gray((32:0)/32))
}
# combine
x <- rbind(x1,x2)
x.mean <- apply(x,2,mean)
# center
x.demean <- x-rep(1,nrow(x))%*%t(x.mean)
x.svd <- svd(x.demean)

# mean image and 8 eigen-images
windows()
par(mfrow = c(3,3),mar=rep(1,4))
image(matrix(x.mean,16),col=gray((32:0)/32))
for(i in 1:8){
  image(matrix(-x.svd$v[,i],16),col=gray((32:0)/32))
}
# scores
plot(x.svd$u[,1],x.svd$u[,2],col = c(rep(2,nrow(x1)),rep(3,nrow(x2))))
# singular values
plot(x.svd$d)
# approximation, take the first 10 principal components
r <- 10
x.reconstruct <- rep(1,nrow(x)) %*% t(x.mean) + 
                 x.svd$u[,1:r] %*% diag(x.svd$d[1:r]) %*% t(x.svd$v[,1:r])
# images of reconstruction
windows()
par(mfrow = c(4,4),mar=rep(1,4))
for(i in 1:8){ #these are 4s
  image(matrix(x.reconstruct[i,],16),col=gray((32:0)/32))
}
for(i in 1:8){ #these are 9s
  image(matrix(x.reconstruct[i+nrow(x1),],16),col=gray((32:0)/32))
}
###########################################################################

