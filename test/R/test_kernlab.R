library(kernlab)
# another example using the iris
data(iris)
test <- sample(1:150,20)
kpc <- kpca(~.,data=iris[-test,-5],kernel="rbfdot",kpar=list(sigma=0.2),features=2)

#print the principal component vectors
pcv(pc)
#plot the data projection on the components
plot(rotated(kpc),col=as.integer(iris[-test,5]))

#embed remaining points
emb <- predict(kpc,iris[test,-5])
points(emb,col=as.integer(iris[test,5]))

pc <- prcomp(iris[-test,-5], scale = TRUE, center = TRUE, retx = TRUE)
plot(pc$x[,1:2],col=as.integer(iris[-test,5]))

