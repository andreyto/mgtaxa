library(grDevices)
library(kernlab)
feat <- read.csv(file="all.csv",head=TRUE,sep=",")
col <- read.csv(file="all.col.csv",head=TRUE,sep=",",as.is=TRUE)
train <- sample(1:nrow(feat),200)
test <- sample(1:nrow(feat),20)
kpc <- kpca(~.,data=feat[train,c(-1,-2)],kernel="rbfdot",kpar=list(sigma=0.2),features=2)
#print the principal component vectors
pcv(kpc)

#plot the data projection on the components
#plot(rotated(kpc),col=as.integer(feat[train,2]))
plot(rotated(kpc),col=col[train,-1])

#embed remaining points
emb <- predict(kpc,feat[test,c(-1,-2)])
points(emb,col=as.integer(feat[test,2]))
#scatterplot3d
#
pc <- prcomp(feat[train,c(-1,-2)], scale = FALSE, center = TRUE, retx = TRUE)
plot(pc$x[,1:2],col=col[train,-1])
