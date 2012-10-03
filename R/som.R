require(klaR)

feat <- read.csv(file="vir.csv",head=TRUE,sep=" ")
library(som)
featsom <- som(feat[,2:ncol(feat)], xdim = 6, ydim = 6)
feat2 <- feat[c(2:ncol(feat),1)]
shardsplot(featsom, data.or = feat2, vertices = FALSE)
opar <- par(xpd = NA)
legend(7.5, 6.1, col = rainbow(3), xjust = 0.5, yjust = 0,legend = levels(factor(feat[,1])), pch = 16, horiz = TRUE)
par(opar)    

