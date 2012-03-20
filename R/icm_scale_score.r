#install.packages(c("RSQLite","popbio","plotrix","hdrcde","Hmisc","evd","NORMT3","biglm","rjson"))

library(plotrix)

#for nclass.Sturges etc which select the histogram breaks
library(grDevices)
histCond<-function(y,x,breaks=NULL) {
    if (is.null(breaks)) breaks = nclass.Sturges(x)
    wh = weighted.hist(x,y,breaks=breaks,freq=FALSE,plot=FALSE)
    h = hist(x,breaks=wh$breaks,freq=FALSE,plot=FALSE)
    h$counts.raw = h$counts
    h$counts = wh$counts/h$counts
    h$density = wh$density/h$density
    h$intensities = h$density
    h
}

library(RSQLite)

drv <- dbDriver("SQLite")
con <- dbConnect(drv, dbname = "bench.comb.sqlite")

dinp = dbGetQuery(con, "select * from samp_acc 
where
not i_lev_per in (80,100)
and
i_samp 
in (select distinct i_samp from samp_acc 
order by random() limit 100000)
")

dinp$is_correct <- dinp$is_correct==TRUE
#dinp$i_lev_per <- ordered(dinp$i_lev_per)
#dinp$score_b <- dinp$score/(dinp$len_samp)^0.5
dinp$score_b <- dinp$score/(dinp$len_samp)
#dinp$score_b <- (dinp$score_b - mean(dinp$score_b))/(dinp$len_samp)^0.5
#dinp$score_b <- dinp$score

#attach(dinp)
#library(lattice)
#histogram(~score_b | is_correct*len_samp)
#bwplot(is_correct ~ score_b | i_lev_per * len_samp)
#detach(dinp)

dfilt = dinp[dinp$len_samp==400 & dinp$i_lev_per==50,]
#dfilt = dinp
attach(dfilt)

#library(popbio)
#logi.hist.plot(score_b,is_correct,boxp=TRUE,rug=TRUE,type="hist",col="gray")


#library(hdrcde)
#is_corr.cde = cde(score_b,is_correct)
#plot(is_corr.cde)
#plot(is_corr.cde,plot.fn="hdr")

gmod <- glm(is_correct ~ score_b,family=binomial)

#plot(score_b,is_correct,xlab="Score",ylab="Probability of being correct")
hc = histCond(is_correct,score_b)
plot(hc$mids,hc$counts,col="red")
#this creates a vectorized function (late binding for 'predict',
#which is then called by 'curve' with different values of x.
#This will get correct range for x from the 'plot' above.
curve(predict(gmod,data.frame(score_b=x),type="resp"),add=TRUE)

library(boot)
plot(hc$mids,logit(hc$counts),col="red")
#plot(hc$mids,hc$counts,col="red",cex=(10*hc$counts.raw/sum(hc$counts.raw)))
curve(logit(predict(gmod,data.frame(score_b=x),type="response")),add=TRUE)
summary(dfilt)
sum(score_b<=1.0)/length(score_b)
#Weibull plot
#library(stats)
#curve(-log(1-ecdf(-score_b)(x)),1,10,log="xy")

library(evd)
#fit generalized extreme value distribution - fits quite well 
#except at the upper quantiles where the observed
#is higher.
#The GEV distribution function with parameters 'loc' = a, 'scale' =
#b and 'shape' = s is
#G(x) = exp[-{1+s(z-a)/b}^(-1/s)]                
#		       
#for 1+s(z-a)/b > 0, where b > 0.  If s = 0 the distribution is
#defined by continuity.  If 1+s(z-a)/b <= 0, the value z is either
#greater than the upper end point (if s < 0), or less than the
#lower end point (if s > 0).
#
#The parametric form of the GEV encompasses that of the Gumbel,
#Frechet and reversed Weibull distributions, which are obtained for
#s = 0, s > 0 and s < 0 respectively. 

#score.evd = fgev(score_b)
#plot(score.evd)
#curve(dgev(x,scale=score.evd$param["scale"],shape=score.evd$param["shape"],loc=score.evd$param["loc"]),from=-2,to=0,col="red")

phymmScoreGenus<-function(score,len_samp)
{
	0.86 - (-score-(1.25*len_samp-550))^(6+12/len_samp)/1.5e17 + 8/len_samp
}


#for erfc
library(NORMT3)

dThisModelWins<-function(x,m0,s0,m1,s1,n) {
	dnorm(x,mean=m0,sd=s0)*(0.5*erfc((x-m1)/(s1*sqrt(2))))^n
}

pCorrectAtScore<-function(x,m1,s1,n) {
	Re((0.5*erfc((x-m1)/(s1*sqrt(2))))^n)
	#y = 1
	#for (index in 1:n) {
	#	y=y*0.5*erfc((x-m1-rnorm(1,0,s1/2))/(s1*sqrt(2)))
	#}
	#Re(y)
}

plotScoreForLenAndLev<-function(dinp)
{
	score.len_samp = 400
	data = dinp[dinp$i_lev_per==30 & dinp$len_samp==score.len_samp,]
	with(data,
	{
		#plot(score,phymmScoreGenus(score,score.len_samp),col="blue",ylim=c(0,1))
		hc = histCond(is_correct,score)
		plot(hc$mids,hc$counts,col="red")
		#points(hc$mids,hc$counts.raw/max(hc$counts.raw),col="green")
		dens = density(score)
		points(dens$x,dens$y/max(dens$y),col="green")
		#points(score,is_correct,col="cyan")
		points(score,pCorrectAtScore(-score_b,3,0.7,100),col="yellow")
	})
}
plotScoreForLenAndLev(dinp)

#plot(-dinp$score_b,pCorrectAtScore(-dinp$score_b,mean(-dinp$score_b),sqrt(var(-dinp$score_b)),3000))

predict.lreg<-function(lreg,x)
{
	b<-coef(lreg)
	return  (b["(Intercept)"] + 		
		b["score_b"]*x$score_b+		
		b["len_samp.pre"]*x$len_samp.pre+		
		b["i_lev_per.pre"]*x$i_lev_per.pre+		
		b["score_b:len_samp.pre"]*x$score_b*x$len_samp.pre)
}

makeAndPlotAccuracyLreg<-function(score.lens,dinp)
{
	#data = data.frame(is_correct=dinp$is_correct,
	#	len_samp=dinp$len_samp,
	#	len_samp.pre=1/log(sqrt(dinp$len_samp)),
	#	score_b=dinp$score_b*log(sqrt(dinp$len_samp)))
	## We  need to deal with the i_lev_per, which is by nature 
	## an ordered factor variable. If we use "ordered(i_lev_per)",
	## glm builds x*lev1,x^2*lev2,x^3*lev3,x^4*lev4.
	## It is not clear yet to me if that model also enforces
	## monotonicity. There is a more flexible way of dealing
	## with factors by using 'C' function that maps
	## factor levels to numbers (see help for C).
	## However, if we just used factor(i_lev_per), the coefficients
	## produced on the current input dataset turn out to be ordered,
	## resulting in a regression that is also isotonic.
	data = data.frame(is_correct=dinp$is_correct,
		i_lev_per=dinp$i_lev_per,
		#i_lev_per.pre=factor(dinp$i_lev_per),
		i_lev_per.pre=dinp$i_lev_per,
		len_samp=dinp$len_samp,
		len_samp.pre=dinp$len_samp^(-1/exp(1)),
		#len_samp.pre=dinp$len_samp^-0.5,
		score_b=dinp$score_b)
	data$i_lev_per.pre[data$i_lev_per==90]=115
	#data = data.frame(is_correct=dinp$is_correct,
	#	len_samp=dinp$len_samp,
	#	len_samp.pre=dinp$len_samp^(-1/exp(1)),
	#	score_b=pCorrectAtScore(-dinp$score_b,mean(-dinp$score_b),sqrt(var(-dinp$score_b)),3000))
	#xlim=c(0.3,4)
	#ylim=c(-10,10)
	xlim=c(-2,0)
	ylim=c(-5,7)
	#xlim=c(0,2)
	#ylim=c(0,2)
	#data = data.frame(is_correct=dinp.len$is_correct,score_b=dinp.len$score_b,dinp.)
    	lreg.ret=NULL
	with(data,{
	score.lreg <- glm(is_correct ~ score_b*len_samp.pre+i_lev_per.pre,
        	family=binomial(link = "logit"))
	assign("lreg.ret",score.lreg,inherits=TRUE)
	print(summary(score.lreg))
	#print(anova(score.lreg))
	#plot(score_b,is_correct,xlab="Score",ylab="Probability of being correct")
	plot.colors = rainbow(length(score.lens))
	#i_lev_per.levels = levels(factor(i_lev_per))
	i_lev_per.levels = c(30,90)
	#par(mfrow=c(2,ceiling(length(i_lev_per.levels)/2)))
	par(mfrow=c(1,2))
	for (score.i_lev_per in i_lev_per.levels)
	{
		i.score.len = 1
		for (score.len in score.lens)
		{
			data.len_fix = data[len_samp==score.len & i_lev_per==score.i_lev_per,]
			with(data.len_fix,{
			#score.lreg <- glm(is_correct ~ score_b,family=binomial)
			#print(score.lreg)
			hc = histCond(is_correct,score_b,breaks=100)
			if (i.score.len==1) plotFunc=plot
			else plotFunc=points
			#plotFunc(hc$mids,hc$counts,col=plot.colors[i.score.len],xlim=xlim,ylim=ylim)
			plotFunc(hc$mids,logit(hc$counts),col=plot.colors[i.score.len],
				xlim=xlim,ylim=ylim,main=paste("Prediction rank level ",score.i_lev_per))
			#points(score_b,logit(pCorrectAtScore(-score_b,2.5,0.4,600)),col=plot.colors[i.score.len],pch="x")
			#this creates a vectorized function (late binding for 'predict',
			#which is then called by 'curve' with different values of x.
			#This will get correct range for x from the 'plot' above.
			n.curve = 101
			#curve(predict(score.lreg,data.frame(score_b=x,len_samp.pre=rep(len_samp.pre[1],n.curve)),type="resp"),n=n.curve,col=plot.colors[i.score.len],add=TRUE)
			curve(logit(predict(score.lreg,
					data.frame(score_b=x,
						len_samp.pre=rep(len_samp.pre[1],n.curve),
						i_lev_per.pre=rep(i_lev_per.pre[1],n.curve)
						),
					type="resp")),
				n=n.curve,col=plot.colors[i.score.len],add=TRUE)
			#points(score_b,predict.lreg(score.lreg,data.len_fix),pch="x",col=plot.colors[i.score.len])
			})
			i.score.len = i.score.len + 1
		}
		legend("topleft", as.character(score.lens), cex=0.8, col=plot.colors, lty=1:3, lwd=2, bty="n");
	}
	})
    lreg.ret
}

score.lens = c(100,400,1000,10000)
lreg.ret = makeAndPlotAccuracyLreg(score.lens,dinp)
print(toJSON(coef(lreg.ret)))

testWith<-function()
{
	ret = NULL
	data=data.frame(a=c(1,2,3),b=c("a","b","c"))
	with(data,{
		assign("ret",1,inherits=TRUE)
	})
	ret
}

