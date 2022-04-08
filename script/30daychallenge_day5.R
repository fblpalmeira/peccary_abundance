tay<-read.csv("tay_pcount.open.csv")
site.covs<-read.csv("ts.sitecovs_500m.csv")
obs1.covs<-read.csv("ts.sample_hour.csv") 
obs2.covs<-read.csv("ts.sample_juliandate.csv") 
obs.covs<-list(hour=obs1.covs, date=obs2.covs);is.na(obs.covs) <-NA
years<-matrix(c('01', '02', '03', '04'), nrow(tay), 4, byrow=T) 

library(unmarked)
upc.tay<-unmarkedFramePCO (y=tay, siteCovs=site.covs, obsCovs=obs.covs, yearlySiteCovs=list(year=years), numPrimary=4)

sc<-siteCovs(upc.tay) 
sc<-scale(sc)
siteCovs(upc.tay)<-sc

oc<-obsCovs(upc.tay) 
oc<-scale(oc)
obsCovs(upc.tay)<-oc

fm1 <- pcountOpen (~1, ~1, ~1, ~1, upc.tay, K=50)
fm2 <- pcountOpen (~1, ~1, ~1, ~hour+date+year-1, upc.tay, K=50)
fm3 <- pcountOpen (~primary, ~1, ~1, ~primary+hour+date+year-1, upc.tay, K=50)
fm4 <- pcountOpen (~secondary, ~1, ~1, ~secondary+hour+date+year-1, upc.tay, K=50)
fm5 <- pcountOpen (~plantation, ~1, ~1, ~plantation+hour+date+year-1, upc.tay, K=50)
fm6 <- pcountOpen (~stream, ~1, ~1, ~stream+hour+date+year-1, upc.tay, K=50)
fm7 <- pcountOpen (~primary+stream, ~1, ~1, ~primary+stream+hour+date+year-1, upc.tay, K=50)
fm8 <- pcountOpen (~secondary+stream, ~1, ~1, ~secondary+stream+hour+date+year-1, upc.tay, K=50)
fm9 <- pcountOpen (~plantation+stream, ~1, ~1, ~plantation+stream+hour+date+year-1, upc.tay, K=50)

names(fm7)

models<-fitList('(fm1)' = fm1,
                '(fm2)' = fm2,
                '(fm3)' = fm3,
                '(fm4)' = fm4, 
                '(fm5)' = fm5,
                '(fm6)' = fm6,
                '(fm7)' = fm7,
                '(fm8)' = fm8,
                '(fm9)' = fm9)

ms<-modSel(models) 
ms

lam.coef <- exp(coef(fm7, type="lambda"))
lam.coef
lam.SE <- exp(SE(fm7, type="lambda"))
lam.SE
lam.confint <-exp(confint(fm7, type="lambda"))
lam.confint

p.coef <- plogis(coef(fm7, type="det"))
p.coef
p.SE <- plogis(SE(fm7, type="det"))
p.SE
p.confint <- plogis(confint(fm7, type="det"))
p.confint

SITE1 <- as.matrix(site.covs[,1])
y.tay <- as.matrix(tay[,1:32])
sd.SITE1 <- sd(c(SITE1), na.rm=TRUE)
mean.SITE1 <- mean(SITE1, na.rm=TRUE)
SITE1 <- (SITE1 - mean.SITE1) / sd.SITE1

SITE2 <- as.matrix(site.covs[,2])
y.tay <- as.matrix(tay[,1:32])
sd.SITE2 <- sd(c(SITE2), na.rm=TRUE)
mean.SITE2 <- mean(SITE2, na.rm=TRUE)
SITE2 <- (SITE2 - mean.SITE2) / sd.SITE2

SITE3 <- as.matrix(site.covs[,3])
y.tay <- as.matrix(tay[,1:32])
sd.SITE3 <- sd(c(SITE3), na.rm=TRUE)
mean.SITE3 <- mean(SITE3, na.rm=TRUE)
SITE3 <- (SITE3 - mean.SITE3) / sd.SITE3

HOUR <- as.matrix(obs1.covs[,1:32])
y.tay <- as.matrix(tay[,1:32])
y.tay[is.na(HOUR) != is.na(y.tay)] <- NA
sd.HOUR <- sd(c(HOUR), na.rm=TRUE)
mean.HOUR <- mean(HOUR, na.rm=TRUE)
HOUR <- (HOUR - mean.HOUR) / sd.HOUR

DATE <- as.matrix(obs2.covs[,1:32])
y.tay <- as.matrix(tay[,1:32])
y.tay[is.na(DATE) != is.na(y.tay)] <- NA
sd.DATE <- sd(c(DATE), na.rm=TRUE)
mean.DATE <- mean(DATE, na.rm=TRUE)
DATE <- (DATE - mean.DATE) / sd.DATE

library(ggplot2)

png(file = "peccary_plot.png", width = 1224, height = 768)

par(mfrow=c(2,3),mar = c(10, 5, 5, 5))

fm1 <- pcountOpen (~primary, ~1, ~1, ~1, upc.tay, K=50)
nd <- data.frame(primary=seq(-1.83, 1.25, length=50))
E.psi <- predict(fm1, type="lambda", newdata=nd, appendData=TRUE)
E.psi$primaryOrig <- E.psi$primary*sd.SITE1 + mean.SITE1
with(E.psi, {
  plot(primaryOrig, Predicted, ylim=c(0,20), type="l", col= "blue",
       xlab="Percent of primary forest",
       ylab="Abundance", cex.lab=2.2, cex.axis=2.2)
  lines(primaryOrig, Predicted+1.96*SE, col=gray(0.7))
  lines(primaryOrig, Predicted-1.96*SE, col=gray(0.7))
})

fm1 <- pcountOpen (~secondary, ~1, ~1, ~1, upc.tay, K=50)
nd <- data.frame(secondary=seq(-.54, 1.93, length=50))
E.psi <- predict(fm1, type="lambda", newdata=nd, appendData=TRUE)
E.psi$secondaryOrig <- E.psi$secondary*sd.SITE2 + mean.SITE2
with(E.psi, {
  plot(secondaryOrig, Predicted, ylim=c(0,20), type="l", col= "red",
       xlab="Percent of secondary forest",
       ylab="Abundance", cex.lab=2.2, cex.axis=2.2)
  lines(secondaryOrig, Predicted+1.96*SE, col=gray(0.7))
  lines(secondaryOrig, Predicted-1.96*SE, col=gray(0.7))
})

fm1 <- pcountOpen (~plantation, ~1, ~1, ~1, upc.tay, K=50)
nd <- data.frame(plantation=seq(-.79, 1.98, length=50))
E.psi <- predict(fm1, type="lambda", newdata=nd, appendData=TRUE)
E.psi$plantationOrig <- E.psi$plantation*sd.SITE3 + mean.SITE3
with(E.psi, {
  plot(plantationOrig, Predicted, ylim=c(0,20), type="l", col= "red",
       xlab="Percent of forest plantation",
       ylab="Abundance", cex.lab=2.2, cex.axis=2.2)
  lines(plantationOrig, Predicted+1.96*SE, col=gray(0.7))
  lines(plantationOrig, Predicted-1.96*SE, col=gray(0.7))
})

fm1 <- pcountOpen (~1, ~1, ~1, ~primary, upc.tay, K=50)
nd <- data.frame(primary=seq(-1.83, 1.25, length=50))
E.p <- predict(fm1, type="det", newdata=nd, appendData=TRUE)
E.p$primaryOrig <- E.p$primary*sd.SITE1 + mean.SITE1
with(E.p, {
  plot(primaryOrig, Predicted, ylim=c(0,.2), type="l", col= "blue", 
       xlab="Percent of primary forest", 
       ylab="Detection", cex.lab=2.2, cex.axis=2.2)
  lines(primaryOrig, Predicted+1.96*SE, col=gray(0.7))
  lines(primaryOrig, Predicted-1.96*SE, col=gray(0.7))
})

fm1 <- pcountOpen (~1, ~1, ~1, ~secondary, upc.tay, K=50)
nd <- data.frame(secondary=seq(-.54,1.93, length=50))
E.p <- predict(fm1, type="det", newdata=nd, appendData=TRUE)
E.p$secondaryOrig <- E.p$secondary*sd.SITE2 + mean.SITE2
with(E.p, {
  plot(secondaryOrig, Predicted, ylim=c(0,.2), type="l", col= "red", 
       xlab="Percent of secondary forest", 
       ylab="Detection", cex.lab=2.2, cex.axis=2.2)
  lines(secondaryOrig, Predicted+1.96*SE, col=gray(0.7))
  lines(secondaryOrig, Predicted-1.96*SE, col=gray(0.7))
})

fm1 <- pcountOpen (~1, ~1, ~1, ~plantation, upc.tay, K=50)
nd <- data.frame(plantation=seq(-.79, 1.98, length=50))
E.p <- predict(fm1, type="det", newdata=nd, appendData=TRUE)
E.p$plantationOrig <- E.p$plantation*sd.SITE3 + mean.SITE3
with(E.p, {
  plot(plantationOrig, Predicted, ylim=c(0,.2), type="l", col= "red", 
       xlab="Percent of forest plantation", 
       ylab="Detection", cex.lab=2.2, cex.axis=2.2)
  lines(plantationOrig, Predicted+1.96*SE, col=gray(0.7))
  lines(plantationOrig, Predicted-1.96*SE, col=gray(0.7))

})

dev.off()

library(magick)
library(magrittr) 

# Call back the plot
plot <- image_read("peccary_plot.png")
plot2<-image_annotate(plot, "Effects of vegetation cover on abundance of the white-lipped peccary in the southern Brazilian Amazonia", color = "blue", size = 25, 
                      location = "10+50", gravity = "north")
plot3<-image_annotate(plot2, "Data source: Palmeira, FBL and Trinca, CTT (unpublished data) | Visualization by @fblpalmeira", color = "gray", size = 20, 
                      location = "10+50", gravity = "southeast")
# And bring in a logo
logo_raw <- image_read("peccary_avatar.png") 
out<-image_composite(plot3,image_scale(logo_raw,"x50"), offset = "+50+50")
image_browse(out)

# And overwrite the plot without a logo
image_write(out, "30daychallenge_day5.png")

