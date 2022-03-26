error
rm(list=ls())

setwd("~/Dropbox/Projects/013_Ant revision/code/")
load("../data/vDateSubplotDepth.rda")
bacdat<-bacdat[c(grep("2011-",bacdat$date), grep("2012-",bacdat$date)),]
head(bacdat)

#Humidity?
bacdat<-bacdat[bacdat$depth>-20,]
bacdat<-bacdat[bacdat$richness<32,]

hum<-as.numeric(as.character(bacdat$meanHumidity))
hm<-with(bacdat, tapply(hum, list(richness, subplot), function(x) quantile(x, 0.5, na.rm=T)))
hl<-with(bacdat, tapply(hum, list(richness, subplot), function(x) quantile(x, 0.1586553, na.rm=T)))
hu<-with(bacdat, tapply(hum, list(richness, subplot), function(x) quantile(x, 0.8413447, na.rm=T)))

agbacdat<-aggregate(x=bacdat$meanTemp, by=list(bacdat$date, bacdat$richness, bacdat$subplot), FUN=function(x) cbind(mean(x), sd(x)))
agbacdat<-data.frame(agbacdat[,1:3], agbacdat[,4][,1], agbacdat[,4][,2])
agbacdat<-cbind(agbacdat, t(matrix(nrow=3, data=unlist(strsplit(as.character(agbacdat[,1]), "-", fixed=T)))))
for(i in 6:8) {
  agbacdat[,i]<-as.numeric(as.character(agbacdat[,i]))
}
colnames(agbacdat)<-c("date", "richness", "subplot", "meantmp", "sdtmp", "Year", "Month", "Day")
dayofyearlookup<-cumsum(c(0,31,28,31,30,31,30,31,31,30,31,30))
agbacdat$dayofyear<-dayofyearlookup[agbacdat$Month]+agbacdat$Day
head(agbacdat)


bacdat<-bacdat[bacdat$depth==1,]

#Make grand average
agbacdatgrand<-agbacdat[agbacdat$richness==1 & agbacdat$subplot=="control",]
agbacdat$grandmeantmp<-agbacdatgrand$meantmp[match(as.character(agbacdat$date), as.character(agbacdatgrand$date))]
#################

#Biomass data
dat2<-read.table("../data/2013MegaE249.dat", header=T, sep=",")
dat2$NumSp[dat2$NumSp==34]<-32
dat2<-dat2[dat2$NumSp<32,]
dat2<-dat2[dat2$Year%in%c(2011, 2012),]





############################################################
####################Make figure
############################################################
pdf("figures/BAC_tmp_figure.pdf", width=8, height=6)
par(mfrow=c(2,2), mar=c(1,4,2,2), oma=c(1,0,0,0))
collst2<-c("firebrick4", "dodgerblue4")
collst2<-adjustcolor(collst2, alpha.f = 0.6)
#FIGURES
#Mean temp
datplot<-subset(agbacdat, Year==2011)
datplot<-unique(datplot[,c("dayofyear", "grandmeantmp")])
datplot<-datplot[order(datplot$dayofyear),]
par(mar=c(c(2, 4.5, 3.4, 2) + 0.1))
collst<-c("black", "darkgrey")
plot(grandmeantmp~dayofyear, data=datplot, type="l", xlab="",
     xlim=c(115,310), ylim=c(0, 33),
     ylab=expression(paste("Mean Temperature, "^"o", "C")), axes=F, lwd=2,
     main="a. Temperature in Control Plots", col=collst[1])
datplot<-subset(agbacdat, Year==2012)
datplot<-unique(datplot[,c("dayofyear", "grandmeantmp")])
datplot<-datplot[order(datplot$dayofyear),]
lines(datplot$dayofyear, datplot$grandmeantmp, lwd=2, col=collst[2])
axis(2)
axis(1, c("May", "June", "July", "Aug", "Sept", "Oct", "Nov"), at=dayofyearlookup[5:11])
box()
legend("bottomleft", c("2011", "2012"), lty=1, lwd=2, col=collst, bty="n")

#Biomass
par(mar=c(c(2, 4.5, 3.4, 4.5) + 0.1))
plot(c(0.05, 0.95), c(0, 600),
     xlab="Sown Plant Species", ylab=expression(paste("Aboveground Biomass, g m"^"-2")), main="b. Plot Biomass and Humidity", col="grey", cex=0.5, axes=F, type="n")
axis(2)
divlst<-c(0.2, 0.5, 0.8)
divlst2<-c(1, 4, 16)
axis(1, c("1", "4", "16"), at=divlst)
meanlst<-numeric(3)
for(i in 1:3) {
  x<-divlst[i]
  y<-dat2[dat2$NumSp==divlst2[i],]$TotalAbovegroundBiomass.g.m2.
  qt<-quantile(y, c(0.1586553, 0.5, 0.8413447))
  segments(x-0.04, qt[1], x-0.04, qt[3], col=1, lwd=2, lend=2)
  points(x-0.04, qt[2], lwd=2, pch=16)
  meanlst[i]<-qt[2]
}
lines(divlst-0.04, meanlst, lwd=2)
par(new=T)
plot(c(0.05,0.95), c(50, 100), xlab="", ylab="", axes=F, type="n")
axis(4, seq(50,100,by=10), at=seq(50,100,by=10))
mtext("% Air Humidity", 4, line=2.5)
mtext("Sown Plant Richness", 1, line=1.8)
segments(divlst, hl[2:4,1], divlst, hu[2:4,1], lwd=2, col=collst2[2], lend=2)
points(divlst, hm[2:4,1], col=collst2[2], lwd=2, pch=16)
lines(divlst, hm[2:4,1], col=collst2[2], lwd=2)
segments(divlst, hl[2:4,1], divlst, hu[2:4,1], lwd=2, col=collst2[2], lend=2)
points(divlst+0.04, hm[2:4,2], col=collst2[1], lwd=2, pch=16)
lines(divlst+0.04, hm[2:4,2], col=collst2[1], lwd=2)
segments(divlst+0.04, hl[2:4,2], divlst+0.04, hu[2:4,2], lwd=2, col=collst2[1], lend=2)
legend("topleft", c("Biomass", "Hum., C", "Hum., H"), lty=1, pch=16, lwd=2, col=c(1, collst2[2], collst2[1]), ncol=3, bty="n", cex=0.8)
box()

#Delta temp
trtlst<-c("control", "high")
titlelst<-c("c. Control Plots", "d. Heated Plots")
legendpos<-c("topright", "bottomright")
collst<-rbind(c("lightblue3", "cornflowerblue", "darkblue"),
              c("lightpink2", "darkorange1", "darkred"))


for(j in 1:2) {
agbacdat$dtemp<-agbacdat$meantmp-agbacdat$grandmeantmp
datplot<-subset(agbacdat, richness==1 & subplot==trtlst[j])
datplot<-aggregate(x=datplot$dtemp, by=list(dayofyear=datplot$dayofyear), FUN=mean)
colnames(datplot)[2]<-"dtemp"
par(mar=c(c(2, 4.5, 3.4, 2) + 0.1))
plot(dtemp~dayofyear, data=datplot, type="l", ylim=c(-4.5, 4), xlim=c(115,310), xlab="",
     ylab=expression(paste("Temperature Difference, "^"o", "C")), axes=F, lwd=1.5,
     main=titlelst[j], col=adjustcolor(collst[j,1], alpha.f = 1), lty=1)
#lines(dtemp~dayofyear, data=datplot, col=collst[j,1], lwd=1.5, lty=4)

axis(2)
axis(1, c("May", "June", "July", "Aug", "Sept", "Oct", "Nov"), at=dayofyearlookup[5:11])
abline(h=0, lty=2, col="darkgrey")
for(i in 1:2) {
  datplot<-subset(agbacdat, richness==c(4, 16, 32)[i]&subplot==trtlst[j])
  datplot<-aggregate(x=datplot$dtemp, by=list(dayofyear=datplot$dayofyear), FUN=mean)
  colnames(datplot)[2]<-"dtemp"
  with(datplot, lines(dayofyear, dtemp, col=adjustcolor(collst[j,i+1], alpha.f = 1), lwd=1.5, lty=1))
  #with(datplot, lines(dayofyear, dtemp, col=collst[j,i+1], lwd=1.5, lty=4-i))
}
box()
legend(legendpos[j], c("1", "4", "16"), col=collst[j,], lwd=1.5, lty=1, bty="n", ncol=3, title="Sown Plant Species")
}
dev.off()

