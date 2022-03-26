error
rm(list=ls())

collst<-rbind(c("lightblue3", "cornflowerblue", "darkblue"),
              c("lightpink2", "darkorange1", "darkred"))
dayofyearlookup<-cumsum(c(0,31,28,31,30,31,30,31,31,30,31,30))

setwd("~/Dropbox/Projects/013_Ant revision/code/")
tempdat = read.csv("../data/vHourMonthSubplot2012feb11.csv")
tempdat = tempdat[tempdat$month%in%c(8,9),]

load("../data/vDateSubplotDepth.rda")
bacdat<-bacdat[c(grep("2011-",bacdat$date), grep("2012-",bacdat$date)),]
bacdat$year = as.numeric(substr(bacdat$date, 1,4))
bacdat$month = as.numeric(substr(bacdat$date, 6,7))
bacdat$day = as.numeric(substr(bacdat$date, 9,10))

# get plant diversity data
plantdat = read.csv("../data/pitfall_ants_bac_140913.csv")
tempdat$Plant_Species = plantdat$Plant_Species[match(tempdat$plot, plantdat$PlotID)]
tempdat = tempdat[tempdat$Plant_Species%in%c(1,4,16),]

agdat = with(tempdat,
     aggregate(meanHourlyMean, by=list(month=month, subplot=subplot, hour=hour, Plant_Species=Plant_Species), function(x) c(mean = mean(x), sd = sd(x), n = length(x))))

### plot average daily trend
xrng = c(0, 23)
yrng = c(12, 35)

pdf("figures/dailytemp_trend.pdf", width=7, height=5, colormodel = "cmyk", useDingbats = FALSE)
  par(mfrow=c(2,3), mar = c(2,2,1,1), oma=c(1,4,1,0))
  for(mnth in 8:9) {
    for(ind in 1:3) {
      i = c(1,4,16)[ind]
      
      subs1 = which(agdat$Plant_Species==i & agdat$subplot=="control" & agdat$month==mnth)
      subs2 = which(agdat$Plant_Species==i & agdat$subplot=="high" & agdat$month==mnth)
      
      subs1 = subs1[order(agdat$hour[subs1])]
      subs2 = subs2[order(agdat$hour[subs2])]
      
      
      plot(xrng, yrng, type = "n", xlab = "", ylab = "", xaxs = "i")
      title(paste(letters[ind+(mnth-8)*3],".", sep=""), line=-0.95, adj = 0.025)
      
      polygon(c(agdat$hour[subs1], rev(agdat$hour[subs1])),
              c(agdat$x[,1][subs1]+agdat$x[,2][subs1],
                rev(agdat$x[,1][subs1]-agdat$x[,2][subs1])),
              border = 1, col=adjustcolor(collst[1,ind], alpha.f = 0.8))
      
      polygon(c(agdat$hour[subs2], rev(agdat$hour[subs2])),
              c(agdat$x[,1][subs2]+agdat$x[,2][subs2],
                rev(agdat$x[,1][subs2]-agdat$x[,2][subs2])),
              border = 1, col=adjustcolor(collst[2,ind], alpha.f = 0.8))
      
      if(ind==1) {
        mtext(c("August", "September")[mnth-7], 2, line=2.2)
      }
      if(mnth==8) {
        mtext(paste("Sown Richness =", i), 3, line=0.2)
      }
      
      
      # identify significance of differences at p = 0.05
      tempdiff = agdat$x[,1][subs1]-agdat$x[,1][subs2][match(agdat$hour[subs1], agdat$hour[subs2])]
      jointsd = sqrt(agdat$x[,2][subs1]^2+agdat$x[,2][subs2][match(agdat$hour[subs1], agdat$hour[subs2])]^2)/sqrt(unique(agdat$x[,3][subs1])[1]-1)
      pv = pnorm(tempdiff, 0, jointsd)
      
      tmp=which(is.finite(pv) & pv < 0.05)
      k = 1
      while(k < length(tmp)) {
        if(is.na(subs1[tmp[k]])) {
          k = k+1
        } else {
          if(is.na(subs1[tmp[k+1]]) || (tmp[k+1]!=(tmp[k]+1))) {
            #abline(v=agdat$hour[subs1][tmp[k]], col=adjustcolor("black", 0.1))
            k = k+1
          } else {
            startk = k
            while((tmp[k+1]==(tmp[k]+1)) & (k < length(tmp))) {
              k = k+1
              endk = k
            }
            polygon(c(agdat$hour[subs1][tmp[c(startk,startk,endk,endk)]]),
                    c(0, 100, 100, 0),
                    col=adjustcolor("black", 0.1), border = NA)
            k = k+1
          }
        }
      }
      
    }
  }
  mtext("Mean Hourly Temperature, °C", side = 2, outer=TRUE, line=2.2)
dev.off()



### plot max temperature over August and September
bacdat$dayofyear = bacdat$day+dayofyearlookup[bacdat$month]

agdat2 = with(bacdat[bacdat$depth == 1,],
             aggregate(maximumTemp, by=list(dayofyear=dayofyear, year = year, subplot=subplot, Plant_Species=richness), function(x) c(mean = mean(x), sd = sd(x), n = length(x))))


pdf("figures/maxtemp.pdf", width=7, height=5, colormodel = "cmyk", useDingbats = FALSE)
  par(mfcol=c(2,3), mar = c(2,2,1,1), oma=c(1,4,1,0))
  xlm = dayofyearlookup[c(8,10)]
  ylm = c(15,55)
  for(ind in 1:3) {
    for(yr in c(2011, 2012)) {
      i = c(1,4,16)[ind]
      subs = which(agdat2$Plant_Species==i & agdat2$year==yr & agdat2$subplot=="control")
      subs = subs[order(agdat2$dayofyear[subs])]
      
      # remove gaps
      steps = diff(agdat2$dayofyear[subs])
      tmp = which(steps > 1)
      for(j in 1:length(tmp)) {
        subs[tmp[j]+0:1] = NA
      }
      plot(agdat2$dayofyear[subs], agdat2$x[,1][subs], type = "l", xlim = xlm, ylim = ylm, axes = F, ylab = "", xlab = "",
           lwd = 1.5, col=adjustcolor(collst[1,ind], alpha.f = 0.8))
      title(paste(letters[ind+(yr-2011)*3],".", sep=""), line=-0.95, adj = 0.025)
      
      subs2 = which(agdat2$Plant_Species==i & agdat2$year==yr & agdat2$subplot=="high")
      subs2 = subs2[order(agdat2$dayofyear[subs2])]
      
      # remove gaps
      steps = diff(agdat2$dayofyear[subs2])
      tmp = which(steps > 1)
      for(j in 1:length(tmp)) {
        subs2[tmp[j]+0:1] = NA
      }
      
      lines(agdat2$dayofyear[subs2], agdat2$x[,1][subs2],
           lwd = 1.5, col=adjustcolor(collst[2,ind], alpha.f = 0.8))
      axis(2)
      axis(1, c("May", "June", "July", "Aug", "Sept", "Oct", "Nov"), at=dayofyearlookup[5:11])
      box()
      
      # identify significance of differences at p = 0.05
      tempdiff = agdat2$x[,1][subs]-agdat2$x[,1][subs2][match(agdat2$dayofyear[subs], agdat2$dayofyear[subs2])]
      jointsd = sqrt(agdat2$x[,2][subs]^2+agdat2$x[,2][subs2][match(agdat2$dayofyear[subs], agdat2$dayofyear[subs2])]^2)/sqrt(unique(agdat2$x[,3][subs])[1]-1)
      pv = pnorm(tempdiff, 0, jointsd)
      
      tmp=which(is.finite(pv) & pv < 0.05)
      k = 1
      while(k < length(tmp)) {
        if(is.na(subs[tmp[k]])) {
          k = k+1
        } else {
          if(is.na(subs[tmp[k+1]]) || (tmp[k+1]!=(tmp[k]+1))) {
            #abline(v=agdat2$dayofyear[subs][tmp[k]], col=adjustcolor("black", 0.1))
            k = k+1
          } else {
            startk = k
            while((tmp[k+1]==(tmp[k]+1)) & (k < length(tmp))) {
              k = k+1
              endk = k
            }
            polygon(c(agdat2$dayofyear[subs][tmp[c(startk,startk,endk,endk)]]),
                    c(0, 100, 100, 0),
                    col=adjustcolor("black", 0.1), border = NA)
            k = k+1
          }
        }
      }
      if(ind==1) {
        mtext(yr, 2, line=2.2)
      }
      if(yr==2011) {
        mtext(paste("Sown Richness =", i), 3, line=0.2)
      }
    }
  }
  mtext("Maximum Daily Temperature, °C", side = 2, outer=TRUE, line=2.2)
dev.off()


### plot max temperature over August and September
bacdat$dayofyear = bacdat$day+dayofyearlookup[bacdat$month]

agdat2 = with(bacdat[bacdat$depth == 1,],
              aggregate(minimumTemp, by=list(dayofyear=dayofyear, year = year, subplot=subplot, Plant_Species=richness), function(x) c(mean = mean(x), sd = sd(x), n = length(x))))


pdf("figures/mintemp.pdf", width=7, height=5, colormodel = "cmyk", useDingbats = FALSE)
  par(mfcol=c(2,3), mar = c(2,2,1,1), oma=c(1,4,1,0))
  xlm = dayofyearlookup[c(8,10)]
  ylm = c(0,25)
  for(ind in 1:3) {
    for(yr in c(2011, 2012)) {
      i = c(1,4,16)[ind]
      subs = which(agdat2$Plant_Species==i & agdat2$year==yr & agdat2$subplot=="control")
      subs = subs[order(agdat2$dayofyear[subs])]
      
      # remove gaps
      steps = diff(agdat2$dayofyear[subs])
      tmp = which(steps > 1)
      for(j in 1:length(tmp)) {
        subs[tmp[j]+0:1] = NA
      }
      plot(agdat2$dayofyear[subs], agdat2$x[,1][subs], type = "l", xlim = xlm, ylim = ylm, axes = F, ylab = "", xlab = "",
           lwd = 1.5, col=adjustcolor(collst[1,ind], alpha.f = 0.8))
      title(paste(letters[ind+(yr-2011)*3],".", sep=""), line=-0.95, adj = 0.025)
      
      subs2 = which(agdat2$Plant_Species==i & agdat2$year==yr & agdat2$subplot=="high")
      subs2 = subs2[order(agdat2$dayofyear[subs2])]
      
      # remove gaps
      steps = diff(agdat2$dayofyear[subs2])
      tmp = which(steps > 1)
      for(j in 1:length(tmp)) {
        subs2[tmp[j]+0:1] = NA
      }
      
      lines(agdat2$dayofyear[subs2], agdat2$x[,1][subs2],
            lwd = 1.5, col=adjustcolor(collst[2,ind], alpha.f = 0.8))
      axis(2)
      axis(1, c("May", "June", "July", "Aug", "Sept", "Oct", "Nov"), at=dayofyearlookup[5:11])
      box()
      
      # identify significance of differences at p = 0.05
      tempdiff = agdat2$x[,1][subs]-agdat2$x[,1][subs2][match(agdat2$dayofyear[subs], agdat2$dayofyear[subs2])]
      jointsd = sqrt(agdat2$x[,2][subs]^2+agdat2$x[,2][subs2][match(agdat2$dayofyear[subs], agdat2$dayofyear[subs2])]^2)/sqrt(unique(agdat2$x[,3][subs])[1]-1)
      pv = pnorm(tempdiff, 0, jointsd)
      
      tmp=which(is.finite(pv) & pv < 0.05)
      k = 1
      while(k < length(tmp)) {
        if(is.na(subs[tmp[k]])) {
          k = k+1
        } else {
          if(is.na(subs[tmp[k+1]]) || (tmp[k+1]!=(tmp[k]+1))) {
            #abline(v=agdat2$dayofyear[subs][tmp[k]], col=adjustcolor("black", 0.1))
            k = k+1
          } else {
            startk = k
            while((tmp[k+1]==(tmp[k]+1)) & (k < length(tmp))) {
              k = k+1
              endk = k
            }
            polygon(c(agdat2$dayofyear[subs][tmp[c(startk,startk,endk,endk)]]),
                    c(-10, 100, 100, -10),
                    col=adjustcolor("black", 0.1), border = NA)
            k = k+1
          }
        }
      }
      if(ind==1) {
        mtext(yr, 2, line=2.2)
      }
      if(yr==2011) {
        mtext(paste("Sown Richness =", i), 3, line=0.2)
      }
    }
  }
  mtext("Minimum Daily Temperature, °C", side = 2, outer=TRUE, line=2.2)
dev.off()

