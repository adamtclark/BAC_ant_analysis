error
rm(list=ls())
setwd("~/Dropbox/Projects/013_Ant revision/code/")

##############################
## load data
##############################
dabac<-read.csv("../data/pitfall_ants_bac_140913.csv", stringsAsFactors = TRUE)
dabac$year<-NA
dabac$year[grep("/11", dabac$smpDate, fixed=T)]<-2011
dabac$year[grep("/2012", dabac$smpDate, fixed=T)]<-2012
dabac$month<-NA
dabac$month[grep("9/", dabac$smpDate, fixed=T)]<-9
dabac$month[grep("8/", dabac$smpDate, fixed=T)]<-8


##############################
## metadata
##############################
plotdat<-unique(dabac[,c("PlotID", "Subplot", "Plant_Species")])
plotdat<-plotdat[order(plotdat$Plant_Species, plotdat$PlotID, plotdat$Subplot),]
timedat<-data.frame(year=c(2011, 2012), month=rep(c(8,9), each=2))
asplst<-sort(unique(dabac$reduced_spCode))

##############################
## make wide data
##############################
#BAC ants, alpha
wantbac<-data.frame(plotdat[rep(1:nrow(plotdat), each=nrow(timedat)),],
                 timedat[rep(1:nrow(timedat), nrow(plotdat)),],
                 richness=NA, activity=NA)
wantbac<-data.frame(wantbac, matrix(nrow=nrow(wantbac), ncol=length(asplst), data=0))
metapos<-1:7
sppos<-8:20
row.names(wantbac)<-NULL
colnames(wantbac)[-metapos]<-as.character(asplst)

windex<-paste(wantbac$PlotID, wantbac$Subplot, wantbac$month,  wantbac$year)
uwindex<-sort(unique(windex))
lindex<-paste(dabac$PlotID, dabac$Subplot, dabac$month, dabac$year)

for(i in 1:length(uwindex)) {
  ps1<-which(lindex==uwindex[i])
  ps2<-which(windex==uwindex[i])
  
  tmp<-tapply(dabac$Count[ps1], dabac$reduced_spCode[ps1], function(x) sum(x, na.rm=T))
  tmp[is.na(tmp)]<-0
  
  wantbac[ps2,-metapos]<-tmp
}
wantbac$richness<-rowSums(wantbac[,-metapos]>0)
wantbac$activity<-rowSums(wantbac[,-metapos])

shannonfun<-function(x) {
  p<-x/sum(x,na.rm=T)
  psz<-which(p!=0)
  if(length(psz)>0) {
    sh<-exp(-sum(p[psz]*log(p[psz])))
  } else {
    sh<-0
  }
  rich<-length(psz)
  return(c(sh, rich))
}

wantbac$shannon<-apply(wantbac[,-metapos], 1, shannonfun)[1,]

require(lme4)
require(lmerTest)

mod<-lmer(log10(shannon+1)~log10(Plant_Species)*as.factor(month)*as.factor(Subplot)+(1|PlotID)+(1|year), wantbac[wantbac$Plant_Species<=16,])
summary(mod)
modstep<-get_model(step(mod))
summary(modstep)

##############################
## rarefaction
##############################
tmp<-wantbac[wantbac$Plant_Species<=16,]
nsmps<-table(unique(tmp[,c("PlotID", "Subplot", "Plant_Species")])[,-1])

niter<-10000
ararf<-NULL

if(FALSE) {
  set.seed(201504)
  for(i in c(sort(unique(tmp$Plant_Species)), 99)) {
    for(j in c(sort(unique(tmp$month)), 3)) {
      if(j==3) {
        if(i==99) {
          #group across plant diversity treatments
          psC<-which(tmp$Subplot=="C")
          psH<-which(tmp$Subplot=="H")
        } else {
          #group by month
          psC<-which(tmp$Plant_Species==i & tmp$Subplot=="C")
          psH<-which(tmp$Plant_Species==i & tmp$Subplot=="H")
        }
      } else {
        if(i==99) {
          #group across plant diversity treatments
          psC<-which(tmp$month==j & tmp$Subplot=="C")
          psH<-which(tmp$month==j & tmp$Subplot=="H")
        } else {
          psC<-which(tmp$Plant_Species==i & tmp$month==j & tmp$Subplot=="C")
          psH<-which(tmp$Plant_Species==i & tmp$month==j & tmp$Subplot=="H")
        }
      }
      
      if(i!=99) {
        nstp<-1
      } else {
        nstp<-1
      }
      
      for(l in 1:nstp) {
        if(i!=99 & (i!=1 & l>9)) {
          richC<-NA
          richH<-NA
          
          abundC<-NA
          abundH<-NA
          
          richCsd<-NA
          richHsd<-NA
          
          abundCsd<-NA
          abundHsd<-NA
        } else {
          shannonCtmp<-numeric(niter)
          shannonHtmp<-numeric(niter)
          
          richCtmp<-numeric(niter)
          richHtmp<-numeric(niter)
          abundCtmp<-numeric(niter)
          abundHtmp<-numeric(niter)
          
          for(k in 1:niter) {
            if(j==3 | i==99) {
              smpCtmp<-sample(unique(tmp$PlotID[psC]),l)
              smpHtmp<-sample(unique(tmp$PlotID[psH]),l)
              
              if(j==3) {
                smpC<-which(tmp$PlotID%in%smpCtmp & tmp$Subplot=="C")
                smpH<-which(tmp$PlotID%in%smpHtmp & tmp$Subplot=="H")
              } else {
                smpC<-which(tmp$PlotID%in%smpCtmp & tmp$Subplot=="C" & tmp$month==j)
                smpH<-which(tmp$PlotID%in%smpHtmp & tmp$Subplot=="H" & tmp$month==j)
              }
            } else {
              smpC<-sample(psC,l)
              smpH<-sample(psH,l)
            }
            
            shannonCtmp[k]<-shannonfun(colSums(tmp[smpC,sppos]))[1]
            shannonHtmp[k]<-shannonfun(colSums(tmp[smpH,sppos]))[1]
            
            richCtmp[k]<-sum(colSums(tmp[smpC,sppos])>0)
            richHtmp[k]<-sum(colSums(tmp[smpH,sppos])>0)
            
            abundCtmp[k]<-sum(colSums(tmp[smpC,sppos]))
            abundHtmp[k]<-sum(colSums(tmp[smpH,sppos]))
          }
          shannonC<-mean(shannonCtmp)
          shannonH<-mean(shannonHtmp)
          
          richC<-mean(richCtmp)
          richH<-mean(richHtmp)
          
          abundC<-mean(abundCtmp)
          abundH<-mean(abundHtmp)
          
          shannonCsd<-sd(shannonCtmp)
          shannonHsd<-sd(shannonHtmp)
          
          richCsd<-sd(richCtmp)
          richHsd<-sd(richHtmp)
          
          abundCsd<-sd(abundCtmp)
          abundHsd<-sd(abundHtmp)
        }
        
        ararf<-rbind(ararf,
                     data.frame(Plant_Species=i, month=j,
                                Sublot=c("C", "H"),
                                smpnum=l,
                                richness=c(richC, richH),
                                richness_sd=c(richCsd, richHsd),
                                activity=c(abundC, abundH),
                                activitysd=c(abundCsd, abundHsd),
                                shannon=c(shannonC, shannonH),
                                shannonsd=c(shannonCsd, shannonHsd)))
      }
    }
    print(i)
  }
  save(list = "ararf", file = "datout/raref_full.rda", version = 2)
} else {
  load("datout/raref_full.rda")
}

#pvalue_comparisons, activity
require(lme4)
p81<-summary(lmer(activity~Subplot+(1|year/PlotID), wantbac[wantbac$month==8 & wantbac$Plant_Species==1,]))$coefficients[2,5]
#*
p91<-summary(lmer(activity~Subplot+(1|year/PlotID), wantbac[wantbac$month==9 & wantbac$Plant_Species==1,]))$coefficients[2,5]

p84<-summary(lmer(activity~Subplot+(1|year/PlotID), wantbac[wantbac$month==8 & wantbac$Plant_Species==4,]))$coefficients[2,5]
p94<-summary(lmer(activity~Subplot+(1|year/PlotID), wantbac[wantbac$month==9 & wantbac$Plant_Species==4,]))$coefficients[2,5]

p816<-summary(lmer(activity~Subplot+(1|year/PlotID), wantbac[wantbac$month==8 & wantbac$Plant_Species==16,]))$coefficients[2,5]
p916<-summary(lmer(activity~Subplot+(1|PlotID), wantbac[wantbac$month==9 & wantbac$Plant_Species==16,]))$coefficients[2,5]

p8a<-summary(lmer(activity~Subplot+(1|year/PlotID), wantbac[wantbac$month==8 & wantbac$Plant_Species<=16,]))$coefficients[2,5]
#*
p9a<-summary(lmer(activity~Subplot+(1|year/PlotID), wantbac[wantbac$month==9 & wantbac$Plant_Species<=16,]))$coefficients[2,5]

#***
p1<-summary(lmer(activity~as.factor(month)+(1|year/PlotID), wantbac[wantbac$Plant_Species==1,]))$coefficients[2,5]
#***
p4<-summary(lmer(activity~as.factor(month)+(1|year/PlotID), wantbac[wantbac$Plant_Species==4,]))$coefficients[2,5]
#***
p16<-summary(lmer(activity~as.factor(month)+(1|year/PlotID), wantbac[wantbac$Plant_Species==16,]))$coefficients[2,5]
#***
pa<-summary(lmer(activity~as.factor(month)+(1|year/PlotID), wantbac[wantbac$Plant_Species<=16,]))$coefficients[2,5]


## plot activity by treatment
collst<-matrix(adjustcolor(cbind(c("lightblue3", "cornflowerblue", "darkblue", "dodgerblue4"),
              c("lightpink2", "darkorange1", "darkred", "firebrick4")), alpha.f = 0.8), nrow=4)

pdf("figures/activity.pdf", width=5, height=4, colormodel = "cmyk", useDingbats = FALSE)
  par(mfrow=c(1,1), mar=c(2,2,2,1), oma=c(2,2,0,0))
  plot(c(0.7,4.3), c(0,250), xlab="Sown Plant Richness", ylab="Ant Activity", type="n", axes=F)
  
  axis(1, at=1:4, c(as.character(c(1,4,16)), "All"))
  axis(2, las=2); box()
  
  abline(h=seq(0, 350, by=50), col=adjustcolor("grey", alpha.f = 0.5))
  abline(v=c(1.5, 2.5, 3.5), col=adjustcolor("grey", alpha.f = 0.5))
  
  mj<-1
  for(j in c(8, 9)) {
    padj<-0.05
    
    n<-1+c(-0.2, 0.2)[mj]
    nsp<-1
    for(i in c(1,4,16,99)) {
      sbsC<-which(ararf$Plant_Species==i & ararf$Sublot=="C" & ararf$month==j & ararf$smpnum==1)
      sbsH<-which(ararf$Plant_Species==i & ararf$Sublot=="H" & ararf$month==j & ararf$smpnum==1)
      
      if(i!=99) {
        nuse<-nsmps[1,nsp]
      } else {
        nuse<-rowSums(nsmps)[1]
      }
      
      segments(n-padj, ararf$activity[sbsC]+ararf$activitysd[sbsC]/sqrt(nuse),
               n-padj, ararf$activity[sbsC]-ararf$activitysd[sbsC]/sqrt(nuse),
               col=collst[4,1], lend=2, lwd=3.5)
      segments(n+padj, ararf$activity[sbsH]+ararf$activitysd[sbsH]/sqrt(nuse),
               n+padj, ararf$activity[sbsH]-ararf$activitysd[sbsH]/sqrt(nuse),
               col=collst[4,2], lend=2, lwd=3.5)
      
      segments(n-padj, ararf$activity[sbsC]+ararf$activitysd[sbsC]/sqrt(nuse)*2,
               n-padj, ararf$activity[sbsC]-ararf$activitysd[sbsC]/sqrt(nuse)*2,
               col=collst[4,1], lend=2, lwd=1.5, lty=mj)
      segments(n+padj, ararf$activity[sbsH]+ararf$activitysd[sbsH]/sqrt(nuse)*2,
               n+padj, ararf$activity[sbsH]-ararf$activitysd[sbsH]/sqrt(nuse)*2,
               col=collst[4,2], lend=2, lwd=1.5, lty=mj)
      
      points(c(n-padj, n+padj), c(ararf$activity[sbsC], ararf$activity[sbsH]), pch=16+(mj-1), col=collst[4,])
      
      n<-n+1
      nsp<-nsp+1
    }
    mj<-mj+1
  }
  legend("topleft", legend = as.character(c("Control", "Heated", "August", "September")), fill = c(collst[4,], NA, NA), border=c(1,1,NA,NA), lty=c(NA,NA,1,2), pch=c(NA,NA,16,17), lwd=1.5, col=c(collst[4,], 1,1), bty="n")
  
  mtext("Sown Species Richness", 1, outer=TRUE, line=0.8)
  mtext("Ant Activity", 2, outer=TRUE, line=0.8)
  
  text(1.2, 75, "*", pos = 3)
  text(4.2, 100, "*", pos = 3)
  
dev.off()










#pvalue_comparisons, shannon
require(lme4)
p81<-summary(lmer(shannon~Subplot+(1|year/PlotID), wantbac[wantbac$month==8 & wantbac$Plant_Species==1,]))$coefficients[2,5]
#*
p91<-summary(lmer(shannon~Subplot+(1|year/PlotID), wantbac[wantbac$month==9 & wantbac$Plant_Species==1,]))$coefficients[2,5]

p84<-summary(lmer(shannon~Subplot+(1|PlotID), wantbac[wantbac$month==8 & wantbac$Plant_Species==4,]))$coefficients[2,5]
p94<-summary(lmer(shannon~Subplot+(1|PlotID), wantbac[wantbac$month==9 & wantbac$Plant_Species==4,]))$coefficients[2,5]

p816<-summary(lmer(shannon~Subplot+(1|PlotID), wantbac[wantbac$month==8 & wantbac$Plant_Species==16,]))$coefficients[2,5]
p916<-summary(lm(shannon~Subplot, wantbac[wantbac$month==9 & wantbac$Plant_Species==16,]))$coefficients[2,4]

p8a<-summary(lmer(shannon~Subplot+(1|year/PlotID), wantbac[wantbac$month==8 & wantbac$Plant_Species<=16,]))$coefficients[2,5]

p9a<-summary(lmer(shannon~Subplot+(1|year/PlotID), wantbac[wantbac$month==9 & wantbac$Plant_Species<=16,]))$coefficients[2,5]

#*
p1<-summary(lmer(shannon~as.factor(month)+(1|year/PlotID), wantbac[wantbac$Plant_Species==1,]))$coefficients[2,5]

p4<-summary(lmer(shannon~as.factor(month)+(1|PlotID), wantbac[wantbac$Plant_Species==4,]))$coefficients[2,5]

p16<-summary(lmer(shannon~as.factor(month)+(1|PlotID), wantbac[wantbac$Plant_Species==16,]))$coefficients[2,5]
#*
pa<-summary(lmer(shannon~as.factor(month)+(1|year/PlotID), wantbac[wantbac$Plant_Species<=16,]))$coefficients[2,5]


## plot activity by treatment
collst<-matrix(adjustcolor(cbind(c("lightblue3", "cornflowerblue", "darkblue", "dodgerblue4"),
                                 c("lightpink2", "darkorange1", "darkred", "firebrick4")), alpha.f = 0.8), nrow=4)

pdf("figures/shannon_richness_overall.pdf", width=5, height=8, colormodel = "cmyk", useDingbats = FALSE)
  par(mfrow=c(2,1), mar=c(2,2,2,1), oma=c(2,2,0,0))
  plot(c(0.7,4.3), c(0,4), xlab="Sown Plant Richness", ylab="Ant Diversity", type="n", axes=F)
  
  axis(1, at=1:4, c(as.character(c(1,4,16)), "All"))
  axis(2, las=2); box()
  
  abline(h=seq(0, 10, by=1), col=adjustcolor("grey", alpha.f = 0.5))
  abline(v=c(1.5, 2.5, 3.5), col=adjustcolor("grey", alpha.f = 0.5))
  
  mj<-1
  for(j in c(8, 9)) {
    padj<-0.05
    
    n<-1+c(-0.2, 0.2)[mj]
    nsp<-1
    for(i in c(1,4,16,99)) {
      sbsC<-which(ararf$Plant_Species==i & ararf$Sublot=="C" & ararf$month==j & ararf$smpnum==1)
      sbsH<-which(ararf$Plant_Species==i & ararf$Sublot=="H" & ararf$month==j & ararf$smpnum==1)
      
      if(i!=99) {
        nuse<-nsmps[1,nsp]
      } else {
        nuse<-rowSums(nsmps)[1]
      }
      
      segments(n-padj, ararf$shannon[sbsC]+ararf$shannonsd[sbsC]/sqrt(nuse),
               n-padj, ararf$shannon[sbsC]-ararf$shannonsd[sbsC]/sqrt(nuse),
               col=collst[4,1], lend=2, lwd=3.5)
      segments(n+padj, ararf$shannon[sbsH]+ararf$shannonsd[sbsH]/sqrt(nuse),
               n+padj, ararf$shannon[sbsH]-ararf$shannonsd[sbsH]/sqrt(nuse),
               col=collst[4,2], lend=2, lwd=3.5)
      
      segments(n-padj, ararf$shannon[sbsC]+ararf$shannonsd[sbsC]/sqrt(nuse)*2,
               n-padj, ararf$shannon[sbsC]-ararf$shannonsd[sbsC]/sqrt(nuse)*2,
               col=collst[4,1], lend=2, lwd=1.5, lty=mj)
      segments(n+padj, ararf$shannon[sbsH]+ararf$shannonsd[sbsH]/sqrt(nuse)*2,
               n+padj, ararf$shannon[sbsH]-ararf$shannonsd[sbsH]/sqrt(nuse)*2,
               col=collst[4,2], lend=2, lwd=1.5, lty=mj)
      
      points(c(n-padj, n+padj), c(ararf$shannon[sbsC], ararf$shannon[sbsH]), pch=16+(mj-1), col=collst[4,])
      
      n<-n+1
      nsp<-nsp+1
    }
    mj<-mj+1
  }
  legend("bottomright", legend = as.character(c("Control", "Heated", "August", "September")), fill = c(collst[4,], NA, NA), border=c(1,1,NA,NA), lty=c(NA,NA,1,2), pch=c(NA,NA,16,17), lwd=1.5, col=c(collst[4,], 1,1), bty="n")
  
  mtext("Sown Species Richness", 1, outer=TRUE, line=0.8)
  mtext(expression(paste("Ant Diversity (e"^H,")")), 2, outer=FALSE, line=1.6)
  
  text(1.2, 2.7, "*", pos = 3)
  








#pvalue_comparisons, richness
require(lme4)
p81<-summary(lmer(richness~Subplot+(1|year/PlotID), wantbac[wantbac$month==8 & wantbac$Plant_Species==1,]))$coefficients[2,5]

p91<-summary(lmer(richness~Subplot+(1|year/PlotID), wantbac[wantbac$month==9 & wantbac$Plant_Species==1,]))$coefficients[2,5]

p84<-summary(lmer(richness~Subplot+(1|year/PlotID), wantbac[wantbac$month==8 & wantbac$Plant_Species==4,]))$coefficients[2,5]
p94<-summary(lmer(richness~Subplot+(1|PlotID), wantbac[wantbac$month==9 & wantbac$Plant_Species==4,]))$coefficients[2,5]

p816<-summary(lmer(richness~Subplot+(1|PlotID), wantbac[wantbac$month==8 & wantbac$Plant_Species==16,]))$coefficients[2,5]
p916<-summary(lm(richness~Subplot, wantbac[wantbac$month==9 & wantbac$Plant_Species==16,]))$coefficients[2,4]

p8a<-summary(lmer(richness~Subplot+(1|year/PlotID), wantbac[wantbac$month==8 & wantbac$Plant_Species<=16,]))$coefficients[2,5]

p9a<-summary(lmer(richness~Subplot+(1|year/PlotID), wantbac[wantbac$month==9 & wantbac$Plant_Species<=16,]))$coefficients[2,5]

#***
p1<-summary(lmer(richness~as.factor(month)+(1|year/PlotID), wantbac[wantbac$Plant_Species==1,]))$coefficients[2,5]
#***
p4<-summary(lmer(richness~as.factor(month)+(1|year/PlotID), wantbac[wantbac$Plant_Species==4,]))$coefficients[2,5]
#***
p16<-summary(lmer(richness~as.factor(month)+(1|PlotID), wantbac[wantbac$Plant_Species==16,]))$coefficients[2,5]
#***
pa<-summary(lmer(richness~as.factor(month)+(1|year/PlotID), wantbac[wantbac$Plant_Species<=16,]))$coefficients[2,5]


## plot activity by treatment
collst<-matrix(adjustcolor(cbind(c("lightblue3", "cornflowerblue", "darkblue", "dodgerblue4"),
                                 c("lightpink2", "darkorange1", "darkred", "firebrick4")), alpha.f = 0.8), nrow=4)

  plot(c(0.7,4.3), c(0,7), xlab="Sown Plant Richness", ylab="Ant Richness", type="n", axes=F)
  
  axis(1, at=1:4, c(as.character(c(1,4,16)), "All"))
  axis(2, las=2); box()
  
  abline(h=seq(0, 10, by=1), col=adjustcolor("grey", alpha.f = 0.5))
  abline(v=c(1.5, 2.5, 3.5), col=adjustcolor("grey", alpha.f = 0.5))
  
  mj<-1
  for(j in c(8, 9)) {
    padj<-0.05
    
    n<-1+c(-0.2, 0.2)[mj]
    nsp<-1
    for(i in c(1,4,16,99)) {
      sbsC<-which(ararf$Plant_Species==i & ararf$Sublot=="C" & ararf$month==j & ararf$smpnum==1)
      sbsH<-which(ararf$Plant_Species==i & ararf$Sublot=="H" & ararf$month==j & ararf$smpnum==1)
      
      if(i!=99) {
        nuse<-nsmps[1,nsp]
      } else {
        nuse<-rowSums(nsmps)[1]
      }
      
      segments(n-padj, ararf$richness[sbsC]+ararf$richness_sd[sbsC]/sqrt(nuse),
               n-padj, ararf$richness[sbsC]-ararf$richness_sd[sbsC]/sqrt(nuse),
               col=collst[4,1], lend=2, lwd=3.5)
      segments(n+padj, ararf$richness[sbsH]+ararf$richness_sd[sbsH]/sqrt(nuse),
               n+padj, ararf$richness[sbsH]-ararf$richness_sd[sbsH]/sqrt(nuse),
               col=collst[4,2], lend=2, lwd=3.5)
      
      segments(n-padj, ararf$richness[sbsC]+ararf$richness_sd[sbsC]/sqrt(nuse)*2,
               n-padj, ararf$richness[sbsC]-ararf$richness_sd[sbsC]/sqrt(nuse)*2,
               col=collst[4,1], lend=2, lwd=1.5, lty=mj)
      segments(n+padj, ararf$richness[sbsH]+ararf$richness_sd[sbsH]/sqrt(nuse)*2,
               n+padj, ararf$richness[sbsH]-ararf$richness_sd[sbsH]/sqrt(nuse)*2,
               col=collst[4,2], lend=2, lwd=1.5, lty=mj)
      
      points(c(n-padj, n+padj), c(ararf$richness[sbsC], ararf$richness[sbsH]), pch=16+(mj-1), col=collst[4,])
      
      n<-n+1
      nsp<-nsp+1
    }
    mj<-mj+1
  }
  
  mtext(expression(paste("Ant Richness")), 2, outer=FALSE, line=1.6)

dev.off()
