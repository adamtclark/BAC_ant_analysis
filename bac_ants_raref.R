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
asplst<-sort(unique(dabac$spCode))
asplst2<-sort(unique(dabac$reduced_spCode))

splookup<-unique(dabac[,c("spCode", "reduced_spCode")])

##############################
## make wide data
##############################
#BAC ants, alpha
wantbac<-data.frame(plotdat[rep(1:nrow(plotdat), each=nrow(timedat)),],
                    timedat[rep(1:nrow(timedat), nrow(plotdat)),],
                    richness=NA, abundance=NA)
wantbac<-data.frame(wantbac, matrix(nrow=nrow(wantbac), ncol=length(asplst), data=0))
metapos<-1:7
sppos<-max(metapos)+(1:length(asplst))
row.names(wantbac)<-NULL
colnames(wantbac)[-metapos]<-as.character(asplst)

windex<-paste(wantbac$PlotID, wantbac$Subplot, wantbac$year, wantbac$month)
uwindex<-sort(unique(windex))
lindex<-paste(dabac$PlotID, dabac$Subplot, dabac$year, dabac$month)

for(i in 1:length(uwindex)) {
  ps1<-which(lindex==uwindex[i])
  ps2<-which(windex==uwindex[i])
  
  tmp<-tapply(dabac$Count[ps1], dabac$spCode[ps1], function(x) sum(x, na.rm=T))
  tmp[is.na(tmp)]<-0
  
  wantbac[ps2,-metapos]<-tmp
}
wantbac$richness<-rowSums(wantbac[,-metapos]>0)
wantbac$abundance<-rowSums(wantbac[,-metapos])

##############################
## rarefaction
##############################
tmp<-wantbac[wantbac$Plant_Species<=16,]
nsmps<-table(unique(tmp[,c("PlotID", "Subplot", "Plant_Species")])[,-1])

niter<-1e4
ararf<-NULL

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


if(FALSE) {
  set.seed(201504)
  
  ararf<-data.frame(Plant_Species=NA, month=NA,
             Sublot=NA,
             smpnum=NA,
             richness=NA,
             shannon=NA,
             abundance=rep(NA, sum(nsmps)*2*niter),
             iter=NA,stringsAsFactors = FALSE)
  
  mm<-1
  for(i in c(sort(unique(tmp$Plant_Species)))) {
    for(j in c(sort(unique(tmp$month)))) {
      psC<-which(tmp$Plant_Species==i & tmp$month==j & tmp$Subplot=="C")
      psH<-which(tmp$Plant_Species==i & tmp$month==j & tmp$Subplot=="H")
      
      pltC<-sort(unique(tmp$PlotID[psC]))
      pltH<-sort(unique(tmp$PlotID[psH]))
      
      #lookup table for faster run below
      matchmatC<-matrix(nrow=length(pltC), ncol=nrow(tmp), data=0)
      matchmatH<-matrix(nrow=length(pltH), ncol=nrow(tmp), data=0)
      
      for(k in 1:nrow(matchmatC)) {
        matchmatC[k,]<-(tmp$PlotID==pltC[k] & tmp$Subplot=="C" & tmp$month==j)
      }
      for(k in 1:nrow(matchmatH)) {
        matchmatH[k,]<-(tmp$PlotID==pltH[k] & tmp$Subplot=="H" & tmp$month==j)
      }
      
      nstp<-min(c(length(pltC), length(pltH)))
      
      for(k in 1:niter) {
        richCtmp<-numeric(nstp)
        richHtmp<-numeric(nstp)
        abundCtmp<-numeric(nstp)
        abundHtmp<-numeric(nstp)
        shannonCtmp<-numeric(nstp)
        shannonHtmp<-numeric(nstp)
        
        for(l in 1:nstp) {
          smpCtmp<-sample(length(pltC),l)
          #use same plots for each comparison
          smpHtmp<-smpCtmp#sample(length(pltH),l)
          
          divC<-shannonfun(colSums(tmp[as.logical(colSums(matchmatC[smpCtmp,,drop=FALSE])),sppos,drop=FALSE]>0))
          divH<-shannonfun(colSums(tmp[as.logical(colSums(matchmatH[smpHtmp,,drop=FALSE])),sppos,drop=FALSE]>0))
          
          richCtmp[l]<-divC[2]
          richHtmp[l]<-divH[2]
          
          shannonCtmp[l]<-divC[1]
          shannonHtmp[l]<-divH[1]
          
          abundCtmp[l]<-sum(tmp[as.logical(colSums(matchmatC[smpCtmp,,drop=FALSE])),sppos])
          abundHtmp[l]<-sum(tmp[as.logical(colSums(matchmatH[smpHtmp,,drop=FALSE])),sppos])
        }
        mmnew<-mm+length(abundCtmp)+length(abundHtmp)
        
        
        
        ararf[mm:(mmnew-1),]<-data.frame(Plant_Species=i, month=j,
                  Sublot=rep(c("C","H"),each=nstp),
                  smpnum=1:nstp,
                  richness=c(richCtmp,richHtmp),
                  shannon=c(shannonCtmp,shannonHtmp),
                  abundance=c(abundCtmp,abundHtmp),
                  iter=k)
        
        mm<-mmnew
        
        if(k/100==floor(k/100)) {
          print(k/niter)
        }
      }
    }
    print(i)
  }
  save(list = "ararf", file = "datout/raref_tmp.rda", version = 2)
} else {
  load("datout/raref_tmp.rda")
}


## Rarefaction
rarf_ag<-with(ararf,
              aggregate(cbind(richness=richness, shannon=shannon, abundance=abundance),
                        list(Plant_Species=Plant_Species, month=month, Sublot=Sublot, smpnum=smpnum),
                        function(x) quantile(x, pnorm(-2:2))))
collst<-cbind(c("lightblue3", "cornflowerblue", "darkblue", "dodgerblue4"),
              c("lightpink2", "darkorange1", "darkred", "firebrick4"))
collst2<-c("dodgerblue4", "firebrick4")
ylm<-c(2,10); xlm<-c(1,14)

pdf("figures/raref.pdf", width=6, height=5, colormodel = "cmyk", useDingbats = FALSE)
  rarf_ag$divmetric<-rarf_ag$shannon
  
  par(mfcol=c(2,2), mar=c(2,2,1,1), oma=c(2.5,2.5,1,0))
  m<-1
  for(k in sort(unique(rarf_ag$month))) {
    for(j in 1:length(unique(rarf_ag$Sublot))) {
      subplot_type = c("C" ,"H")[j]
      
      xlmuse<-xlm
      ylmuse<-ylm
      
      if(m==1) {
        mn<-"August"
      } else if(m==3) {
        mn<-"September"
      } else {
        mn<-""
      }
      plot(1,1, ylim=ylmuse, xlim=xlmuse, xlab="", ylab="", type="n")
      title(mn, line=0.5, xpd=NA)
      abline(h=seq(ylmuse[1],ylmuse[2],by=2), v=seq(xlmuse[1]-1,xlmuse[2],by=2), col=adjustcolor("grey", alpha.f = 0.5))
      title(paste(letters[m], ".", sep=""), outer = FALSE, adj=0.02, line=-1)
      m<-m+1
      
      n<-1
      for(i in sort(unique(rarf_ag$Plant_Species))) {
        ps<-which(rarf_ag$Plant_Species==i & rarf_ag$Sublot==subplot_type & rarf_ag$month==k)
        
        lines(rarf_ag$smpnum[ps], rarf_ag$divmetric[ps,3], lwd=1.5, col=collst[n,(j==2)+1], lty=5-n)
        polygon(c(rarf_ag$smpnum[ps], rev(rarf_ag$smpnum[ps])),
                c(rarf_ag$divmetric[ps,2],rev(rarf_ag$divmetric[ps,4])),
                col=adjustcolor(collst[n,(j==2)+1], alpha.f = 0.2), border=adjustcolor(collst[n,(j=="H")+1], alpha.f = 0.5))
        n<-n+1
      }
      if(m<=3) {
        legend("bottomright", legend = as.character(c(1,4,16)), lty=4:2, lwd=1.5, col=collst[,(j==2)+1], bty="n", title = "Sown Richness")
      }
    }
  }
  
  mtext("Number of Samples", 1, outer=TRUE, line=0.8)
  mtext(expression(paste("Ant Species Diversity (e"^H, ")")), 2, outer=TRUE, line=0.2)
dev.off()



pdf("figures/raref_by_plantdiv.pdf", width=6, height=7.5, colormodel = "cmyk", useDingbats = FALSE)
  rarf_ag$divmetric<-rarf_ag$shannon
  
  par(mfcol=c(3,2), mar=c(2,2,1,1), oma=c(2.5,4,1,0))
  m<-1
  for(k in sort(unique(rarf_ag$month))) {
    n<-1
    for(i in sort(unique(rarf_ag$Plant_Species))) {
      xlmuse<-xlm
      ylmuse<-ylm
      
      if(m==1) {
        mn<-"August"
      } else if(m==4) {
        mn<-"September"
      } else {
        mn<-""
      }
      
      plot(1,1, ylim=ylmuse, xlim=xlmuse, xlab="", ylab="", type="n")
      title(mn, line=0.5, xpd=NA)
      abline(h=seq(ylmuse[1],ylmuse[2],by=2), v=seq(xlmuse[1]-1,xlmuse[2],by=2), col=adjustcolor("grey", alpha.f = 0.5))
      title(paste(letters[m], ".", sep=""), outer = FALSE, adj=0.02, line=-1)
      if(m%in%c(1,2,3)) {
        mtext(side = 2, text = paste("Sown Richness =", i), line = 2.4)
      }
      m<-m+1
      
      
      for(j in 1:length(unique(rarf_ag$Sublot))) {
        subplot_type = c("C" ,"H")[j]
        ps<-which(rarf_ag$Plant_Species==i & rarf_ag$Sublot==subplot_type & rarf_ag$month==k)
        
        lines(rarf_ag$smpnum[ps], rarf_ag$divmetric[ps,3], lwd=1.5, col=collst[n,(j==2)+1], lty=5-n)
        polygon(c(rarf_ag$smpnum[ps], rev(rarf_ag$smpnum[ps])),
                c(rarf_ag$divmetric[ps,2],rev(rarf_ag$divmetric[ps,4])),
                col=adjustcolor(collst[n,(j==2)+1], alpha.f = 0.2), border=adjustcolor(collst[n,(j=="H")+1], alpha.f = 0.5))
      }
      #if(m<=3) {
      #  legend("bottomright", legend = as.character(c(1,4,16)), lty=4:2, lwd=1.5, col=collst[,(j==2)+1], bty="n", title = "Sown Richness")
      #}
      n<-n+1
    }
  }
  
  mtext("Number of Samples", 1, outer=TRUE, line=0.8)
  mtext(expression(paste("Ant Species Diversity (e"^H, ")")), 2, outer=TRUE, line=2.2)
dev.off()

# rarefaction results for richness
pdf("figures/raref_by_plantdiv_rich.pdf", width=6, height=7.5, colormodel = "cmyk", useDingbats = FALSE)
  rarf_ag$divmetric<-rarf_ag$richness
  
  par(mfcol=c(3,2), mar=c(2,2,1,1), oma=c(2.5,4,1,0))
  m<-1
  for(k in sort(unique(rarf_ag$month))) {
    n<-1
    for(i in sort(unique(rarf_ag$Plant_Species))) {
      xlmuse<-xlm
      ylmuse<-c(2, 15)
      
      if(m==1) {
        mn<-"August"
      } else if(m==4) {
        mn<-"September"
      } else {
        mn<-""
      }
      
      plot(1,1, ylim=ylmuse, xlim=xlmuse, xlab="", ylab="", type="n")
      title(mn, line=0.5, xpd=NA)
      abline(h=seq(ylmuse[1],ylmuse[2],by=2), v=seq(xlmuse[1]-1,xlmuse[2],by=2), col=adjustcolor("grey", alpha.f = 0.5))
      title(paste(letters[m], ".", sep=""), outer = FALSE, adj=0.02, line=-1)
      if(m%in%c(1,2,3)) {
        mtext(side = 2, text = paste("Sown Richness =", i), line = 2.4)
      }
      m<-m+1
      
      
      for(j in 1:length(unique(rarf_ag$Sublot))) {
        subplot_type = c("C" ,"H")[j]
        ps<-which(rarf_ag$Plant_Species==i & rarf_ag$Sublot==subplot_type & rarf_ag$month==k)
        
        lines(rarf_ag$smpnum[ps], rarf_ag$divmetric[ps,3], lwd=1.5, col=collst[n,(j==2)+1], lty=5-n)
        polygon(c(rarf_ag$smpnum[ps], rev(rarf_ag$smpnum[ps])),
                c(rarf_ag$divmetric[ps,2],rev(rarf_ag$divmetric[ps,4])),
                col=adjustcolor(collst[n,(j==2)+1], alpha.f = 0.2), border=adjustcolor(collst[n,(j=="H")+1], alpha.f = 0.5))
      }
      #if(m<=3) {
      #  legend("bottomright", legend = as.character(c(1,4,16)), lty=4:2, lwd=1.5, col=collst[,(j==2)+1], bty="n", title = "Sown Richness")
      #}
      n<-n+1
    }
  }
  
  mtext("Number of Samples", 1, outer=TRUE, line=0.8)
  mtext(expression(paste("Ant Species Richess")), 2, outer=TRUE, line=2.2)
dev.off()





ararf$plrich<-(ararf$Plant_Species-mean(c(1,4,16)))
ararf$sbplt<-ararf$Sublot#c("C", "H")[ararf$Sublot]
ararf$mnt<-as.factor(ararf$month)

if(FALSE) {
  div_out<-array(dim=c(length(unique(ararf$smpnum)),
                           length(unique(ararf$iter)),
                           8,2))
  
  for(j in sort(unique(ararf$iter))) {
    for(i in sort(unique(ararf$smpnum))) {
      ps<-which(ararf$smpnum==i & ararf$iter==j)
      if(i<=9) {
        div_out[i,j,,1]<-coef(lm((shannon)~-1+sbplt:mnt+plrich:sbplt:mnt, ararf[ps,]))
        div_out[i,j,,2]<-coef(lm((richness)~-1+sbplt:mnt+plrich:sbplt:mnt, ararf[ps,]))
      } else {
        div_out[i,j,1:4,1]<-coef(lm((shannon)~-1+sbplt:mnt, ararf[ps,]))
        div_out[i,j,1:4,2]<-coef(lm((richness)~-1+sbplt:mnt, ararf[ps,]))
      }
    }
    
    if(j/100 == floor(j/100)) {
      print(j/niter)
    }
  }
  save(list = "div_out", file = "datout/div_out_tmp.rda", version = 2)
} else {
  load("datout/div_out_tmp.rda")
}
ps<-which(ararf$smpnum==1 & ararf$iter==1)
cfnames<-names(coef(lm((richness)~-1+sbplt:mnt+plrich:sbplt:mnt, ararf[ps,])))

ag_div_out<-apply(div_out, c(1,3:4), function(x) quantile(x, pnorm(-1:1), na.rm=T))#cbind(mean(x,na.rm=T)-sd(x,na.rm=T),

## plot coefs
divps<-1
mxda<-9
mxdb<-9

pdf("figures/regression_out.pdf", width=7.5, height=5, colormodel = "cmyk", useDingbats = FALSE)
  #August, Intercept
  par(mfrow=c(2,2), mar=c(1,2,2,3), oma=c(2.5,2,0,1.5))
  plot(c(1,mxda), c(3,8.1), type="n", xlab="Number of Samples", ylab="Mean Shannon Diversity")
  abline(v=seq(0,12,by=2), h=seq(3,8,by=1), col=adjustcolor("grey", alpha.f = 0.5))
  title("August")
  title(paste("a", ".", sep=""), outer = FALSE, adj=0.02, line=-1)
  matlines(1:mxda, ag_div_out[2,1:mxda,1:2,divps], lty=c(1,2), col=collst2[c(1,2)], lwd=2)
  polygon(c(1:mxda, rev(1:mxda)),
          c(ag_div_out[1,1:mxda,1,divps], rev(ag_div_out[3,1:mxda,1,divps])),
          col=adjustcolor(collst2[1], alpha.f = 0.2), border = NA)
  polygon(c(1:mxda, rev(1:mxda)),
          c(ag_div_out[1,1:mxda,2,divps], rev(ag_div_out[3,1:mxda,2,divps])),
          col=adjustcolor(collst2[2], alpha.f = 0.2), border = NA)
  pv<-rowMeans((div_out[,,1,1]-div_out[,,2,1])>0)
  par(new=TRUE)
  plot(1:mxda, pmin(1-pv[1:mxda], pv[1:mxda]), type="b", xlab="", ylab="", axes=F, ylim=c(0.001,0.85), log="y")
  axis(4, las=2, at = c(0.001, 0.01, 0.05, 0.1, 0.5), labels = as.character(c(0.001, 0.01, 0.05, 0.1, 0.5)))
  legend("bottomright", c("Control", "Warmed"), lwd=1.5, lty=1:2, col=collst[4,], bty="n")
  
  #September, Intercept
  plot(c(1,mxda), c(3,8.1), type="n", xlab="Number of Samples", ylab="Mean Shannon Diversity")
  abline(v=seq(0,10,by=2), h=seq(3,8,by=1), col=adjustcolor("grey", alpha.f = 0.5))
  title(paste("b", ".", sep=""), outer = FALSE, adj=0.02, line=-1)
  
  matlines(1:mxda, ag_div_out[2,1:mxda,3:4,divps], lty=c(1,2), col=collst2[c(1,2)], lwd=2)
  title("September")
  polygon(c(1:mxda, rev(1:mxda)),
          c(ag_div_out[1,1:mxda,3,divps], rev(ag_div_out[3,1:mxda,3,divps])),
          col=adjustcolor(collst2[1], alpha.f = 0.2), border = NA)
  polygon(c(1:mxda, rev(1:mxda)),
          c(ag_div_out[1,1:mxda,4,divps], rev(ag_div_out[3,1:mxda,4,divps])),
          col=adjustcolor(collst2[2], alpha.f = 0.2), border = NA)
  pv<-rowMeans((div_out[,,3,1]-div_out[,,4,1])>0)
  par(new=TRUE)
  plot(1:mxda, pmin(1-pv[1:mxda], pv[1:mxda]), type="b", xlab="", ylab="", axes=F, ylim=c(0.001,0.85), log="y")
  axis(4, las=2, at = c(0.001, 0.01, 0.05, 0.1, 0.5), labels = as.character(c(0.001, 0.01, 0.05, 0.1, 0.5)))
  
  #August, Slope
  plot(c(1,mxdb), c(-0.1, 0.15), type="n", xlab="Number of Samples", ylab="Plant Richness Effect")
  abline(v=seq(0,10,by=2), h=seq(-0.1,0.15,by=0.05), col=adjustcolor("grey", alpha.f = 0.5))
  title(paste("c", ".", sep=""), outer = FALSE, adj=0.02, line=-1)
  
  matlines(1:mxdb, ag_div_out[2,1:mxdb,5:6,divps], lty=c(1,2), col=collst2[c(1,2)], lwd=2)
  polygon(c(1:mxdb, rev(1:mxdb)),
          c(ag_div_out[1,1:mxdb,5,divps], rev(ag_div_out[3,1:mxdb,5,divps])),
          col=adjustcolor(collst2[1], alpha.f = 0.2), border = NA)
  polygon(c(1:mxdb, rev(1:mxdb)),
          c(ag_div_out[1,1:mxdb,6,divps], rev(ag_div_out[3,1:mxdb,6,divps])),
          col=adjustcolor(collst2[2], alpha.f = 0.2), border = NA)
  abline(h=0, lty=3)
  pv<-rowMeans((div_out[,,5,1]-div_out[,,6,1])>0)
  par(new=TRUE)
  plot(1:mxdb, pmin(1-pv[1:mxdb], pv[1:mxdb]), type="b", xlab="", ylab="", axes=F, ylim=c(0.001,0.85), log="y")
  axis(4, las=2, at = c(0.001, 0.01, 0.05, 0.1, 0.5), labels = as.character(c(0.001, 0.01, 0.05, 0.1, 0.5)))
  
  
  #September, Slope
  plot(c(1,mxdb), c(-0.1, 0.15), type="n", xlab="Number of Samples", ylab="Plant Richness Effect")
  abline(v=seq(0,10,by=2), h=seq(-0.1,0.15,by=0.05), col=adjustcolor("grey", alpha.f = 0.5))
  title(paste("d", ".", sep=""), outer = FALSE, adj=0.02, line=-1)
  
  matlines(1:mxdb, ag_div_out[2,1:mxdb,7:8,divps], lty=c(1,2), col=collst2[c(1,2)], lwd=2)
  polygon(c(1:mxdb, rev(1:mxdb)),
          c(ag_div_out[1,1:mxdb,7,divps], rev(ag_div_out[3,1:mxdb,7,divps])),
          col=adjustcolor(collst2[1], alpha.f = 0.2), border = NA)
  polygon(c(1:mxdb, rev(1:mxdb)),
          c(ag_div_out[1,1:mxdb,8,divps], rev(ag_div_out[3,1:mxdb,8,divps])),
          col=adjustcolor(collst2[2], alpha.f = 0.2), border = NA)
  abline(h=0, lty=3)
  pv<-rowMeans((div_out[,,7,1]-div_out[,,8,1])>0)
  par(new=TRUE)
  plot(1:mxdb, pmin(1-pv[1:mxdb], pv[1:mxdb]), type="b", xlab="", ylab="", axes=F, ylim=c(0.001,0.85), log="y")
  axis(4, las=2, at = c(0.001, 0.01, 0.05, 0.1, 0.5), labels = as.character(c(0.001, 0.01, 0.05, 0.1, 0.5)))
  
  mtext("Number of Samples", 1, outer=TRUE, line=1.2)
  
  mtext(expression(paste("Mean Ant Diversity (e"^H,")")), 2, outer=TRUE, line=0.2, adj=0.87)
  mtext(expression(paste("Plant Richness Effect")), 2, outer=TRUE, line=0.22, adj = 0.1)
  
  mtext("Pr(|C-H|>0)", 4, outer=TRUE, line=0)
dev.off()



divps<-2
pdf("figures/regression_out_richness.pdf", width=7.5, height=5, colormodel = "cmyk", useDingbats = FALSE)
  #August, Intercept
  par(mfrow=c(2,2), mar=c(1,2,2,3), oma=c(2.5,2,0,1.5))
  plot(c(1,mxda), c(3,12), type="n", xlab="Number of Samples", ylab="Mean Shannon Diversity")
  abline(v=seq(0,12,by=2), h=seq(3,15,by=1), col=adjustcolor("grey", alpha.f = 0.5))
  title("August")
  title(paste("a", ".", sep=""), outer = FALSE, adj=0.02, line=-1)
  matlines(1:mxda, ag_div_out[2,1:mxda,1:2,divps], lty=c(1,2), col=collst2[c(1,2)], lwd=2)
  polygon(c(1:mxda, rev(1:mxda)),
          c(ag_div_out[1,1:mxda,1,divps], rev(ag_div_out[3,1:mxda,1,divps])),
          col=adjustcolor(collst2[1], alpha.f = 0.2), border = NA)
  polygon(c(1:mxda, rev(1:mxda)),
          c(ag_div_out[1,1:mxda,2,divps], rev(ag_div_out[3,1:mxda,2,divps])),
          col=adjustcolor(collst2[2], alpha.f = 0.2), border = NA)
  pv<-rowMeans((div_out[,,1,divps]-div_out[,,2,divps])>0)
  par(new=TRUE)
  plot(1:mxda, pmin(1-pv[1:mxda], pv[1:mxda]), type="b", xlab="", ylab="", axes=F, ylim=c(0.0001,0.85), log="y")
  axis(4, las=2, at = c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.5), labels = as.character(c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.5)))
  legend("bottomright", c("Control", "Warmed"), lwd=1.5, lty=1:2, col=collst[4,], bty="n")
  
  #September, Intercept
  plot(c(1,mxda), c(3,12), type="n", xlab="Number of Samples", ylab="Mean Shannon Diversity")
  abline(v=seq(0,10,by=2), h=seq(3,15,by=1), col=adjustcolor("grey", alpha.f = 0.5))
  title(paste("b", ".", sep=""), outer = FALSE, adj=0.02, line=-1)
  
  matlines(1:mxda, ag_div_out[2,1:mxda,3:4,divps], lty=c(1,2), col=collst2[c(1,2)], lwd=2)
  title("September")
  polygon(c(1:mxda, rev(1:mxda)),
          c(ag_div_out[1,1:mxda,3,divps], rev(ag_div_out[3,1:mxda,3,divps])),
          col=adjustcolor(collst2[1], alpha.f = 0.2), border = NA)
  polygon(c(1:mxda, rev(1:mxda)),
          c(ag_div_out[1,1:mxda,4,divps], rev(ag_div_out[3,1:mxda,4,divps])),
          col=adjustcolor(collst2[2], alpha.f = 0.2), border = NA)
  pv<-rowMeans((div_out[,,3,divps]-div_out[,,4,divps])>0)
  par(new=TRUE)
  plot(1:mxda, pmin(1-pv[1:mxda], pv[1:mxda]), type="b", xlab="", ylab="", axes=F, ylim=c(0.0001,0.85), log="y")
  axis(4, las=2, at = c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.5), labels = as.character(c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.5)))
  
  #August, Slope
  plot(c(1,mxdb), c(-0.1, 0.25), type="n", xlab="Number of Samples", ylab="Plant Richness Effect")
  abline(v=seq(0,10,by=2), h=seq(-0.1,0.25,by=0.05), col=adjustcolor("grey", alpha.f = 0.5))
  title(paste("c", ".", sep=""), outer = FALSE, adj=0.02, line=-1)
  
  matlines(1:mxdb, ag_div_out[2,1:mxdb,5:6,divps], lty=c(1,2), col=collst2[c(1,2)], lwd=2)
  polygon(c(1:mxdb, rev(1:mxdb)),
          c(ag_div_out[1,1:mxdb,5,divps], rev(ag_div_out[3,1:mxdb,5,divps])),
          col=adjustcolor(collst2[1], alpha.f = 0.2), border = NA)
  polygon(c(1:mxdb, rev(1:mxdb)),
          c(ag_div_out[1,1:mxdb,6,divps], rev(ag_div_out[3,1:mxdb,6,divps])),
          col=adjustcolor(collst2[2], alpha.f = 0.2), border = NA)
  abline(h=0, lty=3)
  pv<-rowMeans((div_out[,,5,divps]-div_out[,,6,divps])>0)
  par(new=TRUE)
  plot(1:mxdb, pmin(1-pv[1:mxdb], pv[1:mxdb]), type="b", xlab="", ylab="", axes=F, ylim=c(0.0001,0.85), log="y")
  axis(4, las=2, at = c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.5), labels = as.character(c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.5)))
  
  
  #September, Slope
  plot(c(1,mxdb), c(-0.1, 0.25), type="n", xlab="Number of Samples", ylab="Plant Richness Effect")
  abline(v=seq(0,10,by=2), h=seq(-0.1,0.25,by=0.05), col=adjustcolor("grey", alpha.f = 0.5))
  title(paste("d", ".", sep=""), outer = FALSE, adj=0.02, line=-1)
  
  matlines(1:mxdb, ag_div_out[2,1:mxdb,7:8,divps], lty=c(1,2), col=collst2[c(1,2)], lwd=2)
  polygon(c(1:mxdb, rev(1:mxdb)),
          c(ag_div_out[1,1:mxdb,7,divps], rev(ag_div_out[3,1:mxdb,7,divps])),
          col=adjustcolor(collst2[1], alpha.f = 0.2), border = NA)
  polygon(c(1:mxdb, rev(1:mxdb)),
          c(ag_div_out[1,1:mxdb,8,divps], rev(ag_div_out[3,1:mxdb,8,divps])),
          col=adjustcolor(collst2[2], alpha.f = 0.2), border = NA)
  abline(h=0, lty=3)
  pv<-rowMeans((div_out[,,7,divps]-div_out[,,8,divps])>0, na.rm=T)
  par(new=TRUE)
  plot(1:mxdb, pmax(pmin(1-pv[1:mxdb], pv[1:mxdb]), 1e-4), type="b", xlab="", ylab="", axes=F, ylim=c(0.0001,0.85), log="y")
  axis(4, las=2, at = c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.5), labels = as.character(c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.5)))
  
  mtext("Number of Samples", 1, outer=TRUE, line=1.2)
  
  mtext(expression(paste("Mean Ant Species Richness")), 2, outer=TRUE, line=0.2, adj=0.87)
  mtext(expression(paste("Plant Richness Effect")), 2, outer=TRUE, line=0.22, adj = 0.1)
  
  mtext("Pr(|C-H|>0)", 4, outer=TRUE, line=0)
dev.off()





# Plot raw trends
# Shannon
ylm = c(2,10)
xlm = c(1,16)
xdf = c(-0.2, 0.2)
divps = 1

pdf("figures/piecewise_effects_shannon.pdf", width=8, height=7.5, colormodel = "cmyk", useDingbats = FALSE)
  rarf_ag$divmetric<-rarf_ag$shannon
  
  par(mfcol=c(3,3), mar=c(2,2,1,1), oma=c(2.5,2.5,2.5,0))
  m<-1
  for(k in sort(unique(rarf_ag$month))) {
    xlmuse<-xlm
    ylmuse<-ylm
    
    if(k==8) {
      mn<-"August"
    } else if(k==9) {
      mn<-"September"
    } else {
      mn<-""
    }
    
    
    for(i in 1:9) {
      n<-1
      
      plot(1,1, ylim=ylmuse, xlim=xlmuse, xlab="", ylab="", type="n")
      title(paste("Number of Samples =", i), line=0.5, xpd=NA)
      abline(h=seq(ylmuse[1],ylmuse[2],by=2), v=seq(xlmuse[1]-1,xlmuse[2],by=2), col=adjustcolor("grey", alpha.f = 0.5))
      title(paste(letters[m], ".", sep=""), outer = FALSE, adj=0.02, line=-1)
      m<-m+1
      
      for(j in 1:length(unique(rarf_ag$Sublot))) {
        subplot_type = c("C" ,"H")[j]
        ps<-which(rarf_ag$smpnum ==i & rarf_ag$Sublot==subplot_type & rarf_ag$month==k)
        
        points(rarf_ag$Plant_Species[ps]+xdf[n], rarf_ag$divmetric[ps,3], col=collst[4,(j==2)+1])
        segments(rarf_ag$Plant_Species[ps]+xdf[n], rarf_ag$divmetric[ps,2], rarf_ag$Plant_Species[ps]+xdf[n], rarf_ag$divmetric[ps,4], col=collst[4,(j==2)+1])
        
        n<-n+1
      }
      
      # add regression lines
      b0_centered = ag_div_out[2,i,(k==9)*2+(1:2),divps]
      b1 = ag_div_out[2,i,(k==9)*2+(5:6),divps]
      b0 = b0_centered-mean(c(1,4,16))*b1
      
      abline(a = b0[1], b = b1[1], col = collst[4,1])
      abline(a = b0[2], b = b1[2], col = collst[4,2])
    }
    mtext("Sown Richness", 1, outer=TRUE, line=0.8)
    mtext(expression(paste("Ant Species Diversity (e"^H, ")")), 2, outer=TRUE, line=0.2)
    
    mtext(mn, 3, outer = TRUE, line = 0.8)
  }
dev.off()


# Plot raw trends
# Richness

ylm = c(4,12)
xlm = c(1,16)
xdf = c(-0.2, 0.2)
divps = 2

pdf("figures/piecewise_effects_richness.pdf", width=8, height=7.5, colormodel = "cmyk", useDingbats = FALSE)
  rarf_ag$divmetric<-rarf_ag$richness
  
  par(mfcol=c(3,3), mar=c(2,2,1,1), oma=c(2.5,2.5,2.5,0))
  m<-1
  for(k in sort(unique(rarf_ag$month))) {
    xlmuse<-xlm
    ylmuse<-ylm
    
    if(k==8) {
      mn<-"August"
    } else if(k==9) {
      mn<-"September"
    } else {
      mn<-""
    }
    
    
    for(i in 1:9) {
      n<-1
      
      plot(1,1, ylim=ylmuse, xlim=xlmuse, xlab="", ylab="", type="n")
      title(paste("Number of Samples =", i), line=0.5, xpd=NA)
      abline(h=seq(ylmuse[1],ylmuse[2],by=2), v=seq(xlmuse[1]-1,xlmuse[2],by=2), col=adjustcolor("grey", alpha.f = 0.5))
      title(paste(letters[m], ".", sep=""), outer = FALSE, adj=0.02, line=-1)
      m<-m+1
      
      for(j in 1:length(unique(rarf_ag$Sublot))) {
        subplot_type = c("C" ,"H")[j]
        ps<-which(rarf_ag$smpnum ==i & rarf_ag$Sublot==subplot_type & rarf_ag$month==k)
        
        points(rarf_ag$Plant_Species[ps]+xdf[n], rarf_ag$divmetric[ps,3], col=collst[4,(j==2)+1])
        segments(rarf_ag$Plant_Species[ps]+xdf[n], rarf_ag$divmetric[ps,2], rarf_ag$Plant_Species[ps]+xdf[n], rarf_ag$divmetric[ps,4], col=collst[4,(j==2)+1])
        
        n<-n+1
      }
      
      # add regression lines
      b0_centered = ag_div_out[2,i,(k==9)*2+(1:2),divps]
      b1 = ag_div_out[2,i,(k==9)*2+(5:6),divps]
      b0 = b0_centered-mean(c(1,4,16))*b1
      
      abline(a = b0[1], b = b1[1], col = collst[4,1])
      abline(a = b0[2], b = b1[2], col = collst[4,2])
    }
    mtext("Sown Richness", 1, outer=TRUE, line=0.8)
    mtext(expression(paste("Ant Species Diversity (e"^H, ")")), 2, outer=TRUE, line=0.2)
    
    mtext(mn, 3, outer = TRUE, line = 0.8)
  }
dev.off()




# calculate average incidence in 32 species plots
nplt = length(unique(wantbac$PlotID[wantbac$Plant_Species==32]))

tc = wantbac[wantbac$Plant_Species==32 & wantbac$Subplot=="C",8:22]>0
th = wantbac[wantbac$Plant_Species==32 & wantbac$Subplot=="H",8:22]>0

sum(colSums(tc)>0) - # total ant species in control plots
  sum(colSums(th)>0)  # total ant species in heated plots

(colSums(tc)>0) & (!(colSums(th)>0)) # species in c but not h
(!colSums(tc)>0) & ((colSums(th)>0)) # species in h but not c

