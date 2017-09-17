require(rworldmap)
require(plotrix)
require(maps)
require(viridis)

# Read in data
dat=read.table('data.carrierFreq.rs140950220.tsv',header=T,sep='\t')

colordf=data.frame(matrix(ncol=3,nrow=3))
names(colordf)=c("min","max","color")
colordf$min=c(-0.01,0.0002,0.01)
colordf$max=c(0.0002,0.01, 1.01)
colordf$color=viridis(3)
colordf

dat$color=NA
for ( i in 1:nrow(dat)) {
  f=dat[i,]$carrier_Frq
  for (j in 1:nrow(colordf)) {
    if ( findInterval( f, c(colordf$min[j], colordf$max[j]) )==TRUE) { dat$color[i]=colordf$color[j]; print (paste(f, colordf$min[j], colordf$max[j], sep=',')) } 
    else {j = j + 1}
  }
}

dat[dat$Count==1& dat$MAF<0.01,]$color="#440154FF"

sizedf=data.frame(matrix(ncol=3,nrow=5))
names(sizedf)=c('min','max','r')
sizedf$min=c(10,101,201,501,2001)
sizedf$max=c(100,200,500,2000,50000)
sizedf$r=c(1.2,1.7,2.5,3.4,4.4)
sizedf

dat$size=NA
for ( i in 1:nrow(dat)) {
  f=dat[i,]$N.Obs
  for (j in 1:nrow(sizedf)) {
    if ( findInterval( f, c(sizedf$min[j], sizedf$max[j]) )==TRUE) { dat$size[i]=sizedf$r[j]; print (paste(f, sizedf$min[j], sizedf$max[j], sep=',')) } 
    else {j = j + 1}
  }
}

# World Map 

pdf(paste('worldMap.rs140950220.pdf',sep=''),width=10,height=5)

print(map("world", fill=T, col="grey80", bg="white", ylim=c(-60, 90),border="white",lwd=0.5))
  
for (i in 1:nrow(dat)) { 
  print(dat$Population)
  draw.circle(dat$lon[i], dat$lat[i],radius=dat$size[i],col=dat$color[i] ,border="white", lwd=0.8)
}

legend(150,80, legend=c("absent/rare","> 1 in 1000 ","> 1 in 100"), cex=0.5,box.col="grey30",bg="white",pch=21,pt.cex=c(1,1,1),pt.bg = c(viridis(3),'black'),pt.lwd = 0.5, title="Carrier frequency")
legend(150,55, legend=c("5-50","51-100","101-250","251-1000",">1000"), cex=0.5,box.col="grey30",bg="white",pch=21,pt.cex=sizedf$r*0.55,pt.bg = c('white'),pt.lwd = 0.5, title="Population size")

dev.off()

# Inset - Carribbean 

pdf('inset.rs140950220.pdf',width=5,height=5)
print(map("world", fill=T, col="grey80", bg="white", ylim=c(10, 28),xlim=c(-90,-60),border="white"))
for (i in 1:nrow(dat)) { 
  print(dat$Population)
  draw.circle(dat$lon[i], dat$lat[i],radius=dat$size[i]*0.3,col=dat$color[i] ,border="white", lwd=0.8)
}
dev.off()
