###
#Use Function
###

#Apply propagationPred to top performer (G2_b27_1 from UTT2)

OUTdat <- propagationPred(n=5, ns=15, r=0.90, s=0.45, g=4, date_user=mdy("10/28/2020")-weeks(8), 
                          growth_time_g1=weeks(8), rooting_time_g1=weeks(2), growth_time_g=weeks(8), 
                          rooting_time_g=weeks(3))

OUTdat <- OUTdat[1:10,]
###~~~
#Plot data
###~~~
#pdf("Schedule_sagebrush_propagation.pdf")
#Populate plot
Type <- levels(OUTdat$Type)[c(2,1,3)]
colSeg <- c("black","green","blue")

#Create empty plot
plot(x=c(as.Date(OUTdat$Date_start)), y= c(as.numeric(as.vector(OUTdat$N_plant_start))), 
     xlim= c(min(as.Date(OUTdat$Date_start)), max(as.Date(OUTdat$Date_end))+20), 
     ylim=c(0, max(as.numeric(as.vector(OUTdat$N_plant_start)))+50), type='n', xlab="Time (days)", ylab = "Number of Plantlets") 

#Add rectangle showing cutting period for generation 2 of propagation
rect(xleft = as.Date(OUTdat$Date_start[5])-days(1), xright = as.Date(OUTdat$Date_start[5])+days(1), ybottom = as.numeric(as.vector(OUTdat$N_plant_start[5])), ytop = as.numeric(as.vector(OUTdat$N_plant_end[5])), col=alpha("green", alpha = 0.4))
text(x=as.Date(OUTdat$Date_start[5])+days(8), y=40, srt=90, label=paste("Cutting period" , "(3 days)", sep='\n'), col='black', adj=0, cex=.6)

#Add rectangle showing cutting period for generation 3 of propagation
rect(xleft = as.Date(OUTdat$Date_start[8])-days(4), xright = as.Date(OUTdat$Date_start[8])+days(14), ybottom = as.numeric(as.vector(OUTdat$N_plant_start[8])), ytop = as.numeric(as.vector(OUTdat$N_plant_end[8])), col=alpha("green", alpha = 0.4))
text(x=as.Date(OUTdat$Date_start[8])+days(9), y=900, srt=90, label=paste("Cutting period (18 days): ", as.Date(OUTdat$Date_start[8])-days(4), " to ", as.Date(OUTdat$Date_start[8])+days(14), sep=''), col='black', adj=0, cex=.6)

#Add rectangle showing biomass production period
rect(xleft = as.Date(OUTdat$Date_end[10])-days(4), xright = as.Date(OUTdat$Date_end[10])+days(14), ybottom = as.numeric(as.vector(OUTdat$N_plant_end[10])), ytop = as.numeric(as.vector(OUTdat$N_plant_end[8])), col=alpha("black", alpha = 0.4))
text(x=as.Date(OUTdat$Date_end[10])+days(9), y=1120, srt=90, label=paste("Biomass production (18 days): ", as.Date(OUTdat$Date_end[10])-days(4), " to ", as.Date(OUTdat$Date_end[10])+days(14), sep=''), col='white', adj=0, cex=.6)

#Add arrows for growth media prep. for generation 3 of propagation
arrows(x0=as.Date(OUTdat$Date_start[7])-days(8), x1=as.Date(OUTdat$Date_start[7])-days(4), y0=as.numeric(as.vector(OUTdat$N_plant_start[6]))+30, y1=as.numeric(as.vector(OUTdat$N_plant_start[6]))+30, length=.03,code = 3)
text(x=as.Date(OUTdat$Date_start[7])-days(6), y=as.numeric(as.vector(OUTdat$N_plant_start[6]))+80, srt=40, label=paste("Growth media prep. (4 days): ", paste("starts on ", as.Date(OUTdat$Date_start[7])-days(8), sep=''), sep='\n'), col='black', adj=0, cex=.6)

#Add arrows for growth media prep. for generation 4 of propagation
arrows(x0=as.Date(OUTdat$Date_end[9])-days(24), x1=as.Date(OUTdat$Date_end[9])-days(4), y0=as.numeric(as.vector(OUTdat$N_plant_start[8]))-40, y1=as.numeric(as.vector(OUTdat$N_plant_start[8]))-40, length=.05,code = 3)
text(x=as.Date(OUTdat$Date_end[9])-days(24), y=as.numeric(as.vector(OUTdat$N_plant_start[8]))-150, srt=0, label=paste("Growth media prep. (20 days): ", paste(as.Date(OUTdat$Date_end[9])-days(24), " to ", as.Date(OUTdat$Date_end[9])-days(4), sep=''), sep='\n'), col='black', adj=0, cex=.6)

#Add segments
for(i in 1:length(Type)){
  tmp <- subset(OUTdat, OUTdat$Type == Type[i])
  segments(x0=as.Date(tmp$Date_start), x1=as.Date(tmp$Date_end), y0=as.numeric(as.vector(tmp$N_plant_start)), y1=as.numeric(as.vector(tmp$N_plant_end)), col=colSeg[i], lwd=1.5)
  points(x=as.Date(tmp$Date_start), y=as.numeric(as.vector(tmp$N_plant_start)), col=colSeg[i], pch=16, cex=.3)
  points(x=as.Date(tmp$Date_end), y=as.numeric(as.vector(tmp$N_plant_end)), col=colSeg[i], pch=16, cex=.3)
  if(Type[i] == "Cutting"){
    text(x=as.Date(tmp$Date_start), y=as.numeric(as.vector(tmp$N_plant_end))+20, labels = as.Date(tmp$Date_start), srt=40, adj=0.5, cex=.6)
  }
  if(Type[i] == "Rooting"){
    text(x=as.Date(tmp$Date_end), y=as.numeric(as.vector(tmp$N_plant_end))+8, labels = as.Date(tmp$Date_end), srt=40, adj=0, cex=.6)
  }
}

#Add when plantlets are ready for sequencing
text(x=as.Date(OUTdat$Date_end[10])+2, y=as.numeric(as.vector(OUTdat$N_plant_end[10]))+15, labels = paste("", as.Date(OUTdat$Date_end[10]), sep=" "), srt=90, adj=0, cex=.6)

#Add objective of 800 plantlets
abline(h=800, col='grey')
text(x=as.Date(OUTdat$Date_start[1])-9, y=840, labels="Number of Plantlets for Genome Sequencing (N: 800 = 120 gr.)", adj=0, col="grey", cex=.65)

#Add cushion/safety net to keep line in vitro
arrows(x0=as.Date(OUTdat$Date_end[10]), x1=as.Date(OUTdat$Date_end[10]), y0=800, y1=as.numeric(as.vector(OUTdat$N_plant_end[10])), length=.05,code = 3)
text(x=as.Date(OUTdat$Date_end[10])-5, y=950, srt=0, label=paste(paste("Cushion of ", as.numeric(as.vector(OUTdat$N_plant_end[10]))-800, sep=" "), "plantlets", sep='\n'), col='black', adj=1, cex=.6)

#Add legend
legend("topleft", legend = paste(Type, c(": 8 w. (45% survival)", ": 15 shoot tips per plantlet", ": 3 w. (90% response)"), sep=""), lty=1, col=colSeg, cex=.6)
dev.off()

