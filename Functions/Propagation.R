#In vitro propagation

#Load packages
library(lubridate)
library(scales)

###
##propagationPred: A function to predict propagation trends for clonal lines
###

#Arguments
# - n: number of plantlets to seed experiment
# - ns: number of shoot tips cut per plantlet
# - r: vector of rooting rate (0 to 1) at rooting phase 
# - s: vector of survival rate (0 to 1) at growth phase
# - g: number of generations (for propagation)
# - biobank: vector with propartion (0 to 1) of plantlets biobanked per generation (at end of growth) (e.g: c(0,0,0.4,0))
# - date_user: starting date for experiment (mdy from ludridate or any date format)
# - growth_time_g: vector with number of days/weeks for growth for g
# - biobanking_time: one value (either weeks() or days()) compatible with lubridate describing the duration of biobanking to generated biomass.
# - rooting_time_g: vector with number of days/weeks for rooting for g
# - biomass_plantlet: numerical value of biomass of a single plantlet.

#Output
# A dataframe with 12 columns:
# [1] "Generation"            "Type"                  "Date_start"            "Date_end"             
# [5] "N_plant_start"         "N_plant_end"           "N_vessels"             "Volume_media_litre"   
# [9] "N_Autoclave"           "Time_Autoclave (hrs)"  "Media_prep_time (hrs)" "Time_Cutting (hrs)"  
propagationPred <- function(n, ns, r, s, g, biobank, date_user, growth_time_g, biobanking_time, rooting_time_g, biomass_plantlet){
  ###~~~  
  #Create matrix to store data
  ###~~~
  OUT <- matrix(ncol=13, nrow=4*g)
  colnames(OUT) <- c("Generation","Type","Date_start","Date_end","N_plant_start", "N_plant_end", "N_vessels", "Volume_media_litre", "N_Autoclave", "Time_Autoclave (hrs)","Media_prep_time (hrs)", "Time_Cutting (hrs)", "Total_Biomass (gr)")
  #Populate generations and types
  OUT[,1] <- sort(rep(seq(from=1, to=g, by=1), 4))
  OUT[,2] <- rep(c("Growth","Biobanking","Cutting","Rooting"), g)
  
  #biobankMat <- matrix(ncol=6, nrow=g)
  #colnames(biobankMat) <- c("Generation","Rate","N_plantlets", "Date_start", "Date_end", "Age")
  #Populate generations and rate
  #biobankMat[,1] <- seq(from=1, to=g, by=1)
  #biobankMat[,2] <- biobank
  
  ###~~~
  #Populate matrix
  ###~~~
  for(i in 1:g){
    if(i == 1){
      ##GROWTH
      rowGrowth <- which(OUT[,1] == i & OUT[,2] == "Growth")
      #N plants start/end
      OUT[rowGrowth,5] <- n-(n*biobank[i])
      OUT[rowGrowth,6] <- n-(n*biobank[i])*s[i]
      #Date start
      OUT[rowGrowth,3] <- as.character(date_user)
      #Date end
      OUT[rowGrowth,4]  <- as.character(as.Date(OUT[rowGrowth,3]) + growth_time_g[i])
      #Number_magenta
      OUT[rowGrowth,7] <- n-(n*biobank[i])
      #Volume_media
      OUT[rowGrowth,8] <- ceiling(n*0.1)
      #N_Autoclave (60 boxes per autoclave)
      OUT[rowGrowth,9] <- ceiling(as.numeric(OUT[rowGrowth,7])/60)
      #Time_Autoclave
      OUT[rowGrowth,10] <- as.numeric(OUT[rowGrowth,9])*90/60
      #Media_prep_time: 1 hr for 60 boxes
      OUT[rowGrowth,11] <- ceiling(as.numeric(OUT[rowGrowth,7])/60)
      
      ##BIOBANKING
      rowBiobank <- which(OUT[,1] == i & OUT[,2] == "Biobanking")
      #Only populate columns if biobanking rate is > 0
      if(biobank[i] > 0){
        #N plants start/end
        OUT[rowBiobank,5] <- n*biobank[i]
        OUT[rowBiobank,6] <- floor((n*biobank[i])*s[i])
        #Date start
        OUT[rowBiobank,3] <- as.character(date_user)
        #Date end
        OUT[rowBiobank,4]  <- as.character(as.Date(OUT[rowBiobank,3]) + biobanking_time)
        #Number_magenta
        OUT[rowBiobank,7] <- n*biobank[i]
        #Volume_media
        OUT[rowBiobank,8] <- ceiling(n*0.1)
        #N_Autoclave (60 boxes per autoclave)
        OUT[rowBiobank,9] <- ceiling(as.numeric(OUT[rowBiobank,7])/60)
        #Time_Autoclave
        OUT[rowBiobank,10] <- as.numeric(OUT[rowBiobank,9])*90/60
        #Media_prep_time: 1 hr for 60 boxes
        OUT[rowBiobank,11] <- ceiling(as.numeric(OUT[rowBiobank,7])/60)
        #Total biomass
        OUT[rowBiobank,13] <- as.numeric(OUT[rowBiobank,6])*biomass_plantlet
      }
            
      ##CUTTING
      rowCut <- which(OUT[,1] == i & OUT[,2] == "Cutting")
      #N plants start
      OUT[rowCut,5] <- n
      #N plants end
      OUT[rowCut,6] <- n*ns
      #Date start
      OUT[rowCut,3] <- as.character(as.Date(OUT[rowGrowth,3]) + growth_time_g[i])
      #Date end
      OUT[rowCut,4]  <- as.character(as.Date(OUT[rowGrowth,3]) + growth_time_g[i])
      #Time_Cutting: 5 individuals per hr (max of 10 ind. per day)
      OUT[rowCut,12]  <-  as.numeric(OUT[rowCut,5])/5
      
      ##ROOTING
      rowRoot <- which(OUT[,1] == i & OUT[,2] == "Rooting")
      #N plants start
      OUT[rowRoot,5] <- n*ns
      #N plants end
      OUT[rowRoot,6] <- floor((n*ns)*r[i])
      #Date start
      OUT[rowRoot,3] <- as.character(as.Date(OUT[rowGrowth,3]) + growth_time_g[i])
      #Date end
      OUT[rowRoot,4]  <- as.character(as.Date(OUT[rowRoot,3]) + rooting_time_g[i])
      #Number_plates
      OUT[rowRoot,7] <- ceiling(n*ns/9)
      #Volume_media
      OUT[rowRoot,8] <- ceiling(as.numeric(OUT[rowRoot,7])*0.05)
      #N_Autoclave (1 litre per autoclave)
      OUT[rowRoot,9] <- as.numeric(OUT[rowRoot,8])
      #Time_Autoclave
      OUT[rowRoot,10] <- as.numeric(OUT[rowRoot,8])*120/60
      #Media_prep_time: 90 minutes for 1 liter
      OUT[rowRoot,11] <- as.numeric(OUT[rowRoot,8])*90/60
      
      #Number of plantlets biobanked
      #biobankMat[i,3] <- n*biobank[i]
      #biobankMat[i,4] <- as.character(as.Date(OUT[rowGrowth,5]))
      
      #Biomass
      #OUT[i,8] <- as.numeric(OUT[i,3])*bio
    }else{
      ##GROWTH
      rowGrowth <- which(OUT[,1] == i & OUT[,2] == "Growth")
      #N plants start/end
      OUT[rowGrowth,5] <- as.numeric(OUT[which(OUT[,1] == i-1 & OUT[,2] == "Rooting"),6])-floor(as.numeric(OUT[which(OUT[,1] == i-1 & OUT[,2] == "Rooting"),6])*biobank[i])
      OUT[rowGrowth,6] <- floor(as.numeric(OUT[rowGrowth,5])*s[i])
      #Date start
      OUT[rowGrowth,3] <- as.character(OUT[which(OUT[,1] == i-1 & OUT[,2] == "Rooting"),4])
      #Date end
      OUT[rowGrowth,4]  <- as.character(as.Date(OUT[rowGrowth,3]) + growth_time_g[i])
      #Number_magenta
      OUT[rowGrowth,7] <- OUT[rowGrowth,5]
      #Volume_media
      OUT[rowGrowth,8] <- ceiling(as.numeric(OUT[rowGrowth,5])*0.1)
      #N_Autoclave (60 boxes per autoclave)
      OUT[rowGrowth,9] <- ceiling(as.numeric(OUT[rowGrowth,7])/60)
      #Time_Autoclave
      OUT[rowGrowth,10] <- as.numeric(OUT[rowGrowth,9])*90/60
      #Media_prep_time: 1 hr for 60 boxes
      OUT[rowGrowth,11] <- ceiling(as.numeric(OUT[rowGrowth,7])/60)
      
      ##BIOBANKING
      rowBiobank <- which(OUT[,1] == i & OUT[,2] == "Biobanking")
      #Only populate columns if biobanking rate is > 0
      if(biobank[i] > 0) {
        #N plants start/end
        OUT[rowBiobank,5] <- floor(as.numeric(OUT[which(OUT[,1] == i-1 & OUT[,2] == "Rooting"),6])*biobank[i])
        OUT[rowBiobank,6] <- floor(as.numeric(OUT[rowBiobank,5])*s[i])
        #Date start
        OUT[rowBiobank,3] <- as.character(OUT[which(OUT[,1] == i-1 & OUT[,2] == "Rooting"),4])
        #Date end
        OUT[rowBiobank,4]  <- as.character(as.Date(OUT[rowBiobank,3]) + biobanking_time)
        #Number_magenta
        OUT[rowBiobank,7] <- OUT[rowBiobank,5]
        #Volume_media
        OUT[rowBiobank,8] <- ceiling(n*0.1)
        #N_Autoclave (60 boxes per autoclave)
        OUT[rowBiobank,9] <- ceiling(as.numeric(OUT[rowBiobank,7])/60)
        #Time_Autoclave
        OUT[rowBiobank,10] <- as.numeric(OUT[rowBiobank,9])*90/60
        #Media_prep_time: 1 hr for 60 boxes
        OUT[rowBiobank,11] <- ceiling(as.numeric(OUT[rowBiobank,7])/60)
        #Total biomass
        OUT[rowBiobank,13] <- floor(as.numeric(OUT[rowBiobank,5])*s[i])*biomass_plantlet
      }
          
      ##CUTTING
      rowCut <- which(OUT[,1] == i & OUT[,2] == "Cutting")
      #N plants start
      OUT[rowCut,5] <- OUT[rowGrowth,6]
      #N plants end
      OUT[rowCut,6] <- as.numeric(OUT[rowCut,5])*ns
      #Date start
      OUT[rowCut,3] <- as.character(as.Date(OUT[rowGrowth,3]) + growth_time_g[i])
      #Date end
      OUT[rowCut,4]  <- as.character(as.Date(OUT[rowGrowth,3]) + growth_time_g[i])
      #Time_Cutting: 5 individuals per hr (max of 10 ind. per day)
      OUT[rowCut,12]  <-  as.numeric(OUT[rowCut,5])/5
      
      ##ROOTING
      rowRoot <- which(OUT[,1] == i & OUT[,2] == "Rooting")
      #N plants start
      OUT[rowRoot,5] <- OUT[rowCut,6]
      #N plants end
      OUT[rowRoot,6] <- floor(as.numeric(OUT[rowRoot,5])*r[i])
      #Date start
      OUT[rowRoot,3] <- as.character(as.Date(OUT[rowGrowth,3]) + growth_time_g[i])
      #Date end
      OUT[rowRoot,4]  <- as.character(as.Date(OUT[rowRoot,3]) + rooting_time_g[i])
      #Number_plates
      OUT[rowRoot,7] <- ceiling(as.numeric(OUT[rowCut,6])/9)
      #Volume_media
      OUT[rowRoot,8] <- ceiling(as.numeric(OUT[rowRoot,7])*0.05)
      #N_Autoclave (1 litre per autoclave)
      OUT[rowRoot,9] <- as.numeric(OUT[rowRoot,8])
      #Time_Autoclave
      OUT[rowRoot,10] <- as.numeric(OUT[rowRoot,8])*120/60
      #Media_prep_time: 90 minutes for 1 liter
      OUT[rowRoot,11] <- as.numeric(OUT[rowRoot,8])*90/60
      
    }
  }
  
  ###~~~
  #Convert into a dataframe and tidy for project
  ###~~~
  OUTdat <- as.data.frame(OUT) 
  
  #Rm Biobanking rows with NAs
  OUTdat <- OUTdat[-which(OUTdat$Type == "Biobanking" & is.na(OUTdat$Date_start) == T),]
  
  #biobankMatOUT <- as.data.frame(biobankMat)
  #biobankMatOUT[,5] <- rep()
  return(OUTdat)
  #We are stopping experiment at generation 4 after growth period
  #OUTdat <- OUTdat[1:10,]
  #write.csv(OUTdat, file="Schedule_sagebrush_propagation.csv", row.names = F, quote=F)
}

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




#b24 11/13/2020 10 week old 
# 57 shoot tips
# 53 alive
# 4 weeks in rooting media
# biomass <- 1x 15-week old plantlet = 0.8gr
# How many plantlets to get at biomass with some wiggle room
#Apply propagationPred to top performer (G2_b24_1 from UTT2)

OUTdat <- propagationPred(n=6, ns=9.5, r=rep(0.93,4), s=rep(0.8,4), g=4, biobank=c(0,0,0.9,0), date_user=mdy("11/13/2020")-weeks(10), 
                          growth_time_g=c(weeks(10),weeks(8),weeks(8),weeks(8)), biobanking_time = weeks(15), rooting_time_g=c(weeks(4),weeks(3),weeks(3),weeks(3)), biomass_plantlet = 0.8)

OUTdat <- OUTdat[1:11,]


write.csv(OUTdat, file="Propagation_schedule_G2_b24_1.csv", row.names=F)

###~~~
#Plot data
###~~~
pdf("Schedule_sagebrush_propagation_G2_b24_1.pdf")
#Populate plot
Type <- levels(OUTdat$Type)[c(3,2,4,1)]
colSeg <- c("black","green","blue","orange")

#Create empty plot
plot(x=c(as.Date(OUTdat$Date_start)), y= c(as.numeric(as.vector(OUTdat$N_plant_start))), 
     xlim= c(min(as.Date(OUTdat$Date_start)), max(as.Date(OUTdat$Date_end))+20), 
     ylim=c(0, max(as.numeric(as.vector(OUTdat$N_plant_start)))+50), type='n', xlab="Time (days)", ylab = "Number of Plantlets") 

biomassTot <- 120
biomass_plantlet <- 0.8
#Add objective of N plantlets
abline(h=biomassTot/biomass_plantlet, col='grey')
text(x=as.Date(OUTdat$Date_start[1]), y=biomassTot/biomass_plantlet+10, labels=paste("Number of Plantlets for Project (N: ", biomassTot/biomass_plantlet, " = ", biomassTot, " gr. )", sep=""), adj=0, col="grey", cex=.65)

#Add segments
for(i in 1:length(Type)){
  tmp <- OUTdat[which(OUTdat$Type == Type[i]),]
  segments(x0=as.Date(tmp$Date_start), x1=as.Date(tmp$Date_end), y0=as.numeric(as.vector(tmp$N_plant_start)), y1=as.numeric(as.vector(tmp$N_plant_end)), col=colSeg[i], lwd=1.5)
  #Add points to better see beginning and end of phases
  points(x=as.Date(tmp$Date_start), y=as.numeric(as.vector(tmp$N_plant_start)), col=colSeg[i], pch=16, cex=.8)
  points(x=as.Date(tmp$Date_end), y=as.numeric(as.vector(tmp$N_plant_end)), col=colSeg[i], pch=16, cex=.8)
  
  if(Type[i] == "Cutting"){
    text(x=as.Date(tmp$Date_start), y=as.numeric(as.vector(tmp$N_plant_end))+20, labels = as.Date(tmp$Date_start), srt=45, adj=0.5, cex=.6)
  }
  if(Type[i] == "Rooting"){
    text(x=as.Date(tmp$Date_end), y=as.numeric(as.vector(tmp$N_plant_end))+8, labels = as.Date(tmp$Date_end), srt=45, adj=0, cex=.6)
  }
  if(Type[i] == "Growth"){
    #Start of experiment
    startExp <- which(tmp$Date_start == as.character(min(as.Date(tmp$Date_start))))
    text(x=as.Date(tmp$Date_start)[startExp], y=as.numeric(as.vector(tmp$N_plant_start))[startExp]+8, labels = as.Date(tmp$Date_start)[startExp], srt=90, adj=0, cex=.6)
    
    #End of experiment
    endExp <- which(tmp$Date_end == as.character(max(as.Date(tmp$Date_end))))
    text(x=as.Date(tmp$Date_end)[endExp], y=as.numeric(as.vector(tmp$N_plant_end))[endExp]+8, labels = as.Date(tmp$Date_end)[endExp], srt=90, adj=0, cex=.6)
  }
  if(Type[i] == "Biobanking"){
    #Add segment(s) showing when rooted shoot tips are split between growth and biobanking
    ystart <- as.numeric(as.vector(OUTdat$N_plant_start[which(OUTdat$Date_start == tmp$Date_start & OUTdat$Type == "Growth")]))
    yend <- as.numeric(as.vector(OUTdat$N_plant_end[which(OUTdat$Date_end == as.character(tmp$Date_start) & OUTdat$Type == "Rooting")]))
    segments(x0=as.Date(tmp$Date_start), x1=as.Date(tmp$Date_start), y0=ystart, y1=yend, lty=2) 
    text(x=as.Date(tmp$Date_start)-4, y=ystart, paste("Split rooted shoot tips:", paste(as.character(tmp$Type), " (", tmp$N_plant_start, ")", " & Growth (", ystart, ")", sep="")), adj=0, srt=90, cex=.6)
    
    #End of experiment
    endExp <- which(tmp$Date_end == as.character(max(as.Date(tmp$Date_end))))
    text(x=as.Date(tmp$Date_end)[endExp], y=as.numeric(as.vector(tmp$N_plant_end))[endExp]+8, labels = paste(as.Date(tmp$Date_end)[endExp], " (", as.character(tmp$`Total_Biomass (gr)`), " gr.)", sep=''), srt=90, adj=0, cex=.6)
  }
  
}

#Add legend
legend("topleft", legend = Type, lty=1, col=colSeg, cex=.6)
dev.off()

