#In vitro propagation
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