#2016-09-17 dinclaux@insa-toulouse.fr

#Script to generate a PCA with an output file iscocor
#Run this sript with the template

#########################################################
###       Installing and loading required packages    ###
#########################################################


if (!require("mixOmics")) {
  install.packages("mixOmicst", dependencies = TRUE)
  library(mixOmics)
}

#############################################################
###                     Reading data                      ###
#############################################################

setwd("~/Labo/Données/Script/Isocor/IsocorToPCA")
file = "Manip2_isocor_res.txt" # txt file
file_res = "transfo finale.txt" # after rearrangement of Isocor data
experiment = "manip 1" # experiment name

#############################################################
###                   Rearrangement of data               ###
#############################################################

data <- read.table(file,
           header = T,
           sep = "\t",
           dec = ".",
           fill = TRUE,
           na.strings = c("","NA"))

data <- as.data.frame(data)


Metabo <- read.table("template.txt",
                     sep = "\t")
Metabo <- as.data.frame(Metabo)


data2<- matrix(as.numeric(data$Isotopologue.distribution),
               nrow=100)
Metabo<-cbind(Metabo, 
              data2)

ii=data$Sample[!is.na(data$Sample)]
ii <- as.character(ii)
colnames(Metabo)=c("Metabolite", 
                   "isotopologues", 
                   ii)

#############################################################
###                   Save to file                        ###
#############################################################


write.table(Metabo, file = file_res,
            sep = "\t",
            row.names = F )


#############################################################
###                         PCA                           ###
#############################################################

data3 <-Metabo[,-c(1,2)]

pca.res<-pca(data3,
             ncomp = NULL,
             center = T,
             scale = T)
plotVar(pca.res,
        style = "3d",
        title =  experiment,
        overlap=T,
        rad.in=0,
        plot = T)

