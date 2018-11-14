#2018-11-12 dinclaux@insa-toulouse.fr

#Script to generate Barplot and bias with an output file iscocor 2.X.X


# Reinitialize the session
rm(list=ls(all=TRUE))

#############################################################
###                     Reading data                      ###
#############################################################

setwd("~/Labo/Data/Script/Isocor/IsocorToPlot")
file = "Data_example.tsv" # tsv file
Samples = 2 # 1 or 2 kind of samples
Name1 = "MC"
Name2 = "proteo"
p = 0.499  # 13C enrichment

#########################################################
###       Installing and loading required packages    ###
#########################################################


if (!require("stringr")){
  install.packages("stringr", dependencies = TRUE)
  library(stringr)
}

if (!require("RColorBrewer")){
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

if (!require("gtools")){
  install.packages("gtools", dependencies = TRUE)
  library(gtools)
}


#Theoretical CID
calc_cid_th <- function(n,p){
  # calculate the theoretical cid of a metabolite
  # with n carbon atoms at 13C enrichment p
  x=seq(0,n)
  cid=choose(n,x)*p**x*(1-p)**(n-x)
  return(cid)
}

#function for error bars
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
  #turn off warning messages 
}


#############################################################
###                   Rearrangement of data               ###
#############################################################

data <- read.delim(file, header=TRUE, allowEscapes=FALSE, sep="\t", na.strings="")

data[is.na(data)] <- 0
data <- data[order(data$sample),]



if(Samples == 1){
  nb = nrow(data)
  data <- rbind(data,data)
  Name2 <- Name1
  Samplesplot <- Name1
} 
if(Samples == 2){
  txt <- grep(substring(paste(data[1,1]),0,5), x = as.vector(data[,1]),value = F)
  nb = length(txt) * (txt[2]-txt[1])
  Samplesplot <- c(Name1,Name2)
}


data <- subset(data, select = c(sample,metabolite, isotopologue, isotopologue_fraction))
data$isotopologue <- as.factor(data$isotopologue)
data[,3]<- data.frame(Massif=paste("M",data[,3], sep = ""))

#############################################################
###           Split data by type of samples               ###
#############################################################

X1 <- data[1:nb,]
X2 <- data[(nb+1):nrow(data),]
nb_X1=(sum(table(as.numeric(X1$isotopologue_fraction))))/ (sum(table((unique(X1[,1])))))
nb_X2=(sum(table(as.numeric(X2$isotopologue_fraction))))/ (sum(table((unique(X2[,1])))))

Metabo_X1 <- as.data.frame(X1[1:nb_X1,c(2,3)])
Metabo_X2 <- as.data.frame(X2[1:nb_X2,c(2,3)])
X1<- matrix(as.numeric(as.matrix(X1$isotopologue_fraction)),
            nrow= nb_X1)

X2<- matrix(as.numeric(as.matrix(X2$isotopologue_fraction)),
            nrow= nb_X2)
t=c()
tmp <- setdiff(Metabo_X1[,1],Metabo_X2[,1])



if(length(tmp)>0){
  for(i in tmp){
    t <- rbind(t,subset(Metabo_X1,Metabo_X1[,1] == i ))
  }
  Metabo_X2 <- rbind(Metabo_X2,t)
  X2 <- rbind(X2,matrix(data=0, nrow = nrow(t), ncol = ncol(X2)))
}

#Theoretical CID
Theoretical_CID_X1=c()
for (n in table(Metabo_X1[,1])[unique(Metabo_X1[,1])]) {
  Theoretical_CID_X1 <-c(Theoretical_CID_X1,calc_cid_th(n-1,p))
}
Theoretical_CID_X2=c()
for (n in table(Metabo_X2[,1])[unique(Metabo_X2[,1])]) {
  Theoretical_CID_X2 <-c(Theoretical_CID_X2,calc_cid_th(n-1,p))
}

#colnames for samples
ii=data[1:nb,1][!is.na(data[1:nb,1])]
ii <- as.character(unique(ii))
X1 <- cbind(Metabo_X1,Theoretical_CID_X1,X1)
colnames(X1)=c("Metabolites", 
               "isotopologues",
               "Theoretical_CID",
               ii)
ii=data[nb+1:nrow(data),1][!is.na(data[nb+1:nrow(data),1])]
ii <- as.character(unique(ii))
X2 <- cbind(Metabo_X2,Theoretical_CID_X2,X2)
colnames(X2)=c("Metabolites", 
               "isotopologues", 
               "Theoretical_CID",
               ii)

#assign names
nmu = unique(X1$Metabolites)
assign(paste(Name1),X1)
assign(paste(Name2),X2)
assign(paste("Theoretical_CID_",Name1,sep=""),Theoretical_CID_X1)
assign(paste("Theoretical_CID_",Name2,sep=""),Theoretical_CID_X2)
rm(data,X1,X2,Metabo_X1,Metabo_X2,Theoretical_CID_X1,Theoretical_CID_X2)

#Opening the pdf

file <- gsub(pattern = ".txt" , replacement = "" , str_sub(file,end=-5), fixed = TRUE)
tmp <- paste(file,".pdf", sep = "")
pdf(file = tmp ,onefile= TRUE, width = 15)
par(mfrow=c(1,2))

for (nm in nmu) {
  for (var in Samplesplot){
    tmp <-paste("Theoretical_CID_",var, sep="")
    Theo<- as.matrix(eval(parse(text=tmp)))
    eval(parse(text=paste("meta <-", var,"",sep = "")))  
    
    a=colnames(meta)
    b=str_sub(a[4:(length(a))],end=-4)
    c= c("Theoretical CID",b)
    d=unique(c)
    res=matrix(NA,nrow=nrow(meta), ncol=length(d), dimnames=list(row=NULL,col=d))
    for (i in d){
      e=which(c==i)
      if ((length(e))>1){
        mi=apply(meta[,(e+2)], 1, mean)
      } else {
        mi=meta[,(e+2)]
      }
      res[,i]=mi
    }
    
    
    #bias samples
    
    bias_res <- cbind(meta[,1:2], as.data.frame(apply(meta[,3:(ncol(meta))], 2, FUN=function(x) (x-Theo)*100)))
    
    tmp <- paste(var,file,"_bias_samples.txt", sep = "_")
    write.table(bias_res, file = tmp,
                sep = "\t",
                row.names = FALSE,
                quote = FALSE)
    #bias mean
    
    bias_res <- cbind(meta[,1:2], as.data.frame(apply(res[,1:(ncol(res))], 2, FUN=function(x) (x-Theo)*100)))
    
    tmp <- paste(var,file,"_bias.txt", sep = "_")
    write.table(bias_res, file = tmp,
                sep = "\t",
                row.names = FALSE,
                quote = FALSE)
    
    #############################################################
    ###                   Barplot                             ###
    #############################################################
    
    #colors
    cols <- as.character(brewer.pal(n=9, name = "Paired"))
    #sd barplot
    bias=matrix(NA,nrow=nrow(meta), ncol=length(d), dimnames=list(row=NULL,col=d))
    for (i in d){
      e=which(c==i)
      if ((length(e))>1){
        mi=apply(meta[,(e+2)], 1, sd)
      }else{
        mi=meta[,(e+2)]
      }
      bias[,i]=mi
    }
    
    ind <- meta[,1]==nm
    res <- t(as.matrix(res[ind,1:ncol(res)]))
    res <- res[c(rownames(res)[1],mixedsort(x = rownames(res)[-1])),]
    bias <- t(as.matrix(bias[ind,1:ncol(bias)]))
    bias[1,] <- NA
    bias[c(rownames(bias)[1],mixedsort(x = rownames(bias)[-1])),]
    #Barplot
    barx <- barplot(res,
                    beside = TRUE,
                    main = c(nm,var),
                    legend = rownames(res),
                    names.arg =meta$isotopologues[ind],
                    ylim = c(0,1.1* max(apply(res, 1, max))),
                    col = c("#FFFF00",cols[1:(nrow(res)-1)]),
                    args.legend = list("topright", cex = 0.7,bty = "n"))
    
    suppressWarnings(error.bar(barx, res, bias))
  }
  
}

dev.off()

