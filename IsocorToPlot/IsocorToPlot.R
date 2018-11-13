#2018-12-20 dinclaux@insa-toulouse.fr

#Script to generate Barplot and bias with an output file iscocor

# Reinitialize the session
rm(list=ls(all=TRUE))

#########################################################
###       Installing and loading required packages    ###
#########################################################


if (!require("stringr")) {
  install.packages("stringr", dependencies = TRUE)
  library(stringr)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

if (!require("data.table")) {
  install.packages("data.table", dependencies = TRUE)
  library(data.table)
}

#############################################################
###                     Reading data                      ###
#############################################################

setwd("~/Labo/Data/Script/Isocor/IsocorToPlot")
file = "Exemple.txt" # txt file
Samples = 2 # 1 or 2 kind of samples
Name1 = "MC"
Name2 = "proteo"
p = 0.499  # 13C enrichment


#############################################################
###                   Rearrangement of data               ###
#############################################################

data <- read.table(file,
                   header = TRUE,
                   sep = "\t",
                   dec = ".",
                   fill = TRUE,
                   na.strings = c("",NA))


nC_metab <- t.default(cbind(matrix(c("Alanine",
                                     "Arginine",
                                     "Aspartate",
                                     "Glutamate",
                                     "Glycine",
                                     "Histidine",
                                     "Isoleucine",
                                     "Leucine",
                                     "Lysine",
                                     "Methionine",
                                     "Phenylalanine",
                                     "Proline",
                                     "Serine",
                                     "Threonine",
                                     "Tyrosine",
                                     "Valine",
                                     "Asparagine",
                                     "Glutamine",
                                     "Tryptophan")), matrix(c(4,7,5,6,3,7,7,7,7,6,10,6,4,5,10,6,5,6,12))))
colnames(nC_metab) = nC_metab[1,]
nC_metab = nC_metab[-1,]

for (row in seq(nrow(data))) {
  if (is.na(data[row,2])) {
    data[row,2] <- j
  }else{
    j=data[row,2]
  }
}

detect_nodata <- function(n){
  if (n[3] == "No data."){
    p = nC_metab[n[2]]
    insert<- matrix(data = c(NA,paste(n[2]),"0",NA,NA,"0",NA,NA),nrow =as.numeric(p),ncol = length(data),byrow = TRUE)
    colnames(insert) <- colnames(data)
    data2 <- insert
  }
  else {
    data2 <- n
  }
  return(data2)
}

data <- do.call(rbind,apply(data, 1, detect_nodata))
data <- as.data.frame(data)

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

data[,3]<- data.frame(Massif=paste("M",data[,3], sep = ""))
X1 <- data[1:nb,]
X2 <- data[(nb+1):nrow(data),]
nb_X1=(sum(table(as.numeric(X1$Isotopologue.distribution))))/ (sum(table(as.numeric(unique(X1[,1])))))
nb_X2=(sum(table(as.numeric(X2$Isotopologue.distribution))))/ (sum(table(as.numeric(unique(X2[,1])))))

Metabo_X1 <- as.data.frame(X1[1:nb_X1,c(2,3)])
Metabo_X2 <- as.data.frame(X2[1:nb_X2,c(2,3)])
X1<- matrix(as.numeric(as.matrix(X1$Isotopologue.distribution)),
               nrow= nb_X1)

X2<- matrix(as.numeric(as.matrix(X2$Isotopologue.distribution)),
                nrow= nb_X2)


X1<-cbind(Metabo_X1, 
             X1)
t=c()
tmp <- setdiff(Metabo_X1[,1],Metabo_X2[,1])
X2<-cbind(Metabo_X2, 
              X2)
for(i in tmp){
  t <- rbind(t,subset(Metabo_X1,Metabo_X1[,1] == i ))
}


t <-cbind(t,matrix(data=0,nrow=nrow(t), ncol=(ncol(X2)-2)))
X2 <- rbind(X2,t)


#Theoretical CID
calc_cid_th <- function(n,p){
  # calculate the theoretical cid of a metabolite
  # with n carbon atoms at 13C enrichment p
  x=seq(0,n)
  cid=choose(n,x)*p**x*(1-p)**(n-x)
  return(cid)
}

Metabo_X1=c()
for (n in table(X1[,1])[unique(X1[,1])]) {
  Metabo_X1 <-c(Metabo_X1,calc_cid_th(n-1,p))
}
Metabo_X2=c()
for (n in table(X2[,1])[unique(X2[,1])]) {
  Metabo_X2 <-c(Metabo_X2,calc_cid_th(n-1,p))
}

#colnames for samples
ii=data[1:nb,1][!is.na(data[1:nb,1])]
ii <- as.character(ii)
colnames(X1)=c("Metabolites", 
                  "isotopologues", 
                  ii)
ii=data[nb:nrow(data),1][!is.na(data[nb:nrow(data),1])]
ii <- as.character(ii)
colnames(X2)=c("Metabolites", 
                   "isotopologues", 
                   ii)

X1 <- cbind(X1, Metabo_X1)
setnames(X1, old = c("Metabo_X1"), new = c("Theoretical CID"))
X2 <- cbind(X2, Metabo_X2)
setnames(X2, old = c("Metabo_X2"), new = c("Theoretical CID"))


#############################################################
###                   Barplot                             ###
#############################################################

#function for error bars
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
  #turn off warning messages 
  
  
}
#colors
cols <- as.character(brewer.pal(n=9, name = "Paired"))

#assign names
nmu = unique(X1$Metabolites)
assign(paste(Name1),X1)
assign(paste(Name2),X2)
assign(paste("Metabo_",Name1,sep=""),Metabo_X1)
assign(paste("Metabo_",Name2,sep=""),Metabo_X2)
rm(X1,X2,Metabo_X1,Metabo_X2)



#barplot
file <- gsub(pattern = ".txt" , replacement = "" , file, fixed = TRUE)
tmp <- paste(file,".pdf", sep = "")
pdf(file = tmp ,onefile= TRUE)
par(mfrow=c(1,2))

for (nm in nmu) {
  for (var in Samplesplot){
    tmp <-paste("Metabo_",var, sep="")
    Theo<- as.matrix(eval(parse(text=tmp)))
    eval(parse(text=paste("meta <-", var,"",sep = "")))  
    
    a=colnames(meta)
    b=str_sub(a[3:(length(a)-1)],end=-3)
    c= c("Metabolites","isotopologues",b,"Theoretical CID")
    d=unique(c)
    r=c(unique(str_sub(a[3:(length(a)-1)],end=-4)),"Theoretical CID")
    res=matrix(NA,nrow=nrow(meta), ncol=length(d), dimnames=list(row=NULL,col=d))
    for (i in d){
      e=which(c==i)
      if ((length(e))>1){
        mi=apply(meta[,e], 1, mean)
      }
      else{
        mi=meta[,e]
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
    
    bias_res <- cbind(meta[,1:2], as.data.frame(apply(res[,3:(ncol(res))], 2, FUN=function(x) (x-Theo)*100)))
    
    tmp <- paste(var,file,"_bias.txt", sep = "_")
    write.table(bias_res, file = tmp,
                sep = "\t",
                row.names = FALSE,
                quote = FALSE)
    
    #sd barplot
    bias=matrix(NA,nrow=nrow(meta), ncol=length(d), dimnames=list(row=NULL,col=d))
    for (i in d){
      e=which(c==i)
      if ((length(e))>1){
        mi=apply(meta[,e], 1, sd)
      }else{
        mi=meta[,e]
      }
      bias[,i]=mi
    }
    
    
    ind <- meta[,1]==nm
    res <- t(as.matrix(res[ind,3:ncol(res)]))
    bias <- t(as.matrix(bias[ind,3:ncol(bias)]))
    bias[nrow(bias),(bias[nrow(bias),]<5) ] <- NA
    barx <- barplot(res,
                    beside = TRUE,
                    main = c(nm,var),
                    #legend = r,
                    names.arg =meta$isotopologues[ind],
                    col = c(cols[1:(nrow(res)-1)],"#FFFF00"),
                    args.legend = list("topright", cex = 0.40))
    
    suppressWarnings(error.bar(barx, res, bias))
  }
  
}

dev.off()
