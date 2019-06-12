setwd('/Users/avinashbarnwal/Desktop/Network_Lasso/TCGA/LUAD/Elastic Net III')
library(sqldf)
library(cluster)
library(fpc)
library(factoextra)
library(mclust)
library(survival)
library(ggplot2)
library(survminer)
library(missForest)
library(DMwR)
library(mice)
library(RTCGA.clinical)
library(ggpubr)
#attach(data)


data_input <- function(name){
  file_name=paste(name,".csv",sep="")
  data=read.table(file=file_name,sep=",",header=TRUE,stringsAsFactors = FALSE)
}

data = data_input("text_gene_coeff_IIB_INTERACTION")
data = data_input("text_gene_coeff_III_INTERACTION")
data = data_input("text_gene_coeff_IIA_INTERACTION")
data = data_input("IA_INTERACTION_IMAC_100000")
data = data_input("IB_INTERACTION_100000")
data = data_input("text_gene_coeff_III_0_elastic_net")

boxplot(OS_MONTHS~GENDER,data=data,main="Normalized Survival Months vs Gender")
boxplot(OS_MONTHS~SMOKING,data=data,main="Normalized Survival Months vs Turmor Stage",las=2 ,cex.axis=0.6)
boxplot(OS_MONTHS~Smoking,data=data,main="Normalized Survival Months vs Smoking History")

data$OS_STATUS_NUM=ifelse(data$OS_STATUS=="DECEASED",1,2)
data$predicted = data$GENDER_COEFF*data$GENDER+data$AGE_COEFF*data$NORMALISED_AGE+data$SMOKING_COEFF*data$SMOKING+data$AGE_SMOKING_COEFF*data$AGE_SMOKING+data$AGE_GENDER_COEFF*data$AGE_GENDER+data$GENDER_SMOKING_COEFF*data$GENDER_SMOKING_COEFF+data$CONSTANT


data=data[order(data$OS_MONTHS),]
data$index=seq(1,length(data[,1]))

plot(data$index,data$NORMALIZED_OS_MONTHS,type="l",ylab="Normalized Survival months",xlab="Index",main="STAGE III")
lines(data$predicted,col="red")
#legend("topleft",c("Original","Predicted"),lty=c(1,1),lwd=c(1,1),col=c("black","red"),cex=1,box.lty=0)

plot(data$NORMALIZED_OS_MONTHS,data$predicted,xlab="Actual Normalized Survival Months",ylab="Predicted Survival Months")

rownames(data)=NULL
set.seed(20)
var_explained<-vector()
cluster=seq(1:20)
flag=0
for(i in 1:20){
  flag=flag+1
  coeffCluster <- kmeans(data[, 1:6], cluster[flag], nstart = 20)
  var_explained[flag] <- 100*coeffCluster$betweenss/coeffCluster$totss
}

plot(cluster,var_explained,main="STAGE III/IV")
coeffCluster = kmeans(data[, 1:6],5, nstart = 20)
data$kmeans_cluster = coeffCluster$cluster

plotcluster(data[, 1:6], data$kmeans_cluster,main="K-Means Clustering")

table(data$kmeans_cluster)

##############################FOR III/IV##########################
#1  2  3  4  5 
#14  8 48  9 20 

##############################FOR IB##########################
#1   2   3   4 
#1   1 137   2 
##############################FOR IA##########################
##1  2  3  4  5  6  7  8  9 
##14 15  1  3  1 19  3 22  1 


##############################FOR IIA##########################
#1  2  3  4  5  6  7  8  9 
#8  6  1 17  1  2 12  5  9 

##############################FOR III/IV##########################
#1  2  3  4  5 
#5 47  4 19  2 

###########################FOR IIB#######################
#1  2  3  4  5  6  7  8  9 
#5  7  2 11 12  1  4  5 36 
###########################################Gaussian Mixture Model###################################
X = data[,1:6]
fit = Mclust(X)
fitsummary = summary(fit, parameters = TRUE)
gmm_cluster = fitsummary$classification
data$gmm_cluster = gmm_cluster
#plot(fit)
plot(fit, what = "classification",main="III/IV",cex.lab=0.8)
table(gmm_cluster)

#1  2  3  4  5  6 
#8 20 39 14 11  7 

#Clustering table for IB:
#1  2  3 
#93  9 39 


#Clustering table for IA:
#  1  2 
#13 66 

#Clustering table for IIA:
#1  2  3  4 
#22 30  5  4


#Clustering table for III/IV:
#1  2  3  4  5  6  7  8  9 
#9 15 30  6  2  2  5  6  2 


#Clustering table for IIB:
#1  2  3  4  5  6 
#9 15 48  3  3  5 
##########################################Survival Plot#############################################

###########################################KMEANS################################################
sfit=survfit(Surv(OS_MONTHS,OS_STATUS_NUM)~kmeans_cluster, data=data)
ggsurvplot(sfit,conf.int=F,ggtheme = theme_bw(),font.main = 18,font.x =  10,pval.size=3,font.y = 10,font.tickslab = 10,legend.labs = c("1", "2", "3","4","5"),pval=T,risk.table=T,risk.table.fontsize = 3,risk.table.title ="",risk.table.height=0.5)
#ggpar(p,font.legend = c(10, "bold", "red"))
#smoke = table(data$cluster,data$OS_STATUS)
#plot(sfit)

############################################GMM##################################################
sfit=survfit(Surv(OS_MONTHS,OS_STATUS_NUM)~gmm_cluster, data=data)
ggsurvplot(sfit,conf.int=F,ggtheme = theme_bw(),font.main = 18,font.x =  10,pval.size=3,font.y = 10,font.tickslab = 10,legend.labs = c("1", "2","3","4","5","6"),pval=T,risk.table=T,risk.table.fontsize = 4,risk.table.title ="",risk.table.height=0.5)
#smoke = table(data$cluster,data$OS_STATUS)
#plot(sfit)

#data$OS_STATUS_NUM=ifelse(data$OS_STATUS=="DECEASED",1,0)
#BRCAOV.survInfo =  survivalTCGA(BRCA.clinical, OV.clinical,extract.cols = "admin.disease_code") 
#data$cluster=ifelse(data$cluster==1|data$cluster==2,1,ifelse(data$cluster==3,2,ifelse(data$cluster==4,3,4)))

#######################################CONTRIBUTION ANALYSIS####################################

#data$Age_Contribution = 100*data$NORMALISED_OS_MONTHS/(data$Age_Coeff*data$AGE)

data$Age_Contribution = data$Age_Coeff*data$AGE
data$Gender_Contribution = data$Gender_Coeff*data$Gender
data$Tumor_Contribution = data$TUMOR_STATUS_Coeff*data$Tumor_Stage
data$Smoking_Contribution = data$Smoking_Coeff*data$Smoking

sd_Age      = sd(data$Age_Contribution)
sd_Gender   = sd(data$Gender_Contribution)
sd_Tumor    = sd(data$Tumor_Contribution)
sd_Smoking  = sd(data$Smoking_Contribution)

sd_os_months         = sd(data$NORMALISED_OS_MONTHS)

data$Age_Error       = data$NORMALISED_OS_MONTHS - data$Age_Contribution
sd_Age_error         = sd(data$Age_Error)

data$Gender_Error    = data$NORMALISED_OS_MONTHS - data$Gender_Contribution
sd_Gender_error      = sd(data$Gender_Error)

data$Tumor_Error     = data$NORMALISED_OS_MONTHS - data$Tumor_Contribution
sd_Tumor_error       = sd(data$Tumor_Error)

data$Smoking_Error   = data$NORMALISED_OS_MONTHS - data$Smoking_Contribution
sd_Smoking_error     = sd(data$Smoking_Error)


########################################Age First Selected Variable########################################

data$Error_1              = data$NORMALISED_OS_MONTHS - data$Age_Contribution
sd_Error_1                = sd(data$Error_1)

data$Gender_Error_2       = data$Error_1 - data$Gender_Contribution
sd_Gender_Error_2         = sd(data$Gender_Error_2)


data$Tumor_Error_2        = data$Error_1 - data$Tumor_Contribution
sd_Tumor_Error_2          = sd(data$Tumor_Error_2)

data$Smoking_Error_2      = data$Error_1 - data$Smoking_Contribution
sd_Smoking_Error_2        = sd(data$Smoking_Error_2)

########################################Smoking Second Selected Variable########################################

data$Error_2              = data$NORMALISED_OS_MONTHS - data$Age_Contribution - data$Gender_Contribution
sd_Error_2                = sd(data$Error_2)

data$Gender_Error_3       = data$Error_2 - data$Gender_Contribution
sd_Gender_Error_3         = sd(data$Gender_Error_3)

data$Tumor_Error_3        = data$Error_2 - data$Tumor_Contribution
sd_Tumor_Error_3          = sd(data$Tumor_Error_3)

#########################################Tumor Third Selected Variable##########################################

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}


#########################################################
### B) Reading in data and transform it into matrix format
#########################################################

data_coeff=data[,1:6]
data_coeff=as.matrix(data_coeff)
colnames(data_coeff)=c("GENDER","AGE","SMOKING","AGE_SMOKING","AGE_GENDER","GENDER_SMOKING")

#########################################################
### C) Customizing and plotting the heat map
#########################################################

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,0,length=100),               # for red
               seq(0.01,0.8,length=100),           # for yellow
               seq(0.81,1,length=100))             # for green

# creates a 5 x 5 inch image
#png("../heatmaps_coeff.png",    # create PNG for the heat map        
#    width = 5*300,        # 5 x 300 pixels
#    height = 5*300,
#    res = 300,            # 300 pixels per inch
#    pointsize = 8)        # smaller font size
par(cex.main=.4)
heatmap.2(data_coeff,
          #cellnote = data_coeff,  # same data set for cell labels
          col = rev(rainbow(20*10, start = 0/6, end = 4/6)), 
          scale="none",
          margins=c(5,0), # ("margin.Y", "margin.X")
          trace='none', 
          symkey=FALSE, 
          symbreaks=FALSE, 
          dendrogram='none',
          density.info='histogram', 
          denscol="black",
          keysize=1,
          #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
          key.par=list(mar=c(3.5,0,3,0)),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2.5, 5), lwid=c(1, 10, 1),cexRow=0.4,cexCol=0.6)         # turn off column clustering

dev.off()               # close the PNG device


