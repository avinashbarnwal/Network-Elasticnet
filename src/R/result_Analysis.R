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
#attach(data)

data_input<-function(name)
{
  file_name=paste(name,".csv",sep="")
  data=read.table(file=file_name,sep=",",header=TRUE,stringsAsFactors = FALSE)
  print(colnames(data))
  data$OS_STATUS_NUM=ifelse(data$OS_STATUS=="DECEASED",1,2)
  return(data)
}

data = data_input("text_gene_coeff_III_0_elastic_net")
###########################################CLUSTERING#######################################################
kmeans_cluster_no=function(data){
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
  return(list(var=var_explained,cluster=cluster))
}

kmeans_cluster = kmeans_cluster_no(data)

plot(kmeans_cluster$cluster,kmeans_cluster$var)

coeffCluster <- kmeans(data[, 1:6],7, nstart = 20)
data$cluster = coeffCluster$cluster
clusplot(data[, 1:6], data$cluster,main="K-Means Clustering",labels=FALSE,color=TRUE)
centers= coeffCluster$centers

###########################################Gaussian Mixture Model###################################

X = data[,1:6]
fit = Mclust(X)
fitsummary = summary(fit, parameters = TRUE)
gmm_cluster = fitsummary$classification
data$gmm_cluster = gmm_cluster
plot(fit, what = "classification")
gmm_center = fitsummary$mean

##########################################Survival Plot#############################################

###########################################KMEANS################################################
sfit=survfit(Surv(OS_MONTHS,OS_STATUS_NUM)~kmeans_cluster, data=data)
ggsurvplot(sfit,conf.int=F,ggtheme = theme_bw(),font.main = 18,font.x =  10,pval.size=3,font.y = 10,font.tickslab = 10,legend.labs = c("1", "2", "3","4"),pval=T,risk.table=T,risk.table.fontsize = 3,risk.table.title ="",risk.table.height=0.5)
#ggpar(p,font.legend = c(10, "bold", "red"))
#smoke = table(data$cluster,data$OS_STATUS)
#plot(sfit)

############################################GMM##################################################
sfit=survfit(Surv(OS_MONTHS,OS_STATUS_NUM)~gmm_cluster, data=data)
ggsurvplot(sfit,conf.int=F,ggtheme = theme_bw(),font.main = 18,font.x =  10,pval.size=3,font.y = 10,font.tickslab = 10,legend.labs = c("1", "2","3","4"),pval=T,risk.table=T,risk.table.fontsize = 4,risk.table.title ="",risk.table.height=0.5)

#smoke = table(data$cluster,data$OS_STATUS)
#plot(sfit)

#data$OS_STATUS_NUM=ifelse(data$OS_STATUS=="DECEASED",1,0)
#BRCAOV.survInfo =  survivalTCGA(BRCA.clinical, OV.clinical,extract.cols = "admin.disease_code") 
#data$cluster=ifelse(data$cluster==1|data$cluster==2,1,ifelse(data$cluster==3,2,ifelse(data$cluster==4,3,4)))


##################################IMPUTING DATA########################################
tempData = mice(data,m=5,maxit=50,meth='pmm',seed=500)
data = complete(tempData)
data$OS_STATUS = ifelse(data$OS_STATUS_NUM==1,"DECEASED","LIVING")
data$OS_STATUS=as.factor(OS_STATUS)

#summary(tempData)

data$predicted = data$Age_Coeff*data$AGE + data$Gender_Coeff*data$Gender + data$TUMOR_STATUS_Coeff*data$Tumor_Stage+data$Smoking_Coeff*data$Smoking+data$Constant

data=data[order(data$OS_MONTHS),]
data$index=seq(1,length(data[,1]))

plot(data$index,data$OS_MONTHS,type="l",ylab="Normalized Survival months",xlab="Index")
lines(data$predicted,col="red")
legend("topleft",c("Original","Predicted"),lty=c(1,1),lwd=c(1,1),col=c("black","red"),cex=1,box.lty=0)

data$Tumor_Stage=factor(data$Tumor_Stage,labels=c("Stage I/IA","Stage IB","Stage II/IIA","Stage IIB","Stage III/IIIA","Stage IV"))
data$Gender=factor(data$Gender,labels=c("Male","Female"))


plot(AGE,OS_MONTHS,main="Normalized Survival Months vs Age",ylab="Normalized Survival Months")
boxplot(OS_MONTHS~Gender,data=data,main="Normalized Survival Months vs Gender")
boxplot(OS_MONTHS~Tumor_Stage,data=data,main="Normalized Survival Months vs Turmor Stage",las=2 ,cex.axis=0.6)
boxplot(OS_MONTHS~Smoking,data=data,main="Normalized Survival Months vs Smoking History")

plot(AGE,Age_Coeff,main="Age Coefficients vs Age",ylab="Age Coefficients")
boxplot(Gender_Coeff~Gender,data=data,main="Gender Coefficients vs Gender",ylab="Gender Coefficients")
boxplot(TUMOR_STATUS_Coeff~Tumor_Stage,data=data,main="Tumor Stage Coefficients vs Tumor Stage",las=2,cex.axis=0.6)
boxplot(Smoking_Coeff~Smoking,data=data,main="Smoking Coefficients vs Smoking history")

data_jnci=read.table(file="text_gene_coeff_jnci.csv",sep=",",header=TRUE)
attach(data_jnci)
data_jnci$predicted = data_jnci$Age_Coeff*data_jnci$Age+data_jnci$Gender_Coeff*data_jnci$Gender+data_jnci$Tumor_Coeff*data_jnci$Tumor+data_jnci$Constant
data_jnci=data_jnci[order(data_jnci$OS_Months),]
data_jnci$index=seq(1,length(data_jnci[,1]))

plot(data_jnci$index,data_jnci$OS_Months,type="l",ylab="Normalized Survival months",xlab="Index")
lines(data_jnci$predicted,col="red")
legend("topleft",c("Original","Predicted"),lty=c(1,1),lwd=c(1,1),col=c("black","red"),cex=1,box.lty=0)

data_jnci$Tumor=factor(data_jnci$Tumor,labels=c("Stage I/IA","Stage IB","Stage II/IIA","Stage IIB","Stage III/IIIA","Stage IIIB"))
data_jnci$Gender=factor(data_jnci$Gender,labels=c("Male","Female"))


plot(Age,OS_Months,main="Normalized Survival Months vs Age",ylab="Normalized Survival Months")
boxplot(OS_Months~Gender,data=data_jnci,main="Normalized Survival Months vs Gender")
boxplot(OS_Months~Tumor,data=data_jnci,main="Normalized Survival Months vs Turmor Stage",las=2 ,cex.axis=0.6)

plot(Age,Age_Coeff,data=data_jnci,main="Age Coefficients vs Age",ylab="Age Coefficients")
boxplot(Gender_Coeff~Gender,data=data_jnci,main="Gender Coefficients vs Gender",ylab="Gender Coefficients")
boxplot(Tumor_Coeff~Tumor,data=data_jnci,main="Tumor Stage Coefficients vs Tumor Stage",las=2,cex.axis=0.6)

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

#########################################################
### C) Customizing and plotting the heat map
#########################################################

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,0,length=100),               # for red
               seq(0.01,0.8,length=100),           # for yellow
               seq(0.81,1,length=100))             # for green

par(cex.main=.8)
heatmap.2(data_coeff,
         #cellnote = data_coeff,  # same data set for cell labels
          main = "Coefficients Analysis", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(6,12),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA",cexRow=0.3,cexCol=0.4)            # turn off column clustering

#dev.off()               # close the PNG device


