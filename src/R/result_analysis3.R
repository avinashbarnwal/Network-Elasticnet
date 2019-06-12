setwd('/Users/avinashbarnwal/Desktop/Research/PhD/network_elasticnet/NYSDS-2019')
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
library(mclust)
#attach(data)

path     = getwd()
dataPath = paste(path,"/result/final_result.csv",sep="")
data     = read.table(dataPath,header=TRUE,sep=",")

data[data['TUMOR_STAGE']=="",'TUMOR_STAGE'] = 'Stage IIIA'
data_coeff = data[,2:101]


kmeans_cluster_no=function(data){
  rownames(data)=NULL
  set.seed(20)
  var_explained<-vector()
  cluster=seq(1:20)
  flag=0
  for(i in 1:20){
    flag=flag+1
    coeffCluster <- kmeans(data, cluster[flag], nstart = 20)
    var_explained[flag] <- 100*coeffCluster$betweenss/coeffCluster$totss
    
  }
  return(list(var=var_explained,cluster=cluster))
}


FDR_FISHER_TEST = function(data,row,col){
  matrix_test = matrix(0,nrow=2,ncol=2)
  matrix_test[1,1] = data[row,col]
  matrix_test[2,1] = sum(data[,col]) - data[row,col]
  matrix_test[1,2] = sum(data[row,]) - data[row,col]
  matrix_test[2,2] = sum(data)       - data[row,col]
  result = fisher.test(matrix_test,alternative='two.sided')
  return(result$p.value)
}



kmeans_cluster = kmeans_cluster_no(data_coeff)
plot(kmeans_cluster$cluster,kmeans_cluster$var)
coeffCluster   = kmeans(data_coeff,7, nstart = 100)
data$cluster   = coeffCluster$cluster
data$cluster = ifelse(data$cluster %in% c(1,7),6,data$cluster)

fit = Mclust(data_coeff)
fitsummary = summary(fit, parameters = TRUE)
gmm_cluster = fitsummary$classification
data$gmm_cluster = gmm_cluster
plot(fit, what = "classification")
gmm_center = fitsummary$mean

#clusplot(data[, 1:6], data$cluster,main="K-Means Clustering",labels=FALSE,color=TRUE)
#centers= coeffCluster$centers

data[,'TUMOR_STAGE_GROUPED'] <- ifelse(data[,'TUMOR_STAGE'] %in% c('Stage I','Stage IA','Stage IB'),1,ifelse(data[,'TUMOR_STAGE'] %in% c('Stage IIA','Stage IIB'),2,ifelse(data[,'TUMOR_STAGE'] %in% c('Stage IIIA','Stage IIIB'),3,4)))
data_tumor        = table(data[,'cluster'],data[,'TUMOR_STAGE_GROUPED'])
data_test_smoking = matrix(0,nrow=5,ncol=4)

for(i in 1:5){
  for(j in 1:4){
    data_test_smoking[i,j] =  FDR_FISHER_TEST(data_tumor,i,j)
  }
}


fannyy             = fanny(data_coeff, k=5, metric = "euclidean", memb.exp = 1.2)
data$fanny_cluster = fannyy$clustering


pamy <- pam(data_coeff, 4)
data$pam_cluster   = pamy$clustering
data_tumor         = table(data[,'pam_cluster'],data[,'TUMOR_STAGE_GROUPED'])






