library("sqldf")
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
library(xtable)
library(igraph)
library(gplots)
library(RColorBrewer)

wd =paste('/Users/avinashbarnwal/Desktop/Machine Learning/Network_Lasso/TCGA/100 GENE VAR')
setwd(wd)
set.seed(20)

LUSC_MEAN_OS_MONTHS    = 32.93
LUSC_STD_DEV_OS_MONTHS = 31.57

LUAD_MEAN_OS_MONTHS    = 29.9
LUAD_STD_DEV_OS_MONTHS = 30.19

data_LUAD = read.table("data_phy_top100_LUAD.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)
data_LUSC = read.table("data_phy_top100_LUSC.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)

data_LUSC$OS_STATUS_NUM = ifelse(data_LUSC$OS_STATUS=="DECEASED",1,2)
data_LUAD$OS_STATUS_NUM = ifelse(data_LUAD$OS_STATUS=="DECEASED",1,2)

data_LUAD[,'OS_MONTHS'] = LUAD_STD_DEV_OS_MONTHS*data_LUAD[,'OS_MONTHS']+ LUAD_MEAN_OS_MONTHS
data_LUSC[,'OS_MONTHS'] = LUSC_STD_DEV_OS_MONTHS*data_LUSC[,'OS_MONTHS'] + LUSC_MEAN_OS_MONTHS

Split <- function(string){
  x<-unlist(strsplit(string, "_"))
  return (x[1])
}

data_massage_coeff = function(data_name){
  file_name = paste(data_name,".csv",sep="")
  data      = read.table(file_name,sep=",",header=TRUE)
  data      = data[,0:102]
  colnames  = colnames(data)
  colnames  = colnames[2:length(colnames)]
  colnames_parent = sapply(colnames,Split)
  colnames_parent = c('INDEX',colnames_parent)
  colnames(data) = colnames_parent
  return (data)
}


kmeans_find<-function(data,name){
  set.seed(20)
  var_explained<-vector()
  cluster=seq(1:20)
  flag=0
  x<-data
  for(i in 1:20){
    flag=flag+1
    coeffCluster <- kmeans(x, cluster[flag], nstart = 20)
    var_explained[flag] <- 100*coeffCluster$betweenss/coeffCluster$totss
  }
  
  print(var_explained)
  title_name=paste("Var Explained",name)
  png(filename=title_name)
  plot(cluster,var_explained,main=title_name)
  dev.off()
}

kmeans_plot<-function(data,name,cluster){
  coeffCluster = kmeans(data,cluster, nstart = 20)
  data$kmeans_cluster = coeffCluster$cluster
  cx <- discrproj(data, data$kmeans_cluster, 'dc', 2)$proj
  stage_name=paste("Stage",name,sep=" ")
  title_name=paste("Kmeans Cluster",name)
  png(filename=title_name)
  plot(cx[,1],cx[,2],col=data$kmeans_cluster,xlab='eigen 1',ylab='eigen 2',main=stage_name)
  legend_name=as.character(seq(1,cluster))
  lty_name=rep(1,cluster)
  col_name=seq(1,cluster)
  legend('bottomleft',legend_name,lty=lty_name,col=col_name,cex=0.7)
  dev.off()
  return(data)
}

combination <- function(data,position,n_row){
  row     = data[position,]
  corr    = vector()
  for(i in 1:n_row){
    corr[i] = (1+cor(as.numeric(data[position,]),as.numeric(data[i,])))/2
  }
  corr[position]=0
  return(corr)
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


#################################################LUAD DATASET#########################################

LUAD_TOP100_COEFF             = data_massage_coeff("LUAD_TOP100")
data_LUAD_SELECTED            = data_LUAD[c('index','TUMOR_STAGE','SMOKING','OS_STATUS_NUM','OS_MONTHS','SMOKING_PACK_YEARS')]
colnames(data_LUAD_SELECTED)  = c('INDEX','TUMOR_STAGE','SMOKING','OS_STATUS_NUM','OS_MONTHS','SMOKING_PACK_YEARS')

LUAD_TOP100_COEFF             = merge(LUAD_TOP100_COEFF,data_LUAD_SELECTED,by ='INDEX')
LUAD_TOP100_COEFF_GENE        = LUAD_TOP100_COEFF[2:101]

kmeans_find(LUAD_TOP100_COEFF_GENE,"LUAD")

LUAD_TOP100_COEFF_CLUSTER = kmeans_plot(LUAD_TOP100_COEFF_GENE,"LUAD",6)
LUAD_TOP100_COEFF_CLUSTER$TUMOR_STAGE = LUAD_TOP100_COEFF[,'TUMOR_STAGE']
LUAD_TOP100_COEFF_CLUSTER$INDEX = LUAD_TOP100_COEFF[,'INDEX']
LUAD_TOP100_COEFF_CLUSTER$SMOKING_PACK_YEARS = LUAD_TOP100_COEFF[,'SMOKING_PACK_YEARS']
LUAD_TOP100_COEFF_CLUSTER$SMOKING = LUAD_TOP100_COEFF[,'SMOKING']
LUAD_TOP100_COEFF_CLUSTER$OS_STATUS_NUM = LUAD_TOP100_COEFF[,'OS_STATUS_NUM']
LUAD_TOP100_COEFF_CLUSTER$OS_MONTHS = LUAD_TOP100_COEFF[,'OS_MONTHS']

LUAD_TOP100_COEFF_CLUSTER[,'TUMOR_STAGE_GROUPED'] <- ifelse(LUAD_TOP100_COEFF_CLUSTER[,'TUMOR_STAGE'] %in% c('Stage I','Stage IA','Stage IB'),1,ifelse(LUAD_TOP100_COEFF_CLUSTER[,'TUMOR_STAGE'] %in% c('Stage IIA','Stage IIB'),2,ifelse(LUAD_TOP100_COEFF_CLUSTER[,'TUMOR_STAGE'] %in% c('Stage IIIA','Stage IIIB'),3,4)))
LUAD_TOP100_COEFF_CLUSTER[,'SMOKING_PACK_YEARS_GROUPED'] <- ifelse(LUAD_TOP100_COEFF_CLUSTER[,'SMOKING_PACK_YEARS'] <= 20,1,ifelse(LUAD_TOP100_COEFF_CLUSTER[,'SMOKING_PACK_YEARS'] <= 40,2,ifelse(LUAD_TOP100_COEFF_CLUSTER[,'SMOKING_PACK_YEARS'] <= 50,3,4)))

LUAD_table_tumor = table(LUAD_TOP100_COEFF_CLUSTER[,'kmeans_cluster'],LUAD_TOP100_COEFF_CLUSTER[,'TUMOR_STAGE_GROUPED'])
LUAD_table_smoking_created = table(LUAD_TOP100_COEFF_CLUSTER[,'kmeans_cluster'],LUAD_TOP100_COEFF_CLUSTER[,'SMOKING_PACK_YEARS_GROUPED'])
LUAD_table_smoking = table(LUAD_TOP100_COEFF_CLUSTER[,'kmeans_cluster'],LUAD_TOP100_COEFF_CLUSTER[,'SMOKING'])


LUAD_TEST_SMOKING = matrix(0,nrow=6,ncol=4)

for(i in 1:6){
  for(j in 1:4){
    LUAD_TEST_SMOKING[i,j] =  FDR_FISHER_TEST(LUAD_table_smoking,i,j)
  }
}


LUAD_TEST_TUMOR = matrix(0,nrow=6,ncol=4)
for(i in 1:6){
  for(j in 1:4){
    LUAD_TEST_TUMOR[i,j] =  FDR_FISHER_TEST(LUAD_table_tumor,i,j)
  }
}

CORRELATION_DATA            = LUAD_TOP100_COEFF_GENE
rownames(CORRELATION_DATA)  = unlist(LUAD_TOP100_COEFF_CLUSTER['INDEX'])

kmean_seq = c("1", "2","3","4","5","6")
#png(filename='Kmeans Survival LUAD') 
sfit=survfit(Surv(OS_MONTHS,OS_STATUS_NUM)~kmeans_cluster,data=LUAD_TOP100_COEFF_CLUSTER)
ggsurvplot(sfit,conf.int=F,ggtheme = theme_bw(),font.main = 18,font.x =  10,pval.size=3,font.y = 10,font.tickslab = 10,legend.labs = kmean_seq,pval=T,risk.table=T,risk.table.fontsize = 3,risk.table.title ="",risk.table.height=0.5)
#dev.off()


n_row = length(CORRELATION_DATA[,1])
n_col = n_row
choosen_mat = matrix(0,nrow=n_row,ncol=n_col)
order_mat   = matrix(0,nrow=n_row,ncol=n_col)

for(i in 1:n_row){
  print(i)
  choosen_all            = combination(CORRELATION_DATA,i,n_row)
  order_all              = order(-choosen_all)
  choosen_mat[i,]        = choosen_all
  order_mat[i,]          = order_all
  #temp                   = as.vector(order_mat[i,1:5])  
  #choosen_mat[i,-(temp)] = 0
}

choosen_mat[choosen_mat<0.6] =0

c("gray50", "tomato", "gold")
cluster_factor = as.character(factor(LUAD_TOP100_COEFF_CLUSTER[,'kmeans_cluster'],labels= c("circle", "square", "csquare", "rectangle","crectangle","vrectangle")))

cluster_factor = as.character(factor(LUAD_TOP100_COEFF_CLUSTER[,'kmeans_cluster'],labels= c("gray50", "blue", "red", "tomato","lightsteelblue2","orange")))

smoking_factor = as.character(factor(LUAD_TOP100_COEFF_CLUSTER[,'SMOKING'],labels= c("green", "blue", "red", "black")))
tumor_stage_factor = as.character(factor(LUAD_TOP100_COEFF_CLUSTER[,'TUMOR_STAGE_GROUPED'],labels= c("green", "blue", "red", "black")))


g2 = graph_from_adjacency_matrix(choosen_mat,weight=TRUE,mode="undirected")
fine = 50 # this will adjust the resolving power.
os_months = as.matrix(LUAD_TOP100_COEFF_CLUSTER['OS_MONTHS'])
V(g2)$color = brewer.pal(as.numeric(cut(os_months,breaks=fine)),"Blues")
#V(g2)$shape = cluster_factor
plot(g2,vertex.size=10,vertex.label.cex=0.2)


choosen_mat_sim = choosen_mat
choosen_mat_sim[choosen_mat_sim<0.6] =0


g2 = graph_from_adjacency_matrix(choosen_mat_sim,weight=TRUE,mode="undirected")
V(g2)$color = cluster_factor
V(g2)$label = LUAD_TOP100_COEFF_CLUSTER[,'kmeans_cluster']
#V(g2)$shape = cluster_factor
plot(g2,layout=layout.kamada.kawai,vertex.size=30,vertex.label.dist=1.5,vertex.label.cex=0.5,edge.color="gray70",edge.width=2)



g3 = graph_from_adjacency_matrix(choosen_mat,weight=TRUE,mode="undirected")
V(g3)$color = smoking_factor
V(g3)$shape = cluster_factor
plot(g3,vertex.size=10,vertex.label.cex=0.2)

g4 = graph_from_adjacency_matrix(choosen_mat,weight=TRUE,mode="undirected")
V(g4)$color = tumor_stage_factor
V(g4)$shape = cluster_factor
plot(g4,vertex.size=10,vertex.label.cex=0.2)


#################################################LUSC DATASET#########################################

LUSC_TOP100_COEFF             = data_massage_coeff("LUSC_TOP100")
data_LUSC_SELECTED            = data_LUSC[c('index','TUMOR_STAGE','SMOKING','OS_STATUS_NUM','OS_MONTHS','SMOKING_PACK_YEARS')]

colnames(data_LUSC_SELECTED)  = c('INDEX','TUMOR_STAGE','SMOKING','OS_STATUS_NUM','OS_MONTHS','SMOKING_PACK_YEARS')
LUSC_TOP100_COEFF             = merge(LUSC_TOP100_COEFF,data_LUSC_SELECTED,by ='INDEX')
LUSC_TOP100_COEFF_GENE        = LUSC_TOP100_COEFF[2:101]

kmeans_find(LUSC_TOP100_COEFF_GENE,"LUSC")
LUSC_TOP100_COEFF_CLUSTER = kmeans_plot(LUSC_TOP100_COEFF_GENE,"LUSC",6)
LUSC_TOP100_COEFF_CLUSTER$TUMOR_STAGE = LUSC_TOP100_COEFF[,'TUMOR_STAGE']
LUSC_TOP100_COEFF_CLUSTER$INDEX = LUSC_TOP100_COEFF[,'INDEX']
LUSC_TOP100_COEFF_CLUSTER$SMOKING_PACK_YEARS = LUSC_TOP100_COEFF[,'SMOKING_PACK_YEARS']
LUSC_TOP100_COEFF_CLUSTER$SMOKING = LUSC_TOP100_COEFF[,'SMOKING']
LUSC_TOP100_COEFF_CLUSTER$OS_STATUS_NUM = LUSC_TOP100_COEFF[,'OS_STATUS_NUM']
LUSC_TOP100_COEFF_CLUSTER$OS_MONTHS = LUSC_TOP100_COEFF[,'OS_MONTHS']


LUSC_TOP100_COEFF_CLUSTER[,'TUMOR_STAGE_GROUPED'] <- ifelse(LUSC_TOP100_COEFF_CLUSTER[,'TUMOR_STAGE'] %in% c('Stage I','Stage IA','Stage IB'),1,ifelse(LUSC_TOP100_COEFF_CLUSTER[,'TUMOR_STAGE'] %in% c('Stage IIA','Stage IIB'),2,ifelse(LUSC_TOP100_COEFF_CLUSTER[,'TUMOR_STAGE'] %in% c('Stage IIIA','Stage IIIB'),3,4)))
LUSC_TOP100_COEFF_CLUSTER[,'SMOKING_PACK_YEARS_GROUPED'] <- ifelse(LUSC_TOP100_COEFF_CLUSTER[,'SMOKING_PACK_YEARS'] <= 30.75,1,ifelse(LUSC_TOP100_COEFF_CLUSTER[,'SMOKING_PACK_YEARS'] <= 50,2,ifelse(LUSC_TOP100_COEFF_CLUSTER[,'SMOKING_PACK_YEARS'] <= 63.75,3,4)))


LUSC_table_tumor = table(LUSC_TOP100_COEFF_CLUSTER[,'kmeans_cluster'],LUSC_TOP100_COEFF_CLUSTER[,'TUMOR_STAGE_GROUPED'])
LUSC_table_smoking_created = table(LUSC_TOP100_COEFF_CLUSTER[,'kmeans_cluster'],LUSC_TOP100_COEFF_CLUSTER[,'SMOKING_PACK_YEARS_GROUPED'])
LUSC_table_smoking = table(LUSC_TOP100_COEFF_CLUSTER[,'kmeans_cluster'],LUSC_TOP100_COEFF_CLUSTER[,'SMOKING'])

LUSC_TEST_SMOKING = matrix(0,nrow=6,ncol=4)

for(i in 1:6){
  for(j in 1:4){
    LUSC_TEST_SMOKING[i,j] =  FDR_FISHER_TEST(LUSC_table_smoking,i,j)
  }
}


LUSC_TEST_TUMOR = matrix(0,nrow=6,ncol=4)
for(i in 1:6){
  for(j in 1:4){
    LUSC_TEST_TUMOR[i,j] =  FDR_FISHER_TEST(LUSC_table_tumor,i,j)
  }
}




kmean_seq = c("1", "2","3","4","5","6")
#png(filename='Kmeans Survival LUAD') 
sfit=survfit(Surv(OS_MONTHS,OS_STATUS_NUM)~kmeans_cluster,data=LUSC_TOP100_COEFF_CLUSTER)
ggsurvplot(sfit,conf.int=F,ggtheme = theme_bw(),font.main = 18,font.x =  10,pval.size=3,font.y = 10,font.tickslab = 10,legend.labs = kmean_seq,pval=T,risk.table=T,risk.table.fontsize = 3,risk.table.title ="",risk.table.height=0.5)
#dev.off()



CORRELATION_DATA            = LUSC_TOP100_COEFF_GENE
rownames(CORRELATION_DATA)  = unlist(LUSC_TOP100_COEFF_CLUSTER['INDEX'])

kmean_seq = c("1", "2","3","4","5","6")
#png(filename='Kmeans Survival LUAD') 
sfit=survfit(Surv(OS_MONTHS,OS_STATUS_NUM)~kmeans_cluster,data=LUSC_TOP100_COEFF_CLUSTER)
ggsurvplot(sfit,conf.int=F,ggtheme = theme_bw(),font.main = 18,font.x =  10,pval.size=3,font.y = 10,font.tickslab = 10,legend.labs = kmean_seq,pval=T,risk.table=T,risk.table.fontsize = 3,risk.table.title ="",risk.table.height=0.5)
#dev.off()


n_row = length(CORRELATION_DATA[,1])
n_col = n_row
choosen_mat = matrix(0,nrow=n_row,ncol=n_col)
order_mat   = matrix(0,nrow=n_row,ncol=n_col)

for(i in 1:n_row){
  print(i)
  choosen_all            = combination(CORRELATION_DATA,i,n_row)
  order_all              = order(-choosen_all)
  choosen_mat[i,]        = choosen_all
  #order_mat[i,]          = order_all
  #temp                   = as.vector(order_mat[i,1:5])  
  #choosen_mat[i,-(temp)] = 0
}

hist(choosen_mat)


cluster_factor = as.character(factor(LUSC_TOP100_COEFF_CLUSTER[,'kmeans_cluster'],labels= c("circle", "square", "csquare", "rectangle","crectangle","vrectangle")))
smoking_factor = as.character(factor(LUSC_TOP100_COEFF_CLUSTER[,'SMOKING'],labels= c("green", "blue", "red", "black")))
tumor_stage_factor = as.character(factor(LUSC_TOP100_COEFF_CLUSTER[,'TUMOR_STAGE_GROUPED'],labels= c("green", "blue", "red", "black")))

cluster_factor = as.character(factor(LUSC_TOP100_COEFF_CLUSTER[,'kmeans_cluster'],labels= c("gray50", "blue", "red", "tomato","lightsteelblue2","orange")))
hist(choosen_mat)
choosen_mat_sim = choosen_mat
choosen_mat_sim[choosen_mat_sim<0.6 & choosen_mat_sim>0.4] =0

g2 = graph_from_adjacency_matrix(choosen_mat_sim,weight=TRUE,mode="undirected")
V(g2)$color = cluster_factor
V(g2)$label = LUSC_TOP100_COEFF_CLUSTER[,'kmeans_cluster']
#V(g2)$shape = cluster_factor
plot(g2,layout=layout.kamada.kawai,vertex.size=30,vertex.label.dist=1.5,vertex.label.cex=0.5,edge.color="gray70",edge.width=2)

g2 = graph_from_adjacency_matrix(choosen_mat,weight=TRUE,mode="undirected")
fine = 50 # this will adjust the resolving power.
os_months = as.matrix(LUSC_TOP100_COEFF_CLUSTER['OS_MONTHS'])
V(g2)$color = brewer.pal(as.numeric(cut(os_months,breaks=fine)),"Blues")
V(g2)$shape = cluster_factor
plot(g2,vertex.size=10,vertex.label.cex=0.2)

g3 = graph_from_adjacency_matrix(choosen_mat,weight=TRUE,mode="undirected")
V(g3)$color = smoking_factor
V(g3)$shape = cluster_factor
plot(g3,vertex.size=10,vertex.label.cex=0.2)

g4 = graph_from_adjacency_matrix(choosen_mat,weight=TRUE,mode="undirected")
V(g4)$color = tumor_stage_factor
V(g4)$shape = cluster_factor
plot(g4,vertex.size=10,vertex.label.cex=0.2)


##############################################################

Cluster_data = function(data,cluster){
        data = data[data$kmeans_cluster == cluster,]
        data = data[c(103,1:101)]
        return(data)
}

LUSC_TOP100_1 = Cluster_data(LUSC_TOP100_COEFF_CLUSTER,1)
LUSC_TOP100_2 = Cluster_data(LUSC_TOP100_COEFF_CLUSTER,2)
LUSC_TOP100_3 = Cluster_data(LUSC_TOP100_COEFF_CLUSTER,3)
LUSC_TOP100_4 = Cluster_data(LUSC_TOP100_COEFF_CLUSTER,4)
LUSC_TOP100_5 = Cluster_data(LUSC_TOP100_COEFF_CLUSTER,5)
LUSC_TOP100_6 = Cluster_data(LUSC_TOP100_COEFF_CLUSTER,6)

LUAD_TOP100_1 = Cluster_data(LUAD_TOP100_COEFF_CLUSTER,1)
LUAD_TOP100_2 = Cluster_data(LUAD_TOP100_COEFF_CLUSTER,2)
LUAD_TOP100_3 = Cluster_data(LUAD_TOP100_COEFF_CLUSTER,3)
LUAD_TOP100_4 = Cluster_data(LUAD_TOP100_COEFF_CLUSTER,4)
LUAD_TOP100_5 = Cluster_data(LUAD_TOP100_COEFF_CLUSTER,5)
LUAD_TOP100_6 = Cluster_data(LUAD_TOP100_COEFF_CLUSTER,6)

write.table(LUSC_TOP100_1,"LUSC_TOP100_1.csv",sep=",")
write.table(LUSC_TOP100_2,"LUSC_TOP100_2.csv",sep=",")
write.table(LUSC_TOP100_3,"LUSC_TOP100_3.csv",sep=",")
write.table(LUSC_TOP100_4,"LUSC_TOP100_4.csv",sep=",")
write.table(LUSC_TOP100_5,"LUSC_TOP100_5.csv",sep=",")
write.table(LUSC_TOP100_6,"LUSC_TOP100_6.csv",sep=",")

write.table(LUAD_TOP100_1,"LUAD_TOP100_1.csv",sep=",")
write.table(LUAD_TOP100_2,"LUAD_TOP100_2.csv",sep=",")
write.table(LUAD_TOP100_3,"LUAD_TOP100_3.csv",sep=",")
write.table(LUAD_TOP100_4,"LUAD_TOP100_4.csv",sep=",")
write.table(LUAD_TOP100_5,"LUAD_TOP100_5.csv",sep=",")
write.table(LUAD_TOP100_6,"LUAD_TOP100_6.csv",sep=",")
