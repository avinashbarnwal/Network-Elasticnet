setwd('/Users/avinashbarnwal/Desktop/Research/PhD/Network_Lasso/NYSDS-2019/')

path     = getwd()
dataPath = paste(path,"/result/final_result.csv",sep="")
data     = read.table(dataPath,header=TRUE,sep=",")

#Missing Value Treatment
data[is.na(data['TUMOR_STAGE']),'TUMOR_STAGE'] = 'Stage IIIA'
