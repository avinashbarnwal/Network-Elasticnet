setwd('/Users/avinashbarnwal/Desktop/Research/PhD/network_elasticnet/NYSDS-2019/')
path     = getwd()
dataPath = paste(path,"/result/final_result.csv",sep="")
data     = read.table(dataPath,header=TRUE,sep=",")

#Missing Value Treatment#
data[is.na(data['TUMOR_STAGE']),'TUMOR_STAGE'] = 'Stage IIIA'
data_matrix = t(data.matrix(data[,2:102]))
#plot(data[,2:102])
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
col_breaks <- c(seq(-1,0,length=100),      # for red
               seq(0,0.8,length=100),      # for yellow
               seq(0.81,1,length=100))     # for green
nba_heatmap <- heatmap(data_matrix, Rowv=NA, Colv=NA, col = my_palette, breaks=col_breaks, scale="column")
