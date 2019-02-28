### Check spatial astrocyte genes for expression in neurons

require(ComplexHeatmap)
require(matrixStats)

data1 = read.delim('/nfs/team205/aa16/AllenData/mouse_ALM_2018-06-14_exon+intron_cpm-matrix.csv', sep = ',')
rownames(data1) = data1[,1]
data1 = data1[,2:dim(data1)[2]]
coldata1 = read.delim('/nfs/team205/aa16/AllenData/mouse_ALM_2018-06-14_samples-columns.csv', sep = ',')
rowdata1 = read.delim('/nfs/team205/aa16/AllenData/mouse_ALM_2018-06-14_genes-rows.csv', sep = ',')

data2 = read.delim('/nfs/team205/aa16/AllenData/mouse_VISp_2018-06-14_exon+intron_cpm-matrix.csv', sep = ',')
rownames(data2) = data2[,1]
data2 = data2[,2:dim(data2)[2]]
coldata2 = read.delim('/nfs/team205/aa16/AllenData/mouse_VISp_2018-06-14_samples-columns.csv', sep = ',')
rowdata2 = read.delim('/nfs/team205/aa16/AllenData/mouse_VISp_2018-06-14_genes-rows.csv', sep = ',')

data = cbind(data1,data2)
coldata = rbind(coldata1,coldata2)
rowdata = rowdata1

rm(data1,data2)
rm(coldata1,coldata2)
rm(rowdata1,rowdata2)

celltypes = coldata[,'cluster']
splitted = lapply(as.character(celltypes), function(x) strsplit(x, " "))
firstEntry = unlist(lapply(splitted, function(x) x[[1]][1]))

remove = (celltypes == 'High Intronic VISp L5 Endou' | firstEntry == 'Low' | firstEntry == 'Batch' | firstEntry == 'Doublet') # Remove low quality cells
data = data[,!remove]
coldata = coldata[!remove,]
celltypes = coldata[,'cluster']

data = log(data + 1,2)

#setwd('/Users/aa16/rotation2/spatialTranscriptomics/')
tab = read.table('Holt_astrocytes_spatial_map_Shalev_p56_p_02.csv', sep = ',')
genes = tab[2:dim(tab)[1],6]

rownames(data) = rowdata1

missing_genes = genes[!genes %in% rownames(data)]
genes = genes[genes %in% rownames(data)]

genes = as.character(genes)

subset = (coldata[,'class'] %in% c('Glutamatergic', 'GABAergic'))
data = data[,subset]
coldata = coldata[subset,]

data = as.matrix(data)

layers = c('L1', 'L2/3', 'L4', 'L5', 'L6')
layer_median = matrix(0,length(layers), length(genes))
rownames(layer_median) = layers
colnames(layer_median) = genes
for (i in 1:length(layers)){
  print(i)
  layer_median[layers[i], genes] = rowMedians(data[genes,coldata[,'brain_subregion'] == layers[i]])
}
  
clusters = as.character(unique(coldata[,'cluster']))
cluster_median = matrix(0,length(clusters), length(genes))
rownames(cluster_median) = clusters
colnames(cluster_median) = genes
for (i in 1:length(clusters)){
  print(i)
  cluster_median[clusters[i], genes] = rowMedians(data[genes,coldata[,'cluster'] == clusters[i]])
}

# Display expression of marker genes per layer:
ht1 = Heatmap(layer_median, name = "log2(CPM + 1)",
              column_title = 'Median Expression per Layer of Astrocyte Zonation Genes \n in all Neurons of AML and VISp',
              cluster_columns = TRUE, cluster_rows = FALSE)
pdf(file = 'AstrocyteZonationGenesP02_ExpressionPerLayer.pdf', width = 100, height = 14)
print(ht1)
dev.off()

# Display expression of marker genes per neuron subtype:
ht1 = Heatmap(cluster_median, name = "log2(CPM + 1)",
              column_title = 'Median Expression of Astrocyte Zonation Genes \n in Neuron Subtypes in AML and VISp',
              cluster_columns = TRUE, cluster_rows = FALSE)
pdf(file = 'AstrocyteZonationGenesP02_ExpressionInEachNeuronSubtype.pdf', width = 100, height = 21)
print(ht1)
dev.off()

