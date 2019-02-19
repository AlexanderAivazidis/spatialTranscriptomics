### Extract good marker genes for RNA-Scope from existing single-cell transcriptomic data

require(limma)
require(pvclust)

## Start with Allen Data:

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

remove = (rowSums(data != 0) < 4) # Remove unexpressed genes
data = data[!remove,]
rowdata = rowdata[!remove,]

## Get pairwise differentially expressed genes of each cell type:

data = log(data + 1,2)

save(data, file = "logTPM.RData")

geneList = list()
cellPairList = list()
count = 0
for (t1 in unique(celltypes)){
  for (t2 in unique(celltypes)){
    if (t1 != t2){
      print(c(t1,t2))
      count = count + 1
      cellPairList[[count]] = c(t1,t2)
      twoTypes = celltypes[celltypes %in% c(t1,t2)]
      design = 1*(twoTypes == t1)
      subData = data[,celltypes %in% c(t1,t2)]
      fit <- lmFit(subData, design)
      fit <- eBayes(fit, trend=TRUE)
      tab = topTable(fit, coef=ncol(design), number = 1000)
      tab = tab[tab[,'adj.P.Val'] < 0.01,]
      if (dim(tab)[1] > 0){
        tab = tab[tab[,'logFC'] > 2,]
        if (dim(tab)[1] > 50){
          tab = tab[1:50,]
        }
      }
      geneList[[count]] = rownames(tab)
    }
  }
}

allGenes = unique(unlist(geneList))

data = data[allGenes,]

# Further filtering down to binary expressed genes:

data = 2^(data)-1
binaryData = 1*(data >= 1)

Pmatrix = matrix(0,dim(data)[1],length(unique(celltypes)))
colnames(Pmatrix) = unique(celltypes)
for (i in 1:dim(Pmatrix)[2]){
  print(i)
  Pmatrix[,i] = apply(binaryData[,celltypes == unique(celltypes)[i]],1,function(x) mean(x))
}
rownames(Pmatrix) = rownames(data)

for (count in 1:length(geneList)){
  print(round(count/length(geneList),2))
  geneList[[count]] = geneList[[count]][Pmatrix[geneList[[count]],cellPairList[[count]][1]] > 0.5]
  Pratio = (Pmatrix[geneList[[count]],cellPairList[[count]][1]] - Pmatrix[geneList[[count]],cellPairList[[count]][2]])/
            apply(Pmatrix[geneList[[count]],cellPairList[[count]]],1,function(x) max(x))
  geneList[[count]] = geneList[[count]][Pratio > 0.7]
}

allGenes = unique(unlist(geneList))

save(allGenes, file = 'allGenes.RData')
save(geneList, file = 'geneList.RData')

data = data[allGenes,]
data = log(data+1,2)

## Cluster into taxonomy:

medianMatrix = matrix(0,dim(data)[1],length(unique(celltypes)))
colnames(medianMatrix) = unique(celltypes)
for (i in 1:dim(medianMatrix)[2]){
  print(i)
  medianMatrix[,i] = apply(data[,celltypes == unique(celltypes)[i]],1,function(x) median(x))
}
rownames(medianMatrix) = rownames(data)

corMatrix = cor(medianMatrix)
res = hclust(dist(corMatrix), method = 'average')

save(res, file = 'res.RData')

res =  pvclust(medianMatrix, method.dist="cor", method.hclust="average", nboot=100,
               parallel = TRUE)

## Chose marker genes for each split in the tree:

## Train an SVM on the top 12 DE genes:

## Now run an genetic optimization algorithm to select the best 12 marker genes:


