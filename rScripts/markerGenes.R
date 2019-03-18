loadAllenData = function(cortical_area, normalization = 'exon+intron_cpm', directory = '/nfs/team205/aa16/AllenData/'){
  
  if (!normalization %in% c('exon_rpkm', 'exon+intron_cpm'))
  {
    stop('No valid normalization argument supplied.')
  }
  
  dataList = list()
  coldataList = list()
  for (i in 1:length(cortical_area)){
    dataList[[i]] = read.delim(paste(directory, 'mouse_', cortical_area[i], '_2018-06-14_', normalization, '-matrix.csv', sep = ''), sep = ',')
    rownames(dataList[[i]]) = dataList[[i]][,1]
    dataList[[i]] = dataList[[i]][,2:dim(dataList[[i]])[2]]
    coldataList[[i]] = read.delim(paste(directory, 'mouse_', cortical_area[i], '_2018-06-14_samples-columns.csv', sep = ''), sep = ',')
    rowdata = read.delim(paste(directory, 'mouse_', cortical_area[i], '_2018-06-14_genes-rows.csv', sep = ''), sep = ',')
  }
  
  data = dataList[[1]]
  coldata = coldataList[[1]]
  if (length(dataList) > 1){
    for (i in 1:(length(dataList)-1)){
      data = cbind(data, dataList[[i+1]])
      coldata = rbind(coldata, coldataList[[i+1]])
    }
  }
  
  rm(dataList)
  rm(coldataList)
  
  celltypes = coldata[,'cluster']
  splitted = lapply(as.character(celltypes), function(x) strsplit(x, " "))
  firstEntry = unlist(lapply(splitted, function(x) x[[1]][1]))
  
  remove = (celltypes == 'High Intronic VISp L5 Endou' | firstEntry == 'Low' | firstEntry == 'Batch' | firstEntry == 'Doublet') # Remove low quality cells
  data = data[,!remove]
  coldata = coldata[!remove,]
  
  remove = (rowSums(data != 0) == 0) # Remove unexpressed genes
  data = data[!remove,]
  rowdata = rowdata[!remove,]
  
  return(list(data, coldata, rowdata))
}

getCellTaxonomy = function(data, celltypes, verbose = TRUE){
  
  require(limma)
  require(pvclust)
  
  data = log(data + 1,2)
  if (verbose == TRUE){
    print('Getting pairwise differentially expressed genes...') 
  }
  geneList = list()
  cellPairList = list()
  count = 0
  for (t1 in unique(celltypes)){
    for (t2 in unique(celltypes)){
      if (t1 != t2){
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
  
  if (verbose == TRUE){
    print('Further filtering down to binary expressed genes...') 
  }
  data = 2^(data)-1
  binaryData = 1*(data >= 1)
  
  Pmatrix = matrix(0,dim(data)[1],length(unique(celltypes)))
  colnames(Pmatrix) = unique(celltypes)
  for (i in 1:dim(Pmatrix)[2]){
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
  
  data = data[allGenes,]
  data = log(data+1,2)
  
  ## Cluster into taxonomy:
  
  if (verbose == TRUE){
    print('Clustering into cell taxonomy...') 
  }
  
  medianMatrix = matrix(0,dim(data)[1],length(unique(celltypes)))
  colnames(medianMatrix) = unique(celltypes)
  for (i in 1:dim(medianMatrix)[2]){
    print(i)
    medianMatrix[,i] = apply(data[,celltypes == unique(celltypes)[i]],1,function(x) median(x))
  }
  rownames(medianMatrix) = rownames(data)
  
  corMatrix = cor(medianMatrix)
  res = hclust(dist(corMatrix), method = 'average')
  
  # 
  # res =  pvclust(medianMatrix, method.dist="cor", method.hclust="average", nboot=100,
  #                parallel = TRUE)
  
  return(res)
}

### Make a cell taxonomy tree for each Allen data set seperately:

allCorticalAreas = c('ALM', 'VISp')
# for (cortical_area in allCorticalAreas){
#   resList = loadAllenData(cortical_area)
#   data = resList[[1]]
#   coldata = resList[[2]]
#   rowdata = resList[[3]]
#   cell_taxonomy = getCellTaxonomy(data, coldata[,'cluster'])
#   save(cell_taxonomy, file = paste('/nfs/team205/aa16/AllenData/celltypeTaxonomy', cortical_area, '.RData', sep = ''))
# }

### And for all of them combined:

resList = loadAllenData(allCorticalAreas)
data = resList[[1]]
coldata = resList[[2]]
rowdata = resList[[3]]
cell_taxonomy = getCellTaxonomy(data, coldata[,'cluster'])
save(cell_taxonomy, file = paste('celltypeTaxonomy',
paste(as.character(allCorticalAreas), collapse = '-'), '.RData', sep = ''))

