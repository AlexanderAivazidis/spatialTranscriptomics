### Get cell type taxonomy, then marker genes for the kidney data
require(Matrix)
require(myUtils)
require(mfishtools)
require(matrixStats)

directory = '/nfs/team205/aa16/KidneyData/'
capture_efficiency = 0.145

coldata = read.delim('/nfs/team205/aa16/KidneyData/cellManifestCompressed.tsv')
data = readMM('/nfs/team205/aa16/KidneyData/tableOfCounts.mtx')
sampleNames = read.delim('/nfs/team205/aa16/KidneyData/tableOfCounts_colLabels.tsv')
geneNames = read.delim('/nfs/team205/aa16/KidneyData/tableOfCounts_rowLabels.tsv')

colnames(data) = sampleNames[,2]
rownames(data) = geneNames[,2]

# Split Kidney Data into the different tissues:

dropletID_List = list()
for (i in 1:length(unique(coldata[,'Compartment']))){
  dropletID_List[[i]] = coldata[coldata[,'Compartment'] == unique(coldata[,'Compartment'])[i],'DropletID']
}
names(dropletID_List) = unique(coldata[,'Compartment'])

for (i in c(1,2,3,5,6,7,8)){
  print(names(dropletID_List)[i])
  dense_matrix = as.matrix(data[,colnames(data) %in% dropletID_List[[i]]])
  dense_matrix = round(cpm(dense_matrix),2)
  write.table(dense_matrix, file = paste(directory, 'cpmMatrix_', names(dropletID_List)[i], '.txt', sep = ''), sep = '\t', quote = F)
  rm(dense_matrix) 
}
rm(data)

for (i in c(1,2,3,5,6,7,8)){
  print(names(dropletID_List)[i])
  coldata_subset = coldata[coldata[,'DropletID'] %in% dropletID_List[[i]],]
  write.table(coldata_subset, file = paste(directory, 'coldata_', names(dropletID_List)[i], '.txt', sep = ''), sep = '\t', quote = F)
  rm(coldata_subset) 
}

# Cell type taxonomy:

getCellTaxonomy = function(data, celltypes, verbose = TRUE, type = 'mean'){
  
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
    if (length(geneList[[count]]) > 1){
      Pratio = (Pmatrix[geneList[[count]],cellPairList[[count]][1]] - Pmatrix[geneList[[count]],cellPairList[[count]][2]])/
        apply(Pmatrix[geneList[[count]],cellPairList[[count]]],1,function(x) max(x))
    }else{
      Pratio = (Pmatrix[geneList[[count]],cellPairList[[count]][1]] - Pmatrix[geneList[[count]],cellPairList[[count]][2]])/
        max(Pmatrix[geneList[[count]],cellPairList[[count]]])
    }
    geneList[[count]] = geneList[[count]][Pratio > 0.7]
  }
  
  allGenes = unique(unlist(geneList))
  
  data = data[allGenes,]
  data = log(data+1,2)
  
  ## Cluster into taxonomy:
  
  if (verbose == TRUE){
    print('Clustering into cell taxonomy...') 
  }
  
  if (type == 'median'){
    medianMatrix = matrix(0,dim(data)[1],length(unique(celltypes)))
    colnames(medianMatrix) = unique(celltypes)
    for (i in 1:dim(medianMatrix)[2]){
      print(i)
      medianMatrix[,i] = apply(data[,celltypes == unique(celltypes)[i]],1,function(x) median(x))
    }
    rownames(medianMatrix) = rownames(data)
  }else{
    medianMatrix = matrix(0,dim(data)[1],length(unique(celltypes)))
    colnames(medianMatrix) = unique(celltypes)
    for (i in 1:dim(medianMatrix)[2]){
      print(i)
      medianMatrix[,i] = apply(data[,celltypes == unique(celltypes)[i]],1,function(x) mean(x))
    }
    rownames(medianMatrix) = rownames(data)
  }

  
  corMatrix = cor(medianMatrix)
  res = hclust(dist(corMatrix), method = 'average')
  
  # 
  # res =  pvclust(medianMatrix, method.dist="cor", method.hclust="average", nboot=100,
  #                parallel = TRUE)
  
  return(res)
}

for (i in c(1,2,3,5,6,7,8)){
  print(names(dropletID_List)[i])
  dense_matrix = read.delim(paste(directory, 'cpmMatrix_', names(dropletID_List)[i], '.txt', sep = ''), sep = '\t')
  cell_taxonomy = getCellTaxonomy(dense_matrix, coldata[coldata[,'Compartment'] == names(dropletID_List)[i],'ClusterID'])
  save(cell_taxonomy, file = paste(directory, 'cellTaxonomy_', names(dropletID_List)[i],'.RData', sep = ''))
  rm(dense_matrix)
}

### Now select 36 markers for each dataset:

options(stringsAsFactors = FALSE)

for (i in c(1,2,3,4,5,6,7,8)){
  dense_matrix = read.delim(paste(directory, 'cpmMatrix_', names(dropletID_List)[i], '.txt', sep = ''), sep = '\t')
  colnames(dense_matrix) = substring(colnames(dense_matrix),2)
  coldata = read.delim(paste(directory, 'coldata_', names(dropletID_List)[i], '.txt', sep = ''))
  normDat = log(dense_matrix+1,2)
  
  cl = coldata[,'ClusterID']
  names(cl) = dropletID_List[[i]]
  
  exprThresh = 1
  medianExpr = do.call("cbind", tapply(names(cl), as.factor(cl), function(x) apply(normDat[,x],1, function(y) mean(y)))) 
  propExpr   = do.call("cbind", tapply(names(cl), cl, function(x) rowMeans(normDat[,x]>exprThresh))) 
  rownames(medianExpr) <- rownames(propExpr) <- genes <- rownames(normDat) 
  
  nonZeroMedian = (!rowSums(medianExpr) == 0)
  medianExpr = medianExpr[nonZeroMedian,]
  propExpr = propExpr[nonZeroMedian,]
  
  runGenes <- filterPanelGenes(
    summaryExpr = 2^medianExpr-1,  # medians (could also try means); We enter linear values to match the linear limits below
    propExpr    = propExpr,    # proportions
    geneLengths = NULL,        # vector of gene lengths (not included here)
    numBinaryGenes = 500,      # Number of binary genes (explained below)
    startingGenes = c(),
    minOn     = 5,   # Minimum required expression in highest expressing cell type
    maxOn     = 1000,  # Maximum allowed expression
    fractionOnClusters = 0.5,  # Max fraction of on clusters (described above)
    excludeGenes    = NULL,    # Genes to exclude.  Often sex chromosome or mitochondrial genes would be input here.
    excludeFamilies = c("LOC","Fam","RIK","RPS","RPL","\\-","Gm","Rnf","BC0")) # Avoid LOC markers, in this case
  
  corDist         <- function(x) return(as.dist(1-cor(x)))
  clusterDistance <- as.matrix(corDist(medianExpr[runGenes,]))
  print(dim(clusterDistance))
  
  fishPanel <- buildMappingBasedMarkerPanel(
    mapDat        = normDat[runGenes,],     # Data for optimization
    medianDat     = medianExpr[runGenes,], # Median expression levels of relevant genes in relevant clusters
    clustersF     = cl,                   # Vector of cluster assignments
    panelSize     = 6,                           # Final panel size
    currentPanel  = NULL,            # Starting gene panel
    subSamp       = 10,
    panelMin = 1,# Maximum number of cells per cluster to include in analysis (20-50 is usually best)
    optimize      = "CorrelationDistance",        # CorrelationDistance maximizes the cluster distance as described
    clusterDistance = clusterDistance,            # Cluster distance matrix
    percentSubset = 10                           # Only consider a certain percent of genes each iteration to speed up calculations (in most cases this is not recommeded)
  ) 
  
  save(fishPanel, file = paste('markerGenes/fishPanel_KidneyProject_', names(dropletID_List)[i], '.RData'))
}
