# Poster Figure 1: This figure shows performance on mapping all 117 neuron types in the Allen mouse RNAseq dataset with a panel of 1 to 96 genes.

getWeightMatrix = function(initialWeightMatrix, clustersToMerge, weight = 0){
  initialWeightMatrix[clustersToMerge, clustersToMerge] = initialWeightMatrix[clustersToMerge, clustersToMerge]*weight
  return(initialWeightMatrix)
}

getFishPanelAndFScore = function(mapDat, runGenes, clusters, medianExpr, weightMatrix, panelSize, panelMin = 1, subSamp = 1000,
                                 startingPanel = c(), focusGroup = NA, lossFunction = 'CorrelationDistance'){
  
  require(mfishtools)
  if (is.na(focusGroup)){
    focusGroup = unique(clusters)
  }
  
  corDist         <- function(x) return(as.dist(1-cor(x)))
  clusterDistance <- as.matrix(corDist(medianExpr))
  weightMatrix = weightMatrix[rownames(clusterDistance), colnames(clusterDistance)]
  clusterDistance = clusterDistance * weightMatrix
  
  fishPanel <- buildMappingBasedMarkerPanel(
    mapDat        = mapDat[runGenes,],                # Data for optimization
    medianDat     = medianExpr[runGenes,],            # Median expression levels of relevant genes in relevant clusters
    clustersF     = clusters,                         # Vector of cluster assignments
    panelSize     = panelSize,                        # Final panel size
    currentPanel  = startingPanel,                    # Starting gene panel
    subSamp       = subSamp,                          # Maximum number of cells per cluster to include in analysis (20-50 is usually best)
    panelMin      = panelMin,
    optimize      = lossFunction,                     # CorrelationDistance maximizes the cluster distance as described
    clusterDistance = clusterDistance,                # Cluster distance matrix (potentiall multiplied by weight matrix)
    focusGroup = focusGroup,
    percentSubset = 100                               # Only consider a certain percent of genes each iteration to speed up calculations (in most cases this is not recommeded)
  )
  
  FList = FscoreWithGenes(fishPanel, mapDat, medianExpr, clusters, focusGroup)
  
  return(list(fishPanel, FList))
}

### Main ###

require(AllenData)
require(mfishtools)
require(matrixStats)
dataDirectory = '/nfs/team205/aa16/AllenData/'
savingDirectory = ''

allData = loadAllenData(cortical_area = c('ALM','VISp'), species = 'mouse', normalization = 'exon+intron_cpm', directory = dataDirectory)
data = as.matrix(log(allData[[1]]+1,2))
coldata = allData[[2]]
rowdata = allData[[3]]
rownames(data) = rowdata[,1]
data = data[,table(coldata[,'cluster'])[coldata[,'cluster']] > 1]
coldata = coldata[table(coldata[,'cluster'])[coldata[,'cluster']] > 1,]

# Prepare some data needed for fishPanel generation:
specific_type = coldata[,'cluster']
names(specific_type) = coldata[,1]
exprThresh = 1
meanExpr = do.call("cbind", tapply(names(specific_type), specific_type, function(x) rowMeans(data[,x])))
medianExpr = do.call("cbind", tapply(names(specific_type), specific_type, function(x) rowMedians(data[,x])))
propExpr   = do.call("cbind", tapply(names(specific_type), specific_type, function(x) rowMeans(data[,x]>exprThresh)))
rownames(medianExpr) <- rownames(propExpr) <- genes <- rownames(data)
nonZeroMedian = (!rowSums(medianExpr) == 0)
medianExpr = medianExpr[nonZeroMedian,]
propExpr = propExpr[nonZeroMedian,]
maxGene = 'Rorb'
maxOn = 2^max(medianExpr[maxGene,])-1
minOn = 0.1
focusGroup = unique(specific_type[coldata['class'] == 'Glutamatergic' | coldata['class'] == 'GABAergic'])
weightMatrix = matrix(1,length(unique(specific_type)), length(unique(specific_type)))
colnames(weightMatrix) = rownames(weightMatrix) = unique(specific_type)
lossFunctionList = c('negative F-Score', 'Correlation Distance')
k = 1

runGenes <- filterPanelGenes(
  summaryExpr = 2^medianExpr-1,  # medians (could also try means); We enter linear values to match the linear limits below
  propExpr    = propExpr,    # proportions
  startingGenes  = c(),  # Starting genes 
  numBinaryGenes = 1000,      # Number of binary genes 
  onClusters = as.character(focusGroup),
  minOn     = minOn,   # Minimum required expression in highest expressing cell type
  maxOn     = maxOn,  # Maximum allowed expression
  fractionOnClusters = 0.5,  # Max fraction of on clusters 
  excludeFamilies = c("LOC","Fam","RIK","RPS","RPL","\\-","Gm","Rnf","BC0")) # Avoid LOC markers, in this case

clustersToMerge = unique(!specific_type %in% focusGroup)
newWeightMatrix = getWeightMatrix(weightMatrix, clustersToMerge, weight = 0)

panelSize = 100
results = list(c(),c())

corDist         <- function(x) return(as.dist(1-cor(x)))
clusterDistance <- as.matrix(corDist(medianExpr))
weightMatrix = weightMatrix[rownames(clusterDistance), colnames(clusterDistance)]
clusterDistance = clusterDistance * newWeightMatrix

fishPanel = c()
for (pS in 1:panelSize){
  print(pS)
  
  fishPanel <- buildMappingBasedMarkerPanel(
    mapDat        = data[runGenes,],                # Data for optimization
    medianDat     = medianExpr[runGenes,],            # Median expression levels of relevant genes in relevant clusters
    clustersF     = specific_type,                         # Vector of cluster assignments
    panelSize     = pS,                               # Final panel size
    currentPanel  = fishPanel,                        # Starting gene panel
    subSamp       = 50,                          # Maximum number of cells per cluster to include in analysis (20-50 is usually best)
    panelMin      = 1,
    optimize      = lossFunctionList[[k]],                     # CorrelationDistance maximizes the cluster distance as described
    clusterDistance = clusterDistance,                # Cluster distance matrix (potentiall multiplied by weight matrix)
    focusGroup = focusGroup,
    percentSubset = 100                               # Only consider a certain percent of genes each iteration to speed up calculations (in most cases this is not recommeded)
  )
  
  save(fishPanel, file = paste(savingDirectory, 'markerGenes/allen_markerGenesNeurons_ALM-VISp_100genes_lossFunction', as.character(k), '.RData', sep = ""))
}



