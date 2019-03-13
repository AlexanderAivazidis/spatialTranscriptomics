# Performance on mapping 3 separate branches of glutamatergic neurons with a panel of 12 genes for different weighting matrix (ALM or VISp).
# 3 different weight matrices: (All cell types), non-neuronal cells weight 0, non-glutamatergic neurons weight 0, not in branch neurons weight 0

getWeightMatrix = function(initialWeightMatrix, clustersToMerge, weight = 0){
  initialWeightMatrix[clustersToMerge, clustersToMerge] = initialWeightMatrix[clustersToMerge, clustersToMerge]*weight
  return(initialWeightMatrix)
}

getFishPanel = function(mapDat, runGenes, clusters, medianExpr, weightMatrix, panelSize, panelMin = 1, subSamp = 1000,
                                 startingPanel = c(), focusGroup = NA, lossFunction = 'CorrelationDistance'){
  
  require(mfishtools)
  if (is.na(focusGroup)){
    focusGroup = unique(clusters)
  }
  
  print(focusGroup)
  
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
  
  #FList = FscoreWithGenes(fishPanel, mapDat, medianExpr, clusters, focusGroup)
  
  #return(list(fishPanel, FList))
  return(fishPanel)
}

### Main ###

require(AllenData)
require(mfishtools)
require(matrixStats)
dataDirectory = '/nfs/team205/aa16/AllenData/'
savingDirectory = ''

allData = loadAllenData(cortical_area = 'ALM', species = 'mouse', normalization = 'exon+intron_cpm', directory = dataDirectory)
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
# We load a hierarchical clustering result of the data to define our focusGroup and clusters to merge:
load('res.RData') # Load a hierarchical clustering result of the data
dendrogramOrder = res[[4]][res[[3]]]
focusGroupsList = list(c(dendrogramOrder[3:10], dendrogramOrder[31:38]), dendrogramOrder[11:30], dendrogramOrder[39:57])
weightMatrix = matrix(1,length(unique(specific_type)), length(unique(specific_type)))
colnames(weightMatrix) = rownames(weightMatrix) = unique(specific_type)
lossFunctionList = c('negative F-Score', 'Correlation Distance')

i = 1
### This is the list we will iterate trough to see if focussing on branches closer to our group of interest (focusGroup) will improve performance 
clustersToMergeList = list(c(), dendrogramOrder[(length(dendrogramOrder)-15):length(dendrogramOrder)],
                           dendrogramOrder[c(1:2, 58:length(dendrogramOrder))], dendrogramOrder[!dendrogramOrder %in% focusGroupsList[[i]]])

runGenes <- filterPanelGenes(
  summaryExpr = 2^medianExpr-1,  # medians (could also try means); We enter linear values to match the linear limits below
  propExpr    = propExpr,    # proportions
  startingGenes  = c(),  # Starting genes 
  numBinaryGenes = 1000,      # Number of binary genes 
  onClusters = focusGroupsList[[i]],
  minOn     = minOn,   # Minimum required expression in highest expressing cell type
  maxOn     = maxOn,  # Maximum allowed expression
  fractionOnClusters = 0.5,  # Max fraction of on clusters 
  excludeFamilies = c("LOC","Fam","RIK","RPS","RPL","\\-","Gm","Rnf","BC0")) # Avoid LOC markers, in this case

for (j in 2:length(clustersToMergeList)){
 newWeightMatrix = getWeightMatrix(weightMatrix, clustersToMergeList[[j]][clustersToMergeList[[j]] %in% specific_type], weight = 0)
 
 k = 1
   fishPanel = getFishPanel(data, runGenes, specific_type, medianExpr, newWeightMatrix, panelSize = 12, focusGroup =unique(specific_type[!specific_type %in% clustersToMergeList[[j]][clustersToMergeList[[j]] %in% specific_type]]), 
                                    panelMin = 1, subSamp = 50, startingPanel = c(), lossFunction = lossFunctionList[[k]])
   save(fishPanel, file = paste(savingDirectory, 'markerGenes/allen_markerGenesALM_','focusGroup', as.character(i), 'clustersToMerge', as.character(j),'lossFunction', as.character(k), '.RData', sep = ""))

}


