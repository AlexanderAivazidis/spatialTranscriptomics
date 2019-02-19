# Get marker genes using mnfishtools
# We use the ALM and VISp datasets from the Allen Institute (seperate first, then combined)
# and compare the results to Omers manually selected panel:

omersGenes = c("Cux2", "Fam84b", "Osr1", "Colq", "Batf3", "Slc17a8", "Scnn1a",
               "Rorb", "Foxp2", "Bcl11b", "Rspo1", "Hsd11b1")

corticalArea = 'VISp'
layer = 'allLayers'

library(matrixStats)   # For rowMedians function, which is fast
require(mfishtools)
library(ComplexHeatmap)
library(gplots)

options(stringsAsFactors = FALSE)  # IMPORTANT

data1 = read.delim('/nfs/team205/aa16/AllenData/mouse_VISp_2018-06-14_exon+intron_cpm-matrix.csv', sep = ',')
rownames(data1) = data1[,1]
data1 = data1[,2:dim(data1)[2]]
data1 = log2(data1+1)

coldata1 = read.delim('/nfs/team205/aa16/AllenData/mouse_VISp_2018-06-14_samples-columns.csv', sep = ',')
rowdata1 = read.delim('/nfs/team205/aa16/AllenData/mouse_VISp_2018-06-14_genes-rows.csv', sep = ',')
rownames(data1) = rowdata1[,1]
data1 = as.matrix(data1)

# Remove celltypes that occur only once
data1 = data1[,table(coldata1[,'cluster'])[coldata1[,'cluster']] > 1]
coldata1 = coldata1[table(coldata1[,'cluster'])[coldata1[,'cluster']] > 1,]
# Get low quality clusters:
celltypes = coldata1[,'cluster']
splitted = lapply(as.character(celltypes), function(x) strsplit(x, " "))
firstEntry = unlist(lapply(splitted, function(x) x[[1]][1]))
lowQuality = (firstEntry == 'Low') # Remove low quality cells
data1 = data1[,!lowQuality]
coldata1 = coldata1[!lowQuality,]

for (layer in c('L1', 'L2/3', 'L4', 'L5', 'L6')){
  
  print(layer)
  broad_type = coldata1[,'class']
  names(broad_type) = coldata1[,'sample_name']
  specific_type = coldata1[,'cluster']
  names(specific_type)   = coldata1[,'sample_name']
  allLayers = coldata1[,'brain_subregion']
  names(allLayers) = coldata1[,'sample_name']
  # Now we define a restricted set of cell types that should be the target for mapping
  specificClass = "Glutamatergic" # For this broad class we want to map each subtype
  broadClass = "GABAergic" # For these cells we just want to map the broad class
  specific_type[broad_type != specificClass] = 'Non-neuronal' # All other cells will be called 'Other'
  specific_type[broad_type == broadClass] = "GABA-ergic"
  # Reduce the data to layer specific cells:
  if (layer == 'allLayers'){
    subset = rep(TRUE,length(allLayers))
  }else{
    subset = (allLayers == layer | allLayers == paste(layer,'a', sep = "") | allLayers == paste(layer,'b', sep = ""))
  }
  subsetNames = unique(specific_type[subset])
  
  if (layer == 'L2/3'){layer = 'L2-3'} # Just need this for saving files later on
  
  exprThresh = 1
  meanExpr = do.call("cbind", tapply(names(specific_type), specific_type, function(x) rowMeans(data1[,x])))
  medianExpr = do.call("cbind", tapply(names(specific_type), specific_type, function(x) rowMedians(data1[,x])))
  propExpr   = do.call("cbind", tapply(names(specific_type), specific_type, function(x) rowMeans(data1[,x]>exprThresh)))
  rownames(medianExpr) <- rownames(propExpr) <- genes <- rownames(data1)
  
  nonZeroMedian = (!rowSums(medianExpr) == 0)
  
  runGenes <- filterPanelGenes(
    summaryExpr = 2^medianExpr[nonZeroMedian,subsetNames]-1,  # medians (could also try means); We enter linear values to match the linear limits below
    propExpr    = propExpr[nonZeroMedian,subsetNames],    # proportions
    startingGenes  = c(),  # Starting genes (from above)
    numBinaryGenes = 500,      # Number of binary genes (explained below)
    minOn     = 1,   # Minimum required expression in highest expressing cell type
    maxOn     = 500,  # Maximum allowed expression
    fractionOnClusters = 0.5,  # Max fraction of on clusters (described above)
    excludeGenes    = NULL,    # Genes to exclude. Often sex chromosome or mitochondrial genes would be input here.
    excludeFamilies = c("LOC","Fam","RIK","RPS","RPL","\\-","Gm","Rnf","BC0")) # Avoid LOC markers, in this case
  
  corDist         <- function(x) return(as.dist(1-cor(x)))
  clusterDistance <- as.matrix(corDist(medianExpr[,subsetNames]))
  
  fishPanel <- buildMappingBasedMarkerPanel(
    mapDat        = data1[runGenes,subset],         # Data for optimization
    medianDat     = medianExpr[runGenes,subsetNames], # Median expression levels of relevant genes in relevant clusters
    clustersF     = specific_type[subset],            # Vector of cluster assignments
    panelSize     = 12,                               # Final panel size
    currentPanel  = c(),                              # Starting gene panel
    subSamp       = 20,                               # Maximum number of cells per cluster to include in analysis (20-50 is usually best)
    optimize      = "CorrelationDistance",            # CorrelationDistance maximizes the cluster distance as described
    clusterDistance = clusterDistance,                # Cluster distance matrix
    percentSubset = 100                               # Only consider a certain percent of genes each iteration to speed up calculations (in most cases this is not recommeded)
  )
  
  frac = fractionCorrectWithGenes(fishPanel, data1[, subset], medianExpr[,subsetNames], specific_type[subset],
                                  return=TRUE, plot = FALSE)
  
  fracOmer = fractionCorrectWithGenes(omersGenes, data1[,subset],medianExpr[,subsetNames],specific_type[subset],
                                      return=TRUE, plot = FALSE)
  
  # For plotting put the clusters in the dendrogram order:
  
  corMatrix = cor(meanExpr)
  res = hclust(dist(corMatrix), method = 'average')
  order = res[[4]][res[[4]] %in% subsetNames]
  
  # Display accuracy as function of included genes and Compare to Omer's gene list:
  pdf(file = paste('/nfs/users/nfs_a/aa16/spatialTranscriptomics/figures/allen_', corticalArea, '_', layer, '_overallPerformance.pdf', sep = ""), width = 6, height = 6)
  plot(1:length(frac), frac,
       type = "l", col = "grey", xlab = "Number of genes in panel",
       ylab = "Percent of cells correctly mapped", main = paste("Mapping Accuracy for Allen ", corticalArea, " ", layer ," data", sep = ""), ylim =  c(-10, 100), lwd = 5)
  lines(1:length(fracOmer), fracOmer, col = 'orange', lwd = 5)
  abline(h = (-2:20) * 5, lty = "dotted", col = 'grey')
  abline(h = 0, col = "black", lwd = 2)
  text(1:length(frac), frac+15, fishPanel, srt = 90, cex = 1)
  text(1:length(fracOmer), fracOmer-9, omersGenes, srt = 90, cex = 1, col = 'red')
  legend(1, 95, legend=c("Greedy Algorithm Selection", "Omer's Selection"),
         col=c("grey", "orange"),lwd = 5, cex=0.8)
  dev.off()
  
  # Display expression of marker genes:
  ht1 = Heatmap(meanExpr[fishPanel, order][,order %in% subsetNames], name = "1",
                column_title = '1: Algorithm Selection',
                cluster_columns = FALSE, cluster_rows = FALSE)
  ht2 = Heatmap(meanExpr[omersGenes, order][,order %in% subsetNames], name = "2",
                column_title = '2: Manual (Omer\'s) Selection',
                cluster_columns = FALSE, cluster_rows = FALSE)
  pdf(file = paste('/nfs/users/nfs_a/aa16/spatialTranscriptomics/figures/allen_', corticalArea, '_', layer, '_markerExpression.pdf', sep = ""), width = 12, height = 7)
  print(ht1 + ht2)
  dev.off()
  
  # Display accuracy per cell type:
  pdf(file = paste('/nfs/users/nfs_a/aa16/spatialTranscriptomics/figures/allen_', corticalArea, '_', layer, '_accuracyByType.pdf', sep = ""), width = 12, height = 8)
  old.par = par(mar = c(12,4,2,2))
  fractionCorrectByType(fishPanel, data1[,subset], medianExpr[,subsetNames], specific_type[subset],
                        main = paste('Mapping Accuracy for ', corticalArea, " ", layer, sep = ""), order = order, axisBreak = TRUE,
                        return = FALSE, plot = TRUE)
  
  dev.off()
  par(old.par)
  
  # Display confusion matrix:
  assignedCluster <- suppressWarnings(getTopMatch(corTreeMapping(mapDat = data1[,subset], 
                                                                 medianDat=medianExpr[,subsetNames], genesToMap=fishPanel)))
  assignedCluster[assignedCluster == 'none'] = names(which.max(table(assignedCluster[,1])))
  membConfusionProp  <- getConfusionMatrix(specific_type[subset],assignedCluster[,1],TRUE)
  pdf(file = paste('/nfs/users/nfs_a/aa16/spatialTranscriptomics/figures/allen_', corticalArea, '_', layer, '_confusionMatrix.pdf', sep = ""), width = 12 , height = 8)
  heatmap.2(pmin(membConfusionProp[order,order],0.25),Rowv=NA,Colv=NA,dendrogram = "none",
            main=paste('Confusion Matrix for \n', corticalArea, layer, sep = " "),
            trace = "none", margins = c(12,12))
  dev.off()
  
  # Finally save marker genes
  write.table(fishPanel, file = paste('/nfs/users/nfs_a/aa16/spatialTranscriptomics/markerGenes/smFISH_markerGenes_', corticalArea, layer, '.txt', sep = ""),
              quote = F, col.names = F, row.names = F)
}