### Another version of marker gene selection

library(matrixStats)   # For rowMedians function, which is fast
require(mfishtools)
library(ComplexHeatmap)
library(gplots)
require(dendextend)

cortical_area_List = c('ALM','VISp')
groups_List = c(1,2,3)
split_List = c('dendrogramSplit','randomSplit')
panelMin_List = c(2)
pallette = c('red','blue','green')

directory = '/nfs/team205/aa16/AllenData/'
#directory = '/Users/aa16/rotation2/spatialTranscriptomics/' 
saveDirectory = '/nfs/users/nfs_a/aa16/spatialTranscriptomics/'

load('res.RData') # Load a hierarchical clustering result of the data
setwd(saveDirectory)
dendrogramOrder = res[[4]][res[[3]]]
relevantCells = dendrogramOrder[3:57]
manualGroups = list()
manualGroups[[1]] = list(relevantCells)
manualGroups[[2]] = list(manualGroups[[1]][[1]][1:36], manualGroups[[1]][[1]][37:55])
manualGroups[[3]] = list(c(manualGroups[[1]][[1]][1:8], manualGroups[[1]][[1]][29:36]),
                        manualGroups[[1]][[1]][9:28], manualGroups[[1]][[1]][37:55])

omersGenes = c("Cux2", "Fam84b", "Osr1", "Colq", "Batf3", "Slc17a8", "Scnn1a",
               "Rorb", "Foxp2", "Bcl11b", "Rspo1", "Hsd11b1")

for (cortical_area in cortical_area_List){
  resList = getData(cortical_area)
  data = resList[[1]]
  coldata = resList[[2]]
  rm(resList)
  # Get glutamatergic neurons and VISp and ALM types
  coldata1 = read.delim(paste(directory, 'mouse_', cortical_area,
                             '_2018-06-14_samples-columns.csv', sep = ''), sep = ',')
  glut_subtypes = coldata1[coldata1[,'class'] == 'Glutamatergic','cluster']
  glut_subtypes = glut_subtypes[table(glut_subtypes)[glut_subtypes] > 1]
  glut_subtypes = unique(glut_subtypes)
  relevantCells = relevantCells[relevantCells %in% glut_subtypes]
  
for (panelMin in panelMin_List){
  for (splitType in split_List){
    for (groups in groups_List){
      if (splitType == 'manual'){
        focusGroups = manualGroups[[groups]]
      }else{
        focusGroups = getFocusGroups(relevantCells, groups, splitType,res) 
      }
      
      # Make a plot of the different focus groups:
      colours = rep('grey',length(dendrogramOrder))
      for (i in 1:length(focusGroups)){
        print(i)
        colours[dendrogramOrder %in% focusGroups[[i]]] = pallette[i] 
      }
      colouredDendro(as.dendrogram(res), paste(saveDirectory, 'figures/allen_',
                                               cortical_area, '_', splitType, '_groups', as.character(groups),
                                               '_colouredByGroup.pdf', sep = ""), 'Dendrogram Coloured By Focus Groups', colours)
      
      resList = runMarkerSelection(data, coldata, omersGenes, focusGroups,saveDirectory, cortical_area, splitType,
                                   panelMin = panelMin, minOn = 0.01, maxOn = 3, panelSize = 12)
      save(resList, file = paste('accuracy/allen_accuracy_', cortical_area, '_', splitType, '_groups',
                                 as.character(groups), '.RData', sep = ""))   
    }
  }
} 
}

getData = function(cortical_area){
  options(stringsAsFactors = FALSE)  # IMPORTANT
  
  data1 = read.delim(paste('/nfs/team205/aa16/AllenData/mouse_', cortical_area,
                           '_2018-06-14_exon_rpkm-matrix.csv', sep = ''),
                     sep = ',')
  rownames(data1) = data1[,1]
  data1 = data1[,2:dim(data1)[2]]
  data1 = log2(data1+1)
  
  coldata1 = read.delim(paste('/nfs/team205/aa16/AllenData/mouse_',
                              cortical_area,'_2018-06-14_samples-columns.csv', sep = ''),
                        sep = ',')
  rowdata1 = read.delim(paste('/nfs/team205/aa16/AllenData/mouse_',
                              cortical_area, '_2018-06-14_genes-rows.csv', sep = ''), sep = ',')
  rownames(data1) = rowdata1[,1]
  data1 = as.matrix(data1)
  
  # Remove celltypes that occur only once
  data1 = data1[,table(coldata1[,'cluster'])[coldata1[,'cluster']] > 1]
  coldata1 = coldata1[table(coldata1[,'cluster'])[coldata1[,'cluster']] > 1,]
  # Only keep celltypes that occur in dendrogram
  keep = (coldata1[,'cluster'] %in% res[[4]])
  coldata1 = coldata1[keep,]
  data1 = data1[,keep]
  
  return(list(data1,coldata1))
}

getFocusGroups = function(relevantCells, groups, splitType, res){
  
  focusGroups = list()
  
  # Random Split:
  N = length(relevantCells)
  n = floor(N/groups)
  if (splitType == 'randomSplit'){
    samp = sample(1:N, N)
    for (i in 1:groups){
      if (i == groups){
        focusGroups[[i]] = relevantCells[samp[(n*(i-1)+1):N]]
      }else{
        focusGroups[[i]] = relevantCells[samp[(n*(i-1)+1):(n*i)]]
      }
    }
  } else if (splitType == 'dendrogramSplit'){
    # Split along dendrogram:
    dendogramOrder = res[[4]][res[[3]]] # Get order in the dendrogram
    dendogramOrder = intersect(dendogramOrder, relevantCells)
    for (i in 1:groups){
      n = floor(N/groups)
      if (i == groups){
        focusGroups[[i]] = dendogramOrder[(n*(i-1)+1):N]
      }else{
        focusGroups[[i]] = dendogramOrder[(n*(i-1)+1):(n*i)]
      }
    }
  }
  
  return(focusGroups)
}  

runMarkerSelection = function(data, coldata, omersGenes, focusGroups, directory,
                              cortical_area, split, minOn = 1, maxOn = 250,
                              subSamp = 50, panelMin = 0, panelSize = 12, numBinaryGenes = 250, mergingFactor = 0)
{
  resList = list()
  for (i in 1:length(focusGroups)){
    print(i)
    specific_type = coldata[,'cluster']
    names(specific_type) = coldata[,'sample_name']
    clustersToMerge = unique(specific_type[!specific_type %in% focusGroups[[i]]])
    numberOfCells = sum(specific_type %in% focusGroups[[i]])
    
    exprThresh = 1
    meanExpr = do.call("cbind", tapply(names(specific_type), specific_type, function(x) rowMeans(data[,x])))
    medianExpr = do.call("cbind", tapply(names(specific_type), specific_type, function(x) rowMedians(data[,x])))
    propExpr   = do.call("cbind", tapply(names(specific_type), specific_type, function(x) rowMeans(data[,x]>exprThresh)))
    rownames(medianExpr) <- rownames(propExpr) <- genes <- rownames(data)
    
    nonZeroMedian = (!rowSums(medianExpr) == 0)
    
    runGenes <- filterPanelGenes(
      summaryExpr = 2^medianExpr[nonZeroMedian,]-1,  # medians (could also try means); We enter linear values to match the linear limits below
      propExpr    = propExpr[nonZeroMedian,],    # proportions
      startingGenes  = c(),  # Starting genes (from above)
      numBinaryGenes = numBinaryGenes,      # Number of binary genes (explained below)
      minOn     = minOn,   # Minimum required expression in highest expressing cell type
      maxOn     = maxOn,  # Maximum allowed expression
      fractionOnClusters = 0.5,  # Max fraction of on clusters (described above)
      excludeGenes    = NULL,    # Genes to exclude. Often sex chromosome or mitochondrial genes would be input here.
      excludeFamilies = c("LOC","Fam","RIK","RPS","RPL","\\-","Gm","Rnf","BC0")) # Avoid LOC markers, in this case
    
    corDist         <- function(x) return(as.dist(1-cor(x)))
    clusterDistance <- as.matrix(corDist(medianExpr))
    
    clusterDistance[clustersToMerge, clustersToMerge] = clusterDistance[clustersToMerge, clustersToMerge]*mergingFactor
    
    if (length(focusGroups) == 1){
      pS = panelSize*2
    }else{
      pS = panelSize
    }
    
    fishPanel <- buildMappingBasedMarkerPanel(
      mapDat        = data[runGenes,],         # Data for optimization
      medianDat     = medianExpr[runGenes,], # Median expression levels of relevant genes in relevant clusters
      clustersF     = specific_type,            # Vector of cluster assignments
      panelSize     = pS,                               # Final panel size
      currentPanel  = c(),                              # Starting gene panel
      subSamp       = subSamp,                               # Maximum number of cells per cluster to include in analysis (20-50 is usually best)
      panelMin      = panelMin,
      optimize      = "CorrelationDistance",            # CorrelationDistance maximizes the cluster distance as described
      clusterDistance = clusterDistance,                # Cluster distance matrix
      percentSubset = 100                               # Only consider a certain percent of genes each iteration to speed up calculations (in most cases this is not recommeded)
    )
    
    frac = fractionCorrectWithGenes(fishPanel, data, medianExpr, specific_type,
                                    return=TRUE, plot = FALSE, clustersToMerge = clustersToMerge)
    
    fracOmer = fractionCorrectWithGenes(omersGenes, data, medianExpr,specific_type,
                                        return=TRUE, plot = FALSE, clustersToMerge = clustersToMerge)
    
    # Display accuracy as function of included genes and Compare to Omer's gene list:
    pdf(file = paste(directory, 'figures/allen_', cortical_area, '_', split, '_group',
                     as.character(i),'of', as.character(groups), '_panelMin', as.character(panelMin), '_overallPerformance.pdf', sep = ""), width = 6, height = 6)
    plot(1:length(frac), frac,
         type = "l", col = "grey", xlab = "Number of genes in panel",
         ylab = "Percent of cells correctly mapped", main = paste("Mapping Accuracy for Allen data ", cortical_area,  ' ', split, ' group ',
                                                                  as.character(i),' of ', as.character(groups), sep = ""), ylim =  c(-10, 100), lwd = 5)
    lines(1:length(fracOmer), fracOmer, col = 'orange', lwd = 5)
    abline(h = (-2:20) * 5, lty = "dotted", col = 'grey')
    abline(h = 0, col = "black", lwd = 2)
    text(1:length(frac), frac+15, fishPanel, srt = 90, cex = 1)
    text(1:length(fracOmer), fracOmer-9, omersGenes, srt = 90, cex = 1, col = 'red')
    legend(1, 95, legend=c("Greedy Algorithm Selection", "Omer's Selection"),
           col=c("grey", "orange"),lwd = 5, cex=0.8)
    dev.off()
    
    
    reduceMatrix = function(matrix, clustersToMerge){
      reduced_matrix = cbind(matrix[,!colnames(matrix) %in% clustersToMerge], apply(matrix[,clustersToMerge],1, function(x) mean(x)))
      colnames(reduced_matrix)[length(colnames(reduced_matrix))] = 'Other'
      return(reduced_matrix)
    }
    
    meanExpr_reduced = reduceMatrix(meanExpr, clustersToMerge)
    
    # Display expression of marker genes:
    ht1 = Heatmap(meanExpr_reduced[fishPanel[1:12], ], name = "1",
                  column_title = '1: Algorithm Selection',
                  cluster_columns = FALSE, cluster_rows = FALSE)
    ht2 = Heatmap(meanExpr_reduced[omersGenes, ], name = "2",
                  column_title = '2: Manual (Omer\'s) Selection',
                  cluster_columns = FALSE, cluster_rows = FALSE)
    pdf(file = paste(directory, 'figures/allen_', 
                     cortical_area, '_', split, '_group',
                     as.character(i),'of', as.character(groups), '_panelMin', as.character(panelMin), '_markerExpression.pdf', sep = ""), width = 12, height = 7)
    print(ht1 + ht2)
    dev.off()
    
    names = colnames(medianExpr)
    order = c(names[names %in% clustersToMerge],
              names[!names %in% clustersToMerge])
    
    # Display accuracy per cell type:
    pdf(file = paste(directory, 'figures/allen_',cortical_area, '_', split, '_group',
                     as.character(i),'of', as.character(groups), '_panelMin', as.character(panelMin),'_accuracyByType.pdf', sep = ""), width = 12, height = 8)
    old.par = par(mar = c(12,4,2,2))
    fractionCorrectByType(fishPanel[1:12], data, medianExpr, specific_type,
                          main = paste('Mapping Accuracy for ', cortical_area,  ' ', split, ' group ',
                                       as.character(i),' of ', as.character(groups), sep = ""), axisBreak = TRUE,
                          return = FALSE, plot = TRUE, order = order)
    
    dev.off()
    par(old.par)
    
    # Display confusion matrix:
    assignedCluster <- suppressWarnings(getTopMatch(corTreeMapping(mapDat = data, 
                                                                   medianDat=medianExpr, genesToMap=fishPanel[1:12])))
    assignedCluster[assignedCluster == 'none'] = names(which.max(table(assignedCluster[,1])))
    membConfusionProp  <- getConfusionMatrix(specific_type,assignedCluster[,1],TRUE)
    membConfusionProp_reduced = reduceMatrix(membConfusionProp, clustersToMerge)
    membConfusionProp_reduced = reduceMatrix(t(membConfusionProp_reduced), clustersToMerge)
    pdf(file = paste(directory, 'figures/allen_',
                     cortical_area, '_', split, '_group',
                     as.character(i),'of', as.character(groups),'_panelMin', as.character(panelMin),
                     '_confusionMatrix.pdf', sep = ""), width = 12 , height = 8)
    
    heatmap.2(pmin(membConfusionProp,0.25),Rowv=NA,Colv=NA,dendrogram = "none",
              main=paste('Confusion Matrix for \n', cortical_area,  ' ', split, ' group ',
                         as.character(i),' of ', as.character(groups), sep = " "),
              trace = "none", margins = c(12,12))
    dev.off()
    
    # Finally save marker genes
    write.table(fishPanel, file = paste(directory, 'smFISH_markerGenes_', cortical_area, ' ', split, ' group ',
                                        as.character(i),' of ', as.character(groups),'_panelMin', as.character(panelMin), '.txt', sep = ""), quote = F, col.names = F, row.names = F)
  resList[[i]] = list(frac, numberOfCells)
  }
  return(resList)
}

colouredDendro = function(dendro, file_name, title, colours)
{
  require(dendextend)
  pdf(file = file_name, width = 21.5, height = 12)
  old_par = par(mar = c(13,5,2,0))
  dendro %>% set("leaves_pch", c(19)) %>%  # node point type
    set("leaves_cex", 1.5) %>%  # node point size
    set("labels_cex", 1) %>%  # node point size
    set("leaves_col", colours) %>% #node point color
    plot(main = title)
  par(old_par)
  dev.off()
}
    
### Make a final figures to compare performance of splitting maker genes into groups:

for (cortical_area in cortical_area_List){
  accMatrix = matrix(0,2,3)
  rownames(accMatrix) = split_List
  colnames(accMatrix) = groups_List
  for (split in split_List){
    for (groups in groups_List){
      load(paste('/Users/aa16/rotation2/spatialTranscriptomics/accuracy/allen_accuracy_',
      cortical_area, '_' , split, '_groups', groups, '.RData', sep = ''))
      accMatrix[split,groups] = round(sum(unlist(lapply(1:length(resList), function(x) resList[[x]][[1]][12]*resList[[x]][[2]])))/
        sum(unlist(lapply(1:length(resList), function(x) resList[[x]][[2]]))),1)
  }
  }
}








