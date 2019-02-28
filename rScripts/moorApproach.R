### Moor et al. 2018 approach ###
require(MAST)
require(plotly)
require(readxl)

scData = read.table('/Users/aa16/rotation2/spatialTranscriptomics/HOLTLAB_LnNormData_Cortex.tsv')
bulk = read_excel('/Users/aa16/rotation2/spatialTranscriptomics/SupTable9_layerAstrocyteRNAseq.xlsx', sheet = 1)

## Prepare data:
celltypes = unname(unlist(scData[1,]))
scData = scData[2:dim(scData)[1],]
gNames = rownames(scData)
scData = apply(scData,2,function(x) as.numeric(x))
rownames(scData) = gNames

bulkDE_Genes = unname(unlist(bulk[,'symbol']))
foldChange = unname(unlist(bulk[,'log2FoldChange']))
foldChange = foldChange[bulkDE_Genes %in% gNames]
bulkDE_Genes = bulkDE_Genes[bulkDE_Genes %in% gNames]
topMarkers = bulkDE_Genes[foldChange > 0]
bottomMarkers = bulkDE_Genes[foldChange < 0]
bottomMarkers = bottomMarkers[length(bottomMarkers):1]
print(length(topMarkers))
print(length(bottomMarkers))
subset_size = 0.8
bottomMarkers_subset = bottomMarkers[1:round(subset_size*length(bottomMarkers))] 
topMarkers_subset = topMarkers[1:round(subset_size*length(topMarkers))] 
print(length(topMarkers_subset))
print(length(bottomMarkers_subset))

topMarkers_subset = c('Chrdl1')
bottomMarkers_subset = c('Il33')

## Now sum normalized counts of top and bottom markers in each cell:

orderBySpatialMarkers = function(scData,topMarkers,bottomMarkers, norm = NULL)
{
  bottomSum = unlist(lapply(1:dim(scData)[2], function(x) sum(scData[bottomMarkers_subset,x]))) 
  topSum = unlist(lapply(1:dim(scData)[2], function(x) sum(scData[topMarkers_subset,x])))
  # Get ratio of top and top + bottom markers
  spatialCoordinate = topSum/(bottomSum + topSum)
  # Now order cells by spatialCoordinate
  scData = scData[,order(spatialCoordinate)]
  spatialCoordinate = spatialCoordinate[order(spatialCoordinate)]
  return(scData)
}

scData = exp(scData)-1
scData = orderBySpatialMarkers(scData, topMarkers, bottomMarkers)

# Split cells into 3 equal bins:
bins = c(rep(1,495),rep(3,495)) # Since there are 3 x 330 = 990 cells

## Calculate differential expression between spatial bins:

cdr = colSums(scData != 0)
cDat = data.frame(colnames(scData))
rDat = data.frame(rownames(scData))

# Fit model
scData = log(scData + 1)
scaRaw <- FromMatrix(scData, fData = rDat, cData = cDat)
colData(scaRaw)$cngeneson <- scale(cdr)
colData(scaRaw)$bin = factor(bins)
zlmCond <- zlm(~cngeneson + bin,scaRaw)

# Do LR test for each bin:
DEgenes = list()
DEgenes_Details = list()
noBins = 2
for (x in 1:length(2:noBins)){
  i = (2:noBins)[x]
  print(i)
  summaryCond <- summary(zlmCond, doLRT=paste('bin',i, sep = '')) 
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast==paste('bin',i, sep = '') & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast==paste('bin',i, sep = '') & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdle = fcHurdle[which(fcHurdle[,'fdr'] < 0.05),]
  fcHurdle = fcHurdle[order(fcHurdle[,'fdr']),]
  DEgenes[[x]] = unname(unlist(fcHurdle[,'primerid']))
  DEgenes_Details[[x]] = fcHurdle[,c('primerid', 'Pr(>Chisq)', 'coef', 'fdr')]
}

allGenes = sort(unique(unlist(DEgenes)))

zlmResidDE <- zlm(~cngeneson, scaRaw, hook=deviance_residuals_hook)
residDE <- zlmResidDE@hookOut
residDEMatrix <- do.call(rbind, residDE)

normData1 = residDEMatrix[allGenes,]
heatmap(t(normData1), Rowv = NA, scale = 'none')

normData2 = scData[allGenes,]
heatmap(t(normData2), Rowv = NA, scale = 'row')

### Bin expression:                                                                                                                        ````````````````````````                                                                                                                                                                                                                                                                                                                                                                                                           

size = 10
breaks = seq(1,990,10)





