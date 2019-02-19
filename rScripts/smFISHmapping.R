### Map single-cell transcriptomic data to smFISH data

require(MAST)
require(plotly)
require(xlsx)

smFISH = read.table('/Users/aa16/rotation2/spatialTranscriptomics/20190122_layerAst_P56cor.tsv')
colnames(smFISH) = smFISH[1,]
smFISH = smFISH[2:dim(smFISH)[1],]
options(stringsAsFactors = FALSE)
scData = read.table('/Users/aa16/rotation2/spatialTranscriptomics/HOLTLAB_LnNormData_Cortex.tsv')

smFISH_genes = unique(smFISH[,'genes'])

## Bin expression of smFISH genes into X bins or according to dissected Layer
noBins = 10
breaks = seq(0,1,length.out = noBins+1)
bins = .bincode(smFISH[,'normalisedDepth'], binBoundaries)

# Make a bin x gene matrix (10 x 6) with averaged expression
binMatrix = matrix(0,length(unique(smFISH[,'genes'])),noBins)
colnames(binMatrix) = as.character(1:10)
rownames(binMatrix) = unique(smFISH[,'genes'])
nMatrix = matrix(0,length(unique(smFISH[,'genes'])),noBins)
colnames(nMatrix) = as.character(1:10)
rownames(nMatrix) = unique(smFISH[,'genes'])

for (b in 1:noBins){
  for (g in unique(smFISH[,'genes'])){
    print(c(b,g))
    binMatrix[g,b] = mean(as.numeric(smFISH[smFISH[,'genes'] == g & bins == b,'spotcounts']))
      nMatrix[g,b] = sum(smFISH[,'genes'] == g & bins == b)
  }
}

p1 = plot_ly(z = t(binMatrix), x = unique(smFISH[,'genes']), y = as.character(1:10),type = "heatmap")
p1

# Calculate correlation of each sc-RNA cell to bins to assign bin

# Prepare data:
celltypes = unname(unlist(scData[1,]))
scData = scData[2:dim(scData)[1],]
gNames = rownames(scData)
scData = apply(scData,2,function(x) as.numeric(x))
rownames(scData) = gNames
cdr = colSums(scData != 0)
cDat = data.frame(colnames(scData))
rDat = data.frame(rownames(scData))

# Fit model
scaRaw <- FromMatrix(scData, fData = rDat, cData = cDat)
colData(scaRaw)$cngeneson <- scale(cdr)
colData(scaRaw)$bin = factor(bestFit)
zlmResidDE <- zlm(~cngeneson, scaRaw, hook=deviance_residuals_hook)
residDE <- zlmResidDE@hookOut
residDEMatrix <- do.call(rbind, residDE)
bestFit = rep(0,dim(scData)[2])
corValues = matrix(0, noBins, dim(scData)[2])
for (i in 1:dim(scData)[2]){
  print(i)
  corValues[,i] = unlist(lapply(1:noBins, function(x) 
    cor(binMatrix[,x],residDEMatrix[unique(smFISH[,'genes']),i])))
  bestFit[i] = which.max(corValues[,i])
}

## Find genes differentially expressed between bins with MAST:

# Fit model
zlmCond <- zlm(~cngeneson + bin,scaRaw)

# Do LR test for each bin:
DEgenes = list()
DEgenes_Details = list()
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

# Make various matrix and heatmap for DE genes: 
allGenes = sort(unique(unlist(DEgenes)))
allGenes = allGenes[order(-table(unlist(DEgenes)))]

# Visualize whether FRD < 0.05 or not:
DE_matrix = matrix(0,noBins, length(allGenes))
colnames(DE_matrix) = allGenes
rownames(DE_matrix) = as.character(1:10) 
for (i in 1:(noBins-1)){
  DE_matrix[i+1,] = 1*(unlist(lapply(allGenes, function(x) x %in% DEgenes[[i]]))) 
}
noGenes = 20
p2 = plot_ly(z = DE_matrix[,1:noGenes], x = allGenes[1:noGenes], y = as.character(1:10),type = "heatmap")
p2

# Visualize discrete model component:
Disc_matrix = matrix(0,noBins, length(allGenes))
colnames(Disc_matrix) = allGenes
rownames(Disc_matrix) = as.character(1:10) 
Disc_matrix = zlmCond@coefD[allGenes,c(1,3:11)] 
noGenes = 20
p3 = plot_ly(z = t(Disc_matrix[1:noGenes,]), x = allGenes[1:noGenes], y = as.character(1:10),type = "heatmap")
p3

# Visualize continuous model component:
Cont_matrix = matrix(0,noBins, length(allGenes))
colnames(Cont_matrix) = allGenes
rownames(Cont_matrix) = as.character(1:10) 
Cont_matrix = zlmCond@coefC[allGenes,c(1,3:11)] 
noGenes = 20
p4 = plot_ly(z = t(Cont_matrix[1:noGenes,]), x = allGenes[1:noGenes], y = as.character(1:10),type = "heatmap")
p4
zlmResidDE <- zlm(~bestFit + cdr, data.frame(t(scData)), hook=deviance_residuals_hook)
residDE <- zlmResidDE@hookOut

# Import known DE genes from previous experiment:
bulk = read_excel('/Users/aa16/rotation2/spatialTranscriptomics/SupTable9_layerAstrocyteRNAseq.xlsx', sheet = 1)
smFISH = read_excel('/Users/aa16/rotation2/spatialTranscriptomics/SupTable9_layerAstrocyteRNAseq.xlsx', sheet = 2)
bulkDE_Genes = bulk[,'symbol']
intersect(as.character(unlist(bulkDE_Genes)), allGenes)
smFISHDE_Genes = smFISH[,'symbol']
intersect(as.character(unlist(smFISHDE_Genes)), allGenes)






