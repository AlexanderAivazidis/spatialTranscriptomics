# Check expression of OpenTargets top hits for ALS for expression in Allen data

require(ComplexHeatmap)
require(mfishtools)
require(AllenData)
require(mfishtools)
require(matrixStats)
dataDirectory = '/home/jovyan/allData/AllenData/'
savingDirectory = '/home/jovyan/spatialTranscriptomics/'

dataDirectory = '/nfs/team205/aa16/AllenData/'
savingDirectory = ''

setwd(savingDirectory)

options(stringsAsFactors = FALSE) 

allData = loadAllenData(cortical_area = c('ALM','VISp'), species = 'mouse', normalization = 'exon+intron_cpm', directory = dataDirectory)
data = as.matrix(log(allData[[1]]+1,2))
coldata = allData[[2]]
rowdata = allData[[3]]
rownames(data) = rowdata[,1]
data = data[,table(coldata[,'cluster'])[coldata[,'cluster']] > 1]
coldata = coldata[table(coldata[,'cluster'])[coldata[,'cluster']] > 1,]
specific_type = coldata[,'cluster']
names(specific_type) = coldata[,1]
exprThresh = 1
medianExpr = do.call("cbind", tapply(names(specific_type), specific_type, function(x) rowMedians(data[,x])))
propExpr   = do.call("cbind", tapply(names(specific_type), specific_type, function(x) rowMeans(data[,x]>exprThresh)))
rownames(medianExpr) <- rownames(propExpr) =  rownames(data)

save(medianExpr, file = 'medianExpr_logCPM_ALM-VISp.RData')
save(propExpr, file = 'propExpr_logCPM_ALM-VISp.RData')
load('medianExpr_logCPM_ALM-VISp.RData')
load('res.RData')

options(stringsAsFactors = FALSE)
tab = read.csv('targets_associated_with_motor_neuron_disease (2).csv')
tab = tab[,c(1,4)]
tab = tab[order(-tab[,2]),]
tab = tab[tab[,2] == 1,]

rownames(medianExpr) = toupper(rownames(medianExpr))

# Replace some genes with mouse orthologs:
leftOver = tab[,1][which(!tab[,1] %in% rownames(medianExpr))]
tab[,1][tab[,1] == 'RAB7A'] = 'RAB7'
leftOver = tab[,1][which(!tab[,1] %in% rownames(medianExpr))]
tab[,1][tab[,1] == 'SMN2'] = 'SMN1'
leftOver = tab[,1][which(!tab[,1] %in% rownames(medianExpr))]
tab[,1][tab[,1] == 'MORC2'] = 'MORC2A'
leftOver = tab[,1][which(!tab[,1] %in% rownames(medianExpr))]
tab[,1][tab[,1] == 'CFAP410'] = '1810043G02RIK'
leftOver = tab[,1][which(!tab[,1] %in% rownames(medianExpr))]

# Keep genes that have expression greater x rawcounts in at least 1 cell type:

geneTab = medianExpr[tab[,1],]
geneTab = 2^geneTab-1
cutoff = 1
binaryTab = geneTab > cutoff
binary = rowSums(binaryTab)
print(sum(binary == 0))
geneTab = geneTab[binary > 0,]

geneTab = log(geneTab+1,2)
associationScore = tab[,2]
names(associationScore) = tab[,1]
associationScore = associationScore[binary > 0]

type = gsub("s\\d+_", "", colnames(geneTab))
ha = HeatmapAnnotation(df = data.frame(type = type))

geneTab = geneTab[,res[[4]][res[[3]]]]

hm = Heatmap(geneTab, cluster_columns = FALSE, name = 'log2(CPM + 1)')
pdf(file = 'MNDgenes_medianExpression.pdf', width = 21, height = 14)
print(hm)
dev.off()

geneTab = cbind(associationScore, geneTab)
geneTab = geneTab[sample(1:dim(geneTab)[1], dim(geneTab)[1]),]
# Calculate beta score for each gene and rank by it:
load('propExpr_logCPM_ALM-VISp.RData')
rownames(propExpr) = toupper(rownames(propExpr))
beta = getBetaScore(propExpr, returnScore = TRUE)
beta = beta[rownames(geneTab)]
geneTab = cbind(beta, geneTab)
geneTab = geneTab[order(-geneTab[,1]),]
write.csv(geneTab, file = 'MNDgenes_medianExpression.csv', quote = FALSE, row.names = TRUE)
MNDgenes = rownames(geneTab)
save(MNDgenes, file = 'markerGenes/MND_genes.RData')

# Now same for ALS

options(stringsAsFactors = FALSE)
tab1 = read.csv('targets_associated_with_amyotrophic_lateral_sclerosis.csv')
tab1 = tab1[,c(1,4)]
additionalGene = tab1[which(tab1[,1] == 'SQSTM1'),]
tab1 = tab1[order(-tab1[,2]),]
tab1 = tab1[tab1[,2] > 0,]
tab1 = rbind(tab1[1:50,], additionalGene)

rownames(medianExpr) = toupper(rownames(medianExpr))
tab1[,1] = toupper(tab1[,1])

leftOver = tab1[,1][which(!tab1[,1] %in% rownames(medianExpr))]
tab1[,1][tab1[,1] == "CFAP410"] = '1810043G02RIK'

leftOver = tab1[,1][which(!tab1[,1] %in% rownames(medianExpr))]
tab1[,1][tab1[,1] == "C9ORF72"] = '3110043O21RIK'

leftOver = tab1[,1][which(!tab1[,1] %in% rownames(medianExpr))]
tab1[,1][tab1[,1] == "WASHC5"] = 'E430025E21RIK'

leftOver = tab1[,1][which(!tab1[,1] %in% rownames(medianExpr))]
tab1[,1][tab1[,1] == "KIAA0513"] = '6430548M08RIK'

leftOver = tab1[,1][which(!tab1[,1] %in% rownames(medianExpr))]

print(leftOver)

tab = tab1
tab1 = tab1[!tab1[,1] %in% leftOver,]

# Keep genes that have expression greater x rawcounts in at least 1 cell type:

geneTab = medianExpr[tab1[,1],]
geneTab = 2^geneTab-1
cutoff = 1
binaryTab = geneTab > cutoff
binary = rowSums(binaryTab)
print(sum(binary == 0))
geneTab = geneTab[binary > 0,]

geneTab = log(geneTab+1,2)
associationScore = tab1[,2]
names(associationScore) = tab1[,1]
associationScore = c(associationScore[binary > 0 ], tab[,2][tab[,1] == leftOver])
names(associationScore)[length(associationScore)] = leftOver

geneTab = geneTab[,res[[4]][res[[3]]]]

hm = Heatmap(geneTab, cluster_columns = FALSE, name = 'log2(CPM + 1)')
pdf(file = 'ALSgenes_medianExpression.pdf', width = 21, height = 14)
print(hm)
dev.off()

ALSgenes = rownames(geneTab)
save(ALSgenes, file = 'markerGenes/ALS_genes.RData')
geneTab = rbind(geneTab, rep(NA, dim(geneTab)[2]))
rownames(geneTab)[dim(geneTab)[1]] = leftOver[1]
geneTab = cbind(associationScore, geneTab)
geneTab = geneTab[order(-geneTab[,1]),]
write.csv(geneTab, file = 'ALSgenes_medianExpression.csv', quote = FALSE, row.names = TRUE)

# Now combine these markers with cell type markers and assess classification performance:

load('markerGenes/ALS_genes.RData')
load('markerGenes/MND_genes.RData')
load('markerGenes/allen_markerGenesNeurons_ALM-VISp_allCells_100genes_lossFunctionNew1.RData')
# Get ensemble names:
rowdata = read.table('../allData/AllenData/mouse_ALM_2018-06-14_genes-rows.csv', sep = ',')
entrezIds = rowdata[rowdata[,1] %in% fishPanel,4]
library('biomaRt')
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "entrezgene", attributes= c("ensembl_gene_id",
                                                          "entrezgene"),values=entrezIds,mart= mart)

fullPanel = cbind(fishPanel[match(G_list[,2], entrezIds)], G_list)
extra1 = c('Rorb', 'ENSMUSG00000036192', '107350')
extra2 = c("Nacc2", 'ENSMUSG00000026932', '21577')
fullPanel = rbind(fullPanel, extra1, extra2)
colnames(fullPanel)[1] = 'gene_symbol'
write.csv(geneTab, file = 'ALSgenes_medianExpression.csv', quote = FALSE, row.names = TRUE)

geneTab = medianExpr[toupper(fishPanel),]
geneTab = geneTab[,res[[4]][res[[3]]]]
hm = Heatmap(geneTab, cluster_columns = FALSE, name = 'log2(CPM + 1)')
pdf(file = 'celltypeMarkers_medianExpression.pdf', width = 21, height = 20)
print(hm)
dev.off()
write.csv(geneTab, file = 'celltypeMarkers_medianExpression.csv', quote = FALSE, row.names = TRUE)

fishPanel = c(fishPanel, ALSgenes, MNDgenes[!MNDgenes %in% ALSgenes])
fishPanel = fishPanel[1:150]
rownames(medianExpr) = toupper(rownames(medianExpr))
fishPanel = toupper(fishPanel)
rownames(data) = toupper(rownames(data))
accuracy = fractionCorrectWithGenes(orderedGenes = fishPanel, mapDat = data, medianDat = medianExpr, plot = FALSE, clustersF = specific_type)
print(accuracy)

plot(1:length(accuracy), accuracy/100, pch = 16,
     col = c(rep("red", 88), rep('blue', length(ALSgenes)), rep('green', length(MNDgenes))), xlab = "Number of genes in panel",
     ylab = "Accuracy", main = '', ylim =  c(0, 1), lwd = 5, cex.lab = 1.5, cex.axis = 1.5)
abline(h = c(0,0.2,0.4,0.6,0.8,1), lty = "dotted", col = 'grey', lwd = 3)
abline(h = 0, col = "black", lwd = 2.5)r
legend(, legend=c("cell type markers", "ALS genes", "MND genes"),
       col=c("red","blue", "green"),lwd = 5, cex=1.5, bty = 'n')


