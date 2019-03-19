# Check expression of OpenTargets top hits for ALS for expression in Allen data

require(ComplexHeatmap)

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
medianExpr = do.call("cbind", tapply(names(specific_type), specific_type, function(x) rowMedians(data[,x])))
rownames(medianExpr) <- rownames(data)

save(medianExpr, file = 'medianExpr_logCPM_ALM-VISp.RData')
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
additionalGenes = c('SQSTM1')

geneTab = medianExpr[c(tab[,1], additionalGenes),]
geneTab = 2^geneTab-1
cutoff = 1
binaryTab = geneTab > cutoff
binary = rowSums(binaryTab)
print(sum(binary == 0))
geneTab = geneTab[binary > 0,]

geneTab = log(geneTab+1,2)
associationScore = tab[,2]
names(associationScore) = tab[,1]

type = gsub("s\\d+_", "", colnames(geneTab))
ha = HeatmapAnnotation(df = data.frame(type = type))

geneTab = geneTab[,res[[4]][res[[3]]]]

hm = Heatmap(geneTab, cluster_columns = FALSE, name = 'log2(CPM + 1)')
pdf(file = 'ALSgenes_medianExpression.pdf', width = 21, height = 14)
print(hm)
dev.off()

ALSgenes = rownames(geneTab)
save(ALSgenes, file = 'markerGenes/ALS_genes.RData')
geneTab = geneTab[sample(1:dim(geneTab)[1], dim(geneTab)[1]),]
write.table(geneTab, file = 'ALSgenes_medianExpression.txt', quote = FALSE, row.names = TRUE, col.names = TRUE)

# Now combine these markers with cell type markers and assess classification performance:

load('markerGenes/')



