# Check expression of OpenTargets top hits for ALS for expression in Allen data

require(AllenData)
require(mfishtools)
require(matrixStats)
dataDirectory = '/home/jovyan/allData/AllenData/'
savingDirectory = '/home/jovyan/spatialTranscriptomics/'

dataDirectory = '/nfs/team205/aa16/AllenData/'
savingDirectory = ''
#setwd(savingDirectory)

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
