### Get cell type taxonomy, then marker genes for the kidney data
require(Matrix)

coldata = read.delim('/nfs/team205/aa16/KidneyData/cellManifestCompressed.tsv')
data = readMM('/nfs/team205/aa16/KidneyData/aat1699_DataS1')
sampleNames = read.delim('/nfs/team205/aa16/KidneyData/tableOfCounts_colLabels.tsv')
geneNames = read.delim('/nfs/team205/aa16/KidneyData/tableOfCounts_rowLabels.tsv')

colnames(data) = sampleNames[,4]
rownames(data) = geneNames[,3]

# Cell type taxonomy:

