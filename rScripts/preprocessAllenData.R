### Preprocess Allen Data ###

require(myUtils)

## Get a combined exon+introns normalized count matrix:

cpmAllenData = function(file){
  # Save a csv file of introns + exons, counts-per-million normalized
  data1_e = read.delim(paste('/nfs/team205/aa16/AllenData/', file ,'_exon-matrix.csv', sep = ""), sep = ',')
  rownames(data1_e) = data1_e[,1]
  data1_e = data1_e[,2:dim(data1_e)[2]]
  data1_i = read.delim(paste('/nfs/team205/aa16/AllenData/', file, '_intron-matrix.csv', sep = ""), sep = ',')
  rownames(data1_i) = data1_i[,1]
  data1_i = data1_i[,2:dim(data1_i)[2]]
  data1 = data1_e+data1_i
  rm(data1_e, data1_i)
  cpm = round(cpm(data1),2)
  rm(data1)
  write.csv(cpm, paste('/nfs/team205/aa16/AllenData/', file, '_exon+intron_cpm-matrix.csv', sep = ""), quote = F)
  rm(cpm)
}

cpmAllenData("mouse_VISp_2018-06-14") # Primary Visual Cortex
cpmAllenData("mouse_ALM_2018-06-14")  # Anterior Lateral Motor Area

cpmAllenData("mouse_ACA_2018-10-04")  # Anterior Cingulate Area
cpmAllenData("mouse_MOp_cells_2018-10-04") # Primary Motor Area

## Get rpkm normalization for exon count matrix

# Get length of all genes:

gtf = read.delim('/nfs/team205/aa16/AllenData/rsem_GRCm38.p3.gtf', header = F)
gtf = gtf[gtf[,3] == 'gene',]
lengths = gtf[,5] - gtf[,4]
entrez_ids = unlist(lapply(1:dim(gtf)[1], function(x)
  strsplit(strsplit(as.character(gtf[x,9]), split = ';')[[1]][1], split = ' ')[[1]][2]))
gene_symbols = unlist(lapply(1:dim(gtf)[1], function(x)
  strsplit(strsplit(as.character(gtf[x,9]), split = ';')[[1]][2], split = ' ')[[1]][3]))
tab = cbind(entrez_ids, gene_symbols, lengths)
colnames(tab) = c('entrez_id', 'gene_symbol', 'length')
write.table(tab, file = '/nfs/team205/aa16/AllenData/gene_lengths_Allen_GRCm38.txt', 
            quote = F, row.names = F, col.names = T, sep = '\t')

rpkmAllenData = function(file){
  # Saves an rpkm file for exons
  data1_e = read.delim(paste('/nfs/team205/aa16/AllenData/', file ,'_exon-matrix.csv', sep = ""), sep = ',')
  rownames(data1_e) = data1_e[,1]
  data1_e = data1_e[,2:dim(data1_e)[2]]
  tab = read.delim('/nfs/team205/aa16/AllenData/gene_lengths_Allen_GRCm38.txt')
  lengths = tab[,3]
  names(lengths) = tab[,1]
  lengths = lengths[rownames(data1_e)]
  normData = round(rpkm(data1_e, lengths),2)
  write.csv(normData, paste('/nfs/team205/aa16/AllenData/', file, '_exon_rpkm-matrix.csv', sep = ""), quote = F)
  rm(rpkm)
}

rpkmAllenData("mouse_VISp_2018-06-14") # Primary Visual Cortex
rpkmAllenData("mouse_ALM_2018-06-14")  # Anterior Lateral Motor Area

rpkmAllenData("mouse_ACA_2018-10-04")  # Anterior Cingulate Area
rpkmAllenData("mouse_MOp_cells_2018-10-04") # Primary Motor Area

subsampleAllenData_RPKM = function(file, n = 10){
  # Downsample the RPKM data matrices so that each cell type occurs at most n times
  data = read.delim(paste('/nfs/team205/aa16/AllenData/', file, '_exon_rpkm-matrix.csv', sep = ""), sep = ',')
  rownames(data) = data[,1]
  data = data[,2:dim(data)[2]]
  coldata = read.delim(paste('/nfs/team205/aa16/AllenData/', file, '_samples-columns.csv', sep = ''), sep = ',')
  celltypes = coldata[,'cluster']
  splitted = lapply(as.character(celltypes), function(x) strsplit(x, " "))
  firstEntry = unlist(lapply(splitted, function(x) x[[1]][1]))
  remove = (celltypes == 'High Intronic VISp L5 Endou' | firstEntry == 'Low' | firstEntry == 'Batch' | firstEntry == 'Doublet') # Remove low quality cells
  data = data[,!remove]
  coldata = coldata[!remove,]
  celltypes = coldata[,'cluster']
  keep = c()
  for (i in 1:length(unique(celltypes))){
    s = which(celltypes == unique(celltypes)[i])
    if (length(s) > n){
      s = sample(s,n)
    }
    keep = c(keep,s)
  }
  keep = keep[sample(1:length(keep),length(keep))]
  data = data[,keep]
  coldata = coldata[keep,]
  write.csv(data, paste('/nfs/team205/aa16/AllenData/', file,'subSampledSize', as.character(n), '_exon_rpkm-matrix.csv', sep = ""), quote = F)
  write.csv(coldata, paste('/nfs/team205/aa16/AllenData/', file, 'subSampledSize', as.character(n), '_samples-columns.csv', sep = ''), quote = T)
}

subsampleAllenData_RPKM("mouse_VISp_2018-06-14") # Primary Visual Cortex
subsampleAllenData_RPKM("mouse_ALM_2018-06-14")  # Anterior Lateral Motor Area


