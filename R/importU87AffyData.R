
# open libraries
library("affy")
library("biomaRt")

# set up biomaRt data for probe/symbol conversion
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
filters = listFilters(ensembl)
martAttributes = listAttributes(ensembl)

# change to starting directory
startDir = "/Users/jeddy/Google Drive/GEO Datasets/"
setwd(startDir)
sampleList = read.delim("~/Google Drive/GEO Datasets/sampleList.txt")
sampleFiles = paste(sampleList[,1], ".CEL.gz", sep = "")

# identify compressed RAW files to untar
tarFiles = dir(pattern = "_RAW")

for (i in 1:length(tarFiles)) {
  tarFile_i = tarFiles[i]
  gseDir = sub("_RAW.tar","",tarFile_i)
  
  # untar GSE file to new folder
  untar(tarFile_i, exdir = gseDir)
  
  # read CEL files
  fileList = list.files(gseDir)
  setwd(gseDir)
  
  fileMatch = fileList[gsub("[_|-][[:alnum:]]*", "", fileList) %in% sampleFiles]
  affydata = read.affybatch(fileMatch)
  
  # use MAS5 to determine presence/absence calls
  eset_PMA = mas5calls(affydata)
  
  # convert ExpressionSignature to data frame
  X = data.frame(exprs(eset_PMA))
  
  # convert affy IDs to gene symbols
  affyids = data.frame(row.names(X))
  platform = annotation(eset_PMA)
  martPlatform = martAttributes[grep(platform, gsub("_", "", martAttributes[,1]), 
                                 value = FALSE),1]
  
  affyMapping = getBM(attributes = c(martPlatform, "entrezgene"), 
                      filters = martPlatform, values = affyids, mart = ensembl)
  
  # return to original directory
  setwd(startDir)

  # write expression data and ID conversion info to text files
  gseFile = paste(gseDir, "_", length(fileMatch), ".txt", sep = "")
  gseMap = paste(gseDir, "_map.txt", sep = "")

  write.table(X, gseFile, sep = "\t", eol = "\n", append = FALSE)
  write.table(affyMapping, gseMap, sep = "\t", eol = "\n", row.names = FALSE, append = FALSE)

  # remove new folder and untarred files
  unlink(gseDir, recursive = "TRUE")
}

