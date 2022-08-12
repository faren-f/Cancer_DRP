rm(list=ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

#Method5)
library(PharmacoGx)
availablePSets()
gCSI = downloadPSet("gCSI_2019")

#Method1)
library(ArrayExpress)
library(Biobase)
library(oligo)
GEOD = read.delim("Raw_data/GDSC_Raw_data/Cellline_info/GE/A-GEOD-13667_comments.txt")
adf = read.delim("Raw_data/GDSC_Raw_data/Cellline_info/GE/A-GEOD-13667.adf.txt")
idf = read.delim("Raw_data/GDSC_Raw_data/Cellline_info/GE/E-MTAB-3610.idf.txt")
SDRF = read.delim("Raw_data/GDSC_Raw_data/Cellline_info/GE/E-MTAB-3610.sdrf.txt")

rownames(SDRF) = SDRF$Array.Data.File
SDRF = AnnotatedDataFrame(SDRF)

data = read.celfiles(filenames = file.path("Raw_data/GDSC_Raw_data/Cellline_info/GE/E-MTAB-3610/", 
                                               SDRF$Array.Data.File),verbose = TRUE)

read.celfiles(filenames = file.path("Raw_data/GDSC_Raw_data/Cellline_info/GE/E-MTAB-3610/"),
              phenoData, featureData,
              experimentData, protocolData, notes, verbose=TRUE, sampleNames,
              rm.mask=FALSE, rm.outliers=FALSE, rm.extra=FALSE, checkType=TRUE)

stopifnot(validObject(data))
#saveRDS(raw_data,"Processed_data/Step12/GDSC_Raw_data.rds")
data = readRDS("Processed_data/Step12/GDSC_Raw_data.rds")
phenoData = pData(data)

f = fData(data)
expr = exprs(data)
expr = log2(expr + 1)


#------------------------------------------------------------------------------
# Background corection & Normalization & Summerization
data_rma = oligo::rma(data, target = core, normalize = TRUE)

boxplot(expr[,1:20], names=NA, cex=.1, outline = FALSE, main="raw data")

expr_rma = exprs(data_rma)
boxplot(expr_rma[,1:20], names=NA, cex=.1, outline = FALSE, main="raw data")


mean_expr = apply(expr,2,mean)
hist(mean_expr,100)
oligo::boxplot(data, target = "core",
               main = "Boxplot of log2-intensitites for the raw data")

f= fData(data)

#--------------------------------------------------------------------------------

#Method2)
library(affy)
celpath = "Raw_data/GDSC_Raw_data/Cellline_info/GE/E-MTAB-3610/"
data = ReadAffy(celfile.path=celpath)


#Method3)
library(oligo)
celpath = "Raw_data/GDSC_Raw_data/Cellline_info/GE/E-MTAB-3610/"
list = list.files(celpath,full.names=TRUE)
data = read.celfiles(list)

expr = exprs(data)
ph = data@phenoData
feat = data@featureData
exp = data@experimentData
#e = cdfName(data)
f = featureNames(data)
length(featureNames(data))
p = pData(data)

length(probeNames(data))

prob = probeNames(data)
cdfName(data)
data.qc = qc(data)



#Method4)
library(ArrayExpress)
library(Biobase)
library(oligo)
celpath = "Raw_data/GDSC_Raw_data/Cellline_info/GE/E-MTAB-3610/"
AE = getAE("E-MTAB-3610", path ="Raw_data/GDSC_Raw_data/Cellline_info/GE2/" , 
      type = "full", extract = TRUE, local = FALSE)

b = ae2bioc(AE, dataCols = NULL, drop = TRUE)

ArrayExpress("E-MTAB-3610", path = "Raw_data/GDSC_Raw_data/Cellline_info/GE/",
             save = FALSE, dataCols = NULL, drop = TRUE)




