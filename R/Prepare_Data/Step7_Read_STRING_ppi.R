rm(list=(ls))

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/ppi_Raw_data/")

#Read Full graph
#1) SRTING website(Download>>Homo sapiens>>protein network data (full network, scored links between proteins))
# Data from SRTING website is new version and compeleted than data that exist in STRINGdb package

protein_links = read.delim2("9606.protein.links.v11.5.txt", header = TRUE, sep = " ")
protein_info = read.delim2("9606.protein.info.v11.5.txt", header = TRUE, sep = "\t",quote="")

#2) STRINGdb package;
# library(STRINGdb)
# library(igraph)
# 
# string_db = STRINGdb$new( version="11.5", species=9606, score_threshold=200, 
#                           input_directory="", protocol="http")
# full_graph = string_db$get_graph()
# edge_list = as_edgelist(full_graph)










#######change ids
# ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 
# head(listFilters(ensembl), 3)             ## filters
# myFilter <- "chromosome_name"
# substr(filterOptions(myFilter, ensembl), 1, 50) ## return values
# myValues <- c("21", "22")
# head(listAttributes(ensembl), 3)          ## attributes
# myAttributes <- c("ensembl_gene_id","chromosome_name")
# 
# ## assemble and query the mart
# res <- getBM(attributes =  myAttributes, filters =  myFilter,
#              values =  myValues, mart = ensembl)




