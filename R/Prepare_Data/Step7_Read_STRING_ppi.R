rm(list=ls())

library(biomaRt)
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/ppi_Raw_data/")

# Read_ppi ----------------------------------------------------------------
#1) SRTING website(Download>>Homo sapiens>>protein network data (full network, scored links between proteins))
# Data from SRTING website is new version and compeleted than data that exist in STRINGdb package

protein_links = read.delim2("9606.protein.links.v11.5.txt", header = TRUE, sep = " ")
protein_info = read.delim2("9606.protein.info.v11.5.txt", header = TRUE, sep = "\t",quote="")

#Reduce ppi links to select links with scores more than e.g. 900
#protein_links_reduced = protein_links[protein_links$score>900,]

# Seperating 9606 from first part of the ENSPs in protein_links
col1 = strsplit(protein_links$protein1, "[.]")
colname_1 = sapply(col1, '[[', 2)
col2 = strsplit(protein_links$protein2, "[.]")
colname_2 = sapply(col2, '[[', 2)
protein_links$protein1 = colname_1
protein_links$protein2 = colname_2

rm(col1)
rm(col2)
rm(colname_1)
rm(colname_2)

###############################3
#2) STRINGdb package;
# library(STRINGdb)
# library(igraph)
# 
# string_db = STRINGdb$new( version="11.5", species=9606, score_threshold=200, 
#                           input_directory="", protocol="http")
# full_graph = string_db$get_graph()
# edge_list = as_edgelist(full_graph)


# Convert ENSP to ENSG,... ------------------------------------------------

## Using biomaRt package
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
myFilter = "ensembl_peptide_id"
myValues_1 = protein_links$protein1
myAttributes = c("ensembl_peptide_id", "ensembl_gene_id","hgnc_symbol")

## assemble and query the mart (for the first column of protein_links)
res1 = getBM(attributes =  myAttributes, filters =  myFilter,
             values =  myValues_1, mart = ensembl)

## assemble and query the mart (for the second column of protein_links)
myValues_2 = protein_links$protein2 
res2 = getBM(attributes =  myAttributes, filters =  myFilter,
            values =  myValues_2, mart = ensembl)
#saveRDS(res1,"res1.rds")
#saveRDS(res2,"res2.rds")
res1 = readRDS("res1.rds")
res2 = readRDS("res2.rds")
 

v = unique(protein_links$protein1)
f = unique(protein_links$protein2)
length(intersect(protein_links$protein1,protein_links$protein1))





p = merge(protein_links, res1, by.x="protein1", by.y="ensembl_peptide_id", all=FALSE)
p = merge(protein_links, res2, by.x="protein2", by.y="ensembl_peptide_id", all=FALSE)





#######change ids
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
peptide_id <- ppi_peptide_id$node1
genes = getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id", "ensembl_gene_id","hgnc_symbol"),
              values=peptide_id , mart= mart)
colnames(genes) = c("peptide_id_1","gene_id_1","gene_name_1")
ppi = merge(ppi, genes, by.x="node1", by.y="peptide_id_1", all=FALSE)

peptide_id <- ppi_peptide_id$node2
genes = getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","ensembl_gene_id","hgnc_symbol"),
              values=peptide_id , mart= mart)

colnames(genes) = c("peptide_id_2","gene_id_2","gene_name_2")
ppi = merge(ppi, genes, by.x="node2", by.y="peptide_id_2", all=FALSE)

ppi = ppi[,c(4,6,5,7)]
colnames(ppi) = c("gene_id_1","gene_id_2","gene_name_1","gene_neme_2")
ppi = ppi[!ppi$gene_id_1%in%"",]
ppi = ppi[!ppi$gene_id_2%in%"",]

saveRDS(ppi, "ppi_gene_900.rds")


