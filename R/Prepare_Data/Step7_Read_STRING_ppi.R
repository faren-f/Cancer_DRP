
#                     Created on Shu Jul 31 2022

#                     @author: Farzaneh Firoozbakht

# Discription: 
# This script reads scored protein_links that are downloaded from STRING network 
# (https://string-db.org/cgi/download?sessionId=bI3MYlSwhzew&species_text=Homo+sapiens)
# that contain protein_ids. 
# Then we use biomart to find conversion table including protein_id, gene_id, gene symbol 
#and Entrez_id, and then we merge protein_id from STRING with protein_id of conversion table to find
# gene_id, gene_symbol and Entrez_id of protein_links that are coming from STRING.

rm(list=ls())

library(biomaRt)
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

# Read_ppi ----------------------------------------------------------------
#1) SRTING website(Download>>Homo sapiens>>protein network data (full network, scored links between proteins))
# Data from SRTING website is new version and compeleted than data that exist in STRINGdb package

protein_links = read.delim2("Raw_data/ppi_Raw_data/9606.protein.links.v11.5.txt", 
                            header = TRUE, sep = " ")
protein_info = read.delim2("Raw_data/ppi_Raw_data/9606.protein.info.v11.5.txt", 
                           header = TRUE, sep = "\t",quote="")

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

###############################
#2) STRINGdb package;
# library(STRINGdb)
# library(igraph)
# 
# string_db = STRINGdb$new( version="11.5", species=9606, score_threshold=200, 
#                           input_directory="", protocol="http")
# full_graph = string_db$get_graph()
# edge_list = as_edgelist(full_graph)


# Convert ENSP to ENSG & Gene symbol Using biomaRt package ----------------------------------------------
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
myFilter = "ensembl_peptide_id"
myValues = unique(protein_links$protein1)
myAttributes = c("ensembl_peptide_id", "ensembl_gene_id","hgnc_symbol","entrezgene_id")

## assemble and query the mart (for the first column of protein_links)
# conv_table = getBM(attributes =  myAttributes, filters =  myFilter,
#                    values =  myValues, mart = ensembl)

saveRDS(conv_table,"Processed_data/Step7/biomart_conversion_table.rds")
conv_table = readRDS("Processed_data/Step7/biomart_conversion_table.rds")

#deplicated_conv_table = which(duplicated(conv_table, incomparables=FALSE, fromLast=FALSE, by=key(conv_table)))

# Remove columns that are compeletly same in conv_table 
conv_table = unique(conv_table, incomparables=FALSE, fromLast=FALSE, by=key(conv_table))
# convert blank cells in the conv_table to "NA"
conv_table = mutate_all(conv_table, na_if,"")


# Making ppi_EdgeList_compelete -------------------------------------------------

merge1 = merge(x=protein_links,y=conv_table,by.x="protein1",
               by.y ="ensembl_peptide_id", all=FALSE)

merge_all = merge(x=merge1,y=conv_table,by.x="protein2",
               by.y ="ensembl_peptide_id", all=FALSE)

merge_all = merge_all[,c(1,2,7,4,8,5,9,6,3)]

colnames(merge_all) = c("protein_id1","protein_id2",
                        "gene_id1","gene_id2",
                        "gene_symbol1","gene_symbol2",
                        "entrez_id1","entrez_id2",
                        "ppi_score")

saveRDS(merge_all, "Processed_data/Step7/ppi_EdgeList_compelete.rds")

