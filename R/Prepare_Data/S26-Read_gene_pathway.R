rm(list=ls())

library(reactome.db)

conv_table = readRDS("Processed_data/S7/biomart_conversion_table.rds")

gene2path = as.list(reactomeEXTID2PATHID)
path2gene = as.list(reactomePATHID2EXTID)

pw = gene2path[["331"]]

gene_pathways = list()
all_genes = c()
for( i in pw){
  gene_pathways[[i]] = path2gene[[i]]
  all_genes = c(all_genes, path2gene[[i]])
}

all_genes = data.frame(unique(all_genes))
all_genes = conv_table[which(conv_table$entrezgene_id %in% all_genes[,1]),c(3,4)]

saveRDS(all_genes,"Processed_data/S26/cisplatin_gene_pathways.rds")
saveRDS(gene_pathways,"Processed_data/S26/cisplatin_gene_pathways_lists.rds")


