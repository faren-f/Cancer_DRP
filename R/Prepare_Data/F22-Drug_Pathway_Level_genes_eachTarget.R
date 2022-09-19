#                    Created on Wed Sep 08 15:11 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives drug name and gives all genes that exists 
# in all the pathway of its targets in each level 
Drug_Pathway_gene_set_eachTarget = function(drug, level){
  
  source("F21-Drug_Pathway_Level_genes.R")
  pathway_gene_set_all = Drug_Pathway_gene_set(drug, level)
  pathway_gene_set_all_Targets = list(all = pathway_gene_set_all)
  
  source("F20-Drug_Pathway_Level_eachTarget_Reactome.R")
  PW_tab_extend_all_targets = Drug_Pathway_Level_eachTarget(drug)
  #pathway_gene_set_all_Targets = list()
  for (i in names(PW_tab_extend_all_targets)){
    PW_tab_extend = data.frame(PW_tab_extend_all_targets[[i]])
    if((max(PW_tab_extend$level)+1) > level){
      
      PW_level_i = PW_tab_extend[PW_tab_extend$level == level,2]
      
      conv_table = readRDS("Processed_data/S7/biomart_conversion_table.rds")
      path2gene = as.list(reactomePATHID2EXTID)
      
      pw = c()
      pathway_gene_set_entrez = c()
      
      for(j in PW_level_i){
        entrez_j = path2gene[[j]]
        pathway_gene_set_entrez = c(pathway_gene_set_entrez,entrez_j)
      }
      
      pathway_gene_set_entrez = unique(pathway_gene_set_entrez)
      pathway_gene_set = conv_table[which(conv_table$entrezgene_id %in% pathway_gene_set_entrez),c(3,4)]
      pathway_gene_set_all_Targets[i] = list(pathway_gene_set[,1])
    }else{
      pathway_gene_set = c()
      pathway_gene_set_all_Targets[i] = list(pathway_gene_set[,1])
    }
  }
  
  return(pathway_gene_set_all_Targets)
}











