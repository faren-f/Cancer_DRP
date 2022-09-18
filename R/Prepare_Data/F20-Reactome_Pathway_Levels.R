rm(list=ls())

library(reactome.db)
conv_table = readRDS("Processed_data/S7/biomart_conversion_table.rds")
drugs = c("doxorubicin", "temozolomide", "bicalutamide", "docetaxel", 
          "etoposide", "capecitabine", "vinorelbine", "paclitaxel", 
          "leucovorin", "cyclophosphamide", "methotrexate", "epirubicin", 
          "dacarbazine", "anastrozole", "sorafenib", "ifosfamide", "gemcitabine", 
          "tamoxifen", "vincristine", "irinotecan","cisplatin","oxaliplatin",
          "carboplatin", "vinblastine")
drug = "vinorelbine"
#Reactome_Pathway_Level = function(drug){
  
  Drug_Pathways = readRDS("Processed_data/S26/Drug_Pathways.rds")
  path2gene = as.list(reactomePATHID2EXTID)
  
  Drug_Pathways_i = Drug_Pathways[[drug]]
  Drug_Pathways_i = unlist(Drug_Pathways_i)
  
  entrez = list()
  len_genes = c()
  for(j in Drug_Pathways_i){
    entrez[[j]] = list(path2gene[[j]])
    len_genes = c(len_genes, length(unlist(entrez[[j]])))
  }

  ##### Level 1
  ind = c()
  len = len_genes
  L1 = c()
  ind_L1 = c()
  while (sum(len)){
    index_max = which.max(len)
    ind_L1 = c(ind_L1, index_max)
    l1 = Drug_Pathways_i[index_max]
    L1 = c(L1, l1)
    
    ind = rbind(ind, cbind("L1", l1, index_max, len[index_max]))

    len_intersect = c()
    ind_rest = which(len != 0)
      for(k in ind_rest){
        len_intersect = c(len_intersect, length(intersect(unlist(entrez[[index_max]]), 
                                                          unlist(entrez[[Drug_Pathways_i[[k]]]]))))
      }
    len_s = ifelse(len[ind_rest]-len_intersect==0,0,len[ind_rest])
    zeros = ind_rest[len_s == 0]
    len[zeros]=0
  }
  
  if(length(len_genes)>nrow(ind)){
    ##### Level 2
    len = len_genes
    len[ind_L1]=0
    L2 = c()
    ind_L2 = c()
    ind_rest = which(len != 0)
    
    while (sum(len)){
      index_max = which.max(len)
      ind_L2 = c(ind_L2, index_max)
      l2 = Drug_Pathways_i[index_max]
      L2 = c(L2, l2)
      
      ind = rbind(ind, cbind("L2", l2, index_max, len[index_max]))
      
      len_intersect = c()
      for(k in ind_rest){
        len_intersect = c(len_intersect, length(intersect(unlist(entrez[[index_max]]), 
                                                          unlist(entrez[[Drug_Pathways_i[[k]]]]))))
      }
      len_s = ifelse(len[ind_rest]-len_intersect==0,0,len[ind_rest])
      zeros = ind_rest[len_s == 0]
      len[zeros]=0
      ind_rest = which(len != 0)
    }
  }
  
  if(length(len_genes)>nrow(ind)){
    ##### Level 3
    len = len_genes
    len[ind_L1]=0
    len[ind_L2]=0
    L3 = c()
    ind_L3 = c()
    ind_rest = which(len != 0)
    
    while (sum(len)){
      index_max = which.max(len)
      ind_L3 = c(ind_L3, index_max)
      l3 = Drug_Pathways_i[index_max]
      L3 = c(L3, l3)
      
      ind = rbind(ind, cbind("L3", l3, index_max, len[index_max]))

      len_intersect = c()
      for(k in ind_rest){
        len_intersect = c(len_intersect, length(intersect(unlist(entrez[[index_max]]), 
                                                          unlist(entrez[[Drug_Pathways_i[[k]]]]))))
      }
      len_s = ifelse(len[ind_rest]-len_intersect==0,0,len[ind_rest])
      zeros = ind_rest[len_s == 0]
      len[zeros]=0
      ind_rest = which(len != 0)
    }
  }
  
  if(length(len_genes)>nrow(ind)){
    ##### Level 4
    len = len_genes
    len[ind_L1]=0
    len[ind_L2]=0
    len[ind_L3]=0
    L4 = c()
    ind_L4 = c()
    ind_rest = which(len != 0)
    
    while (sum(len)){
      index_max = which.max(len)
      ind_L4 = c(ind_L4, index_max)
      l4 = Drug_Pathways_i[index_max]
      L4 = c(L4, l4)
      
      ind = rbind(ind, cbind("L4", l4, index_max, len[index_max]))
      
      len_intersect = c()
      for(k in ind_rest){
        len_intersect = c(len_intersect, length(intersect(unlist(entrez[[index_max]]), 
                                                          unlist(entrez[[Drug_Pathways_i[[k]]]]))))
      }
      len_s = ifelse(len[ind_rest]-len_intersect==0,0,len[ind_rest])
      zeros = ind_rest[len_s == 0]
      len[zeros]=0
      ind_rest = which(len != 0)
    }
  }  
  
  if(length(len_genes)>nrow(ind)){
    ##### Level 5
    len = len_genes
    len[ind_L1]=0
    len[ind_L2]=0
    len[ind_L3]=0
    len[ind_L4]=0
    L5 = c()
    ind_L5 = c()
    ind_rest = which(len != 0)
    
    while (sum(len)){
      index_max = which.max(len)
      ind_L5 = c(ind_L5, index_max)
      l5 = Drug_Pathways_i[index_max]
      L5 = c(L5, l5)
      
      ind = rbind(ind, cbind("L5", l5, index_max, len[index_max]))
      
      len_intersect = c()
      for(k in ind_rest){
        len_intersect = c(len_intersect, length(intersect(unlist(entrez[[index_max]]), 
                                                          unlist(entrez[[Drug_Pathways_i[[k]]]]))))
      }
      len_s = ifelse(len[ind_rest]-len_intersect==0,0,len[ind_rest])
      zeros = ind_rest[len_s == 0]
      len[zeros]=0
      ind_rest = which(len != 0)
    }
  }
  
  if(length(len_genes)>nrow(ind)){
    ##### Level 6
    len = len_genes
    len[ind_L1]=0
    len[ind_L2]=0
    len[ind_L3]=0
    len[ind_L4]=0
    len[ind_L5]=0
    
    L6 = c()
    ind_L6 = c()
    ind_rest = which(len != 0)
    
    while (sum(len)){
      index_max = which.max(len)
      ind_L6 = c(ind_L6, index_max)
      
      l6 = Drug_Pathways_i[index_max]
      L6 = c(L6, l6)
      
      ind = rbind(ind, cbind("L6", l6, index_max, len[index_max]))
      
      len_intersect = c()
      for(k in ind_rest){
        len_intersect = c(len_intersect, length(intersect(unlist(entrez[[index_max]]), 
                                                          unlist(entrez[[Drug_Pathways_i[[k]]]]))))
      }
      len_s = ifelse(len[ind_rest]-len_intersect==0,0,len[ind_rest])
      zeros = ind_rest[len_s == 0]
      len[zeros]=0
      ind_rest = which(len != 0)
    }
  }
  
  if(length(len_genes)>nrow(ind)){
    ##### Level 7
    len = len_genes
    len[ind_L1]=0
    len[ind_L2]=0
    len[ind_L3]=0
    len[ind_L4]=0
    len[ind_L5]=0
    len[ind_L6]=0
    
    L7 = c()
    ind_L7 = c()
    ind_rest = which(len != 0)
    
    while (sum(len)){
      index_max = which.max(len)
      ind_L7 = c(ind_L7, index_max)
      l7 = Drug_Pathways_i[index_max]
      L7 = c(L7, l7)
      
      ind = rbind(ind, cbind("L7", l7, index_max, len[index_max]))
      
      len_intersect = c()
      for(k in ind_rest){
        len_intersect = c(len_intersect, length(intersect(unlist(entrez[[index_max]]), 
                                                          unlist(entrez[[Drug_Pathways_i[[k]]]]))))
      }
      len_s = ifelse(len[ind_rest]-len_intersect==0,0,len[ind_rest])
      zeros = ind_rest[len_s == 0]
      len[zeros]=0
      ind_rest = which(len != 0)
    }
  }
  
  if(length(len_genes)>nrow(ind)){
    ##### Level 8
    len = len_genes
    len[ind_L1]=0
    len[ind_L2]=0
    len[ind_L3]=0
    len[ind_L4]=0
    len[ind_L5]=0
    len[ind_L6]=0
    len[ind_L7]=0
    
    L8 = c()
    ind_L8 = c()
    ind_rest = which(len != 0)
    
    while (sum(len)){
      index_max = which.max(len)
      ind_L8 = c(ind_L8, index_max)
      l8 = Drug_Pathways_i[index_max]
      L8 = c(L8, l8)
      
      ind = rbind(ind, cbind("L8", l8, index_max, len[index_max]))
      
      len_intersect = c()
      for(k in ind_rest){
        len_intersect = c(len_intersect, length(intersect(unlist(entrez[[index_max]]), 
                                                          unlist(entrez[[Drug_Pathways_i[[k]]]]))))
      }
      len_s = ifelse(len[ind_rest]-len_intersect==0,0,len[ind_rest])
      zeros = ind_rest[len_s == 0]
      len[zeros]=0
      ind_rest = which(len != 0)
    }
  }
  
  
  if(length(len_genes)>nrow(ind)){
    ##### Level 9
    len = len_genes
    len[ind_L1]=0
    len[ind_L2]=0
    len[ind_L3]=0
    len[ind_L4]=0
    len[ind_L5]=0
    len[ind_L6]=0
    len[ind_L7]=0
    len[ind_L8]=0
    
    L9 = c()
    ind_L9 = c()
    ind_rest = which(len != 0)
    
    while (sum(len)){
      index_max = which.max(len)
      ind_L9 = c(ind_L9, index_max)
      l9 = Drug_Pathways_i[index_max]
      L9 = c(L9, l9)
      
      ind = rbind(ind, cbind("L9", l9, index_max, len[index_max]))
      
      len_intersect = c()
      for(k in ind_rest){
        len_intersect = c(len_intersect, length(intersect(unlist(entrez[[index_max]]), 
                                                          unlist(entrez[[Drug_Pathways_i[[k]]]]))))
      }
      len_s = ifelse(len[ind_rest]-len_intersect==0,0,len[ind_rest])
      zeros = ind_rest[len_s == 0]
      len[zeros]=0
      ind_rest = which(len != 0)
    }
  }
  
  
  if(length(len_genes)>nrow(ind)){
    ##### Level 10
    len = len_genes
    len[ind_L1]=0
    len[ind_L2]=0
    len[ind_L3]=0
    len[ind_L4]=0
    len[ind_L5]=0
    len[ind_L6]=0
    len[ind_L7]=0
    len[ind_L8]=0
    len[ind_L9]=0
    
    L10 = c()
    ind_L10 = c()
    ind_rest = which(len != 0)
    
    while (sum(len)){
      index_max = which.max(len)
      ind_L10 = c(ind_L10, index_max)
      l10 = Drug_Pathways_i[index_max]
      L10 = c(L10, l10)
      
      ind = rbind(ind, cbind("L10", l10, index_max, len[index_max]))
      
      len_intersect = c()
      for(k in ind_rest){
        len_intersect = c(len_intersect, length(intersect(unlist(entrez[[index_max]]), 
                                                          unlist(entrez[[Drug_Pathways_i[[k]]]]))))
      }
      len_s = ifelse(len[ind_rest]-len_intersect==0,0,len[ind_rest])
      zeros = ind_rest[len_s == 0]
      len[zeros]=0
      ind_rest = which(len != 0)
    }
  }
  
#  return(ind)
#}






