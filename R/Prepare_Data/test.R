rm(list=ls())

library(reactome.db)
conv_table = readRDS("Processed_data/S7/biomart_conversion_table.rds")
drugs = c("doxorubicin", "temozolomide", "bicalutamide", "docetaxel", 
          "etoposide", "capecitabine", "vinorelbine", "paclitaxel", 
          "leucovorin", "cyclophosphamide", "methotrexate", "epirubicin", 
          "dacarbazine", "anastrozole", "sorafenib", "ifosfamide", "gemcitabine", 
          "tamoxifen", "vincristine", "irinotecan","cisplatin","oxaliplatin",
          "carboplatin", "vinblastine")
drug = "doxorubicin"
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
L2 = c()
L3 = c()
L4 = c()
L5 = c()
L6 = c()
L7 = c()
L8 = c()
L9 = c()
L10 = c()
ind_L1 = c()
ind_L2 = c()
ind_L3 = c()
ind_L4 = c()
ind_L5 = c()
ind_L6 = c()
ind_L7 = c()
ind_L8 = c()
ind_L9 = c()
ind_L10 = c()


while (sum(len)){
  index_max1 = which.max(len)
  ind_L1 = c(ind_L1, index_max1)
  l1 = Drug_Pathways_i[index_max1]
  L1 = c(L1, l1)
  
  ind = rbind(ind, cbind("L1", l1, index_max1, len[index_max1]))
  
  len_intersect = c()
  ind_rest = which(len != 0)
  for(k in ind_rest){
    len_intersect = c(len_intersect, length(intersect(unlist(entrez[[index_max1]]), 
                                                      unlist(entrez[[Drug_Pathways_i[[k]]]]))))
  }
  len_s1 = ifelse(len[ind_rest]-len_intersect==0,0,len[ind_rest])
  zeros1 = ind_rest[len_s1 == 0]

  len2 = len
  len2[!((1:length(len)) %in% zeros1)]=0
  len2[index_max1]=0
  
  if(sum(len2)!=0){
    index_max2 = which.max(len2)
    ind_L2 = c(ind_L2, index_max2)
    l2 = Drug_Pathways_i[index_max2]
    L2 = c(L2, l2)
    
    ind = rbind(ind, cbind("L2", l2, index_max2, len[index_max2]))
    
    len_intersect = c()
    ind_rest = which(len2 != 0)
    for(k in ind_rest){
      len_intersect = c(len_intersect, length(intersect(unlist(entrez[[index_max2]]), 
                                                        unlist(entrez[[Drug_Pathways_i[[k]]]]))))
    }
    len_s = ifelse(len[ind_rest]-len_intersect==0,0,len[ind_rest])
    zeros2 = ind_rest[len_s == 0]

    len3 = len
    len3[!((1:length(len)) %in% zeros2)]=0
    len3[index_max2]=0
    
    if(sum(len3)!=0){
      index_max3 = which.max(len3)
      ind_L3 = c(ind_L3, index_max3)
      l3 = Drug_Pathways_i[index_max3]
      L3 = c(L3, l3)
      
      ind = rbind(ind, cbind("L3", l3, index_max3, len[index_max3]))
      
      len_intersect = c()
      ind_rest = which(len3 != 0)
      for(k in ind_rest){
        len_intersect = c(len_intersect, length(intersect(unlist(entrez[[index_max3]]), 
                                                          unlist(entrez[[Drug_Pathways_i[[k]]]]))))
      }
      len_s = ifelse(len[ind_rest]-len_intersect==0,0,len[ind_rest])
      
      zeros = ind_rest[len_s == 0]
      if(length(zeros)==1){
        len[zeros]=0
        next
      }
      len_zeros = len[zeros]
      
      len4 = len
      len4[!((1:length(len)) %in% zeros)]=0
      len4[index_max3]=0
      
      if(sum(len4)!=0){    
        index_max4 = which.max(len4)
        ind_L4 = c(ind_L4, index_max4)
        l4 = Drug_Pathways_i[index_max4]
        L4 = c(L4, l4)
            
        ind = rbind(ind, cbind("L4", l4, index_max4, len[index_max4]))
        
        len_intersect = c()
        ind_rest = which(len4 != 0)
        for(k in ind_rest){
          len_intersect = c(len_intersect, length(intersect(unlist(entrez[[index_max4]]), 
                                                            unlist(entrez[[Drug_Pathways_i[[k]]]]))))
        }
        len_s = ifelse(len[ind_rest]-len_intersect==0,0,len[ind_rest])
        
        zeros = ind_rest[len_s == 0]
        len_zeros = len[zeros]
        
        len5 = len
        len5[!((1:length(len)) %in% zeros)]=0
        len5[index_max4]=0
        
        if(sum(len5)!=0){
          index_max5 = which.max(len5)
          ind_L5 = c(ind_L5, index_max5)
          l5 = Drug_Pathways_i[index_max5]
          L5 = c(L5, l5)
              
          ind = rbind(ind, cbind("L5", l5, index_max5, len[index_max5]))
          
          len_intersect = c()
          ind_rest = which(len5 != 0)
          for(k in ind_rest){
            len_intersect = c(len_intersect, length(intersect(unlist(entrez[[index_max5]]), 
                                                              unlist(entrez[[Drug_Pathways_i[[k]]]]))))
          }
          len_s = ifelse(len[ind_rest]-len_intersect==0,0,len[ind_rest])
          
          zeros = ind_rest[len_s == 0]
          len_zeros = len[zeros]
          
          len6 = len
          len6[!((1:length(len)) %in% zeros)]=0
          len6[index_max5]=0
          if(sum(len6)!=0){
            index_max6 = which.max(len6)
            ind_L6 = c(ind_L6, index_max6)
            l6 = Drug_Pathways_i[index_max6]
            L6 = c(L6, l6)
            
            ind = rbind(ind, cbind("L6", l6, index_max6, len[index_max6]))
            
            len_intersect = c()
            ind_rest = which(len6 != 0)
            for(k in ind_rest){
              len_intersect = c(len_intersect, length(intersect(unlist(entrez[[index_max6]]), 
                                                                unlist(entrez[[Drug_Pathways_i[[k]]]]))))
            }
            len_s = ifelse(len[ind_rest]-len_intersect==0,0,len[ind_rest])
            
            zeros = ind_rest[len_s == 0]
            len_zeros = len[zeros]
            
            len7 = len
            len7[!((1:length(len)) %in% zeros)]=0
            len7[index_max6]=0
            if(sum(len7)!=0){
              index_max7 = which.max(len7)
              ind_L7 = c(ind_L7, index_max7)
              l7 = Drug_Pathways_i[index_max7]
              L7 = c(L7, l7)
              
              ind = rbind(ind, cbind("L7", l7, index_max7, len[index_max7]))
              
              len_intersect = c()
              ind_rest = which(len7 != 0)
              for(k in ind_rest){
                len_intersect = c(len_intersect, length(intersect(unlist(entrez[[index_max7]]), 
                                                                  unlist(entrez[[Drug_Pathways_i[[k]]]]))))
              }
              len_s = ifelse(len[ind_rest]-len_intersect==0,0,len[ind_rest])
              
              zeros = ind_rest[len_s == 0]
              len_zeros = len[zeros]
              
              len8 = len
              len8[!((1:length(len)) %in% zeros)]=0
              len8[index_max7]=0
              if(sum(len8)!=0){
                index_max8 = which.max(len8)
                ind_L8 = c(ind_L8, index_max8)
                l8 = Drug_Pathways_i[index_max8]
                L8 = c(L8, l8)
                
                ind = rbind(ind, cbind("L8", l8, index_max8, len[index_max8]))
                
                len_intersect = c()
                ind_rest = which(len8 != 0)
                for(k in ind_rest){
                  len_intersect = c(len_intersect, length(intersect(unlist(entrez[[index_max8]]), 
                                                                    unlist(entrez[[Drug_Pathways_i[[k]]]]))))
                }
                len_s = ifelse(len[ind_rest]-len_intersect==0,0,len[ind_rest])
                
                zeros = ind_rest[len_s == 0]
                len_zeros = len[zeros]
                
                len9 = len
                len9[!((1:length(len)) %in% zeros)]=0
                len9[index_max8]=0
                if(sum(len9)!=0){
                  index_max9 = which.max(len9)
                  ind_L9 = c(ind_L9, index_max9)
                  l9 = Drug_Pathways_i[index_max9]
                  L9 = c(L9, l9)
                  
                  ind = rbind(ind, cbind("L9", l9, index_max9, len[index_max9]))
                  
                  len_intersect = c()
                  ind_rest = which(len9 != 0)
                  for(k in ind_rest){
                    len_intersect = c(len_intersect, length(intersect(unlist(entrez[[index_max9]]), 
                                                                      unlist(entrez[[Drug_Pathways_i[[k]]]]))))
                  }
                  len_s = ifelse(len[ind_rest]-len_intersect==0,0,len[ind_rest])
                  
                  zeros = ind_rest[len_s == 0]
                  len_zeros = len[zeros]
                  
                  len10 = len
                  len10[!((1:length(len)) %in% zeros)]=0
                  len10[index_max9]=0
                  if(sum(len10)!=0){
                    print("ooooooh!")
            
            
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  #if(sum(len_s1)){
    len[zeros1]=0
  #}
}




#  return(ind)
#}






