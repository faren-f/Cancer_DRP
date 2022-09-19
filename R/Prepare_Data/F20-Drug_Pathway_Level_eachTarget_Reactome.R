#                    Created on Wed Sep 18 10:11 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives drug name and gives the level of all its 
# drug target pathways 

library(reactome.db)
Drug_Pathway_Level_eachTarget = function(drug){

  setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
  Drug_Pathways = readRDS("Processed_data/S26/Drug_Pathways.rds")
  path2gene = as.list(reactomePATHID2EXTID)
  PW_tab_extend_all_targets = list()
  for(target in names(Drug_Pathways[[drug]])){
  
    Drug_Pathways_i = Drug_Pathways[[drug]][target]
    Drug_Pathways_i = unlist(Drug_Pathways_i[!duplicated(Drug_Pathways_i)])
  
    if(length(Drug_Pathways_i)){
      entrez = list()
      len = c()
      for(j in Drug_Pathways_i){
        entrez[[j]] = list(path2gene[[j]])
        len = c(len, length(unlist(entrez[[j]])))
      }
      
      entrez_extra = c()
      for(i in 1:length(entrez)){
        for(j in i:length(entrez)){
          I = intersect(unlist(entrez[[i]]),unlist(entrez[[j]]))
          if(2*length(I) == length(unlist(entrez[[i]]))+length(unlist(entrez[[j]])) & i!=j){
            entrez_extra = c(entrez_extra, names(entrez)[j])
          }
        }
      }
      
      length(entrez_extra)
      length(unique(entrez_extra))
      entrez_extra = unique(entrez_extra)
      l = !(names(entrez) %in% entrez_extra)
      entrez = entrez[names(entrez)[l]]
      len = len[l]
      Drug_Pathways_i = Drug_Pathways_i[!(Drug_Pathways_i%in% entrez_extra)]
    
    
    
      ##Function
      inter_pathway = function(pw_max, len, entrez){
        
        len_intersect = c()
        for(k in names(len)){
          len_intersect = c(len_intersect,length(intersect(unlist(entrez[[pw_max]]),unlist(entrez[[k]]))))
        }
        len_s1 = ifelse(len-len_intersect==0,0,len)
        inner_pw = len[len_s1 == 0]
        inner_pw = inner_pw[names(inner_pw) != pw_max]
        return(inner_pw)
      }
      #######
      names(len) = Drug_Pathways_i
    
      PW = c()
      while (sum(len)){
        print(length(len))
        pw_max = names(which.max(len))
        PW = rbind(PW, cbind(1, pw_max, len[pw_max]))
        inner_pw = inter_pathway(pw_max, len, entrez)
        
        i=2
        while(sum(inner_pw)){
          pw_max = names(which.max(inner_pw))
          PW = rbind(PW, cbind(i, pw_max, len[pw_max]))
          inner_pw = inter_pathway(pw_max, inner_pw, entrez)
          i = i+1
        }
        len = len[names(len) != PW[nrow(PW),2]]
      }
      
      colnames(PW) = c("level","pw_max", "No_genes")
      PW = data.frame(PW)
    
      PW_tab = PW[,c(1,2)]
      PW_tab = PW_tab[!duplicated(PW$pw_max),]
      rownames(PW_tab) = PW_tab$pw_max
      
      t = table(PW$pw_max)
      PW_tab$num = rep(0,nrow(PW_tab))
      PW_tab[names(t),3] = t
      
      PW_tab$level = as.numeric(PW_tab$level)
      max_level = max(PW_tab$level)
      
      PW_tab_extend = PW_tab
      #PW_tab_extend = PW_tab[,-3]
      
      for (i in 1:nrow(PW_tab)){
        
        if ((PW_tab$num[i]==1) & (PW_tab$level[i]!=max_level)){
          level = (PW_tab$level[i]+1) : max_level
          t = data.frame(level = level, pw_max = rep(PW_tab$pw_max[i], length(level)),
                         num = rep(PW_tab$num[i], length(level)))
          PW_tab_extend = rbind(PW_tab_extend, t)
        }
      }
    }else{
      PW_tab_extend = c()
    }
    PW_tab_extend_all_targets[[target]] = list(PW_tab_extend)
    
  }
  return(PW_tab_extend_all_targets)
}

