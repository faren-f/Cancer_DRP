#                    Created on Wed Aug 19 14:27 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives gene expresion data to reduce the dimension of 
# genes by finding the transcription factors.


Drug_Targets = function(X){
  
  Drug_Targets = readRDS("Processed_data/S1/drug_targets.rds")
  rownames(Drug_Targets) = Drug_Targets$name  
  
  DTs = list()
  for(i in rownames(Drug_Targets)){  
    Drug_Targets_i = strsplit(Drug_Targets[i,2],", ")
    DT = intersect(colnames(X),Drug_Targets_i[[1]])
    
    if (length(DT)<1){
    DTs[[i]]= NA
    next
    }
    DTs[[i]]= DT
  }
  
  return(DTs)
}
  
  
  