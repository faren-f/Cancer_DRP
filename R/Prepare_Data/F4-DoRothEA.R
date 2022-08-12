#                    Created on Wed Aug 11 11:14 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives gene expresion data to reduce the dimension of 
# genes by finding the transcription factors.


DoRothEA = function(X){
  X = t(X)
  data(dorothea_hs, package = "dorothea")
  
  TF = dorothea::run_viper(X, dorothea_hs, options =  list(method = "scale", minsize = 4,
                                            eset.filter = FALSE, cores = 1,
                                            verbose = FALSE))
  
  
  return(TF)
}




