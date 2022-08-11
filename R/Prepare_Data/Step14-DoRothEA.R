#                    Created on Wed Aug 11 11:14 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives X data to .......



# use example gene expression matrix from bcellViper package
library(bcellViper)
data(bcellViper, package = "bcellViper")
# acessing (human) dorothea regulons
# for mouse regulons: data(dorothea_mm, package = "dorothea")
data(dorothea_hs, package = "dorothea")
# run viper
tf_activities <- run_viper(dset, dorothea_hs,
                           options =  list(method = "scale", minsize = 4,
                                           eset.filter = FALSE, cores = 1,
                                           verbose = FALSE))