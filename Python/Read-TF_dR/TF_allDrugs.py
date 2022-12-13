import decoupler as dc
import os
os.chdir("Desktop/Cancer_DRP/Python/Read-TF_dR")
import pandas as pd

dorothea = dc.get_dorothea(organism='human', levels=['A','B','C'])


path_GE = '~/Desktop/Cancer_DRP/R/Prepare_Data/Processed_data/S1/expresion_matrix.csv'

GE = pd.read_csv(path_GE, sep = ',', header = 0, index_col=0)
GE.head()
GE.shape

GE_norm = (GE - GE.mean())/GE.std()

tf_acts_gsea = dc.run_gsea(mat = GE_norm, net = dorothea, source = 'source', 
                           target = 'target', min_n = 5)

tf_acts_gsea_2 = tf_acts_gsea[1]
#tf_acts_gsea_2_PR1.iloc[1,].hist(bins=3)


tf_acts_gsea_2.to_csv("~/Desktop/Cancer_DRP/R/Prepare_Data/Processed_from_Python/TF(gsea2)_PRISM/TF(gsea2)_PRISM.csv", 
                      sep=',', columns=None,header=True)
os.system('say "your program has finished"')
