import decoupler as dc
import os
os.chdir("Desktop/Codes/Cancer_DRP/Python/decouple_R")

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from numpy import savetxt

# all the methods in decouple r that we can get TFs from them
dc.show_methods()

path = 'Raw_data/expresion_matrix.csv'
GE = pd.read_csv(path, sep = ',')
GE.head()
#sns.heatmap(GE, cmap='viridis')
#plt.show()

dorothea = dc.get_dorothea(organism='human', levels=['A','B','C'])

# 1) get TFs using mlm method
#tf_acts, tf_pvals = dc.run_mlm(mat=GE, net=dorothea, source='source', 
                                    #target='target', weight='weight', 
                                    #verbose=True, min_n=5)

#sns.heatmap(tf_acts, cmap='viridis')
#plt.show()

# 2) get TFs using gsea method

tf_acts_gsea = dc.run_gsea(mat=GE, net=dorothea, source='source', 
                                    target='target', min_n=5)
tf_acts_gsea_1 = tf_acts_gsea[0]
tf_acts_gsea_2 = tf_acts_gsea[1]
tf_acts_gsea_3 = tf_acts_gsea[2]



# 3) get TFs using gsva method
tf_acts_gsva = dc.run_gsva(mat=GE, net=dorothea, source='source', 
                                    target='target', min_n=5)
#sns.heatmap(tf_acts_gsva, cmap='viridis')
#plt.show()

tf_acts_mlm = dc.run_mlm(mat=GE, net=dorothea, source='source', 
                                    target='target', min_n=5)
tf_acts_mlm_1 = tf_acts_mlm[0]
tf_acts_mlm_2 = tf_acts_mlm[1]


tf_acts_ora = dc.run_ora(mat=GE, net=dorothea, source='source', 
                                    target='target', min_n=5)

tf_acts_ora_1 = tf_acts_ora[0]
tf_acts_ora_2 = tf_acts_ora[1]



tf_acts_udt = dc.run_udt(mat=GE, net=dorothea, source='source', 
                                    target='target', min_n=5)


tf_acts_wmean  = dc.run_wmean(mat=GE, net=dorothea, source='source', 
                                    target='target', min_n=5)

tf_acts_wmean_1 = tf_acts_wmean[0]
tf_acts_wmean_2 = tf_acts_wmean[1]
tf_acts_wmean_3 = tf_acts_wmean[2]
tf_acts_wmean_4 = tf_acts_wmean[3]



tf_acts_viper = dc.run_viper(mat=GE, net=dorothea, source='source', 
                                    target='target', min_n=5)

tf_acts_viper_1 = tf_acts_viper[0]
tf_acts_viper_2 = tf_acts_viper[1]

tf_acts_ulm = dc.run_ulm(mat=GE, net=dorothea, source='source', 
                                    target='target', min_n=5)

tf_acts_ulm_1 = tf_acts_ulm[0]
tf_acts_ulm_2 = tf_acts_ulm[1]

tf_acts_udt = dc.run_udt(mat=GE, net=dorothea, source='source', 
                                    target='target', min_n=5)



savetxt('Saved_data/node_attr_cellline_TF_.csv', tf_acts_, delimiter=',')







