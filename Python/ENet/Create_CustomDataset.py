import torch
import os
import numpy as np
import pandas as pd
from torch.utils.data import Dataset
from sklearn.preprocessing import StandardScaler

# Need to override __init__, __len__, __getitem__
# as per datasets requirement

class CustomDataset_oneDrug(Dataset):
    def __init__(self, root):
        
        path = os.path.join(root, 'expresion_matrix.csv')
        X = pd.read_csv(path, header = 0, index_col=0, sep = ',' )
        
        path = os.path.join(root, 'sensitivity_matrix.csv')
        Y = pd.read_csv(path, header= None, sep = ',')
        
        self.X_samples = X.index  
        
        X = np.array(X)
        Y = np.array(Y)
        Y = Y[:,324]
        Y = Y[np.where(~np.isnan(Y))].copy()
        X = X[np.where(~np.isnan(Y))[0],:].copy()
        Y = np.expand_dims(Y, axis = 1)
  
        
        # Normalization
    
        ss = StandardScaler()
        X = ss.fit_transform(X)
        Y = ss.fit_transform(Y)
        
        
        #Min-Max normalization    
        #Y = np.divide(np.subtract(Y,np.max(Y)), np.subtract(np.min(Y),np.max(Y)))
        
        #Y = np.divide(Y, np.absolute(np.max(Y)))          #or
        #Y = Y/np.abs(Y).max() 

        self.X = torch.tensor(X, dtype=torch.float)
        self.Y = torch.tensor(Y, dtype=torch.float)

    def __len__(self):
        return len(self.Y)

    def __getitem__(self, idx):
        X = self.X[idx]
        Y = self.Y[idx]
        data = {"Input": X, "Output": Y}
        return data
    
    
#Data = CustomDataset_oneDrug(root= "Raw_Data")
#print(Data)   
    
    
    
