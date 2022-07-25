import torch
import os
from torch_geometric.data import Data
import pandas as pd
import numpy as np


class MyGraphDataset():
    def __init__(self, root):
                
        path = os.path.join(root, 'node_attr.csv')
        node_attr = pd.read_csv(path, header = None, sep = ',' )
        
        path = os.path.join(root, 'edge_index.csv')
        edge_index = pd.read_csv(path, header = None, sep= ',')
        
        #For node classification is required
        path = os.path.join(root, 'y.csv')
        y = pd.read_csv(path, header = None)
        
        #path = os.path.join(root, 'edge_attrs.csv')
        #edge_attrs = pd.read_csv(path, header=None, sep=',')
        
        node_attr = np.array(node_attr)
        self.node_attr = torch.tensor(node_attr, dtype=torch.float)
        
        edge_index = np.array(edge_index)
        self.edge_index = torch.tensor(edge_index, dtype=torch.long)
        
        y = np.array(y)
        self.y = torch.tensor(y, dtype=torch.float)
        
        #edge_attrs = np.array(edge_attrs)
        #self.edge_attrs = torch.tensor(edge_attrs, dtype=torch.float)
        
        self.N_Nodes = len(self.y)
        self.num_node_features = node_attr.shape[1]
        
        self.data = self.create_graph()
    
        
    def create_graph(self):
       
        data = Data(x = self.node_attr, edge_index = self.edge_index, 
                     y = np.squeeze(self.y))
                    
        return (data)    
    

