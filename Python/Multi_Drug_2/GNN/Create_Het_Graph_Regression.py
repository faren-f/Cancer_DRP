import torch
import os
import pandas as pd
import numpy as np
from torch_geometric.data import HeteroData


class MyGraphDataset():
    def __init__(self, root):
                
        path = os.path.join(root, 'node_attr_cellline_TF_gsea_row1.csv')
        node_attr_cellline = pd.read_csv(path, header = None, sep = ',' )
        
        path = os.path.join(root, 'edge_index_cellline.csv')
        edge_index_cellline = pd.read_csv(path, header = None, sep= ',')
        
        path = os.path.join(root, 'node_attr_drug.csv')
        node_attr_drug = pd.read_csv(path, header = None, sep = ',' )
        
        path = os.path.join(root, 'edge_index_drug.csv')
        edge_index_drug = pd.read_csv(path, header = None, sep= ',')
        
        path = os.path.join(root, 'edge_index_cellline_drug.csv')
        edge_index_cellline_drug = pd.read_csv(path, header = None, sep= ',')
        
        path = os.path.join(root, 'edge_lable_cellline_drug.csv')
        edge_lable_cellline_drug = pd.read_csv(path, header = None, sep= ',')
        
        
        
        node_attr_cellline = np.array(node_attr_cellline)
        self.node_attr_cellline = torch.tensor(node_attr_cellline, dtype=torch.float)
        
        edge_index_cellline = np.array(edge_index_cellline)
        self.edge_index_cellline = torch.tensor(edge_index_cellline, dtype=torch.long)
        
        node_attr_drug = np.array(node_attr_drug)
        self.node_attr_drug = torch.tensor(node_attr_drug, dtype=torch.float)
        
        edge_index_drug = np.array(edge_index_drug)
        self.edge_index_drug = torch.tensor(edge_index_drug, dtype=torch.long)
        
        edge_index_cellline_drug = np.array(edge_index_cellline_drug)
        self.edge_index_cellline_drug = torch.tensor(edge_index_cellline_drug, dtype=torch.long)
        
        edge_lable_cellline_drug = np.array(edge_lable_cellline_drug)
        self.edge_lable_cellline_drug = torch.tensor(edge_lable_cellline_drug, dtype=torch.float)
        
        
        
        self.num_nodes = node_attr_cellline.shape[0]
        self.num_nodes_drug = node_attr_drug.shape[0]

        
        self.num_node_features_cellline = node_attr_cellline.shape[1]
        self.num_node_features_drug = node_attr_drug.shape[1]
        self.data = self.create_graph()
    
        
    def create_graph(self):
        
        
        data = HeteroData()

        data['cellline'].x = self.node_attr_cellline    # [num_papers, num_features_paper]
        data['drug'].x = self.node_attr_drug            # [num_authors, num_features_author]

        data['cellline', 'sen', 'drug'].edge_index =  self.edge_index_cellline_drug # [2, num_edges_cites]
        # data['cellline', 'sen_NA', 'drug'].edge_index =  self.edge_index_cellline_drug_NA # [2, num_edges_cites]

        data['cellline', 'sen', 'drug'].edge_label =  torch.squeeze(self.edge_lable_cellline_drug) # [num_edges_cites]



        data['cellline', 'sim_GE', 'cellline'].edge_index = self.edge_index_cellline # [2, num_edges_writes]
        data['drug', 'sim_FP', 'drug'].edge_index = self.edge_index_drug # [2, num_edges_affiliated]

        return (data)    
    




