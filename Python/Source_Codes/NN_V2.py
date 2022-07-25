#Version 2: use NN.sequential to make the expression more concise


import torch
from torch import nn
class myNetwork(nn.Module):
    def __init__(self):
        super().__init__()
        
        self.layers = nn.Sequential(
            nn.Linear(2, 16),
            nn.ReLU(),
            nn.Linear(16, 1),
            nn.Sigmoid()
        )
    
    def forward(self, x):
        x = self.layers(x)
        return x
model = myNetwork()
print(model)


