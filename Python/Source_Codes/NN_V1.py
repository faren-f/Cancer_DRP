#Version 1: use the NN module


from torch import nn
class myNetwork(nn.Module):
    def __init__(self):
        super().__init__()
        
        # Hidden layer with linear transformation
        self.fc1 = nn.Linear(2, 16)
        self.fc2 = nn.Linear(16, 1)
    
    def forward(self, x):
        # Pass the input through all the layers
        x = self.fc1(x)
        x = self.relu(x)
        x = self.fc2(x)
        x = self.softmax(x)
        return x
model = myNetwork()
print(model)



