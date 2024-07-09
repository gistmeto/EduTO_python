import shutil, os
import numpy as np

from DeepONet import DeepONet
from Post import post_losses, post_R_line, post_fields

######################### User inputs ###############################
# Layer setting
LayerDepth_Trunk = 4          # Number of layers at Trunk networks
LayerWidth_Trunk = 128        # Number of units at each Trunk hidden layers
LayerDepth_Branch = 2         # Number of layers at Branch networks
LayerWidth_Branch = 128       # Number of units at each Branch hidden layers

LayerWidth_BT = 100   # Number of output units of Trunk and branch networks
actfunc = "tanh"    

# Training setting - Optimization 
max_epoch = int(1e3)
lr = 4e-3               # Learning rate
b1 = 0.98               # b1 at Adam optimizer
e_print = 100            # Printing frequency (epoch)

# Training setting - Data
num_data1 = 3           # Number of total data for Volume
num_data2 = 5           # Number of total data for Load location

#   Layer setting
layers_Trunk = [1]
for j in range(LayerDepth_Trunk):
    layers_Trunk.append(LayerWidth_Trunk)
layers_Trunk.append(LayerWidth_BT)

layers_Branch = [1]
for j in range(LayerDepth_Branch):
    layers_Branch.append(LayerWidth_Branch)
layers_Branch.append(LayerWidth_BT)

Case_range = np.arange(15)

######################### Network Training ###############################
#   Network initialization
DeepONet = DeepONet(layers_Trunk, layers_Branch, actfunc)

## Training
DeepONet.load_data()
DeepONet.init_optimizer(lr, b1)
(loss_list) = DeepONet.train(int(max_epoch),e_print)
print('Training End')

##plot results
post_losses(loss_list)
post_fields(DeepONet,Case_range)
post_R_line(DeepONet,Case_range)

print('Postprocessing End')

del DeepONet