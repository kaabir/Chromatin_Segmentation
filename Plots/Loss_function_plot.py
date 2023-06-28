# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 16:15:30 2023

@author: kaabi
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)

# Fixing random state for reproducibility
loss = [0.8475, 0.1729, 0.1134, 0.1155, 0.1073, 0.1043, 0.1024, 0.1005, 0.0966, 0.0933]
        #,0.0338,0.0334,0.0329,0.0329,0.0322,0.0336,0.0344,0.0344,0.0336,0.0330,0.0333]
Epochs = np.arange(0, 100, 10)           # white noise 2

# Two signals with a coherent part at 10 Hz and a random part

fig, axs = plt.subplots(1, 1)
axs.plot(Epochs,loss, color='g',linestyle='-.',linewidth=1.5, label='Train')

#axs[0].set_xlim(0, 2)
axs.set_xlabel('Actin Training Epochs')
axs.set_ylabel('Actin Loss')
#axs.grid(True)

axs.set_title('Actin Loss Funtion Over Training ') # (Non-Shrinking)
axs.spines[['right', 'top']].set_visible(False)  
axs.grid(which='minor', color='#d6d6d6', linestyle=':')
axs.grid(which='major', color='#949494', linestyle='--')
#axs.xaxis.set_minor_locator(AutoMinorLocator(5))
#axs.yaxis.set_minor_locator(AutoMinorLocator(5))
plt.legend()
plt.savefig('Loss_Function_actin_.png', dpi=400, bbox_inches="tight", pad_inches=0.1)
  
plt.show()