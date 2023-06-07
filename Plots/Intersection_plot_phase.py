# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 13:26:30 2023

@author: kaabi
"""

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib.colors as mcolors
from decimal import Decimal

# Generate sample data

Q_lift = np.array([0, 2.5E-12, 5.0E-12, 7.5E-12, 1.0E-11])
trap_widths = np.array([7.961211776256691e-12,6.579813183257682e-12,5.852041130083344e-12,
                   5.484782218167517e-12,5.244229884635557e-12,5.11731509027492e-12,
                   5.0090269779292335e-12,4.871276850953692e-12,4.759999786105732e-12,
                   4.636811579419479e-12,4.553836439162763e-12,4.490225281887017e-12,
                   4.473323101008935e-12,4.563046835361546e-12,4.756943888324993e-12,
                   5.2617654250817256e-12,6.5596210951035404e-12])  # Critical flow rates

# Define shear regions
max_shear_region = (Q_lift >= 0.2)
min_shear_region = (Q_lift <= 0.3)
shear_free_region = (~max_shear_region & ~min_shear_region)

# Create phase diagram plot
plt.figure(figsize=(8, 6))

# Plot maximum shear region
plt.fill_betweenx(Q_lift, trap_widths[0], trap_widths[-1], where=max_shear_region, color='red', alpha=0.3,
                 label='Maximum Shear')

# Plot minimum shear region
plt.fill_betweenx(Q_lift, trap_widths[0], trap_widths[-1], where=min_shear_region, color='blue', alpha=0.3,
                 label='Minimum Shear')

# Plot shear-free region
plt.fill_betweenx(Q_lift, trap_widths[0], trap_widths[-1], where=shear_free_region, color='green', alpha=0.3,
                 label='Shear-Free')

# Plot trap width vs. Q_lift
for trap_width in trap_widths:
    plt.plot(trap_width * np.ones_like(Q_lift), Q_lift, color='black', linestyle='--', linewidth=0.5)

plt.xlabel('Trap Width')
plt.ylabel('Q_lift (Critical Flow Rate)')
plt.title('Phase Diagram of Critical Flow Rates for Tissue Ejection')
plt.legend()
plt.grid(True)
plt.show()