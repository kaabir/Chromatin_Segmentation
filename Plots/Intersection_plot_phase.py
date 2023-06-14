# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 13:26:30 2023

@author: kaabir
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
                   5.2617654250817256e-12,6.5596210951035404e-12]) # for 38um 7.289681258817985e-10 # Critical flow rates

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

#####################

# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 09:37:13 2023

@author: Kaabir

1. Flow rate is converted to flow velocity by dividing the cross section 
area by the intermediate layer that connects to the wells

2. The extreme ends of the plot is where Lift force is intercepting far from the simulated 
flow rate are excluded from Sediment regions

3. The graph is plotted for daimeter of well - 43um default then 55um, 70um and 91.5um
Actually >70um is difficult to represent in the same graph

"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)

# Data points for phase 1 (intersecting forces)
## 43um well
y_phase_val_43 = np.array([1.0904065715365985e-11,7.961211776256691e-12, 6.579813183257682e-12, 5.852041130083344e-12,
                   5.484782218167517e-12, 5.244229884635557e-12, 5.11731509027492e-12,
                   5.0090269779292335e-12, 4.871276850953692e-12, 4.759999786105732e-12,
                   4.636811579419479e-12, 4.553836439162763e-12, 4.490225281887017e-12,
                   4.473323101008935e-12, 4.563046835361546e-12, 4.756943888324993e-12,
                   5.2617654250817256e-12, 6.5596210951035404e-12])  # Critical flow rates

with np.printoptions(precision=18):
    y_phase_43 = np.multiply(y_phase_val_43,1e+18) # m to um conversion
    y_phase_43 = np.divide(y_phase_43,13890) # cross section area 300*46.3
    
#x_phase1 = np.linspace(min(y_phase1), max(y_phase1), len(y_phase1))  # Corresponding x-values
# 0.1*d_p,0.05*d_p,0.95*d_p
x_phase43 = np.linspace(0.05, 0.95, 18)    

## 50um well
y_phase_val_50 = np.array([2.3974129833890082e-11,1.8457016221285147e-11,1.580104307448339e-11,1.4322859702398651e-11,
                           1.3213177488803221e-11,1.2497676786841492e-11,1.1735126653383476e-11,1.1119979222314083e-11,
                           1.0343372524155945e-11,9.746111293026835e-12,9.05882523891098e-12,8.476189525565274e-12,
                           7.965848990423126e-12,7.553598136178166e-12,7.2506509547708726e-12,7.182628811161164e-12,
                           7.740686997891078e-12])  # Critical flow rates

with np.printoptions(precision=18):
    y_phase_50 = np.multiply(y_phase_val_50,1e+18) # m to um conversion
    y_phase_50 = np.divide(y_phase_50,13890) # cross section area 300*46.3

x_phase50 = np.linspace(0.1, 0.9, 17)

## 70um well
y_phase_val_70 = np.array([2.2977536136354447e-09,1.0447833546021277e-09,6.916636809467948e-10,5.01023405055832e-10,
                           3.7131232342670273e-10,2.884966941614419e-10,2.2905050944275414e-10,1.7986035734834205e-10,
                           1.457862249119165e-10,1.1987328015436576e-10,9.808994057320384e-11,8.278784147730222e-11,
                           7.155493712406736e-11,6.73136658865651e-11,2.001558664806889e-08])  # Critical flow rates

with np.printoptions(precision=18):
    y_phase_70 = np.multiply(y_phase_val_70,1e+18) # m to um conversion
    y_phase_70 = np.divide(y_phase_70,13890) # cross section area 300*46.3
    
x_phase70 = np.linspace(0.25, 0.95, 15)



# 0.0025, 0.0225

# Fill regions to distinguish Sedimentation and Lift regimes
#Filling points are from 0.1 till 0.9 i.e 4um to 36um 
## 43um well
y_fill_sediment1 = np.array([7.961211776256691e-12, 6.579813183257682e-12, 5.852041130083344e-12,
                   5.484782218167517e-12, 5.244229884635557e-12, 5.11731509027492e-12,
                   5.0090269779292335e-12, 4.871276850953692e-12, 4.759999786105732e-12,
                   4.636811579419479e-12, 4.553836439162763e-12, 4.490225281887017e-12,
                   4.473323101008935e-12, 4.563046835361546e-12, 4.756943888324993e-12,
                   5.2617654250817256e-12, 6.5596210951035404e-12])


#y_fill_sediment = np.divide(y_fill_sediment,1389)

with np.printoptions(precision=18):
    y_fill_sediment = np.multiply(y_fill_sediment1,1e+18) # m to um conversion
    y_fill_sediment = np.divide(y_fill_sediment,13890)

x_fill_sediment = np.linspace(0.1, 0.90, 17)

# Base range for the phase diagram
x_base_range = np.linspace(0, 1, 100)
y_base_range = np.linspace(min(y_phase_43), 1000, 100)

# Plotting 
fig, ax = plt.subplots(figsize=(6, 6))

# Plotting the phase diagram
# 43um
plt.plot(x_phase43, y_phase_43, color='#FF6D60', marker='o', label=r'$43 \mu m \ well $')
# 50um
plt.plot(x_phase50, y_phase_50, color='#c260ff', marker='v', label=r'$50 \mu m \ well $')

# 70um
plt.plot(x_phase70, y_phase_70, color='#ff6085', marker='*', label=r'$70 \mu m \ well $')

plt.fill_between(x_base_range, y_base_range, color='#F3E99F', alpha=0.3, label=r'$ Lift \ (F_g  \ <  \ F_{lift}) $')
plt.fill_between(x_fill_sediment, y_fill_sediment, color='#98D8AA', alpha=0.9, label=r'$ Sediment \ (F_g > F_{lift} )$')
plt.xlabel(r'$Ratio  \ of \ d/w $')
plt.ylabel(r'$Flow \ Velocity \ (\mu m/s)$' )
plt.title('Transition Phase Diagram')

# Place legends on the right side
plt.legend(loc='upper right', bbox_to_anchor=(1.45, 1))
#plt.legend()
plt.grid(True)

# Remove the thin borders near the axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

ax.spines[['right', 'top']].set_visible(False)  
#ax.grid(which='minor', color='#a8a8a8', linestyle=':')
ax.grid(which='major', color='#7a7a7a', linestyle='--')
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))

# Adjust the visible range of the plot
plt.xlim(0, max(x_base_range))
plt.ylim(250, max(y_base_range))
#plt.savefig('Transition Phase Diagram.png', dpi=500, bbox_inches="tight", pad_inches=0.2)
plt.show()

#########
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures

# Given data
y_known = np.array([7.961211776256691e-12, 6.579813183257682e-12, 5.852041130083344e-12,
                   5.484782218167517e-12, 5.244229884635557e-12, 5.11731509027492e-12,
                   5.0090269779292335e-12, 4.871276850953692e-12, 4.759999786105732e-12,
                   4.636811579419479e-12, 4.553836439162763e-12, 4.490225281887017e-12,
                   4.473323101008935e-12, 4.563046835361546e-12, 4.756943888324993e-12,
                   5.2617654250817256e-12, 6.5596210951035404e-12])  # Known y-values
x_known = np.linspace(0.1, 0.90, 17)  # Known x-values

# Reshape the data to 2D arrays
x_known = x_known.reshape(-1, 1)
y_known = y_known.reshape(-1, 1)

# Create PolynomialFeatures object
poly_features = PolynomialFeatures(degree=2)

# Transform the input features to include polynomial terms
x_poly = poly_features.fit_transform(x_known)

# Fit the linear regression model
model = LinearRegression()
model.fit(x_poly, y_known)

# Generate a sequence of x values for prediction
x_pred = np.linspace(0.0, 1.0, 100).reshape(-1, 1)

# Transform the prediction features to include polynomial terms
x_pred_poly = poly_features.transform(x_pred)

# Make predictions using the fitted model
y_pred = model.predict(x_pred_poly)

# Get the statistical values
r_sq = model.score(x_poly, y_known)
intercept = model.intercept_
coefficients = model.coef_

print("R-squared: ", r_sq)
print("Intercept: ", intercept)
print("Coefficients: ", coefficients)

# Plot the known data points and the predicted curve
import matplotlib.pyplot as plt

plt.scatter(x_known, y_known, color='red', label='Known Data')
plt.plot(x_pred, y_pred, color='blue', label='Polynomial Regression')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Polynomial Regression')
plt.legend()
plt.show()
