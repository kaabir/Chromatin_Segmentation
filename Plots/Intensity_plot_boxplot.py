# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 14:40:02 2023

@author: kaabi
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)

df_actin = pd.read_excel('E:/Quantified/eNUC Staining/230406 OF1 D1 Actin 10uMVCA Gactin RH3K9me2 FRH3K27me2/Result/Export_Actin_Excel.xlsx')
df_ctrl = pd.read_excel('E:/Quantified/eNUC Staining/230406 OF1 D1 Actin 10uMVCA Gactin RH3K9me2 FRH3K27me2/Result/Export_Ctrl_Excel.xlsx')
df_actin_ = pd.DataFrame(df_actin, columns= ['actin_intens_H3k27me3','actin_intens_H3K9me2'])
df_actin_np  = np.array(df_actin_).astype(float)

df_ctrl_ = pd.DataFrame(df_ctrl, columns= ['ctrl_intens_H3k27me3','ctrl_intens_H3K9me2'])
df_ctrl_np  = np.array(df_ctrl_).astype(float)

fig, ax = plt.subplots(figsize=(6, 6))
#plt.figure(figsize=(6,5))
# option 1, specify props dictionaries
median_color = '#c9c9c9'
c = "#FF6347" # Ctrl -Red
plt.boxplot(df_ctrl_np, positions=[1,3], patch_artist=True,
            boxprops=dict(facecolor=c),#, color=c),
            #capprops=dict(color=c),
            #whiskerprops=dict(color=c),
            #flierprops=dict(color=c, markeredgecolor=c),
            medianprops=dict(color=median_color),
            showmeans=True, meanline=True
            )

c2 = "#008080" # Actin - Teal
plt.boxplot(df_actin_np, positions=[2,4,], patch_artist=True,
            boxprops=dict(facecolor=c2),#, color=c),
            #capprops=dict(color=c),
            #whiskerprops=dict(color=c),
            #flierprops=dict(color=c, markeredgecolor=c),
            medianprops=dict(color=median_color),
            showmeans=True, meanline=True
            )
# =============================================================================
# # option 2, set all colors individually
# c2 = "blue"
# box1 = plt.boxplot(df_actin_np, positions=[1.75,3.75,5.75], patch_artist=True)
# for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
#         plt.setp(box1[item], color=c2)
# plt.setp(box1["boxes"], facecolor=c2)
# plt.setp(box1["fliers"], markeredgecolor=c2)
# =============================================================================
ax.set_title('Fixed E-Nuc. H3K- Markers (04/06/23)') # (Non-Shrinking)
categories = ['',r'$H3k9me2$',r'$H3k9me2$', r'$H3k27me3$',  r'$H3k27me3$']
ax.set_xticks(range(len(categories)), categories)
ax.set_ylabel(r'$Intensity \  (A.U) $') #(\mu m^{3})
#ax.set_xlabel('(Ctrl vs Actin)')          
ax.spines[['right', 'top']].set_visible(False)  
ax.grid(which='minor', color='#d6d6d6', linestyle=':')
ax.grid(which='major', color='#949494', linestyle='--')
ax.xaxis.set_minor_locator(AutoMinorLocator(8))
ax.yaxis.set_minor_locator(AutoMinorLocator(8))

def define_box_properties(color_code, label):
	# use plot function to draw a small line to name the legend.
	ax.plot([], c=color_code, label=label, linewidth=3)
	ax.legend(bbox_to_anchor=(1.22, 1.0))
# setting colors for each groups
define_box_properties('#FF6347', 'Ctrl') 
define_box_properties('#008080', 'Actin')

#plt.xticks([1,2,3], [1,2,3])
plt.savefig('E:/Quantified/eNUC Staining/230406 OF1 D1 Actin 10uMVCA Gactin RH3K9me2 FRH3K27me2/Result/Mean_Intensity_01.png', dpi=400, bbox_inches="tight", pad_inches=0.2)
plt.show()


