# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 16:24:41 2023

@author: kaabir
"""
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np

# =============================================================================
# ACTIN CHANNEL
#
# Chromocenter Volume Implementation
# =============================================================================
df_actin = pd.read_excel('C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Export_VCA_Excel.xlsx')

Chromo_volume = pd.DataFrame(df_actin, columns= ['VCA_chromo_volume_0','VCA_chromo_volume_1','VCA_chromo_volume_2'])
Chromo_volume_np = np.array(Chromo_volume)

coordinates = [
    dict(
        y=Chromo_volume_np[i,:],
        x= np.arange(0, len(Chromo_volume_np[0]),1))
        #x=len(Chromo_volume_np[0]))
    for i in range(len(Chromo_volume_np))
]

fig, ax = plt.subplots(figsize=(6, 4))
for c in coordinates:
    ax.plot(c['x'], c['y'],marker="o", alpha=0.9, linestyle='dashed',linewidth=0.65)
    #ax.plot(x, y, marker="o", color='#960056', alpha=0.8, linestyle='dashed', label='Nuclei') 
    ax.set_title('Chromocenter Volume (VCA)')
    categories = ['Actin 0','Actin 15', 'Actin 30']
    ax.set_xticks(range(len(categories)), categories)
    ax.set_ylabel('Volume ' r'$\mu m^{3}$')
    ax.set_xlabel('Actin time (mins)')

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend()
    
ax.grid()
#plt.savefig('Chromocenter_Volume1.png', dpi=300, bbox_inches="tight", pad_inches=0.0)
plt.show()      

# =============================================================================
# Chromocenter Count Implementation
# =============================================================================

Chromo_count = pd.DataFrame(df_actin, columns= ['VCA_chromo_count_0','VCA_chromo_count_1','VCA_chromo_count_2'])
Chromo_count_np = np.array(Chromo_count)

coordinates = [
    dict(
        y=Chromo_count_np[i,:],
        x= np.arange(0, len(Chromo_volume_np[0]),1))
    for i in range(len(Chromo_volume_np))
]

fig, ax = plt.subplots(figsize=(6, 4))
for c in coordinates:
    ax.plot(c['x'], c['y'],marker="o", alpha=0.9, linestyle='dashed',linewidth=0.65)
    #ax.plot(x, y, marker="o", color='#960056', alpha=0.8, linestyle='dashed', label='Nuclei') 
    ax.set_title('Chromocenter Count (VCA)')
    categories = ['Actin 0','Actin 15', 'Actin 30']
    ax.set_xticks(range(len(categories)), categories)
    ax.set_ylabel('Counts')
    ax.set_xlabel('Actin time (mins)')

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend()

ax.grid()
#plt.savefig('Chromocenter_count.png', dpi=300, bbox_inches="tight", pad_inches=0.0)

plt.show()   

# =============================================================================
# Nuclei Volume Implementation
# =============================================================================

Nuclei_volume = pd.DataFrame(df_actin, columns= ['nuclei_vca_0','nuclei_vca_1','nuclei_vca_2'])
Nuclei_volume_np = np.array(Nuclei_volume)

coordinates = [
    dict(
        y=Nuclei_volume_np[i,:],
        x= np.arange(0, len(Chromo_volume_np[0]),1))
    for i in range(len(Chromo_volume_np))
]

fig, ax = plt.subplots(figsize=(6, 4))
for c in coordinates:
    ax.plot(c['x'], c['y'],marker="o", alpha=0.9, linestyle='dashed',linewidth=0.65)
    #ax.plot(x, y, marker="o", color='#960056', alpha=0.8, linestyle='dashed', label='Nuclei') 
    ax.set_title('Nulcei Volume (VCA)')
    categories = ['Actin 0','Actin 15', 'Actin 30']
    ax.set_xticks(range(len(categories)), categories)
    ax.set_ylabel('Volume ' r'$\mu m^{3}$')
    ax.set_xlabel('Actin time (mins)')

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend()

ax.grid()
#plt.savefig('Nuclei_volume.png', dpi=300, bbox_inches="tight", pad_inches=0.0)
plt.show()    


# =============================================================================
# Control CHANNEL
#
# Chromocenter Volume Implementation
# =============================================================================
df_Ctrl = pd.read_excel('C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Export_Ctrl_Excel.xlsx',index_col=None)

Chromo_volume = pd.DataFrame(df_Ctrl, columns= ['Ctrl_chromo_volume_0','Ctrl_chromo_volume_1','Ctrl_chromo_volume_2'])
Chromo_volume_np = np.array(Chromo_volume)

coordinates = [
    dict(
        y=Chromo_volume_np[i,:],
        x= np.arange(0, len(Chromo_volume_np[0]),1))
    for i in range(len(Chromo_volume_np))
]

fig, ax = plt.subplots(figsize=(6, 4))
for c in coordinates:
    ax.plot(c['x'], c['y'],marker="o", alpha=0.9, linestyle='dashed',linewidth=0.65)
    #ax.plot(x, y, marker="o", color='#960056', alpha=0.8, linestyle='dashed', label='Nuclei') 
    ax.set_title('Chromocenter Volume (Ctrl)')
    categories = ['Actin 0','Actin 15', 'Actin 30']
    ax.set_xticks(range(len(categories)), categories)
    ax.set_ylabel('Volume ' r'$\mu m^{3}$')
    ax.set_xlabel('Actin time (mins)')

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend()

ax.grid()
#plt.savefig('Chromo_volume_ctrl.png', dpi=300, bbox_inches="tight", pad_inches=0.0)
plt.show()    

# =============================================================================
# Chromocenter Count Implementation
# =============================================================================

Chromo_count = pd.DataFrame(df_Ctrl, columns= ['Ctrl_chromo_count_0','Ctrl_chromo_count_1','Ctrl_chromo_count_2'])
Chromo_count_np = np.array(Chromo_count)

coordinates = [
    dict(
        y=Chromo_count_np[i,:],
        x= np.arange(0, len(Chromo_volume_np[0]),1))
    for i in range(len(Chromo_volume_np))
]

fig, ax = plt.subplots(figsize=(6, 4))
for c in coordinates:
    ax.plot(c['x'], c['y'],marker="o", alpha=0.9, linestyle='dashed',linewidth=0.65)
    #ax.plot(x, y, marker="o", color='#960056', alpha=0.8, linestyle='dashed', label='Nuclei') 
    ax.set_title('Chromocenter Count (Ctrl)')
    categories = ['Actin 0','Actin 15', 'Actin 30']
    ax.set_xticks(range(len(categories)), categories)
    ax.set_ylabel('$Counts$')
    ax.set_xlabel('Actin time (mins)')

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend()

ax.grid()
#plt.savefig('Chromo_count_ctrl.png', dpi=300, bbox_inches="tight", pad_inches=0.0)
plt.show()   

# =============================================================================
# Nuclei Volume Implementation
# =============================================================================

Nuclei_volume = pd.DataFrame(df_Ctrl, columns= ['nuclei_ctrl_0','nuclei_ctrl_1','nuclei_ctrl_2'])
Nuclei_volume_np = np.array(Nuclei_volume)

coordinates = [
    dict(
        y=Nuclei_volume_np[i,:],
        x= np.arange(0, len(Chromo_volume_np[0]),1))
    for i in range(len(Chromo_volume_np))
]

fig, ax = plt.subplots(figsize=(6, 4))
for c in coordinates:
    ax.plot(c['x'], c['y'],marker="o", alpha=0.9, linestyle='dashed',linewidth=0.65)
    #ax.plot(x, y, marker="o", color='#960056', alpha=0.8, linestyle='dashed', label='Nuclei') 
    ax.set_title('Nulcei Volume (Ctrl)')
    categories = ['Actin 0','Actin 15', 'Actin 30']
    ax.set_xticks(range(len(categories)), categories)
    ax.set_ylabel('$Volume$ ' r'$\mu m^{3}$')
    ax.set_xlabel('Actin time (mins)')

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend()

ax.grid()
#plt.savefig('Nuclei_volume_ctrl.png', dpi=300, bbox_inches="tight", pad_inches=0.0)

plt.show() 

# =============================================================================
# # =============================================================================
# # Actin Coverage
# # =============================================================================
# y = [24.3,31.0]
# y1 = [31.1,26.5]
# x = range(len(y))
# 
# fig, ax = plt.subplots(figsize=(6, 4))
# 
# ax.plot(x, y, marker="o", color='#960056', alpha=0.8, linestyle='dashed', label='Nuclei 1')
# ax.plot(x,y1, marker="s",color='#0b5509',alpha=0.8, linestyle='dashed', label='Nuclei 2')
# #ax.margins(0.15, 0.15)   
# ax.set_title('Actin Coverage')
# for index in range(len(x)):
#   ax.text(x[index], y[index], y[index], size=12)
#   ax.text(x[index], y1[index], y1[index], size=12)
# 
# categories = ['Actin 30', 'Actin 60']
# ax.set_xticks(range(len(y)), categories)
# ax.set_ylabel('Coverage (%)')
# ax.set_xlabel('Actin time (mins)')
# 
# A = [0,0,24,31]
# A2 = [0,0,31,26]
# for i in range(len(x)):
#     if i > 1:
#         plt.annotate(A[i], (x[i], y[i]), textcoords='offset pixels', size=12)
#         plt.annotate(A2[i], (x[i], y1[i]), textcoords='offset pixels', size=12)
# 
# plt.gca().spines['top'].set_visible(False)
# plt.gca().spines['right'].set_visible(False)
# plt.legend()
# ax.grid()
# #plt.savefig('Actin_Coverage.png', dpi=300, bbox_inches="tight", pad_inches=0.0)
# 
# plt.show()
# =============================================================================
