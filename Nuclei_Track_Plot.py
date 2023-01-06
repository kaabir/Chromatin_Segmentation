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
# Actin Coverage added Mapped with Nuclei Volume ratio to chromocenter volume
# =============================================================================
df_actin = pd.read_excel('C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Export_VCA_Excel.xlsx')

Chromo_volume = pd.DataFrame(df_actin, columns= ['VCA_chromo_volume_0','VCA_chromo_volume_1','VCA_chromo_volume_2'])
Chromo_volume_np = np.array(Chromo_volume)

Nuclei_volume_c = pd.DataFrame(df_actin, columns= ['nuclei_vca_2','nuclei_vca_1','nuclei_vca_0'])
Nuclei_volume_c_np = np.array(Nuclei_volume_c)

# Ratio of Nuclei Volume Column 2/1 (T60/T30) et 2/0 (T60/T0)
ratio_2 = np.divide(Nuclei_volume_c_np[:,2],Nuclei_volume_c_np[:,0])
#ratio_2 = np.insert(ratio_2, 0, 0., axis=0) # Insert zero initally to skip annotation on graph
ratio_2 =np.round(ratio_2, 2)
ratio_2 =list(map(str, ratio_2))

ratio_1 = np.divide(Nuclei_volume_c_np[:,1],Nuclei_volume_c_np[:,0])
#ratio_1 = np.insert(ratio_1, 0, 0., axis=0)
ratio_1 =np.round(ratio_1, 2)
ratio_1 =list(map(str, ratio_1))

Coverage_per = pd.DataFrame(df_actin, columns= ['actin_area_60'])
Coverage_per_np = np.array(Coverage_per)

FileNum = pd.DataFrame(df_actin, columns= ['File_Name_VCA'])
FileNum_np = np.array(FileNum)

# get annotation coordinates
ac_2_Y=Chromo_volume_np[:,2] 
ac_1_Y=Chromo_volume_np[:,1]
ac_1_Y =list(map(float, ac_1_Y))
ac_2_X = range(len(ac_2_Y))
ac_1_X = range(len(ac_1_Y))
ac_1_X = list(map(float, ac_1_X))

coordinates = [
    dict(
        y=Chromo_volume_np[i,:], # Comparing rows of timestamps
        x= np.arange(0, len(Chromo_volume_np[0]),1),legend =FileNum_np[i],
        act_Cov=Coverage_per_np[i] )
        #x=len(Chromo_volume_np[0])
        
    for i in range(len(Chromo_volume_np))
]

fig, ax = plt.subplots(figsize=(6, 6))
ax2 = ax.twinx()
for c in coordinates:
    ax.plot(c['x'], c['y'],marker="o", alpha=0.9, linestyle='dashed',linewidth=0.65, label=c['legend'])

    ax.set_title('Chromocenter Volume (VCA)')
    categories = ['Actin 0','Actin 15', 'Actin 30']
    ax.set_xticks(range(len(categories)), categories)
    ax.set_ylabel('Volume ' r'($\mu m^{3}$)')
    ax.set_xlabel('Actin time (mins)')
    #ax.annotate(c['A2'], xy=(0.5, 0.5), xycoords=ax.transAxes)
    
    
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend()
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5),title="File:",fancybox=True,shadow=True )
    
    patches = [ plt.plot([],[], marker='o', ms=2, ls="",c='k', mec=None, 
            label="{:.2f}".format(float(Coverage_per_np[i])) )[0]  for i in range(len(Coverage_per_np)) ]    

    ax2.legend(handles=patches,loc='center left', bbox_to_anchor=(1.25, 0.5),title="Actin (%)",fancybox=True,shadow=True)
    ax2.axes.get_yaxis().set_visible(False)
    
for i in range(len(ac_2_Y)):
    #print(ratio_1[i])
    ax.annotate(ratio_1[i], xy=(1, ac_1_Y[i]), textcoords='offset pixels')
    ax.annotate(ratio_2[i], xy=(2, ac_2_Y[i]), textcoords='offset pixels')
    #ax.annotate(ratio_1[i], xy=(ac_1_X[i], ac_1_Y[i]), xycoords='data')
    
ax.annotate("Nuc. V. T60:T0",
            xy=(2, 2.8), xycoords='data',
            xytext=(2, 3), textcoords='data',
            arrowprops=dict(arrowstyle="-",
                            connectionstyle="bar,angle=180,fraction=-0.2"))
ax.annotate("Nuc. V. T30:T0",
            xy=(1, 2.6), xycoords='data',
            xytext=(1, 2.8), textcoords='data',
            arrowprops=dict(arrowstyle="-",
                            connectionstyle="bar,angle=180,fraction=-0.2"))
t = ax.text(
    2.65, 2.40, "Experiment Number", ha="center", va="center", size=11,
    bbox=dict(boxstyle="rarrow,pad=0.3", fc="white", ec="k", lw=1))
bb = t.get_bbox_patch()
bb.set_boxstyle("round", pad=0.6)

ax.set_ylim([0, 3])
ax.spines[['right', 'top']].set_visible(False)  
ax.grid()
plt.savefig('Chromocenter_Volume_ACTIN.png', dpi=400, bbox_inches="tight", pad_inches=0.2)
plt.show()    