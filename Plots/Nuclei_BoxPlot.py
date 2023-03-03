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
# To numpy array
Chromo_volume_np_actin = np.array(Chromo_volume)

fig, ax = plt.subplots(figsize=(6, 4))
categories = ['Actin 0', 'Actin 15', 'Actin 30']
# rectangular box plot
bplot = ax.boxplot(Chromo_volume_np_actin,
                     vert=True,  # vertical box alignment
                     patch_artist=True,  # fill with color
                     labels=categories,showmeans=True, meanline=True)

for median in bplot['medians']:
    median.set_color('black')

# fill with colors
colors = ['#fa878b', '#98fa9a', '#98f3fa','#e098fa']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)

# Labels 
ax.set_title('Chromocenter Volume (Actin)')
ax.set_ylabel('Volume ' r'$\mu m^{3}$')
ax.set_xlabel('Actin time (mins)')

# To get Stats but boxplot dicitionary does the same
# Dictionary in case of normal plot
coordinates = [
    dict(
        y=Chromo_volume_np_actin[:,i],
        x=np.arange(0, len(Chromo_volume_np_actin[0]),1))
    for i in range(len(Chromo_volume_np_actin[0]))
]

for c in coordinates:
    mean = np.mean(c['y'])
    Mean = np.round(mean)
    var = np.var(c['y'])
    std = np.std(c['y'])
    
# Graph Fixes
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.legend()
ax.grid()
#plt.savefig('Chromocenter_Volume_actin.png', dpi=300, bbox_inches="tight", pad_inches=0.0)
plt.show()    


# =============================================================================
# # Array of 8x4 length and prints ROWS
# print(len(Chromo_volume_np))
# # Array of 8x4 length and prints COLUMNS
# print(len(Chromo_volume_np[1]))
# =============================================================================
# =============================================================================
# Chromocenter Count Implementation
# =============================================================================

Chromo_count = pd.DataFrame(df_actin, columns= ['VCA_chromo_count_0','VCA_chromo_count_1','VCA_chromo_count_2'])
Chromo_count_np_actin = np.array(Chromo_count)

fig, ax = plt.subplots(figsize=(6, 4))
categories = ['Actin 0', 'Actin 15', 'Actin 30']
# rectangular box plot
bplot = ax.boxplot(Chromo_count_np_actin,
                     vert=True,  # vertical box alignment
                     patch_artist=True,  # fill with color
                     labels=categories,showmeans=True, meanline=True)

for median in bplot['medians']:
    median.set_color('black')
    
# fill with colors
colors = ['#fa878b', '#98fa9a', '#98f3fa','#e098fa']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
    
# Labels 
ax.set_title('Chromocenter Count (Actin)')
ax.set_ylabel('Counts')
ax.set_xlabel('Actin time (mins)')
# Graph Fixes
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.legend()
ax.grid()
#plt.savefig('Chromocenter_count_actin.png', dpi=300, bbox_inches="tight", pad_inches=0.0)

plt.show()   

# =============================================================================
# Nuclei Volume Implementation
# =============================================================================

Nuclei_volume = pd.DataFrame(df_actin, columns= ['nuclei_vca_0','nuclei_vca_1','nuclei_vca_2'])
Nuclei_volume_np_actin = np.array(Nuclei_volume)

fig, ax = plt.subplots(figsize=(6, 4))
categories = ['Actin 0', 'Actin 15', 'Actin 30']
# rectangular box plot
bplot = ax.boxplot(Nuclei_volume_np_actin,
                     vert=True,  # vertical box alignment
                     patch_artist=True,  # fill with color
                     labels=categories,showmeans=True, meanline=True)

for median in bplot['medians']:
    median.set_color('black')
    
# fill with colors
colors = ['#fa878b', '#98fa9a', '#98f3fa','#e098fa']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
    
# Labels 
ax.set_title('Nuclei Volume (Actin)')
ax.set_ylabel('Volume ' r'$\mu m^{3}$')
ax.set_xlabel('Actin time (mins)')
# Graph Fixes
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.legend()
ax.grid()
#plt.savefig('Nuclei_volume_actin.png', dpi=300, bbox_inches="tight", pad_inches=0.0)
plt.show()    


# =============================================================================
# Control CHANNEL
#
# Chromocenter Volume Implementation
# =============================================================================
df_Ctrl = pd.read_excel('C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Export_Ctrl_Excel.xlsx',index_col=None)

Chromo_volume = pd.DataFrame(df_Ctrl, columns= ['Ctrl_chromo_volume_0','Ctrl_chromo_volume_1','Ctrl_chromo_volume_2'])
Chromo_volume_np_ctrl = np.array(Chromo_volume)

fig, ax = plt.subplots(figsize=(6, 4))
categories = ['Actin 0', 'Actin 15', 'Actin 30']
# rectangular box plot
bplot = ax.boxplot(Chromo_volume_np_ctrl,
                     vert=True,  # vertical box alignment
                     patch_artist=True,  # fill with color
                     labels=categories,showmeans=True, meanline=True)

for median in bplot['medians']:
    median.set_color('black')

# fill with colors
colors = ['#e8792e', '#22e084', '#22bae0','#872fad']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)

# Labels 
ax.set_title('Chromocenter Volume (Ctrl)')
ax.set_ylabel('Volume ' r'$\mu m^{3}$')
ax.set_xlabel('Actin time (mins)')
    
# Graph Fixes
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
Chromo_count_np_ctrl = np.array(Chromo_count)

fig, ax = plt.subplots(figsize=(6, 4))
categories = ['Actin 0', 'Actin 15', 'Actin 30']
# rectangular box plot
bplot = ax.boxplot(Chromo_count_np_ctrl,
                     vert=True,  # vertical box alignment
                     patch_artist=True,  # fill with color
                     labels=categories,showmeans=True, meanline=True)

for median in bplot['medians']:
    median.set_color('black')
    
# fill with colors
colors = ['#e8792e', '#22e084', '#22bae0','#872fad']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
    
# Labels 
ax.set_title('Chromocenter Count (Ctrl)')
ax.set_ylabel('Counts')
ax.set_xlabel('Actin time (mins)')
# Graph Fixes
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
Nuclei_volume_np_ctrl = np.array(Nuclei_volume)

fig, ax = plt.subplots(figsize=(6, 4))
categories = ['Actin 0', 'Actin 15', 'Actin 30']
# rectangular box plot
bplot = ax.boxplot(Nuclei_volume_np_ctrl,
                     vert=True,  # vertical box alignment
                     patch_artist=True,  # fill with color
                     labels=categories,showmeans=True, meanline=True)

for median in bplot['medians']:
    median.set_color('black')
    
# fill with colors
colors = ['#e8792e', '#22e084', '#22bae0','#872fad']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
    
# Labels 
ax.set_title('Nuclei Volume (Ctrl)')
ax.set_ylabel('Volume ' r'$\mu m^{3}$')
ax.set_xlabel('Actin time (mins)')
# Graph Fixes
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.legend()
ax.grid()
#plt.savefig('Nuclei_volume_ctrl.png', dpi=300, bbox_inches="tight", pad_inches=0.0)

plt.show() 