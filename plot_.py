# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 14:52:50 2022

@author: kaabir
"""
import pandas as pd
import glob
import os
import re
import csv
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cbook as cbook
from pathlib import Path

actin_files =[]
chromo_files =[]
nuclei_files = []
Actin_area = []        
Nuclei_area = []
Ar_Df = []

Imaris_ctrl = {"File_Name_Ctrl":[],"Ctrl_chromo_volume_0":[],"Ctrl_chromo_volume_1":[],"Ctrl_chromo_volume_2":[],"Ctrl_chromo_volume_3":[],
           "Ctrl_chromo_count_0":[],"Ctrl_chromo_count_1":[],"Ctrl_chromo_count_2":[],"Ctrl_chromo_count_3":[],
           "nuclei_ctrl_0":[],"nuclei_ctrl_1":[],"nuclei_ctrl_2":[],"nuclei_ctrl_3":[],
          "sphericity_ctrl_0":[],"sphericity_ctrl_1":[],"sphericity_ctrl_2":[],"sphericity_ctrl_3":[]}

Imaris_vca = {"File_Name_VCA":[],"VCA_chromo_volume_1":[],"VCA_chromo_volume_2":[],"VCA_chromo_volume_3":[],
           "VCA_chromo_count_1":[],"VCA_chromo_count_2":[],"VCA_chromo_count_3":[],
          "nuclei_vca_1":[],"nuclei_vca_2":[],"nuclei_vca_3":[],
           "actin_area_60":[],"sphericity_vca_0":[],"sphericity_vca_1":[],"sphericity_vca_2":[],"sphericity_vca_3":[]}

# =============================================================================
#VCA Actin 60min
# =============================================================================
actin_files =[]
chromo_files =[]
nuclei_files = []
Actin_area = []        
Nuclei_area = []

# Actin after  30min  - value assigned is 3
path = 'C:/Users/kaabi/Documents/Nuceli_Data/20221116 E-Nuc VCA_Ctrl/Actin/Output/T30/'
extension = 'csv'
os.chdir(path)
# Ctrl
path_ = glob.glob('*/*.{}'.format(extension))

for f_name in path_:
    if f_name.startswith('(Chromo)'):
        chromo_files.append(f_name)
    elif f_name.startswith('(NUC04)'):
        nuclei_files.append(f_name)
    elif f_name.startswith('(Actin)'):
        actin_files.append(f_name)
    else:
        print('Some random file in the folder')
for chromoread in chromo_files:
     # Delimiter
     filename_delimiter = ','
     # The max column count a line in the file could have
     largest_column_count = 0
     # Loop the data lines
     with open(chromoread) as temp_f:
     # get No of columns in each line
         col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
     # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
         column_names = [i for i in range(0, max(col_count))]
     # Read csv
         df = pd.read_csv(chromoread, delimiter=",", names=column_names)
         df.columns = df.iloc[2]
         df = df[3:]
         df = df.set_index(df.columns[0])
         
         #Get File Name
         File_Name=chromoread[48:]
         Imaris_vca["File_Name_VCA"].append(File_Name)
         ###############         
         
         Imaris_vca["VCA_chromo_volume_3"].append(df.loc['Volume','Sum'])
         Imaris_vca["VCA_chromo_count_3"].append(df.loc['Volume','Count'])
         
for nucread in nuclei_files:
     # Delimiter
     filename_delimiter = ','
     # The max column count a line in the file could have
     largest_column_count = 0
     # Loop the data lines
     with open(nucread) as temp_f:
     # get No of columns in each line
         col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
     # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
         column_names = [i for i in range(0, max(col_count))]
     # Read csv
         df1 = pd.read_csv(nucread, delimiter=",", names=column_names)
         #print(df1)
         df1.columns = df1.iloc[2]
         df1 = df1[3:]
         df1 = df1.set_index(df1.columns[0])
         Imaris_vca["nuclei_vca_3"].append(df1.loc['Volume','Sum'])
         Nuclei_area.append(df1.loc['Volume','Sum'])
         Imaris_vca["sphericity_vca_3"].append(df1.loc['Sphericity','Mean']) 
         
         
# For Actin Channel
for actread in actin_files:
     # Delimiter
     filename_delimiter = ','
     # The max column count a line in the file could have
     largest_column_count = 0
     # Loop the data lines
     with open(actread) as temp_f:
     # get No of columns in each line
         col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
     # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
         column_names = [i for i in range(0, max(col_count))]
     # Read csv
         df2 = pd.read_csv(actread, delimiter=",", names=column_names)
         df2.columns = df2.iloc[2]
         df2 = df2[3:]
         df2 = df2.set_index(df2.columns[0])

         Actin_area.append(df2.loc['Area','Sum'])
         
#Actin_area_divided = [float(x) for x in Actin_area]
Actin_area_divided = list(map(float, Actin_area))
value = 2
Actin_area_divided = [x / value for x in Actin_area_divided] 
   
# Actin coverage in perentage       
numerator = []
for number in Actin_area_divided:
    numerator.append(number * 100)

Nuclei_area_conv = list(map(float,Nuclei_area))

Actin_coverage_per = []
# get last index for the lists for iteration
end_index = len(numerator)

for i in range(end_index):
    Actin_coverage_per = (numerator[i]/Nuclei_area_conv[i])
    Imaris_vca["actin_area_60"].append(Actin_coverage_per)       


# =============================================================================
# Actin after  30min  - value assigned is 2
actin_files =[]
chromo_files =[]
nuclei_files = []
Actin_area = []        
Nuclei_area = []

path = 'C:/Users/kaabi/Documents/Nuceli_Data/20221116 E-Nuc VCA_Ctrl/Actin/Output/T15/'
extension = 'csv'
os.chdir(path)
# Ctrl
path_ = glob.glob('*/*.{}'.format(extension))

for f_name in path_:
    if f_name.startswith('(Chromo)'):
        chromo_files.append(f_name)
    elif f_name.startswith('(NUC04)'):
        nuclei_files.append(f_name)
    else:
        print('Some random file in the folder')
for chromoread in chromo_files:
     # Delimiter
     filename_delimiter = ','
     # The max column count a line in the file could have
     largest_column_count = 0
     # Loop the data lines
     with open(chromoread) as temp_f:
     # get No of columns in each line
         col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
     # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
         column_names = [i for i in range(0, max(col_count))]
     # Read csv
         df = pd.read_csv(chromoread, delimiter=",", names=column_names)
         df.columns = df.iloc[2]
         df = df[3:]
         df = df.set_index(df.columns[0])
         
         Imaris_vca["VCA_chromo_volume_2"].append(df.loc['Volume','Sum'])
         Imaris_vca["VCA_chromo_count_2"].append(df.loc['Volume','Count'])
         
for nucread in nuclei_files:
     # Delimiter
     filename_delimiter = ','
     # The max column count a line in the file could have
     largest_column_count = 0
     # Loop the data lines
     with open(nucread) as temp_f:
     # get No of columns in each line
         col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
     # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
         column_names = [i for i in range(0, max(col_count))]
     # Read csv
         df1 = pd.read_csv(nucread, delimiter=",", names=column_names)
         #print(df1)
         df1.columns = df1.iloc[2]
         df1 = df1[3:]
         df1 = df1.set_index(df1.columns[0])
         Imaris_vca["nuclei_vca_2"].append(df1.loc['Volume','Sum'])
         Nuclei_area.append(df1.loc['Volume','Sum'])
         Imaris_vca["sphericity_vca_2"].append(df1.loc['Sphericity','Mean']) 
         

# =============================================================================
#VCA 15min - value 1 and no actin here
# =============================================================================
actin_files =[]
chromo_files =[]
nuclei_files = []
Actin_area = []        
Nuclei_area = []

path = 'C:/Users/kaabi/Documents/Nuceli_Data/20221116 E-Nuc VCA_Ctrl/Actin/Output/T0/'
extension = 'csv'
os.chdir(path)
# Ctrl
path_ = glob.glob('*/*.{}'.format(extension))

for f_name in path_:
    if f_name.startswith('(Chromo)'):
        chromo_files.append(f_name)
    elif f_name.startswith('(NUC04)'):
        nuclei_files.append(f_name)
    elif f_name.startswith('(Actin)'):
        actin_files.append(f_name)
    else:
        print('Some random file in the folder')
for chromoread in chromo_files:
     # Delimiter
     filename_delimiter = ','
     # The max column count a line in the file could have
     largest_column_count = 0
     # Loop the data lines
     with open(chromoread) as temp_f:
     # get No of columns in each line
         col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
     # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
         column_names = [i for i in range(0, max(col_count))]
     # Read csv
         df = pd.read_csv(chromoread, delimiter=",", names=column_names)
         df.columns = df.iloc[2]
         df = df[3:]
         df = df.set_index(df.columns[0])
         
         Imaris_vca["VCA_chromo_volume_1"].append(df.loc['Volume','Sum'])
         Imaris_vca["VCA_chromo_count_1"].append(df.loc['Volume','Count'])
         
for nucread in nuclei_files:
     # Delimiter
     filename_delimiter = ','
     # The max column count a line in the file could have
     largest_column_count = 0
     # Loop the data lines
     with open(nucread) as temp_f:
     # get No of columns in each line
         col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
     # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
         column_names = [i for i in range(0, max(col_count))]
     # Read csv
         df1 = pd.read_csv(nucread, delimiter=",", names=column_names)
         #print(df1)
         df1.columns = df1.iloc[2]
         df1 = df1[3:]
         df1 = df1.set_index(df1.columns[0])
         Imaris_vca["nuclei_vca_1"].append(df1.loc['Volume','Sum'])
         Nuclei_area.append(df1.loc['Volume','Sum'])
         Imaris_vca["sphericity_vca_1"].append(df1.loc['Sphericity','Mean']) 
         



# =============================================================================
# Control Positions Only Nuceli and Chromocenter 60min - value = 3
# Ctrl_0min
# Ctrl_15min
# Ctrl_30min
# Ctrl_60min
# =============================================================================
actin_files =[]
chromo_files =[]
nuclei_files = []
Actin_area = []        
Nuclei_area = []

path = 'C:/Users/kaabi/Documents/Nuceli_Data/20221116 E-Nuc VCA_Ctrl/Ctrl/Output/T30/'
extension = 'csv'
os.chdir(path)
# Ctrl
path_ = glob.glob('*/*.{}'.format(extension))

for f_name in path_:
    if f_name.startswith('(Chromo)'):
        chromo_files.append(f_name)
    elif f_name.startswith('(NUC04)'):
        nuclei_files.append(f_name)
    else:
        print('Some random file in the folder')
for chromoread in chromo_files:
     # Delimiter
     filename_delimiter = ','
     # The max column count a line in the file could have
     largest_column_count = 0
     # Loop the data lines
     with open(chromoread) as temp_f:
     # get No of columns in each line
         col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
     # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
         column_names = [i for i in range(0, max(col_count))]
     # Read csv
         df = pd.read_csv(chromoread, delimiter=",", names=column_names)
         df.columns = df.iloc[2]
         df = df[3:]
         df = df.set_index(df.columns[0])
         
         Imaris_ctrl["Ctrl_chromo_volume_3"].append(df.loc['Volume','Sum'])
         Imaris_ctrl["Ctrl_chromo_count_3"].append(df.loc['Volume','Count'])
         
for nucread in nuclei_files:
     # Delimiter
     filename_delimiter = ','
     # The max column count a line in the file could have
     largest_column_count = 0
     # Loop the data lines
     with open(nucread) as temp_f:
     # get No of columns in each line
         col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
     # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
         column_names = [i for i in range(0, max(col_count))]
     # Read csv
         df1 = pd.read_csv(nucread, delimiter=",", names=column_names)
         #print(df1)
         df1.columns = df1.iloc[2]
         df1 = df1[3:]
         df1 = df1.set_index(df1.columns[0])
         #Get File Name
         File_Name_ctrl=nucread[48:]
         Imaris_ctrl["File_Name_Ctrl"].append(File_Name_ctrl)
         ###############               
         
         Imaris_ctrl["nuclei_ctrl_3"].append(df1.loc['Volume','Sum'])
         Nuclei_area.append(df1.loc['Volume','Sum'])
         Imaris_ctrl["sphericity_ctrl_3"].append(df1.loc['Sphericity','Mean'])  
# =============================================================================
# Control Positions Only Nuceli and Chromocenter 30min - value = 2
# =============================================================================
actin_files =[]
chromo_files =[]
nuclei_files = []
Actin_area = []        
Nuclei_area = []

path = 'C:/Users/kaabi/Documents/Nuceli_Data/20221116 E-Nuc VCA_Ctrl/Ctrl/Output/T15/'
extension = 'csv'
os.chdir(path)
# Ctrl
path_ = glob.glob('*/*.{}'.format(extension))

for f_name in path_:
    if f_name.startswith('(Chromo)'):
        chromo_files.append(f_name)
    elif f_name.startswith('(NUC04)'):
        nuclei_files.append(f_name)
    else:
        print('Some random file in the folder')
for chromoread in chromo_files:
     # Delimiter
     filename_delimiter = ','
     # The max column count a line in the file could have
     largest_column_count = 0
     # Loop the data lines
     with open(chromoread) as temp_f:
     # get No of columns in each line
         col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
     # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
         column_names = [i for i in range(0, max(col_count))]
     # Read csv
         df = pd.read_csv(chromoread, delimiter=",", names=column_names)
         df.columns = df.iloc[2]
         df = df[3:]
         df = df.set_index(df.columns[0])
         
         Imaris_ctrl["Ctrl_chromo_volume_2"].append(df.loc['Volume','Sum'])
         Imaris_ctrl["Ctrl_chromo_count_2"].append(df.loc['Volume','Count'])
         
for nucread in nuclei_files:
     # Delimiter
     filename_delimiter = ','
     # The max column count a line in the file could have
     largest_column_count = 0
     # Loop the data lines
     with open(nucread) as temp_f:
     # get No of columns in each line
         col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
     # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
         column_names = [i for i in range(0, max(col_count))]
     # Read csv
         df1 = pd.read_csv(nucread, delimiter=",", names=column_names)
         #print(df1)
         df1.columns = df1.iloc[2]
         df1 = df1[3:]
         df1 = df1.set_index(df1.columns[0])
         Imaris_ctrl["nuclei_ctrl_2"].append(df1.loc['Volume','Sum'])
         Nuclei_area.append(df1.loc['Volume','Sum'])
         Imaris_ctrl["sphericity_ctrl_2"].append(df1.loc['Sphericity','Mean'])  
# =============================================================================
# Control Positions Only Nuceli and Chromocenter 15min - value = 1
# =============================================================================
actin_files =[]
chromo_files =[]
nuclei_files = []
Actin_area = []        
Nuclei_area = []

path = 'C:/Users/kaabi/Documents/Nuceli_Data/20221116 E-Nuc VCA_Ctrl/Ctrl/Output/T0/'
extension = 'csv'
os.chdir(path)
# Ctrl
path_ = glob.glob('*/*.{}'.format(extension))

for f_name in path_:
    if f_name.startswith('(Chromo)'):
        chromo_files.append(f_name)
    elif f_name.startswith('(NUC04)'):
        nuclei_files.append(f_name)
    else:
        print('Some random file in the folder')
for chromoread in chromo_files:
     # Delimiter
     filename_delimiter = ','
     # The max column count a line in the file could have
     largest_column_count = 0
     # Loop the data lines
     with open(chromoread) as temp_f:
     # get No of columns in each line
         col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
     # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
         column_names = [i for i in range(0, max(col_count))]
     # Read csv
         df = pd.read_csv(chromoread, delimiter=",", names=column_names)
         df.columns = df.iloc[2]
         df = df[3:]
         df = df.set_index(df.columns[0])
         
         Imaris_ctrl["Ctrl_chromo_volume_1"].append(df.loc['Volume','Sum'])
         Imaris_ctrl["Ctrl_chromo_count_1"].append(df.loc['Volume','Count'])
         
for nucread in nuclei_files:
     # Delimiter
     filename_delimiter = ','
     # The max column count a line in the file could have
     largest_column_count = 0
     # Loop the data lines
     with open(nucread) as temp_f:
     # get No of columns in each line
         col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
     # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
         column_names = [i for i in range(0, max(col_count))]
     # Read csv
         df1 = pd.read_csv(nucread, delimiter=",", names=column_names)
         #print(df1)
         df1.columns = df1.iloc[2]
         df1 = df1[3:]
         df1 = df1.set_index(df1.columns[0])
         Imaris_ctrl["nuclei_ctrl_1"].append(df1.loc['Volume','Sum'])
         Nuclei_area.append(df1.loc['Volume','Sum'])
         Imaris_ctrl["sphericity_ctrl_1"].append(df1.loc['Sphericity','Mean'])  
         
# =============================================================================
# =============================================================================

Chromo_volume = pd.DataFrame(Imaris_vca, columns= ['VCA_chromo_volume_1','VCA_chromo_volume_2','VCA_chromo_volume_3'])
Chromo_volume_np = np.array(Chromo_volume)
# Export Data Ctrl

with open ('C:/Users/kaabi/Documents/Nuceli_Data/20221116 E-Nuc VCA_Ctrl/export_ctrl.csv','w', newline='') as csv_file:
    csv.writer(csv_file).writerows([k, *v] for k,v in Imaris_ctrl.items())
    
df4 = pd.DataFrame.from_dict(Imaris_ctrl, orient='index')
df4 = df4.transpose()
pd.DataFrame(df4).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/20221116 E-Nuc VCA_Ctrl/export_ctrl.xlsx')

# Export Data actin
with open ('C:/Users/kaabi/Documents/Nuceli_Data/20221116 E-Nuc VCA_Ctrl/export_vca.csv','w', newline='') as csv_file:
    csv.writer(csv_file).writerows([k, *v] for k,v in Imaris_vca.items())
    
df5 = pd.DataFrame.from_dict(Imaris_vca, orient='index')
df5 = df5.transpose()
pd.DataFrame(df5).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/20221116 E-Nuc VCA_Ctrl/export_actin.xlsx')

Actin_L_30 = df5[df5['actin_area_60'] < 30]
Actin_L_30.to_csv(index=False, path_or_buf='C:/Users/kaabi/Documents/Nuceli_Data/20221116 E-Nuc VCA_Ctrl/Actin_L_30.csv')
Actin_B_30_to_70 = df5[df5['actin_area_60'].between(30, 70, inclusive="both")]
Actin_B_30_to_70.to_csv(index=False, path_or_buf='C:/Users/kaabi/Documents/Nuceli_Data/20221116 E-Nuc VCA_Ctrl/Actin_B_30_to_70.csv')
Actin_O_70 = df5[df5['actin_area_60'] > 70]
Actin_O_70.to_csv(index=False, path_or_buf='C:/Users/kaabi/Documents/Nuceli_Data/20221116 E-Nuc VCA_Ctrl/Actin_O_70.csv')



coordinates = [
    dict(
        y=Chromo_volume_np[i,:],
        x= np.arange(0, 3,1))
    for i in range(len(Chromo_volume_np))
]


fig, ax = plt.subplots(figsize=(6, 3))
for c in coordinates:
    ax.plot(c['x'], c['y'],marker="o", alpha=0.8, linestyle='dashed')
    #ax.plot(x, y, marker="o", color='#960056', alpha=0.8, linestyle='dashed', label='Nuclei') 
    ax.set_title('Chromocenter Volume (Nuclei1)')
    categories = ['Actin 0', 'Actin 15', 'Actin 30']
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
# import re
# 
# Imaris_vca = {"VCA_DNAdots_Intensity":[], "VCA_DNAdots_Cnt":[],"VCA_DNAdots_Volume":[],
#               "VCA_DNAd_Intensity":[], "VCA_Nuclei_Volume":[], "VCA_DNAd_Sphericity":[],
#               "VCA_Actin_Cov_pER":[]}
# # https://medium.com/factory-mind/regex-tutorial-a-simple-cheatsheet-by-examples-649dc1c3f285
# 
# text='(Chromo) 221116 RGC E-Nuc 9_Statistics\(Chromo) 221116 RGC E-Nuc 20_Average.csv'
# txt_op_re=[str(s) for s in re.findall(r'(Chromo)*\[0-9]', text)]
# txt_op_re1=[str(s) for s in re.findall(r'[0-9]', text)] 
# i = len(txt_op_re1)-1
# someVar=txt_op_re1[i]
# Imaris_vca["VCA_DNAdots_Intensity"].append(someVar)
# =============================================================================
  
         
         #Get File Name
         #File_Name=chromoread.join([n for n in chromoread if n.isdigit()])
         #txt_op_re=[str(s) for s in re.findall(r'\b\w+\b', chromoread)]
         
         #txt = chromoread
         #txt_op= [int(s) for s in txt if s.isdigit()]
         
         #Imaris_vca["File_Name_VCA"].append(File_Name)
# =============================================================================    