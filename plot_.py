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

Imaris = {"Ctrl_chromo_volume_0":[],"Ctrl_chromo_volume_1":[],"Ctrl_chromo_volume_2":[],"Ctrl_chromo_volume_3":[],
           "Ctrl_chromo_count_0":[],"Ctrl_chromo_count_1":[],"Ctrl_chromo_count_2":[],"Ctrl_chromo_count_3":[],
           "VCA_chromo_volume_0":[],"VCA_chromo_volume_1":[],"VCA_chromo_volume_2":[],"VCA_chromo_volume_3":[],
           "VCA_chromo_count_0":[],"VCA_chromo_count_1":[],"VCA_chromo_count_2":[],"VCA_chromo_count_3":[],
           "nuclei_ctrl_0":[],"nuclei_ctrl_1":[],"nuclei_ctrl_2":[],"nuclei_ctrl_3":[], 
           "nuclei_vca_0":[],"nuclei_vca_1":[],"nuclei_vca_2":[],"nuclei_vca_3":[],
           "actin_area_30":[],"actin_area_60":[],
          "sphericity_ctrl_0":[],"sphericity_ctrl_1":[],"sphericity_ctrl_2":[],"sphericity_ctrl_3":[],
          "sphericity_vca_0":[],"sphericity_vca_1":[],"sphericity_vca_2":[],"sphericity_vca_3":[]}
# =============================================================================
#Actin
#Actin_1h
#Ctrl
# =============================================================================

# =============================================================================
#VCA Actin 60min
# =============================================================================
actin_files =[]
chromo_files =[]
nuclei_files = []
Actin_area = []        
Nuclei_area = []

# Actin after  3min  - value assigned is 3
path = 'C:/Users/kaabi/Documents/Nuceli_Data/20221006_Enculeation_RGC/Actin/'
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
         
         Imaris["VCA_chromo_volume_3"].append(df.loc['Volume','Sum'])
         Imaris["VCA_chromo_count_3"].append(df.loc['Volume','Count'])
         
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
         Imaris["nuclei_vca_3"].append(df1.loc['Volume','Sum'])
         Nuclei_area.append(df1.loc['Volume','Sum'])
         Imaris["sphericity_vca_3"].append(df1.loc['Sphericity','Mean']) 
         
         
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
    Imaris["actin_area_60"].append(Actin_coverage_per)       


# =============================================================================
# Control Positions Only Nuceli and Chromocenter 60min - value = 3
# Ctrl_0min
# Ctrl_15min
# Ctrl_30min
# Ctrl
# =============================================================================
actin_files =[]
chromo_files =[]
nuclei_files = []
Actin_area = []        
Nuclei_area = []

path = 'C:/Users/kaabi/Documents/Nuceli_Data/20221006_Enculeation_RGC/Ctrl/'
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
         
         Imaris["Ctrl_chromo_volume_3"].append(df.loc['Volume','Sum'])
         Imaris["Ctrl_chromo_count_3"].append(df.loc['Volume','Count'])
         
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
         Imaris["nuclei_ctrl_3"].append(df1.loc['Volume','Sum'])
         Nuclei_area.append(df1.loc['Volume','Sum'])
         Imaris["sphericity_ctrl_3"].append(df1.loc['Sphericity','Mean'])  

# =============================================================================
# Export Data
with open ('C:/Users/kaabi/Documents/Nuceli_Data/20221006_Enculeation_RGC/export.csv','w', newline='') as csv_file:
    csv.writer(csv_file).writerows([k, *v] for k,v in Imaris.items())
    
df4 = pd.DataFrame.from_dict(Imaris, orient='index')
df4 = df4.transpose()
pd.DataFrame(df4).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/20221006_Enculeation_RGC/export.xlsx')
# Export Data separately for Actin coverage percentage
Actin_L_30 = df4[df4['actin_area_60'] < 30]
Actin_L_30.to_csv(index=False, path_or_buf='C:/Users/kaabi/Documents/Nuceli_Data/20221006_Enculeation_RGC/Actin_L_30.csv')
Actin_B_30_to_70 = df4[df4['actin_area_60'].between(30, 70, inclusive="neither")]
Actin_B_30_to_70.to_csv(index=False, path_or_buf='C:/Users/kaabi/Documents/Nuceli_Data/20221006_Enculeation_RGC/Actin_B_30_to_70.csv')
Actin_O_70 = df4[df4['actin_area_60'] >= 70]
Actin_O_70.to_csv(index=False, path_or_buf='C:/Users/kaabi/Documents/Nuceli_Data/20221006_Enculeation_RGC/Actin_O_70.csv')

#pd.DataFrame(df4).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/20221006_Enculeation_RGC/export.xlsx')



# for filename in path_ctrl_chromo:
#     ##### plot U1
# 
#  	fig0= plt.figure(0)
#  	plt.plot(df["Area"],df["Sum"] ,'-')
#  	plt.xlabel("Nuclei")
#  	plt.ylabel("Chromocenter")
#  	plt.legend(loc="center left")
#  	plt.title("Nuclei vs Chromocenter")
#  	#fig0.savefig("comparison.png", dpi=400)
#plt.show()
# =============================================================================
# =============================================================================
# np.random.seed(8520)
# Per = range(len(Actin_coverage_per))
# N = 13
# colors = np.random.rand(N)
# fig, ax = plt.subplots(figsize=[5,4])
# plt.scatter(Per,Actin_coverage_per,c=colors, alpha=1)
# #ax.plot(Actin_coverage_per,'.', ls='-.', linewidth=1, label='Actin Coverage Percentage')
# #ax.plot(Chromo_number_ctrl,'.', ls='-.', linewidth=1,label='Hoescht Dots Volume')
# 
# 
# 
# plt.legend()
# plt.xlabel(r'$Nuclei$ $Count$')
# plt.ylabel(r'$Actin$ $Coverage$ $(percent)$')
# 
# #plt.savefig('1.png', dpi = 300, bbox_inches = 'tight')
# plt.show()
# 
# =============================================================================
