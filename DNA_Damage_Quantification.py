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

# =============================================================================
# actin_files =[]
# DNAdots_files =[]
# DNAinten_files = []
# Actin_area = []        
# Nuclei_area = []
# Ar_Df = []
# =============================================================================
Imaris_vca = {"VCA_DNAdots_Intensity":[], "VCA_DNAdots_Cnt":[],"VCA_DNAdots_Volume":[],
              "VCA_DNAd_Intensity":[], "VCA_Nuclei_Volume":[], "VCA_DNAd_Sphericity":[],
              "VCA_Actin_Area":[]}

Imaris_ctrl = {"Ctrl_DNAdots_Intensity":[], "Ctrl_DNAdots_Cnt":[],"Ctrl_DNAdots_Volume":[],
               "Ctrl_DNAd_Intensity":[], "Ctrl_Nuclei_Volume":[], "Ctrl_DNAd_Sphericity":[]}


# =============================================================================
# With Actin -AvecActin
# =============================================================================
actin_files =[]
DNAdots_files =[]
DNAinten_files = []
Actin_area = []        
Nuclei_area = []

# Actin after  3min  - value assigned is 3
path = 'C:/Users/kaabi/Documents/Nuceli_Data/20221028_DNA_Damage_RgH2AX/Output/AvecActin/'
extension = 'csv'
os.chdir(path)
# Ctrl
path_ = glob.glob('*/*.{}'.format(extension))

for f_name in path_:
    if f_name.startswith('(DNAdots)'):
        DNAdots_files.append(f_name)
    elif f_name.startswith('(DNAinten)'):
        DNAinten_files.append(f_name)
    elif f_name.startswith('(Actin)'):
        actin_files.append(f_name)        
    else:
        print('Some random file in the folder')
for DNAdots in DNAdots_files:
     # Delimiter
     filename_delimiter = ','
     # The max column count a line in the file could have
     largest_column_count = 0
     # Loop the data lines
     with open(DNAdots) as temp_f:
     # get No of columns in each line
         col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
     # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
         column_names = [i for i in range(0, max(col_count))]
     # Read csv
         df = pd.read_csv(DNAdots, delimiter=",", names=column_names)
         df.columns = df.iloc[2]
         df = df[3:]
         df = df.set_index(df.columns[0])
         
         # When multiple entry of same paramter is encountered it will take the last one
         # it is trciky to know the exact channel we used in imaris (actin or nucleus)
         # maintain consistent workflow to create channels, i.e nucleus -> actin
         someVar = []
         df_temp = []
         df_temp.append(df.loc['Intensity Mean','Mean'])
         for row in df_temp:
             for i in range(len(row)):
                 if i == len(row)-1:
                     someVar=row[i]
                     Imaris_vca["VCA_DNAdots_Intensity"].append(someVar)

         Imaris_vca["VCA_DNAdots_Volume"].append(df.loc['Volume','Sum'])
         Imaris_vca["VCA_DNAdots_Cnt"].append(df.loc['Volume','Count'])
         
for nucread in DNAinten_files:
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
         
         someVar2 = []
         df_temp2 = []
         df_temp2.append(df1.loc['Intensity Mean','Mean'])
         for row in df_temp2:
             for i in range(len(row)):
                 if i == len(row)-1:
                     someVar2=row[i]
                     Imaris_vca["VCA_DNAd_Intensity"].append(someVar2)
         
         
         Imaris_vca["VCA_Nuclei_Volume"].append(df1.loc['Volume','Sum'])
         Nuclei_area.append(df1.loc['Volume','Sum'])
         Imaris_vca["VCA_DNAd_Sphericity"].append(df1.loc['Sphericity','Mean']) 
         
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
    Imaris_vca["VCA_Actin_Area"].append(Actin_coverage_per)   
    
# =============================================================================
# Ctrl -SansActin
# =============================================================================

actin_files =[]
DNAdots_files =[]
DNAinten_files = []
Actin_area = []        
Nuclei_area = []

# Actin after  3min  - value assigned is 3
path = 'C:/Users/kaabi/Documents/Nuceli_Data/20221028_DNA_Damage_RgH2AX/Output/SansActin/'
extension = 'csv'
os.chdir(path)
# Ctrl
path_ = glob.glob('*/*.{}'.format(extension))
         
for f_name in path_:
    if f_name.startswith('(DNAdots)'):
        DNAdots_files.append(f_name)
    elif f_name.startswith('(DNAinten)'):
        DNAinten_files.append(f_name)
    elif f_name.startswith('(Actin)'):
        actin_files.append(f_name)          
    else:
        print('Some random file in the folder')
for DNAdots in DNAdots_files:
     # Delimiter
     filename_delimiter = ','
     # The max column count a line in the file could have
     largest_column_count = 0
     # Loop the data lines
     with open(DNAdots) as temp_f:
     # get No of columns in each line
         col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
     # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
         column_names = [i for i in range(0, max(col_count))]
     # Read csv
         df = pd.read_csv(DNAdots, delimiter=",", names=column_names)
         df.columns = df.iloc[2]
         df = df[3:]
         df = df.set_index(df.columns[0])
         
         # When multiple entry of same paramter is encountered it will take the last one
         # it is trciky to know the exact channel we used in imaris (actin or nucleus)
         # maintain consistent workflow to create channels, i.e nucleus -> actin
         someVar = []
         df_temp = []
         df_temp.append(df.loc['Intensity Mean','Mean'])
         for row in df_temp:
             for i in range(len(row)):
                 if i == len(row)-1:
                     someVar=row[i]
                     Imaris_ctrl["Ctrl_DNAdots_Intensity"].append(someVar)

         Imaris_ctrl["Ctrl_DNAdots_Volume"].append(df.loc['Volume','Sum'])
         Imaris_ctrl["Ctrl_DNAdots_Cnt"].append(df.loc['Volume','Count'])
         
for nucread in DNAinten_files:
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
         
         someVar2 = []
         df_temp2 = []
         df_temp2.append(df1.loc['Intensity Mean','Mean'])
         for row in df_temp2:
             for i in range(len(row)):
                 if i == len(row)-1:
                     someVar2=row[i]
                     Imaris_ctrl["Ctrl_DNAd_Intensity"].append(someVar2)
         
         
         Imaris_ctrl["Ctrl_Nuclei_Volume"].append(df1.loc['Volume','Sum'])
         Nuclei_area.append(df1.loc['Volume','Sum'])
         Imaris_ctrl["Ctrl_DNAd_Sphericity"].append(df1.loc['Sphericity','Mean']) 

# =============================================================================
# Export Data Ctrl
with open ('C:/Users/kaabi/Documents/Nuceli_Data/20221028_DNA_Damage_RgH2AX/Output/Output/export_ctrl.csv','w', newline='') as csv_file:
    csv.writer(csv_file).writerows([k, *v] for k,v in Imaris_ctrl.items())
    
df4 = pd.DataFrame.from_dict(Imaris_ctrl, orient='index')
df4 = df4.transpose()
pd.DataFrame(df4).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/20221028_DNA_Damage_RgH2AX/Output/Output/export_ctrl.xlsx')

# Export Data actin
with open ('C:/Users/kaabi/Documents/Nuceli_Data/20221028_DNA_Damage_RgH2AX/Output/Output/export_vca.csv','w', newline='') as csv_file:
    csv.writer(csv_file).writerows([k, *v] for k,v in Imaris_vca.items())
    
df5 = pd.DataFrame.from_dict(Imaris_vca, orient='index')
df5 = df5.transpose()
pd.DataFrame(df5).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/20221028_DNA_Damage_RgH2AX/Output/Output/export_vca.xlsx')

#Actin coverage less than 30%, 30-70% et >70%
Actin_L_30 = df5[df5['VCA_Actin_Area'] < 30]
Actin_L_30.to_csv(index=False, path_or_buf='C:/Users/kaabi/Documents/Nuceli_Data/20221028_DNA_Damage_RgH2AX/Output/Output/Actin_L_30.csv')

Actin_B_30_to_70 = df5[df5['VCA_Actin_Area'].between(30, 70, inclusive="both")]
Actin_B_30_to_70.to_csv(index=False, path_or_buf='C:/Users/kaabi/Documents/Nuceli_Data/20221028_DNA_Damage_RgH2AX/Output/Output/Actin_B_30_to_70.csv')

Actin_O_70 = df5[df5['VCA_Actin_Area'] > 70]
Actin_O_70.to_csv(index=False, path_or_buf='C:/Users/kaabi/Documents/Nuceli_Data/20221028_DNA_Damage_RgH2AX/Output/Output/Actin_O_70.csv')


# =============================================================================
# Plotting


"""
Imaris_vca = {"VCA_DNAdots_Intensity":[], "VCA_DNAdots_Cnt":[],"VCA_DNAdots_Volume":[],
              "VCA_DNAd_Intensity":[], "VCA_Nuclei_Volume":[], "VCA_DNAd_Sphericity":[]}

Imaris_ctrl = {"Ctrl_DNAdots_Intensity":[], "Ctrl_DNAdots_Cnt":[],"Ctrl_DNAdots_Volume":[],
               "Ctrl_DNAd_Intensity":[], "Ctrl_Nuclei_Volume":[], "Ctrl_DNAd_Sphericity":[]}
"""

df_sns = pd.concat([df4, df5], axis=1)


import seaborn as sns
sns.set_theme()

plt.figure()

#sns.pairplot(df_sns['VCA_DNAdots_Intensity', 'VCA_DNAdots_Cnt', 'Ctrl_DNAdots_Intensity', 'Ctrl_DNAdots_Cnt'], hue = "Outcome", markers=["o", "s"])

plt.show()

#g = sns.catplot(x=df_sns, y="vals", hue='cols', data=df_sns, kind='point')
