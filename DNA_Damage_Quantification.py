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
              "VCA_DNAd_Intensity":[], "VCA_Nuclei_Volume":[], "VCA_DNAd_Sphericity":[]}

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
