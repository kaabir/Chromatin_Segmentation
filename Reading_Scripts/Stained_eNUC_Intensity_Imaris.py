# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 21:55:14 2023
@author: Kaabir
"""
import pandas as pd
import glob
import os
import csv
from pathlib import Path

"""
Output
├── Ctrl
│   ├── T0
│   │   ├── Nuclei├── Volume/Sphericity/
│   │   └── Chromo└── Count/Volume/
│   ├── T15...
│   └── T30...
└── VCA
    ├── T0
    │   ├── Nuclei├── Volume/Sphericity/
    │   ├── Chromo├── Count/Volume/
    │   └── Actin └── Area/
    ├── T15...
    └── T30...
"""

Imaris_vca = {"File_Name_VCA":[],"VCA_H3K27me3_intensity":[],"VCA_H3K9me2_intensity":[],
           "VCA_nucleus_intensity":[],"VCA_chromo_count_0":[],"nuclei_vca_0":[],
          "nuclei_vca_1":[],"actin_area_30":[], "actin_area_60":[],
          "sphericity_vca_0":[],"sphericity_vca_1":[]}

Imaris_ctrl = {"File_Name_Ctrl":[],"Ctrl_H3K27me3_intensity":[],"Ctrl_H3K9me2_intensity":[],
           "Ctrl_nucleus_intensity":[],"Ctrl_chromo_count_0":[],"Ctrl_chromo_count_1":[],
           "nuclei_ctrl_0":[],"nuclei_ctrl_1":[],"nuclei_ctrl_2":[],"nuclei_ctrl_3":[],
          "sphericity_ctrl_0":[],"sphericity_ctrl_1":[]}
            
# Type- Actin, Ctrl  && Time - T0,T15,T30,T60
# ===========================================

# Sorting Files
actin_files =[]
h3K9me2_files =[]
h3K27me3_files =[]
nuclei_files = [] 

folder_path = 'C:/Users/kaabi/Documents/Nuceli_Data/Enucleation/230404_Imaris/'

def folder_Scan(Type):
    # Actin after  30min  - value assigned is 3
    directory = folder_path
    os.chdir(directory)
    global_path = glob.glob('*')
    extension = 'csv'
    for g_name in global_path:
        # Scan Actin/Ctrl first
        if g_name.startswith(Type):
                os.chdir(g_name)
                path_local = glob.glob('*/*.{}'.format(extension))
                for f_name in path_local:
                    #print(f_name)
                    if f_name.startswith('(H3K9me2)'):
                        h3K9me2_files.append(f_name)
                    elif f_name.startswith('(H3K27me3)'):
                        h3K27me3_files.append(f_name)                       
                    elif f_name.startswith('(Nuclei)'):
                        nuclei_files.append(f_name)
                    elif f_name.startswith('(Actin)'):
                        actin_files.append(f_name)
                    else:
                        print('Some random file or folders in the directory',f_name)
                                         
# Reading CSV
def read_CSV(fil_Nam):
    filename_delimiter = ','
    largest_column_count = 0

    with open(fil_Nam) as temp_fil:
        # Get the maximum column count from any row
        for line in temp_fil:
            column_count = len(line.split(filename_delimiter))
            largest_column_count = max(largest_column_count, column_count)
        
        # Read the CSV file without specifying column names
        df = pd.read_csv(fil_Nam, delimiter=filename_delimiter, header=None, skiprows=3, index_col=0)

    # Assign the values of the first row as column names
    df.columns = df.iloc[0]

    # Drop the first row since it's now used as column names
    df = df[1:]  
    
    return df    


def process_string(input_string, marker = None):
    # Find the last occurrence of "(Chromo)"
    
    directory = os.path.splitext(input_string)[0]
    path = os.path.basename(directory )  
    identifier = path.split(marker)[1]
    
    return identifier

# ===============
# Actin Channel
# ===============
#### T0
##
actin_files0 = folder_Scan('Actin') # Scan

# Nuclei Read
nuclei_data = {}
Nuclei_area_T0 = []

for nucleiread in nuclei_files:
    
    file_name = process_string(nucleiread, marker ="(Nuclei) ")

    #file_name = nucleiread.split("(Nucleus)")[1].split(".")[0].strip()
    
    Imaris_vca["File_Name_VCA"].append(file_name) 
    
    if file_name in Imaris_vca["File_Name_VCA"]:
        CT_nuclei_TO = read_CSV(nucleiread)
        nuclei_volume = CT_nuclei_TO.loc['Volume','Sum']
        nuclei_intensity = CT_nuclei_TO.loc['Intensity Mean'].iloc[3]['Mean']
        nuclei_sphericity = CT_nuclei_TO.loc['Sphericity','Mean']
        nuclei_data[file_name] = (nuclei_volume, nuclei_intensity, nuclei_sphericity)

    else:
        # If the file_name is not found in the 'stats_Actin' DataFrame, set both values in the tuple to 'None'
        nuclei_data[file_name] = ('None', 'None', 'None')

# After processing all files, append the values to the 'stats_Actin' DataFrame
Imaris_vca["nuclei_vca_0"] = [nuclei_data.get(file_name, ('None', 'None','None'))[0] for file_name in Imaris_vca["File_Name_VCA"]]
Imaris_vca["VCA_nucleus_intensity"] = [nuclei_data.get(file_name, ('None', 'None','None'))[1] for file_name in Imaris_vca["File_Name_VCA"]]
Imaris_vca["sphericity_vca_0"] = [nuclei_data.get(file_name, ('None', 'None','None'))[2] for file_name in Imaris_vca["File_Name_VCA"]]
Nuclei_area_T0 = [nuclei_data.get(file_name, ('None', 'None','None'))[0] for file_name in Imaris_vca["File_Name_VCA"]]

# H3k27me3

h3K27me3_data = {}

for h3K27me3read in h3K27me3_files:
    
    file_name = process_string(h3K27me3read, marker = "(H3K27me3) ")
    #file_name = h3K27me3read.split("(H3K27me3)")[1].split(".")[0].strip()
    
    if file_name in Imaris_vca["File_Name_VCA"]:

        AC_h3K27me3 = read_CSV(h3K27me3read)
        h3K27me3_intensity = AC_h3K27me3.loc['Intensity Mean'].iloc[1]['Mean']
        #nuclei_sphericity = CT_h3K27me3.loc['Sphericity','Mean']
        h3K27me3_data[file_name] = (h3K27me3_intensity)

    else:
        # If the file_name is not found in the 'stats_Actin' DataFrame, set both values in the tuple to 'None'
        h3K27me3_data[file_name] = ('None')

# After processing all files, append the values to the 'stats_Actin' DataFrame
Imaris_vca["VCA_H3K27me3_intensity"] = [h3K27me3_data.get(file_name, ('None')) for file_name in Imaris_vca["File_Name_VCA"]]

# H3K9me2

h3K9me2_data = {}

for h3K9me2read in h3K9me2_files:
    
    file_name = process_string(h3K9me2read, marker = "(H3K9me2) ")
    #file_name = h3K9me2read.split("(H3K9me2)")[1].split(".")[0].strip()
    
    if file_name in Imaris_vca["File_Name_VCA"]:

        AC_h3K9me2 = read_CSV(h3K9me2read)
        h3K9me2_intensity = AC_h3K9me2.loc['Intensity Mean'].iloc[2]['Mean']
        #nuclei_sphericity = CT_h3K27me3.loc['Sphericity','Mean']
        h3K9me2_data[file_name] = (h3K9me2_intensity)

    else:
        # If the file_name is not found in the 'stats_Actin' DataFrame, set both values in the tuple to 'None'
        h3K9me2_data[file_name] = ('None')

# After processing all files, append the values to the 'stats_Actin' DataFrame
Imaris_vca["VCA_H3K9me2_intensity"] = [h3K9me2_data.get(file_name, ('None')) for file_name in Imaris_vca["File_Name_VCA"]]

actin_data = {}

Actin_area_T0 = [] 

for actinread in actin_files:
    
    file_name = process_string(actinread, marker = "(Actin) ")
    #file_name = actinread.split("(Actin)_")[0].split(".")[0].strip()

    if file_name in Imaris_vca["File_Name_VCA"]:
   
        AC_actin_TO = read_CSV(actinread)
        # Store the 'Actin Coverage' value for the current file in the dictionary
        actin_data[file_name] = AC_actin_TO.loc['Area','Sum']

        
        # print(file_name, df.loc[0, 'Actin Coverage'])
    else:
        # If the file_name is not found in stats_Actin["File_Name_Actin"], set it to 'None'
        actin_data[file_name] = 'None'

# After processing all files, append the values to the 'stats_Actin' DataFrame
Actin_area_T0 = [actin_data.get(file_name, 'None') for file_name in Imaris_vca["File_Name_VCA"]]  

# Actin_area = preprocess_numeric_data(Actin_area_T0)
# Nuclei_area = preprocess_numeric_data(Nuclei_area_T0)
          
def actin_Coverage(nuc_area, actin_area):
    if actin_area == "None":
        Imaris_vca["actin_area_30"].append(0)
        return None
    return (float(actin_area) / 2) / float(nuc_area) * 100

coverage_dict = {}

for idx, (nuc, actin) in enumerate(zip(Nuclei_area_T0, Actin_area_T0)):
    coverage = actin_Coverage(nuc, actin)
    if coverage is not None:
        Imaris_vca["actin_area_30"].append(coverage)
      
# # ================
# # Control Channel
# # ================
# #### T0
# ##
actin_files.clear()  
h3K9me2_files.clear()  
h3K27me3_files.clear()  
nuclei_files.clear()  

ctrl_files0 = folder_Scan('Ctrl') # Scan

# Nuclei Read
nuclei_data = {}

for nucleiread in nuclei_files:
    
    file_name = process_string(nucleiread, marker ="(Nuclei) ")
    Imaris_ctrl["File_Name_Ctrl"].append(file_name) 
    
    if file_name in Imaris_ctrl["File_Name_Ctrl"]:
        CT_nuclei_TO = read_CSV(nucleiread)
        nuclei_volume = CT_nuclei_TO.loc['Volume','Sum']
        nuclei_intensity = CT_nuclei_TO.loc['Intensity Mean'].iloc[3]['Mean']
        nuclei_sphericity = CT_nuclei_TO.loc['Sphericity','Mean']
        nuclei_data[file_name] = (nuclei_volume, nuclei_intensity, nuclei_sphericity)

    else:
        # If the file_name is not found in the 'stats_Actin' DataFrame, set both values in the tuple to 'None'
        nuclei_data[file_name] = ('None', 'None', 'None')

# After processing all files, append the values to the 'stats_Actin' DataFrame
Imaris_ctrl["nuclei_ctrl_0"] = [nuclei_data.get(file_name, ('None', 'None','None'))[0] for file_name in Imaris_ctrl["File_Name_Ctrl"]]
Imaris_ctrl["Ctrl_nucleus_intensity"] = [nuclei_data.get(file_name, ('None', 'None','None'))[1] for file_name in Imaris_ctrl["File_Name_Ctrl"]]
Imaris_ctrl["sphericity_ctrl_0"] = [nuclei_data.get(file_name, ('None', 'None','None'))[2] for file_name in Imaris_ctrl["File_Name_Ctrl"]]

# H3k27me3

h3K27me3_data = {}

for h3K27me3read in h3K27me3_files:
    
    file_name = process_string(h3K27me3read, marker ="(H3K27me3) ")
    
    if file_name in Imaris_ctrl["File_Name_Ctrl"]:
        CT_h3K27me3 = read_CSV(h3K27me3read)
        h3K27me3_intensity = CT_h3K27me3.loc['Intensity Mean'].iloc[1]['Mean']
        #nuclei_sphericity = CT_h3K27me3.loc['Sphericity','Mean']
        h3K27me3_data[file_name] = (h3K27me3_intensity)

    else:
        # If the file_name is not found in the 'stats_Actin' DataFrame, set both values in the tuple to 'None'
        h3K27me3_data[file_name] = ('None')

# After processing all files, append the values to the 'stats_Actin' DataFrame
Imaris_ctrl["Ctrl_H3K27me3_intensity"] = [h3K27me3_data.get(file_name, ('None')) for file_name in Imaris_ctrl["File_Name_Ctrl"]]

# H3K9me2

h3K9me2_data = {}

for h3K9me2read in h3K9me2_files:
    
    file_name = process_string(h3K9me2read, marker ="(H3K9me2) ")
    
    if file_name in Imaris_ctrl["File_Name_Ctrl"]:
        CT_h3K9me2 = read_CSV(h3K9me2read)
        h3K9me2_intensity = CT_h3K9me2.loc['Intensity Mean'].iloc[2]['Mean']
        #nuclei_sphericity = CT_h3K27me3.loc['Sphericity','Mean']
        h3K9me2_data[file_name] = (h3K9me2_intensity)

    else:
        # If the file_name is not found in the 'stats_Actin' DataFrame, set both values in the tuple to 'None'
        h3K9me2_data[file_name] = ('None')

# After processing all files, append the values to the 'stats_Actin' DataFrame
Imaris_ctrl["Ctrl_H3K9me2_intensity"] = [h3K9me2_data.get(file_name, ('None')) for file_name in Imaris_ctrl["File_Name_Ctrl"]]

##
# Export Ctrl Data 
Result_folder = os.path.join(folder_path, 'Result')
if not os.path.exists(Result_folder):
    os.makedirs(Result_folder)
    
with open (Result_folder + '/Export_Ctrl.csv','w', newline='') as csv_file:
    csv.writer(csv_file).writerows([k, *v] for k,v in Imaris_ctrl.items())
    
df_ctrl = pd.DataFrame.from_dict(Imaris_ctrl, orient='index')
df_ctrl = df_ctrl.transpose()

colsCtrl = ['nuclei_ctrl_0', 'nuclei_ctrl_1', 'nuclei_ctrl_2','nuclei_ctrl_3']
df_ctrl[colsCtrl] = df_ctrl[colsCtrl].apply(pd.to_numeric, errors='coerce', axis=1)
df_ctrl['Ratio1'] = df_ctrl['nuclei_ctrl_0']/df_ctrl['nuclei_ctrl_1']
df_ctrl['Ratio2'] = df_ctrl['nuclei_ctrl_0']/df_ctrl['nuclei_ctrl_2']

# Sort based on Nuclei shrinkage >1.2  or <1.2 if inflates
Ctrl_O_R2 = df_ctrl[df_ctrl['Ratio2'] > 1.2]
pd.DataFrame(Ctrl_O_R2).to_excel(Result_folder + '/Export_Ctrl_O_R2.xlsx')
Ctrl_L_R2 = df_ctrl[df_ctrl['Ratio2'] < 1.2]
pd.DataFrame(Ctrl_L_R2).to_excel(Result_folder + '/Export_Ctrl_L_R2.xlsx')

pd.DataFrame(df_ctrl).to_excel(Result_folder + '/Export_Ctrl_Excel.xlsx')

# Export Actin Data
with open (Result_folder + '/Export_Actin.csv','w', newline='') as csv_file:
    csv.writer(csv_file).writerows([k, *v] for k,v in Imaris_vca.items())

df_actin = pd.DataFrame.from_dict(Imaris_vca, orient='index')
df_actin = df_actin.transpose()

# colsActin = ['nuclei_vca_0', 'nuclei_vca_1', 'nuclei_vca_2','nuclei_vca_3']
# df_actin[colsActin] = df_actin[colsActin].apply(pd.to_numeric, errors='coerce', axis=1)
# df_actin['Ratio1'] = df_actin['nuclei_vca_0']/df_actin['nuclei_vca_1']
# df_actin['Ratio2'] = df_actin['nuclei_vca_0']/df_actin['nuclei_vca_2']

pd.DataFrame(df_actin).to_excel(Result_folder + '/Export_Actin_Excel.xlsx')
# Sort Based on Coverage Percentage
Actin_L_30 = df_actin[df_actin['actin_area_30'] < 30]
#Actin_L_30.to_csv(index=False, path_or_buf='C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Actin_Low_30.csv')
# Actin_L_30_R2 = Actin_L_30[Actin_L_30['Ratio2'] > 1.2]
# pd.DataFrame(Actin_L_30_R2).to_excel(Result_folder + '/Export_Actin_Low_30_Over_SHR_1p2_R2.xlsx')
# Actin_L_30_R1 = Actin_L_30[Actin_L_30['Ratio2'] < 1.2]
# pd.DataFrame(Actin_L_30_R1).to_excel(Result_folder + '/Export_Actin_Low_30_Less_noShrink_1p2_R2.xlsx')

pd.DataFrame(Actin_L_30).to_excel(Result_folder + 'Actin_Low_30.xlsx')
Actin_B_30_to_70 = df_actin[df_actin['actin_area_60'].between(30, 70, inclusive="both")]
#Actin_B_30_to_70.to_csv(index=False, path_or_buf='C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Actin_Bt_30_to_70.csv')
Actin_O_70 = df_actin[df_actin['actin_area_60'] > 70]
#Actin_O_70.to_csv(index=False, path_or_buf='C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Actin_Ovr_70.csv')
Actin_O_30 = df_actin[df_actin['actin_area_30'] > 30]
Actin_O_30.to_csv(index=False, path_or_buf=Result_folder + '/Actin_Ovr_30.csv')

# # Sort based on Nuclei shrinkage >1.2  or <1.2 if inflates
# Actin_O_30_O_R2 = Actin_O_30[Actin_O_30['Ratio2'] > 1.2]
# pd.DataFrame(Actin_O_30_O_R2).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Export_Actin_Ovr_30_SHR_1p2_R2.xlsx')
# Actin_O_30_L_R2 = Actin_O_30[Actin_O_30['Ratio2'] < 1.2]
# pd.DataFrame(Actin_O_30_L_R2).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Export_Actin_Ovr_30_Less_noShrink_1p2_R2.xlsx')