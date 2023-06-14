"""
Created on Wed Jan  4 21:55:14 2023
@author: Kaabir
"""
import pandas as pd
import glob
import os
import csv
from pathlib import Path
import numpy as np

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

Imaris_vca = {"File_Name_VCA":[],"VCA_chromo_volume_0":[],"VCA_chromo_volume_1":[],"VCA_chromo_volume_2":[],"VCA_chromo_volume_3":[],
           "VCA_chromo_count_0":[],"VCA_chromo_count_1":[],"VCA_chromo_count_2":[],"VCA_chromo_count_3":[],
          "nuclei_vca_0":[],"nuclei_vca_1":[],"nuclei_vca_2":[],"nuclei_vca_3":[],
           "actin_area_30":[],"actin_area_60":[],"sphericity_vca_0":[],"sphericity_vca_1":[],"sphericity_vca_2":[],"sphericity_vca_3":[]}


Imaris_ctrl = {"File_Name_Ctrl":[],"Ctrl_chromo_volume_0":[],"Ctrl_chromo_volume_1":[],"Ctrl_chromo_volume_2":[],"Ctrl_chromo_volume_3":[],
           "Ctrl_chromo_count_0":[],"Ctrl_chromo_count_1":[],"Ctrl_chromo_count_2":[],"Ctrl_chromo_count_3":[],
           "nuclei_ctrl_0":[],"nuclei_ctrl_1":[],"nuclei_ctrl_2":[],"nuclei_ctrl_3":[],
          "sphericity_ctrl_0":[],"sphericity_ctrl_1":[],"sphericity_ctrl_2":[],"sphericity_ctrl_3":[]}
            
# Type- Actin, Ctrl  && Time - T0,T15,T30,T60
# ===========================================

# Read the Tif/CZI File
def folder_scan(directory, marker):
    # Actin after 30min - value assigned is 3
    get_files = []
    extension = ".xlsx"  # ".czi"
    for f_name in os.listdir(directory):
        #print(f_name)
        if f_name.startswith(marker) and f_name.endswith(extension):
            get_files.append(os.path.join(directory, f_name))
    return get_files

directory = 'D:/Nuceli_Data/Fixed Nulcei Staining/230406 OF1 D1 Actin 10uMVCA Gactin RH3K9me2 FRH3K27me2/Result/Actin/'

# H3K9me2 files
h3K9me2_files = folder_scan(directory,'(H3K9me2)') # 

for temp_h3K9me2 in h3K9me2_files:
    df_ctrl_2 = pd.read_excel(temp_h3K9me2)
    Imaris_vca["VCA_chromo_count_0"].append(df_ctrl_2.loc[0,'mean'])
h3K9me2_files.clear() 

# # Sorting Files

                        
def actin_Coverage(Nuc_Area,Actin_Area):
           # Pd to np et divide
    Nuclei_area_den = list(map(float, Nuc_Area))
    Act_area_to_Flt = list(map(float, Actin_Area))
    value = 2
    Area_Div = [x / value for x in Act_area_to_Flt] 
    # Actin coverage in perentage       
    Actin_area_num = []
    for i in Area_Div :
        Actin_area_num.append(i * 100)
    #Actin_coverage_per = []
    # get last index for the lists for iteration
    end_index = len(Actin_area_num)
    for i in range(end_index):
        Actin_coverage_per = (Actin_area_num[i]/Nuclei_area_den[i])
    return Actin_coverage_per
    #Actin_coverage_per.clear()
    


def get_Filename(fil_Str):
    #filename = Path(file_path).stem
    File_Name=fil_Str[50:]
    index = File_Name.find('Average') # Find endswith key to locate the image number
    index = index - 3 
    File_Name = File_Name[index:index+2]
    File_Name = File_Name.split(',') 
    File_Name = File_Name[0] #extracts the first field
    return File_Name

