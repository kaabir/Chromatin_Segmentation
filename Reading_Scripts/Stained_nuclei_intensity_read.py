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
│      ├── Nuclei├── Volume/Sphericity/
│      └── Chromo└── Count/Volume/
│   
│  
└── VCA
    ├── T0
       ├── Nuclei├── Volume/Sphericity/
       ├── Chromo├── Count/Volume/
       └── Actin └── Area/

    
    Hoescht B H3K9me2-A568 H3k27me3-A647
    
"""

region_Prop_ctrl = {"ctrl_intens_H3k27me3":[],"ctrl_intens_H3K9me2":[],"ctrl_intens_Hoechst":[]}
region_Prop_actin = {"actin_intens_H3k27me3":[],"actin_intens_H3K9me2":[],"actin_intens_Hoechst":[]}


# Imaris_ctrl = {"File_Name_Ctrl":[],"Ctrl_chromo_volume_0":[],"Ctrl_chromo_volume_1":[],"Ctrl_chromo_volume_2":[],"Ctrl_chromo_volume_3":[],
#            "Ctrl_chromo_count_0":[],"Ctrl_chromo_count_1":[],"Ctrl_chromo_count_2":[],"Ctrl_chromo_count_3":[],
#            "nuclei_ctrl_0":[],"nuclei_ctrl_1":[],"nuclei_ctrl_2":[],"nuclei_ctrl_3":[],
#           "sphericity_ctrl_0":[],"sphericity_ctrl_1":[],"sphericity_ctrl_2":[],"sphericity_ctrl_3":[]}
            
# Type- Actin, Ctrl  && Time - T0,T15,T30,T60
# ===========================================

# Read the Tif/CZI File
def folder_scan(directory, extension, marker=None):
    # folder_scan(chromo_folder, ".xlsx", marker="(Chromo)")
    get_files = []
    extension = extension  # ".xlsx"  # ".czi"
    for f_name in os.listdir(directory):
        if marker is not None:
            if f_name.find(marker) != -1 and f_name.endswith(extension):
                get_files.append(os.path.join(directory, f_name))
        else:
            if f_name.endswith(extension):
                get_files.append(os.path.join(directory, f_name))
    return get_files
        
def matchingFiles(sortedCells, sourceCells,sortedMarker=None, sourceMarker=None):
    # Here I have filenames which are segmented labels of individual cells
    # I want to sort the RGC images with their nucleus and choromo markers
    # Note: not all nucleus with the segmentation would have chromocenter 
    # To sort these and match correct files this function matches the files
    # Example Use - ChromoFiles =  matchingFiles(rgc_Scan, chromo_Scan,sortedMarker=None, sourceMarker="(Chromo)_")
    # Example Use -NucleusFiles =  matchingFiles(ChromoFiles, nucleus_Scan, sortedMarker="(Chromo)_", sourceMarker="(Nucleus)_")
    
    identifiers = []

    matching_files = []
    
    if sortedMarker is None:
        for temp_files in sortedCells:
            file_name = os.path.basename(temp_files)
            identifier = os.path.splitext(file_name)[0]
            identifiers.append(identifier)
    else:
        for temp_files in sortedCells:
            file_name = temp_files.split(sortedMarker)[1].split(".")[0].strip() #
            identifier = file_name
            identifiers.append(identifier)
            
    for identifier in identifiers:
        #matching_files = []
        # Find matching files in the source folder
        for source_file in sourceCells:
            file_name = source_file.split(sourceMarker)[1].split(".")[0].strip() #

            if identifier == file_name:
                matching_files.append(source_file)
                #print(matching_files)
                
    return matching_files   

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

# Ctrl

folder_path = 'C:/Users/kaabi/Documents/Nuceli_Data/Fixed Nulcei Staining/230406 OF1 D1 Actin 10uMVCA Gactin RH3K9me2 FRH3K27me2/'


directory = folder_path + 'Ctrl/'

# H3K9me2 files
h3K9me2_files = folder_scan(directory,"(H3K9me2)") # 

for temp_h3K9me2 in h3K9me2_files:
    df = pd.read_excel(temp_h3K9me2)
    region_Prop_ctrl["ctrl_intens_H3K9me2"].append(df.loc[0,'mean'])
h3K9me2_files.clear() 

# H3k27me3 files
h3k27me3_files = folder_scan(directory,"(H3K27me3)") #

for temp_h3k27me3 in h3k27me3_files:
    df = pd.read_excel(temp_h3k27me3)
    region_Prop_ctrl["ctrl_intens_H3k27me3"].append(df.loc[0,'mean'])
h3k27me3_files.clear() 

# Hoechst files
nuclei_Files = folder_scan(directory,"(Nucleus)Nucleus")

for temp_nuclei in nuclei_Files:
    df = pd.read_excel(temp_nuclei)
    region_Prop_ctrl["ctrl_intens_Hoechst"].append(df.loc[0,'mean'])
nuclei_Files.clear() 

df_ctrl = pd.DataFrame.from_dict(region_Prop_ctrl, orient='index')
df_ctrl = df_ctrl.transpose()
pd.DataFrame(df_ctrl).to_excel(folder_path + 'Export_Ctrl_Excel.xlsx')


# Actin
directory = folder_path + 'Actin/'


h3K9me2_files = folder_scan(directory,"(H3K9me2)") # 

for temp_h3K9me2 in h3K9me2_files:
    df = pd.read_excel(temp_h3K9me2)
    region_Prop_actin["actin_intens_H3K9me2"].append(df.loc[0,'mean'])
h3K9me2_files.clear() 

# H3k27me3 files
h3k27me3_files = folder_scan(directory,"(H3K27me3)") #

for temp_h3k27me3 in h3k27me3_files:
    df = pd.read_excel(temp_h3k27me3)
    region_Prop_actin["actin_intens_H3k27me3"].append(df.loc[0,'mean'])
h3k27me3_files.clear() 

# Hoechst files
nuclei_Files = folder_scan(directory,"(Nucleus)Nucleus")

for temp_nuclei in nuclei_Files:
    df = pd.read_excel(temp_nuclei)
    region_Prop_actin["actin_intens_Hoechst"].append(df.loc[0,'mean'])
nuclei_Files.clear() 

df_actin = pd.DataFrame.from_dict(region_Prop_actin, orient='index')
df_actin = df_actin.transpose()
pd.DataFrame(df_actin).to_excel(folder_path + 'Export_Actin_Excel.xlsx')
