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

stats_Actin = {"File_Name_Actin":[],"Actin_chromo_volume_T0":[],"Actin_chromo_volume_T1":[],"Actin_chromo_volume_T2":[],
                      "Actin_chromo_intensity_T0":[],"Actin_chromo_intensity_T1":[],"Actin_chromo_intensity_T2":[],
          "Nuclei_volume_actin_T0":[],"Nuclei_volume_actin_T1":[],"Nuclei_volume_actin_T2":[],
		  "Nuclei_intensity_actin_T0":[],"Nuclei_intensity_actin_T1":[],"Nuclei_intensity_actin_T2":[],
		  "Actin_coverage_T30":[],"Actin_coverage_T60":[],"sphericity_Actin_T0":[],"sphericity_Actin_T1":[],"sphericity_Actin_T2":[]}

stats_Ctrl = {"File_Name_Ctrl":[],"Ctrl_chromo_volume_T0":[],"Ctrl_chromo_volume_T1":[],"Ctrl_chromo_volume_T2":[],
           "Ctrl_chromo_intensity_T0":[],"Ctrl_chromo_intensity_T1":[],"Ctrl_chromo_intensity_T2":[],
           "Nuclei_volume_ctrl_T0":[],"Nuclei_volume_ctrl_T1":[],"Nuclei_volume_ctrl_T2":[],
		   "Nuclei_intensity_ctrl_T0":[],"Nuclei_intensity_ctrl_T1":[],"Nuclei_intensity_ctrl_T2":[],
          "Sphericity_ctrl_T0":[],"Sphericity_ctrl_T1":[],"Sphericity_ctrl_T2":[]}
            
# Type- Actin, Ctrl  && Time - T0, T15, T30, T60

def matchingFiles(sortedCells, sourceCells, sortedMarker=None, sourceMarker=None):   
    identifiers = [] # For my sorted/analyzed cell type/ less in number
    # Here I have filenames which are segmented labels of individual nuclei marker of hoescht
    # I want to sort the cell type (RGC, multi, flower, halo) images with their nucleus and choromo markers
    # Note: not all nucleus with the segmentation would have chromocenter 
    # To sort these and match correct files this function matches the files
    # Example Use - ChromoFiles =  matchingFiles(rgc_Scan, chromo_Scan,sortedMarker=None, sourceMarker="(Chromo)_")
    # Example Use -NucleusFiles =  matchingFiles(ChromoFiles, nucleus_Scan, sortedMarker="(Chromo)_", sourceMarker="(Nucleus)_")
        
    matching_files = []

    # First I collect files from my cell types that are Radial glial or multi ciliated etc.
    # In the first case the files are tif files and not markers
    # I can extract the filename simply 
    # Then I need to find the chormocenters that were segmented for this cell type	
    for temp_files in sortedCells:
        
        if sortedMarker in temp_files:
            identifier = temp_files.split(sortedMarker)[1].split(".")[0].strip()
            identifier = identifier.replace("T0", "T30")  # Update T0 to T30
            # This is needed when I have chromocenter files with (Chromo)_ in front of each filename
            # I remove it and just extract the filename  
            identifiers.append(identifier)
    # Now I iterate over the chromocenter files I have received
    # Here I would always be comparing to my markers so I have to remove again the
    # Markers like chormo, fop, foxJ etc
    # Then I check if the filname match to read only those files in sequence
    for identifier in identifiers:
        # Find matching files in the source folder
        for source_file in sourceCells:
            if sourceMarker and sourceMarker in source_file:
                file_name = source_file.split(sourceMarker)[1].split(".")[0].strip()
                if identifier == file_name:
                    matching_files.append(source_file)

    return matching_files


# Sorting Files

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

# Main Folder
folder_Path = "E:/Quantified/eNUC Live/230120 OF1 D0 Enuc Ctrl5_actin5/"
###########################
#   Actin    Reading      #
###########################
# Time T0
directory = folder_Path + 'T0/Result/Actin/'

ChromoFilesT0 = folder_scan(directory, ".xlsx", marker="(Chromo)_")
NucleusFilesT0 = folder_scan(directory, ".xlsx", marker="(Nucleus)_")

ChromoFilesT0 = sorted(ChromoFilesT0)
NucleusFilesT0 = sorted(NucleusFilesT0)

# Time T30
directory = folder_Path + 'T30/Result/Actin/'

ChromoFilesT1 = folder_scan(directory, ".xlsx", marker="(Chromo)_")
NucleusFilesT1 = folder_scan(directory, ".xlsx", marker="(Nucleus)_")
ActinFilesT1 = folder_scan(directory, ".xlsx", marker="(Actin)_")

ChromoFilesT1 = sorted(ChromoFilesT1)
NucleusFilesT1 = sorted(NucleusFilesT1)
ActinFilesT1 = sorted(ActinFilesT1)

# Actin T0
for nucleusread in NucleusFilesT0:
    df = pd.read_excel(nucleusread)
    file_name = nucleusread.split("(Nucleus)_")[1].split(".")[0].strip()
    stats_Actin["File_Name_Actin"].append(file_name)
    stats_Actin["Nuclei_volume_actin_T0"].append(df.loc[0,'Nucleus Area'])
    stats_Actin["Nuclei_intensity_actin_T0"].append(df.iloc[0, df.columns.get_loc('mean_intensity')])
    
# Actin T30             
for nucleusread in NucleusFilesT1:
    file_name = nucleusread.split("(Nucleus)_")[1].split(".")[0].strip()
    file_name = file_name.replace("T30", "T0")
    
    if file_name in stats_Actin["File_Name_Actin"]:
        df = pd.read_excel(nucleusread)
        stats_Actin["Nuclei_volume_actin_T1"].append(df.loc[0,'Nucleus Area'])
        stats_Actin["Nuclei_intensity_actin_T1"].append(df.iloc[0, df.columns.get_loc('mean_intensity')])
    # else:
    #     stats_Actin["File_Name_Actin"].append(file_name)
    #     stats_Actin["Nuclei_volume_actin_T1"].append(None)
    #     stats_Actin["Nuclei_intensity_actin_T1"].append(None)
             
# Chromocenter Read T0
for chromo_read in ChromoFilesT0:
    df = pd.read_excel(chromo_read)        
    stats_Actin["Actin_chromo_volume_T0"].append(df.loc[0, 'Chromocenter Area'])
    stats_Actin["Actin_chromo_intensity_T0"].append(df.iloc[0, df.columns.get_loc('mean_intensity')])     
# Chromocenter Read T30
for chromo_read in ChromoFilesT1:
    file_name = chromo_read.split("(Chromo)_")[1].split(".")[0].strip()
    file_name = file_name.replace("T30", "T0") 
    if file_name in stats_Actin["File_Name_Actin"]:
        df = pd.read_excel(chromo_read)        
        stats_Actin["Actin_chromo_volume_T1"].append(df.loc[0, 'Chromocenter Area'])
        stats_Actin["Actin_chromo_intensity_T1"].append(df.iloc[0, df.columns.get_loc('mean_intensity')])
    # else:
    #     stats_Actin["Actin_chromo_volume_T1"].append(None)
    #     stats_Actin["Actin_chromo_intensity_T1"].append(None)   

# Actin Read T30
for actinread in ActinFilesT1:
    file_name = actinread.split("(Actin)_")[1].split(".")[0].strip()
    file_name = file_name.replace("T30", "T0") 
    if file_name in stats_Actin["File_Name_Actin"]:
        df = pd.read_excel(actinread)        
        stats_Actin["Actin_coverage_T30"].append(df.loc[0,'Actin Coverage']) 
    # else:        
    #     stats_Actin["Actin_coverage_T30"].append(None) 
        
df_Actin = pd.DataFrame.from_dict(stats_Actin, orient='index')
df_Actin = df_Actin.transpose()         

###########################
#    Ctrl    Reading      #
###########################

directory = folder_Path + 'T0/Result/Ctrl/'

ChromoFilesT0 = folder_scan(directory, ".xlsx", marker="(Chromo)_")
NucleusFilesT0 = folder_scan(directory, ".xlsx", marker="(Nucleus)_")

directory = folder_Path + 'T30/Result/Ctrl/'

ChromoFilesT1 = folder_scan(directory, ".xlsx", marker="(Chromo)_")
NucleusFilesT1 = folder_scan(directory, ".xlsx", marker="(Nucleus)_")

# Ctrl Nucleus Read T0
for nucleusread in NucleusFilesT0:
    df = pd.read_excel(nucleusread)
    file_name = nucleusread.split("(Nucleus)_")[1].split(".")[0].strip()
    stats_Ctrl["File_Name_Ctrl"].append(file_name)
    stats_Ctrl["Nuclei_volume_ctrl_T0"].append(df.loc[0,'Nucleus Area'])
    stats_Ctrl["Nuclei_intensity_ctrl_T0"].append(df.iloc[0, df.columns.get_loc('mean_intensity')])

# Ctrl Nucleus Read T30        
for nucleusread in NucleusFilesT1:
    file_name = nucleusread.split("(Nucleus)_")[1].split(".")[0].strip()
    file_name = file_name.replace("T30", "T0")
    if file_name in stats_Ctrl["File_Name_Ctrl"]:    
        df = pd.read_excel(nucleusread)
        #file_name = nucleusread.split("(Nucleus)_")[1].split(".")[0].strip()
        #stats_Actin["File_Name_Actin"].append(file_name)
        stats_Ctrl["Nuclei_volume_ctrl_T1"].append(df.loc[0,'Nucleus Area'])
        stats_Ctrl["Nuclei_intensity_ctrl_T1"].append(df.iloc[0, df.columns.get_loc('mean_intensity')])
        # else:
           
# Chromocenter Read T0
for temp_read in ChromoFilesT0:
    df = pd.read_excel(temp_read)        
    stats_Ctrl["Ctrl_chromo_volume_T0"].append(df.loc[0, 'Chromocenter Area'])
    stats_Ctrl["Ctrl_chromo_intensity_T0"].append(df.iloc[0, df.columns.get_loc('mean_intensity')])     

# Chromocenter Read T30
for chromo_read in ChromoFilesT1:
    file_name = chromo_read.split("(Chromo)_")[1].split(".")[0].strip()
    file_name = file_name.replace("T30", "T0")   
    if file_name in stats_Ctrl["File_Name_Ctrl"]:    
        df = pd.read_excel(chromo_read)        
        stats_Ctrl["Ctrl_chromo_volume_T1"].append(df.loc[0, 'Chromocenter Area'])
        stats_Ctrl["Ctrl_chromo_intensity_T1"].append(df.iloc[0, df.columns.get_loc('mean_intensity')])             
     
df_Ctrl = pd.DataFrame.from_dict(stats_Ctrl, orient='index')
df_Ctrl = df_Ctrl.transpose()  
# =============================================================================
# # Filter out the incorrectly Segmented Chromo on T0
# filterChromo = df_Actin['Actin_chromo_volume_T1'] > 25
# df_ActinF = df_Actin.drop(df_Actin.index[filterChromo])
# =============================================================================


# =============================================================================
# filterChromo = df_Ctrl['Ctrl_chromo_volume_T0'] > 25
# df_CtrlF = df_Ctrl.drop(df_Ctrl.index[filterChromo])
# =============================================================================
