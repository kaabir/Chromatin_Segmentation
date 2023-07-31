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
import re

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

stats_Actin = {"File_Number_Actin":[], "File_Name_Actin":[],"Actin_chromo_volume_T0":[],"Actin_chromo_volume_T1":[],"Actin_chromo_volume_T2":[],
                      "Actin_chromo_intensity_T0":[],"Actin_chromo_intensity_T1":[],"Actin_chromo_intensity_T2":[],
          "Nuclei_volume_actin_T0":[],"H3K27me3_Intensity_actin_T0":[],"H3K9me2_Intensity_actin_T0":[],"Nuclei_volume_actin_T1":[],
          "Nuclei_volume_actin_T2":[],
		  "Nuclei_intensity_actin_T0":[],"Nuclei_intensity_actin_T1":[],"Nuclei_intensity_actin_T2":[],
		  "Actin_coverage_T30":[],"Actin_coverage_T60":[],"sphericity_Actin_T0":[],"sphericity_Actin_T1":[],"sphericity_Actin_T2":[]}

stats_Ctrl = {"File_Number_Ctrl":[],"File_Name_Ctrl":[],"Ctrl_chromo_volume_T0":[],"Ctrl_chromo_volume_T1":[],"Ctrl_chromo_volume_T2":[],
           "Ctrl_chromo_intensity_T0":[],"Ctrl_chromo_intensity_T1":[],"Ctrl_chromo_intensity_T2":[],
           "Nuclei_volume_ctrl_T0":[],"H3K27me3_Intensity_ctrl_T0":[],"H3K9me2_Intensity_ctrl_T0":[],"Nuclei_volume_ctrl_T1":[],"Nuclei_volume_ctrl_T2":[],
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

def folder_scan(directory, extension, marker=None, Key=str()):
    # folder_scan(chromo_folder, ".xlsx", marker="(Chromo)")
    get_files = []
    extension = extension  # ".xlsx"  # ".czi"
    final_dir = directory + Key
    
    try:
        for f_name in os.listdir(final_dir):
            if marker is not None:
                if f_name.find(marker) != -1 and f_name.endswith(extension):
                    get_files.append(os.path.join(final_dir, f_name))
            else:
                if f_name.endswith(extension):
                    get_files.append(os.path.join(final_dir, f_name))
        return get_files              

    except FileNotFoundError:
        print(">>> This Folder Not Present", Key)
        print(final_dir)

# Main Folder
folder_Path = "E:/Quantified/eNUC Staining/230406 OF1 D1 Actin 10uMVCA Gactin RH3K9me2 FRH3K27me2/Result/"
###########################
#   Actin    Reading      #
###########################
# Time T0
directory = folder_Path + 'Actin/'

ChromoFiles = folder_scan(directory, ".xlsx", marker="(Chromo)_", Key='Chromo/')
NucleusFiles = folder_scan(directory, ".xlsx", marker="(Nucleus)_", Key='Nucleus/')
H3K27me3Files = folder_scan(directory, ".xlsx", marker="(H3K27me3)_", Key='H3K27me3/')
H3K9me2Files = folder_scan(directory, ".xlsx", marker="(H3K9me2)_", Key='H3K9me2/')
ActinFiles = folder_scan(directory, ".xlsx", marker="(Actin)_", Key='Actin/')

for nucleusread in NucleusFiles:
    df = pd.read_excel(nucleusread)
    file_name = nucleusread.split("(Nucleus)_")[1].split(".")[0].strip()
    stats_Actin["File_Name_Actin"].append(file_name)      
    stats_Actin["Nuclei_volume_actin_T0"].append(df.loc[0,'Nucleus Area'])
    stats_Actin["Nuclei_intensity_actin_T0"].append(df.iloc[0, df.columns.get_loc('mean')])

# # # Just saving File Number for Plots  
# # # .group(0) allows us to access and extract the matched substring directly as a string.  
# # for digitread in stats_Actin["File_Name_Actin"]:   
# #     match = re.search(r' #\d+_\d+', digitread)
# #     if match:
# #         number = match.group(0)
# #         stats_Actin["File_Number_Actin"].append(number)  

# # .group(0) allows us to access and extract the matched substring directly as a string.
# stats_Actin["File_Number_Actin"] = [re.search(r' #\d+_\d+', digitread).group(0) for digitread in stats_Actin["File_Name_Actin"] if re.search(r' #\d+_\d+', digitread)]        

h3k27me3_data= {}

for h3k27me3read in H3K27me3Files:
    file_name = h3k27me3read.split("(H3K27me3)_")[1].split(".")[0].strip()

    if file_name in stats_Actin["File_Name_Actin"]:
        df = pd.read_excel(h3k27me3read)
        # Store the 'Actin Coverage' value for the current file in the dictionary
        h3k27me3_data[file_name] = df.iloc[0, df.columns.get_loc('mean_intensity')]
        # print(file_name, df.loc[0, 'Actin Coverage'])
    else:
        # If the file_name is not found in stats_Actin["File_Name_Actin"], set it to 'None'
        h3k27me3_data[file_name] = 'None'

# After processing all files, append the values to the 'stats_Actin' DataFrame
stats_Actin["H3K27me3_Intensity_actin_T0"] = [h3k27me3_data.get(file_name, 'None') for file_name in stats_Actin["File_Name_Actin"]]        

h3k9me2_data= {}

for h3k9me2read in H3K9me2Files:
    file_name = h3k9me2read.split("(H3K9me2)_")[1].split(".")[0].strip()

    if file_name in stats_Actin["File_Name_Actin"]:
        df = pd.read_excel(h3k9me2read)
        # Store the 'Actin Coverage' value for the current file in the dictionary
        h3k9me2_data[file_name] = df.iloc[0, df.columns.get_loc('mean_intensity')]
        # print(file_name, df.loc[0, 'Actin Coverage'])
    else:
        # If the file_name is not found in stats_Actin["File_Name_Actin"], set it to 'None'
        h3k9me2_data[file_name] = 'None'

# After processing all files, append the values to the 'stats_Actin' DataFrame
stats_Actin["H3K9me2_Intensity_actin_T0"] = [h3k27me3_data.get(file_name, 'None') for file_name in stats_Actin["File_Name_Actin"]]        


# Chromocenter 
chromo_data = {}

for chromo_read in ChromoFiles:
    file_name = chromo_read.split("(Chromo)_")[1].split(".")[0].strip()

    if file_name in stats_Actin["File_Name_Actin"]:
        df = pd.read_excel(chromo_read)
        noValues = df['mean_intensity'].count()

        # Check if there are more than 1 values in 'Chromocenter Area' column
        if noValues > 1:
            # Store the mean intensity and Chromocenter Area values for the current file in the dictionary as a tuple
            chromo_data[file_name] = (df['mean_intensity'].mean(), df.loc[0, 'Chromocenter Area'])
        else:
            # If there is only one value or no value, set both values in the tuple to the single value in the DataFrame
            mean_intensity_value = df['mean_intensity'].iloc[0] if noValues == 1 else 'None'
            chromocenter_area_value = df['Chromocenter Area'].iloc[0] if noValues == 1 else 'None'
            chromo_data[file_name] = (mean_intensity_value, chromocenter_area_value)
    else:
        # If the file_name is not found in the 'stats_Actin' DataFrame, set both values in the tuple to 'None'
        chromo_data[file_name] = ('None', 'None')

# After processing all files, append the values to the 'stats_Actin' DataFrame
stats_Actin["Actin_chromo_intensity_T0"] = [chromo_data.get(file_name, ('None', 'None'))[0] for file_name in stats_Actin["File_Name_Actin"]]
stats_Actin["Actin_chromo_volume_T0"] = [chromo_data.get(file_name, ('None', 'None'))[1] for file_name in stats_Actin["File_Name_Actin"]]

# Initialize a dictionary to store the data for each file
actin_data = {}

for actinread in ActinFiles:
    file_name = actinread.split("(Actin)_")[1].split(".")[0].strip()

    if file_name in stats_Actin["File_Name_Actin"]:
        df = pd.read_excel(actinread)
        # Store the 'Actin Coverage' value for the current file in the dictionary
        actin_data[file_name] = df.loc[0, 'Actin Coverage']
        # print(file_name, df.loc[0, 'Actin Coverage'])
    else:
        # If the file_name is not found in stats_Actin["File_Name_Actin"], set it to 'None'
        actin_data[file_name] = 'None'

# After processing all files, append the values to the 'stats_Actin' DataFrame
stats_Actin["Actin_coverage_T30"] = [actin_data.get(file_name, 'None') for file_name in stats_Actin["File_Name_Actin"]]        

df_Actin = pd.DataFrame.from_dict(stats_Actin, orient='index')
df_Actin = df_Actin.transpose()    

# filterChromo = (df_Actin['Nuclei_volume_actin_T0'] < 120)
# df_Actin.loc[filterChromo, ['Nuclei_volume_actin_T0']] = 0

filterChromo = df_Actin['Nuclei_volume_actin_T0'] < 120
df_Actin = df_Actin.drop(df_Actin.index[filterChromo])

pd.DataFrame(df_Actin).to_excel(folder_Path + '/Export_Actin_Excel.xlsx')  

# ###########################
# #    Ctrl    Reading      #
# ###########################

directory = folder_Path + 'Ctrl/'

ChromoFiles = folder_scan(directory, ".xlsx", marker="(Chromo)_", Key='Chromo/')
NucleusFiles = folder_scan(directory, ".xlsx", marker="(Nucleus)_", Key='Nucleus/')
H3K27me3Files = folder_scan(directory, ".xlsx", marker="(H3K27me3)_", Key='H3K27me3/')
H3K9me2Files = folder_scan(directory, ".xlsx", marker="(H3K9me2)_", Key='H3K9me2/')

# Ctrl Nucleus Read T0
for nucleusread in NucleusFiles:
    df = pd.read_excel(nucleusread)
    file_name = nucleusread.split("(Nucleus)_")[1].split(".")[0].strip()
    stats_Ctrl["File_Name_Ctrl"].append(file_name)
    stats_Ctrl["Nuclei_volume_ctrl_T0"].append(df.loc[0,'Nucleus Area'])
    stats_Ctrl["Nuclei_intensity_ctrl_T0"].append(df.iloc[0, df.columns.get_loc('mean')])

# # Just saving File Number for Plots    
# for digitread in stats_Ctrl["File_Name_Ctrl"]:   
#     match = re.search(r' #\d+_\d+', digitread)
#     if match:
#         number = match.group(0)
#         stats_Ctrl["File_Number_Ctrl"].append(number) 
# Chromocenter 
chromo_data = {}

for chromo_read in ChromoFiles:
    file_name = chromo_read.split("(Chromo)_")[1].split(".")[0].strip()

    if file_name in stats_Ctrl["File_Name_Ctrl"]:
        df = pd.read_excel(chromo_read)
        noValues = df['mean_intensity'].count()

        # Check if there are more than 1 values in 'Chromocenter Area' column
        if noValues > 1:
            # Store the mean intensity and Chromocenter Area values for the current file in the dictionary as a tuple
            chromo_data[file_name] = (df['mean_intensity'].mean(), df.loc[0, 'Chromocenter Area'])
        else:
            # If there is only one value or no value, set both values in the tuple to the single value in the DataFrame
            mean_intensity_value = df['mean_intensity'].iloc[0] if noValues == 1 else 'None'
            chromocenter_area_value = df['Chromocenter Area'].iloc[0] if noValues == 1 else 'None'
            chromo_data[file_name] = (mean_intensity_value, chromocenter_area_value)
    else:
        # If the file_name is not found in the 'stats_Actin' DataFrame, set both values in the tuple to 'None'
        chromo_data[file_name] = ('None', 'None')

# After processing all files, append the values to the 'stats_Actin' DataFrame
stats_Actin["Ctrl_chromo_intensity_T0"] = [chromo_data.get(file_name, ('None', 'None'))[0] for file_name in stats_Ctrl["File_Name_Ctrl"]]
stats_Actin["Ctrl_chromo_volume_T0"] = [chromo_data.get(file_name, ('None', 'None'))[1] for file_name in stats_Ctrl["File_Name_Ctrl"]]

# for chromo_read in ChromoFiles:
#     df = pd.read_excel(chromo_read)
#     noValues = df['mean_intensity'].count()
#     # Check if there are more than 1 values in 'Chromocenter Area' column
#     if noValues > 1:
#         # Calculate the average of 'Chromocenter Area' values and append it to the list
#         stats_Ctrl["Ctrl_chromo_volume_T0"].append(df.loc[0, 'Chromocenter Area'])
#         stats_Ctrl["Ctrl_chromo_intensity_T0"].append(df['mean_intensity'].mean())
#     else:
#         # If there is only one value, directly append it to the list         
#         stats_Ctrl["Ctrl_chromo_volume_T0"].append(df.loc[0, 'Chromocenter Area'])
#         stats_Ctrl["Ctrl_chromo_intensity_T0"].append(df.iloc[0, df.columns.get_loc('mean_intensity')])

h3k27me3_data= {}

for h3k27me3read in H3K27me3Files:
    file_name = h3k27me3read.split("(H3K27me3)_")[1].split(".")[0].strip()

    if file_name in stats_Ctrl["File_Name_Ctrl"]:
        df = pd.read_excel(h3k27me3read)
        # Store the 'Actin Coverage' value for the current file in the dictionary
        h3k27me3_data[file_name] = df.iloc[0, df.columns.get_loc('mean_intensity')]
        # print(file_name, df.loc[0, 'Actin Coverage'])
    else:
        # If the file_name is not found in stats_Actin["File_Name_Actin"], set it to 'None'
        h3k27me3_data[file_name] = 'None'

# After processing all files, append the values to the 'stats_Actin' DataFrame
stats_Ctrl["H3K27me3_Intensity_ctrl_T0"] = [h3k27me3_data.get(file_name, 'None') for file_name in stats_Ctrl["File_Name_Ctrl"]]        

h3k9me2_data= {}

for h3k9me2read in H3K9me2Files:
    file_name = h3k9me2read.split("(H3K9me2)_")[1].split(".")[0].strip()

    if file_name in stats_Ctrl["File_Name_Ctrl"]:
        df = pd.read_excel(h3k9me2read)
        # Store the 'Actin Coverage' value for the current file in the dictionary
        h3k9me2_data[file_name] = df.iloc[0, df.columns.get_loc('mean_intensity')]
        # print(file_name, df.loc[0, 'Actin Coverage'])
    else:
        # If the file_name is not found in stats_Actin["File_Name_Actin"], set it to 'None'
        h3k9me2_data[file_name] = 'None'

# After processing all files, append the values to the 'stats_Actin' DataFrame
stats_Ctrl["H3K9me2_Intensity_ctrl_T0"] = [h3k27me3_data.get(file_name, 'None') for file_name in stats_Ctrl["File_Name_Ctrl"]]        

df_Ctrl = pd.DataFrame.from_dict(stats_Ctrl, orient='index')
df_Ctrl = df_Ctrl.transpose()  

# filterChromo = (df_Ctrl['Nuclei_volume_ctrl_T0'] < 120)
# #df_Ctrl.loc[filterChromo, ['Nuclei_volume_ctrl_T0']] = 0
# df_Ctrl.loc[filterChromo, ['Nuclei_volume_ctrl_T0']] = 0

filterChromo = df_Ctrl['Nuclei_volume_ctrl_T0'] < 120
df_Ctrl = df_Ctrl.drop(df_Ctrl.index[filterChromo])

pd.DataFrame(df_Ctrl).to_excel(folder_Path + '/Export_Ctrl_Excel.xlsx')   
# =============================================================================
# # Filter out the incorrectly Segmented Chromo on T0
# filterChromo = df_Actin['Actin_chromo_volume_T1'] > 25
# df_ActinF = df_Actin.drop(df_Actin.index[filterChromo])
# =============================================================================


# =============================================================================
# filterChromo = df_Ctrl['Ctrl_chromo_volume_T0'] > 25
# df_CtrlF = df_Ctrl.drop(df_Ctrl.index[filterChromo])
# =============================================================================
