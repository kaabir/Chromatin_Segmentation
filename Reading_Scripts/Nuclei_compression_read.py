# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 11:10:42 2023

@author: kaabi
"""

import os
import pandas as pd
import numpy as np


# Main Folder
folder_Path = "E:/Quantified/Agarose Compression/20221103 IF OF1 D1 H3K9me3 FOP FOXJ1/Result/Overnight/Agarose"
# Cell Type Folder
rgc_Folder = folder_Path + "/01-RGC/"
halo_Folder = folder_Path + "/02-Halo/"
flower_Folder = folder_Path + "/03-Flower/"
multi_Folder = folder_Path + "/04-Multi/"
# Cell Marker Folder
chromo_Folder = folder_Path + "/Chromo/"
nucleus_Folder = folder_Path + "/Nucleus/"
h3k27me3_Folder = folder_Path + "/H3K27me3/"
foxj_Folder = folder_Path + "/FoxJ/"
fop_Folder = folder_Path + "/FOP/"

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
    # Here I have filenames which are segmented labels of individual nuclei marker of hoescht
    # I want to sort the cell type (RGC, multi, flower, halo) images with their nucleus and choromo markers
    # Note: not all nucleus with the segmentation would have chromocenter 
    # To sort these and match correct files this function matches the files
    # Example Use - ChromoFiles =  matchingFiles(rgc_Scan, chromo_Scan,sortedMarker=None, sourceMarker="(Chromo)_")
    # Example Use -NucleusFiles =  matchingFiles(ChromoFiles, nucleus_Scan, sortedMarker="(Chromo)_", sourceMarker="(Nucleus)_")
    
    identifiers = [] # For my sorted/analyzed cell type/ less in number

    matching_files = []
    
    # First I collect files from my cell types that are Radial glial or multi ciliated etc.
    # In the first case the files are tif files and not markers
    # I can extract the filename simply 
    # Then I need to find the chormocenters that were segmented for this cell type
    if sortedMarker is None:
        for temp_files in sortedCells:
            file_name = os.path.basename(temp_files)
            identifier = os.path.splitext(file_name)[0]
            identifiers.append(identifier)
    # This is needed when I have chromocenter files with (Chromo)_ in front of each filename
    # I remove it and just extract the filename        
    else:
        for temp_files in sortedCells:
            file_name = temp_files.split(sortedMarker)[1].split(".")[0].strip() #
            identifier = file_name
            identifiers.append(identifier)
    # Now I iterate over the chromocenter files I have received
    # Here I would always be comparing to my markers so I have to remove again the
    # Markers like chormo, fop, foxJ etc
    # Then I check if the filname match to read only those files in sequence        
    for identifier in identifiers:
        #matching_files = []
        # Find matching files in the source folder
        for source_file in sourceCells:
            file_name = source_file.split(sourceMarker)[1].split(".")[0].strip() #

            if identifier == file_name:
                matching_files.append(source_file)
                #print(matching_files)                
    return matching_files   

# Scan Folder for all relevant files
chromo_Scan = folder_scan(chromo_Folder, ".xlsx", marker="(Chromo)_")
nucleus_Scan = folder_scan(nucleus_Folder, ".xlsx", marker="(Nucleus)_")
h3k27me3_Scan = folder_scan(h3k27me3_Folder, ".xlsx", marker="(H3K27me3)_")
foxj_Scan = folder_scan(foxj_Folder, ".xlsx", marker="(FOXJ)_")
fop_Scan = folder_scan(fop_Folder, ".xlsx", marker="(FOP)_")

###########################
#       RGC    Reading    #
###########################
# Scan Cell Type Folders
rgc_Scan = folder_scan(rgc_Folder,".tif")

# Find matching files
ChromoFiles =  matchingFiles(rgc_Scan, chromo_Scan,sortedMarker=None, sourceMarker="(Chromo)_")

NucleusFiles =  matchingFiles(ChromoFiles, nucleus_Scan, sortedMarker="(Chromo)_", sourceMarker="(Nucleus)_")
H3k27me3Files = matchingFiles(ChromoFiles, h3k27me3_Scan, sortedMarker="(Chromo)_", sourceMarker="(H3K27me3)_")
FoxjFiles = matchingFiles(ChromoFiles, foxj_Scan, sortedMarker="(Chromo)_", sourceMarker="(FOXJ)_")
FopFiles = matchingFiles(ChromoFiles, fop_Scan, sortedMarker="(Chromo)_", sourceMarker="(FOP)_")

# Dictionary to store the properties
# region_Prop_ctrl = {"ctrl_intens_H3k27me3":[],"ctrl_intens_FOP":[],"ctrl_intens_Nucleus":[],
#                     "ctrl_intens_Chromo":[]}
region_Prop_agar = {"agar_intens_H3k27me3":[],"agar_intens_FOP":[],"agar_intens_Nucleus":[],
                    "agar_intens_Chromo":[]}

###########################
#      RGC    Reading    #
###########################

for temp_read in NucleusFiles:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_Nucleus"].append(df.loc[0,'mean'])
    
for temp_read in H3k27me3Files:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_H3k27me3"].append(df.loc[0,'mean'])    
    
for temp_read in FopFiles:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_FOP"].append(df.loc[0,'mean'])    
    
for temp_read in ChromoFiles:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_Chromo"].append(df.loc[0,'mean']) 

df_agar = pd.DataFrame.from_dict(region_Prop_agar, orient='index')
df_agar = df_agar.transpose()
pd.DataFrame(df_agar).to_excel(folder_Path + '/RGC_Agar_Excel.xlsx')

region_Prop_agar.clear()

###########################
#  Halo_Folder    Reading #
###########################
# Scan Cell Type Folders
halo_Scan = folder_scan(halo_Folder,".tif")

# Find matching files
ChromoFiles =  matchingFiles(halo_Scan, chromo_Scan,sortedMarker=None, sourceMarker="(Chromo)_")

NucleusFiles =  matchingFiles(ChromoFiles, nucleus_Scan, sortedMarker="(Chromo)_", sourceMarker="(Nucleus)_")
H3k27me3Files = matchingFiles(ChromoFiles, h3k27me3_Scan, sortedMarker="(Chromo)_", sourceMarker="(H3K27me3)_")
FoxjFiles = matchingFiles(ChromoFiles, foxj_Scan, sortedMarker="(Chromo)_", sourceMarker="(FOXJ)_")
FopFiles = matchingFiles(ChromoFiles, fop_Scan, sortedMarker="(Chromo)_", sourceMarker="(FOP)_")

# Dictionary to store the properties
# region_Prop_ctrl = {"ctrl_intens_H3k27me3":[],"ctrl_intens_FOP":[],"ctrl_intens_Nucleus":[],
#                     "ctrl_intens_Chromo":[]}
region_Prop_agar = {"agar_intens_H3k27me3":[],"agar_intens_FOP":[],"agar_intens_Nucleus":[],
                    "agar_intens_Chromo":[]}

###########################
#  Halo   Reading #
###########################

for temp_read in NucleusFiles:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_Nucleus"].append(df.loc[0,'mean'])
    
for temp_read in H3k27me3Files:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_H3k27me3"].append(df.loc[0,'mean'])    
    
for temp_read in FopFiles:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_FOP"].append(df.loc[0,'mean'])    
    
for temp_read in ChromoFiles:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_Chromo"].append(df.loc[0,'mean']) 

df_agar = pd.DataFrame.from_dict(region_Prop_agar, orient='index')
df_agar = df_agar.transpose()
    
pd.DataFrame(df_agar).to_excel(folder_Path + '/Halo_Agarose_Excel.xlsx') 

region_Prop_agar.clear()

###########################
#  Flower_Folder  Reading #
###########################
# Scan Cell Type Folders
flower_Scan = folder_scan(flower_Folder,".tif")

# Find matching files
ChromoFiles =  matchingFiles(flower_Scan, chromo_Scan,sortedMarker=None, sourceMarker="(Chromo)_")

NucleusFiles =  matchingFiles(ChromoFiles, nucleus_Scan, sortedMarker="(Chromo)_", sourceMarker="(Nucleus)_")
H3k27me3Files = matchingFiles(ChromoFiles, h3k27me3_Scan, sortedMarker="(Chromo)_", sourceMarker="(H3K27me3)_")
FoxjFiles = matchingFiles(ChromoFiles, foxj_Scan, sortedMarker="(Chromo)_", sourceMarker="(FOXJ)_")
FopFiles = matchingFiles(ChromoFiles, fop_Scan, sortedMarker="(Chromo)_", sourceMarker="(FOP)_")

# Dictionary to store the properties
# region_Prop_ctrl = {"ctrl_intens_H3k27me3":[],"ctrl_intens_FOP":[],"ctrl_intens_Nucleus":[],
#                     "ctrl_intens_Chromo":[]}
region_Prop_agar = {"agar_intens_H3k27me3":[],"agar_intens_FOP":[],"agar_intens_Nucleus":[],
                    "agar_intens_Chromo":[]}

###########################
#  Flower    Reading #
###########################

for temp_read in NucleusFiles:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_Nucleus"].append(df.loc[0,'mean'])
    
for temp_read in H3k27me3Files:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_H3k27me3"].append(df.loc[0,'mean'])    
    
for temp_read in FopFiles:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_FOP"].append(df.loc[0,'mean'])    
    
for temp_read in ChromoFiles:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_Chromo"].append(df.loc[0,'mean']) 

df_agar = pd.DataFrame.from_dict(region_Prop_agar, orient='index')
df_agar = df_agar.transpose()
    
pd.DataFrame(df_agar).to_excel(folder_Path + '/Flower_Agarose_Excel.xlsx') 

region_Prop_agar.clear()

###########################
#  Multi_Folder  Reading #
###########################
# Scan Cell Type Folders
multi_Scan = folder_scan(multi_Folder,".tif")

# Find matching files
ChromoFiles =  matchingFiles(multi_Scan, chromo_Scan,sortedMarker=None, sourceMarker="(Chromo)_")

NucleusFiles =  matchingFiles(ChromoFiles, nucleus_Scan, sortedMarker="(Chromo)_", sourceMarker="(Nucleus)_")
H3k27me3Files = matchingFiles(ChromoFiles, h3k27me3_Scan, sortedMarker="(Chromo)_", sourceMarker="(H3K27me3)_")
FoxjFiles = matchingFiles(ChromoFiles, foxj_Scan, sortedMarker="(Chromo)_", sourceMarker="(FOXJ)_")
FopFiles = matchingFiles(ChromoFiles, fop_Scan, sortedMarker="(Chromo)_", sourceMarker="(FOP)_")

# Dictionary to store the properties
# region_Prop_ctrl = {"ctrl_intens_H3k27me3":[],"ctrl_intens_FOP":[],"ctrl_intens_Nucleus":[],
#                     "ctrl_intens_Chromo":[]}
region_Prop_agar = {"agar_intens_H3k27me3":[],"agar_intens_FOP":[],"agar_intens_Nucleus":[],
                    "agar_intens_Chromo":[]}

###########################
#  Multi_Folder  Reading #
###########################

for temp_read in NucleusFiles:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_Nucleus"].append(df.loc[0,'mean'])
    
for temp_read in H3k27me3Files:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_H3k27me3"].append(df.loc[0,'mean'])    
    
for temp_read in FopFiles:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_FOP"].append(df.loc[0,'mean'])    
    
for temp_read in ChromoFiles:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_Chromo"].append(df.loc[0,'mean']) 

df_agar = pd.DataFrame.from_dict(region_Prop_agar, orient='index')
df_agar = df_agar.transpose()
    
pd.DataFrame(df_agar).to_excel(folder_Path + '/Multi_Agarose_Excel.xlsx') 

region_Prop_agar.clear()
