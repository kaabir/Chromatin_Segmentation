# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 11:10:42 2023

@author: kaabir
"""

import os
import pandas as pd
import numpy as np


# Dictionary to store the properties
region_Prop_ctrl = {"ctrl_intens_H3k27me3":[],"ctrl_intens_FOP":[],"ctrl_intens_Nucleus":[],
                    "ctrl_intens_Chromo":[]}
region_Prop_agar = {"agar_intens_H3k27me3":[],"agar_intens_FOP":[],"agar_intens_Nucleus":[],
                    "agar_intens_Chromo":[]}

# Main Folder
folder_Path = "E:/Quantified/Agarose Compression/20221103 IF OF1 D1 H3K9me3 FOP FOXJ1/Result/1H/"


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
    
    for temp_files in sortedCells:
        if sortedMarker is None:
            file_name = os.path.basename(temp_files)
            file_name = os.path.splitext(file_name)[0]
    # This is needed when I have chromocenter files with (Chromo)_ in front of each filename
    # I remove it and just extract the filename        
        else:
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

###########################
#   Agarose    Reading    #
###########################

# Cell Type Folder
rgc_Folder_A = folder_Path + "Agarose/01-RGC/"
halo_Folder_A = folder_Path + "Agarose/02-Halo/"
flower_Folder_A = folder_Path + "Agarose/03-Flower/"
multi_Folder_A = folder_Path + "Agarose/04-Multi/"
# Cell Marker Folder
chromo_Folder_A = folder_Path + "Agarose/Chromo/"
nucleus_Folder_A = folder_Path + "Agarose/Nucleus/"
h3k27me3_Folder_A = folder_Path + "Agarose/H3K27me3/"
foxj_Folder_A = folder_Path + "Agarose/FoxJ/"
fop_Folder_A = folder_Path + "Agarose/FOP/"

###########################
#   Agarose    Reading    #
###########################

# Scan Folder for all relevant files
chromo_Scan_A = folder_scan(chromo_Folder_A, ".xlsx", marker="(Chromo)_")
nucleus_Scan_A = folder_scan(nucleus_Folder_A, ".xlsx", marker="(Nucleus)_")
h3k27me3_Scan_A = folder_scan(h3k27me3_Folder_A, ".xlsx", marker="(H3K27me3)_")
foxj_Scan_A = folder_scan(foxj_Folder_A, ".xlsx", marker="(FOXJ)_")
fop_Scan_A = folder_scan(fop_Folder_A, ".xlsx", marker="(FOP)_")

# Scan Cell Type Folders
rgc_Scan_A = folder_scan(rgc_Folder_A,".tif")

# Find matching files
ChromoFiles_A =  matchingFiles(rgc_Scan_A, chromo_Scan_A,sortedMarker=None, sourceMarker="(Chromo)_")

NucleusFiles_A =  matchingFiles(ChromoFiles_A, nucleus_Scan_A, sortedMarker="(Chromo)_", sourceMarker="(Nucleus)_")
H3k27me3Files_A = matchingFiles(ChromoFiles_A, h3k27me3_Scan_A, sortedMarker="(Chromo)_", sourceMarker="(H3K27me3)_")
FoxjFiles_A = matchingFiles(ChromoFiles_A, foxj_Scan_A, sortedMarker="(Chromo)_", sourceMarker="(FOXJ)_")
FopFiles_A = matchingFiles(ChromoFiles_A, fop_Scan_A, sortedMarker="(Chromo)_", sourceMarker="(FOP)_")


###########################
#   Agarose    Reading    #
###########################

for temp_read in NucleusFiles_A:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_Nucleus"].append(df.loc[0,'mean'])
    
for temp_read in H3k27me3Files_A:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_H3k27me3"].append(df.loc[0,'mean'])    
    
for temp_read in FopFiles_A:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_FOP"].append(df.loc[0,'mean'])    
    
for temp_read in ChromoFiles_A:
    df = pd.read_excel(temp_read)
    region_Prop_agar["agar_intens_Chromo"].append(df.loc[0,'mean']) 

df_agar = pd.DataFrame.from_dict(region_Prop_agar, orient='index')
df_agar = df_agar.transpose()
df_agar['Ratio_chromo_nucle'] = df_agar['agar_intens_Chromo']/df_agar['agar_intens_Nucleus']

pd.DataFrame(df_agar).to_excel(folder_Path + '/Export_Agar_Excel.xlsx')

###########################
#      Ctrl    Reading    #
###########################

# Cell Type Folder
rgc_Folder = folder_Path + "Ctrl/01-RGC/"
halo_Folder = folder_Path + "Ctrl/02-Halo/"
flower_Folder = folder_Path + "Ctrl/03-Flower/"
multi_Folder = folder_Path + "Ctrl/04-Multi/"
# Cell Marker Folder
chromo_Folder = folder_Path + "Ctrl/Chromo/"
nucleus_Folder = folder_Path + "Ctrl/Nucleus/"
h3k27me3_Folder = folder_Path + "Ctrl/H3K27me3/"
foxj_Folder = folder_Path + "Ctrl/FoxJ/"
fop_Folder = folder_Path + "Ctrl/FOP/"

# Scan Folder for all relevant files
chromo_Scan = folder_scan(chromo_Folder, ".xlsx", marker="(Chromo)_")
nucleus_Scan = folder_scan(nucleus_Folder, ".xlsx", marker="(Nucleus)_")
h3k27me3_Scan = folder_scan(h3k27me3_Folder, ".xlsx", marker="(H3K27me3)_")
foxj_Scan = folder_scan(foxj_Folder, ".xlsx", marker="(FOXJ)_")
fop_Scan = folder_scan(fop_Folder, ".xlsx", marker="(FOP)_")

# Scan Cell Type Folders
rgc_Scan = folder_scan(rgc_Folder,".tif")

# Find matching files
ChromoFiles =  matchingFiles(rgc_Scan, chromo_Scan,sortedMarker=None, sourceMarker="(Chromo)_")

NucleusFiles =  matchingFiles(ChromoFiles, nucleus_Scan, sortedMarker="(Chromo)_", sourceMarker="(Nucleus)_")
H3k27me3Files = matchingFiles(ChromoFiles, h3k27me3_Scan, sortedMarker="(Chromo)_", sourceMarker="(H3K27me3)_")
FoxjFiles = matchingFiles(ChromoFiles, foxj_Scan, sortedMarker="(Chromo)_", sourceMarker="(FOXJ)_")
FopFiles = matchingFiles(ChromoFiles, fop_Scan, sortedMarker="(Chromo)_", sourceMarker="(FOP)_")

###########################
#      Ctrl    Reading    #
###########################

for temp_read in NucleusFiles:
    df = pd.read_excel(temp_read)
    region_Prop_ctrl["ctrl_intens_Nucleus"].append(df.loc[0,'mean'])
    
for temp_read in H3k27me3Files:
    df = pd.read_excel(temp_read)
    region_Prop_ctrl["ctrl_intens_H3k27me3"].append(df.loc[0,'mean'])    
    
for temp_read in FopFiles:
    df = pd.read_excel(temp_read)
    region_Prop_ctrl["ctrl_intens_FOP"].append(df.loc[0,'mean'])    
    
for temp_read in ChromoFiles:
    df = pd.read_excel(temp_read)
    region_Prop_ctrl["ctrl_intens_Chromo"].append(df.loc[0,'mean']) 

df_ctrl = pd.DataFrame.from_dict(region_Prop_ctrl, orient='index')
df_ctrl = df_ctrl.transpose()
df_ctrl['Ratio_chromo_nucle'] = df_ctrl['ctrl_intens_Chromo']/df_ctrl['ctrl_intens_Nucleus']

pd.DataFrame(df_ctrl).to_excel(folder_Path + '/Export_Ctrl_Excel.xlsx')
