# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 11:10:42 2023

@author: kaabir
"""

import os
import pandas as pd
import numpy as np


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
        
def matchingFiles(sortedCells, sourceCells, sortedMarker=None, sourceMarker=None):
    identifiers = set()
    matching_files_dict = {}

    # First, collect identifiers from sortedCells
    if sortedMarker is None:
        identifiers = {os.path.splitext(os.path.basename(file))[0] for file in sortedCells}
    else:
        identifiers = {file.split(sortedMarker)[1].split(".")[0].strip() for file in sortedCells}

    # Next, create a dictionary to store the matching files for each identifier
    for identifier in identifiers:
        matching_files_dict[identifier] = []

    # Now, iterate over the sourceCells and find matching files for each identifier
    for source_file in sourceCells:
        file_name = os.path.basename(source_file)
        if sourceMarker is not None:
            file_name = file_name.split(sourceMarker)[1].split(".")[0].strip()

        if file_name in matching_files_dict:
            matching_files_dict[file_name].append(source_file)

    # Finally, flatten the matching_files_dict into a list of matching files
    matching_files = [file for files in matching_files_dict.values() for file in files]

    return matching_files

def process_string(input_string, marker = None):
    # Find the last occurrence of "(Chromo)"or "(Actin)" or "(NUC04)"
    # 
    # print(input_string)
    directory = os.path.splitext(input_string)[0]
    path = os.path.basename(directory )  
    identifier = path.split(marker)[1]
    # print(identifier)
    
    return identifier

def ProcessData(folder_path = None, FolderType= None):
	Imaris_Type = {}

	# Cell Marker Folder
	chromo_Folder = folder_Path + "/Chromo/"
	nucleus_Folder = folder_Path + "/Nucleus/"
	h3k27me3_Folder = folder_Path + "/H3K27me3/"
	foxj_Folder = folder_Path + "/FoxJ/"
	fop_Folder = folder_Path + "/FOP/"

	# Scan Folder for all relevant files
	chromo_Scan = folder_scan(chromo_Folder, ".xlsx", marker="(Chromo)_")
	nucleus_Scan = folder_scan(nucleus_Folder, ".xlsx", marker="(Nucleus)_")
	h3k27me3_Scan = folder_scan(h3k27me3_Folder, ".xlsx", marker="(H3K27me3)_")
	foxj_Scan = folder_scan(foxj_Folder, ".xlsx", marker="(FOXJ)_")
	fop_Scan = folder_scan(fop_Folder, ".xlsx", marker="(FOP)_")

	# Scan Cell Type Folders
	rgc_Scan = folder_scan(FolderType,".tif")

	# Find matching files
	ChromoFiles =  matchingFiles(rgc_Scan, chromo_Scan,sortedMarker=None, sourceMarker="(Chromo)_")

	NucleusFiles =  matchingFiles(ChromoFiles, nucleus_Scan, sortedMarker="(Chromo)_", sourceMarker="(Nucleus)_")
	H3k27me3Files = matchingFiles(ChromoFiles, h3k27me3_Scan, sortedMarker="(Chromo)_", sourceMarker="(H3K27me3)_")
	FoxjFiles = matchingFiles(ChromoFiles, foxj_Scan, sortedMarker="(Chromo)_", sourceMarker="(FOXJ)_")
	FopFiles = matchingFiles(ChromoFiles, fop_Scan, sortedMarker="(Chromo)_", sourceMarker="(FOP)_")

	for nucleiread in NucleusFiles:
		file_name = process_string(nucleiread, marker="(Nucleus)")
		if file_name not in Imaris_Type:
			Imaris_Type[file_name] = {}
		CT_nuclei_TO = pd.read_excel(nucleiread)
			#print('Nuclei Actin', file_name)
			# Handle normal files
		Imaris_Type[file_name]['Nucleus Mean Intensity'] = CT_nuclei_TO.loc[0,'mean']
		Imaris_Type[file_name]['Nucleus Median Intensity'] = CT_nuclei_TO.loc[0,'median']
		
	# Process chromo Imaris_Type
	for chromoread in ChromoFiles:
		file_name = process_string(chromoread, marker="(Chromo)")
		if file_name not in Imaris_Type:
			Imaris_Type[file_name] = {}
		CT_chromo_TO = pd.read_excel(chromoread)		
			#print('Chromo Actin', file_name)
			# Handle normal files
		Imaris_Type[file_name]['Chromo Mean Intensity'] = CT_chromo_TO.loc[0,'mean']
		Imaris_Type[file_name]['Chromo Median Intensity'] = CT_chromo_TO.loc[0,'median']

	# Process H3K27me3 Imaris_Type
	for h3K27me3read in H3k27me3Files:
		file_name = process_string(h3K27me3read, marker="(H3K27me3)")
		if file_name not in Imaris_Type:
			Imaris_Type[file_name] = {}
		CT_h3K27me3_TO = pd.read_excel(h3K27me3read)

		Imaris_Type[file_name]['H3K27me3 Mean Intensity'] = CT_h3K27me3_TO.loc[0,'mean']

	# Process FOXJ Imaris_Type
	for foxjread in FoxjFiles:
		file_name = process_string(foxjread, marker="(FOXJ)")
		if file_name not in Imaris_Type:
			Imaris_Type[file_name] = {}
		CT_foxj_TO = pd.read_excel(foxjread)
			# Handle normal files
		Imaris_Type[file_name]['FOXJ Intensity'] = CT_foxj_TO.loc[0,'mean']

	for fopread in FopFiles:
		file_name = process_string(fopread, marker="(FOP)")
		if file_name not in Imaris_Type:
			Imaris_Type[file_name] = {}
		CT_fop_TO = pd.read_excel(fopread)
			# Handle normal files
		Imaris_Type[file_name]['FOP Intensity'] = CT_fop_TO.loc[0,'mean']

            
	df = pd.DataFrame.from_dict(Imaris_Type, orient='index')
	df["R"] = df["Chromo Median Intensity"]/df["Nucleus Median Intensity"]
	return df

###########################
#           Agarose       #
###########################

# Main Folder
folder_Path = "E:/Quantified/Agarose Compression/20221103 IF OF1 D1 H3K9me3 FOP FOXJ1/Result/Overnight/Agarose"
Result_folder = 'E:/Quantified/Agarose Compression/20221103 IF OF1 D1 H3K9me3 FOP FOXJ1/Result/'

rgc_Folder = folder_Path + "/01-RGC/"
halo_Folder = folder_Path + "/02-Halo/"
flower_Folder = folder_Path + "/03-Flower/"
multi_Folder = folder_Path + "/04-Multi/"

rgc_Reading = ProcessData(folder_Path,FolderType = rgc_Folder)

rgc_Reading.to_excel(Result_folder + '/Export_RGC_Agarose_Excel.xlsx')  
halo_Reading = ProcessData(folder_Path,FolderType = halo_Folder)

halo_Reading.to_excel(Result_folder + '/Export_Halo_Agarose_Excel.xlsx')  
flower_Reading = ProcessData(folder_Path,FolderType = flower_Folder)

flower_Reading.to_excel(Result_folder + '/Export_Flower_Agarose_Excel.xlsx')  
multi_Reading = ProcessData(folder_Path,FolderType = multi_Folder)

multi_Reading.to_excel(Result_folder + '/Export_Multi_Agarose_Excel.xlsx')  

###########################
#           Ctrl          #
###########################
folder_Path = "E:/Quantified/Agarose Compression/20221103 IF OF1 D1 H3K9me3 FOP FOXJ1/Result/Overnight/Ctrl"
rgc_Folder = folder_Path + "/01-RGC/"
halo_Folder = folder_Path + "/02-Halo/"
flower_Folder = folder_Path + "/03-Flower/"
multi_Folder = folder_Path + "/04-Multi/"

rgc_Reading = ProcessData(folder_Path,FolderType = rgc_Folder)
rgc_Reading.to_excel(Result_folder + '/Export_RGC_Ctrl_Excel.xlsx')  
halo_Reading = ProcessData(folder_Path,FolderType = halo_Folder)
halo_Reading.to_excel(Result_folder + '/Export_Halo_Ctrl_Excel.xlsx')  
flower_Reading = ProcessData(folder_Path,FolderType = flower_Folder)
flower_Reading.to_excel(Result_folder + '/Export_Flower_Ctrl_Excel.xlsx')  
multi_Reading = ProcessData(folder_Path,FolderType = multi_Folder)
multi_Reading.to_excel(Result_folder + '/Export_Multi_Ctrl_Excel.xlsx')
