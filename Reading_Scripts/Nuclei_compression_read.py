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
	
	###########################
    #       RGC    Reading    #
    ###########################

	# Cell Marker Folder
	chromo_Folder = folder_Path + "/Chromo/"
	nucleus_Folder = folder_Path + "/Nucleus/"
	phosRB_Folder = folder_Path + "/PhosRB/"
	foxj_Folder = folder_Path + "/FoxJ/"
	fop_Folder = folder_Path + "/FOP/"

	###########################
	#           Agarose       #
	###########################
	# Scan Folder for all relevant files
	chromo_Scan = folder_scan(chromo_Folder, ".xlsx", marker="(Chromo)_")
	nucleus_Scan = folder_scan(nucleus_Folder, ".xlsx", marker="(Nucleus)_")
	phosRB_Scan = folder_scan(phosRB_Folder, ".xlsx", marker="(PhosRB)_")
	foxj_Scan = folder_scan(foxj_Folder, ".xlsx", marker="(FOXJ)_")
	fop_Scan = folder_scan(fop_Folder, ".xlsx", marker="(FOP)_")

	###########################
	#       RGC    Reading    #
	###########################
	# Scan Cell Type Folders
	rgc_Scan = folder_scan(FolderType,".tif")

	# Find matching files
	ChromoFiles =  matchingFiles(rgc_Scan, chromo_Scan,sortedMarker=None, sourceMarker="(Chromo)_")

	NucleusFiles =  matchingFiles(ChromoFiles, nucleus_Scan, sortedMarker="(Chromo)_", sourceMarker="(Nucleus)_")
	PhosRBFiles = matchingFiles(ChromoFiles, phosRB_Scan, sortedMarker="(Chromo)_", sourceMarker="(PhosRB)_")
	FoxjFiles = matchingFiles(ChromoFiles, foxj_Scan, sortedMarker="(Chromo)_", sourceMarker="(FOXJ)_")
	FopFiles = matchingFiles(ChromoFiles, fop_Scan, sortedMarker="(Chromo)_", sourceMarker="(FOP)_")

	for nucleiread in NucleusFiles:
		file_name = process_string(nucleiread, marker="(Nucleus)_")
		if file_name not in Imaris_Type:
			Imaris_Type[file_name] = {}
		CT_nuclei_TO = pd.read_excel(nucleiread)
			#print('Nuclei Actin', file_name)
			# Handle normal files
		Imaris_Type[file_name]['Nucleus Mean Intensity'] = CT_nuclei_TO.loc[0,'mean']
		Imaris_Type[file_name]['Nucleus Median Intensity'] = CT_nuclei_TO.loc[0,'median']
		Imaris_Type[file_name]['Nucleus Area'] = CT_nuclei_TO.loc[0,'Nucleus Area']
        
	# Process chromo Imaris_Type
	for chromoread in ChromoFiles:
		file_name = process_string(chromoread, marker="(Chromo)_")
		if file_name not in Imaris_Type:
			Imaris_Type[file_name] = {}
		CT_chromo_TO = pd.read_excel(chromoread)		
			#print('Chromo Actin', file_name)
			# Handle normal files
		Imaris_Type[file_name]['Chromo Mean Intensity'] = CT_chromo_TO.loc[0,'mean']
		Imaris_Type[file_name]['Chromo Median Intensity'] = CT_chromo_TO.loc[0,'median']
		Imaris_Type[file_name]['Chromocenter Area'] = CT_chromo_TO.loc[0,'Chromocenter Area']
		Imaris_Type[file_name]['Chromocenter Spots'] = CT_chromo_TO.loc[0,'Chromocenter Spots']
        
	# Process H3K27me3 Imaris_Type
	for phosrbread in PhosRBFiles:
		file_name = process_string(phosrbread, marker="(PhosRB)_")
		if file_name not in Imaris_Type:
			Imaris_Type[file_name] = {}
		CT_PhosRB_TO = pd.read_excel(phosrbread)

		Imaris_Type[file_name]['PhosRB Mean Intensity'] = CT_PhosRB_TO.loc[0,'mean']

	# Process FOXJ Imaris_Type
	for foxjread in FoxjFiles:
		file_name = process_string(foxjread, marker="(FOXJ)_")
		if file_name not in Imaris_Type:
			Imaris_Type[file_name] = {}
		CT_foxj_TO = pd.read_excel(foxjread)
			# Handle normal files
		Imaris_Type[file_name]['FOXJ Intensity'] = CT_foxj_TO.loc[0,'mean']

	for fopread in FopFiles:
		file_name = process_string(fopread, marker="(FOP)_")
		if file_name not in Imaris_Type:
			Imaris_Type[file_name] = {}
		CT_fop_TO = pd.read_excel(fopread)
			# Handle normal files
		Imaris_Type[file_name]['FOP Intensity'] = CT_fop_TO.loc[0,'mean']

    # dimensionless_ratio = (number_of_bright_spots * bright_spot_intensity) / (circle_area * circle_intensity)
	df = pd.DataFrame.from_dict(Imaris_Type, orient='index')
	# df["R"] = (df["Chromocenter Spots"] * df["Chromo Mean Intensity"])/(df["Nucleus Area"] * df["Nucleus Mean Intensity"])
	df["R Mean"] = (df["Chromo Mean Intensity"]/ df["Nucleus Mean Intensity"])
	df["R Median"] = (df["Chromo Median Intensity"]/ df["Nucleus Median Intensity"])
	# df["Chromocenter Occupancy"] = (df["Chromocenter Area"] * df["Chromocenter Spots"]/ df["Nucleus Area"]) * 100
	df["Chromocenter Occupancy"] = (df["Chromocenter Area"] * df["Chromocenter Spots"])
	
	return df

###########################
#           Agarose       #
###########################
#Result_folder = 'E:/Quantified/Agarose Compression/PhosphoRb Experiment/04-231220 OF1 D6  Pad2perc3R Gprb Rdsred FRfop/Result/Pad/'


# folder_Path = "E:/Quantified/Agarose Compression/20221103 IF OF1 D1 H3K9me3 FOP FOXJ1/Test/Result_Ctrl"

folder_Path = "E:/Quantified/Agarose Compression/PhosphoRb Experiment/04-231220 OF1 D6  Pad2perc3R Gprb Rdsred FRfop/Result/Pad/"

if not os.path.exists(folder_Path + '/Result'):
    os.makedirs(folder_Path + '/Result')
        
    Result_folder = folder_Path + '/Result/' 

# Main Folder

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
folder_Path = "E:/Quantified/Agarose Compression/PhosphoRb Experiment/04-231220 OF1 D6  Pad2perc3R Gprb Rdsred FRfop/Result/Ctrl"

if not os.path.exists(folder_Path + '/Result'):
    os.makedirs(folder_Path + '/Result')
        
    Result_folder = folder_Path + '/Result/' 

rgc_Folder = folder_Path + "/01-RGC/"
halo_Folder = folder_Path + "/02-Halo/"
flower_Folder = folder_Path + "/03-Flower/"
multi_Folder = folder_Path + "/04-Multi/"

rgc_Reading = ProcessData(folder_Path,FolderType = rgc_Folder)
rgc_Reading.to_excel(Result_folder + '/Export_RGC_Ctrl_Excel.xlsx' , index=True)  
halo_Reading = ProcessData(folder_Path,FolderType = halo_Folder)
halo_Reading.to_excel(Result_folder + '/Export_Halo_Ctrl_Excel.xlsx', index=True)  
flower_Reading = ProcessData(folder_Path,FolderType = flower_Folder)
flower_Reading.to_excel(Result_folder + '/Export_Flower_Ctrl_Excel.xlsx', index=True)  
multi_Reading = ProcessData(folder_Path,FolderType = multi_Folder)
multi_Reading.to_excel(Result_folder + '/Export_Multi_Ctrl_Excel.xlsx', index=True)
