# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 11:10:42 2023

@author: kaabi
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
        
# def matchingFiles(sortedCells, sourceCells,sortedMarker=None, sourceMarker=None):
#     # Here I have filenames which are segmented labels of individual nuclei marker of hoescht
#     # I want to sort the cell type (RGC, multi, flower, halo) images with their nucleus and choromo markers
#     # Note: not all nucleus with the segmentation would have chromocenter 
#     # To sort these and match correct files this function matches the files
#     # Example Use - ChromoFiles =  matchingFiles(rgc_Scan, chromo_Scan,sortedMarker=None, sourceMarker="(Chromo)_")
#     # Example Use -NucleusFiles =  matchingFiles(ChromoFiles, nucleus_Scan, sortedMarker="(Chromo)_", sourceMarker="(Nucleus)_")
    
#     identifiers = [] # For my sorted/analyzed cell type/ less in number

#     matching_files = []
    
#     # First I collect files from my cell types that are Radial glial or multi ciliated etc.
#     # In the first case the files are tif files and not markers
#     # I can extract the filename simply 
#     # Then I need to find the chormocenters that were segmented for this cell type
#     if sortedMarker is None:
#         for temp_files in sortedCells:
#             file_name = os.path.basename(temp_files)
#             identifier = os.path.splitext(file_name)[0]
#             identifiers.append(identifier)
#     # This is needed when I have chromocenter files with (Chromo)_ in front of each filename
#     # I remove it and just extract the filename        
#     else:
#         for temp_files in sortedCells:
#             file_name = temp_files.split(sortedMarker)[1].split(".")[0].strip() #
#             identifier = file_name
#             identifiers.append(identifier)
#     # Now I iterate over the chromocenter files I have received
#     # Here I would always be comparing to my markers so I have to remove again the
#     # Markers like chormo, fop, foxJ etc
#     # Then I check if the filname match to read only those files in sequence        
#     for identifier in identifiers:
#         #matching_files = []
#         # Find matching files in the source folder
#         for source_file in sourceCells:
#             file_name = source_file.split(sourceMarker)[1].split(".")[0].strip() #

#             if identifier == file_name:
#                 matching_files.append(source_file)
#                 #print(matching_files)                
#     return matching_files   

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

###########################
#           Agarose       #
###########################
# Scan Folder for all relevant files
chromo_Scan = folder_scan(chromo_Folder, ".xlsx", marker="(Chromo)_")
nucleus_Scan = folder_scan(nucleus_Folder, ".xlsx", marker="(Nucleus)_")
h3k27me3_Scan = folder_scan(h3k27me3_Folder, ".xlsx", marker="(H3K27me3)_")
foxj_Scan = folder_scan(foxj_Folder, ".xlsx", marker="(FOXJ)_")
fop_Scan = folder_scan(fop_Folder, ".xlsx", marker="(FOP)_")

###########################
#       RGC    Reading    #
###########################
def exportStatsAgarose(folder=None, label = None):
    
	# Find matching files
	NucleusFiles =  matchingFiles(folder, nucleus_Scan,sortedMarker=None, sourceMarker="(Nucleus)_")

	ChromoFiles =  matchingFiles(NucleusFiles, chromo_Scan, sortedMarker="(Nucleus)_", sourceMarker="(Chromo)_")
	H3k27me3Files = matchingFiles(NucleusFiles, h3k27me3_Scan, sortedMarker="(Nucleus)_", sourceMarker="(H3K27me3)_")
	FoxjFiles = matchingFiles(NucleusFiles, foxj_Scan, sortedMarker="(Nucleus)_", sourceMarker="(FOXJ)_")
	FopFiles = matchingFiles(NucleusFiles, fop_Scan, sortedMarker="(Nucleus)_", sourceMarker="(FOP)_")

	# Dictionary to store the properties
	# region_Prop_ctrl = {"ctrl_intens_H3k27me3":[],"ctrl_intens_FOP":[],"ctrl_intens_Nucleus":[],
	#                     "ctrl_intens_Chromo":[]}
	region_Prop_agar = {"agar_File_Name":[],"agar_intens_H3k27me3":[],"agar_intens_FOP":[],
						"agar_volume_Chromo":[],"agar_intens_mean_Chromo":[], "agar_intens_median_Chromo":[],
						"agar_intens_mean_Nucleus":[], "agar_intens_median_Nucleus":[] ,"R":[]}

	###########################
	#  Multi_Folder  Reading #
	###########################

	for temp_read in NucleusFiles:
		df = pd.read_excel(temp_read)
		file_name = temp_read.split("(Nucleus)_")[1].split(".")[0].strip()
		region_Prop_agar["agar_File_Name"].append(file_name)    
		region_Prop_agar["agar_intens_mean_Nucleus"].append(df.loc[0,'mean'])
		region_Prop_agar["agar_intens_median_Nucleus"].append(df.loc[0,'median'])

	h3k27me3_data= {}

	for temp_read in H3k27me3Files:
		file_name = temp_read.split("(H3K27me3)_")[1].split(".")[0].strip()

		if file_name in region_Prop_agar["agar_File_Name"]:
			df = pd.read_excel(temp_read)
			# Store the value in the temporary dictionary
			h3k27me3_data[file_name] = df.iloc[0, df.columns.get_loc('mean')]
		else:
			# If the file_name is not found in stats_Actin["File_Name_Actin"], set it to 'None'
			h3k27me3_data[file_name] = 'None'

	# After processing all files, append the values to the 'stats_Actin' DataFrame
	region_Prop_agar["agar_intens_H3k27me3"] = [h3k27me3_data.get(file_name, 'None') for file_name in region_Prop_agar["agar_File_Name"]]        
	 
	fopfiles_data= {}

	for temp_read in FopFiles:
		file_name = temp_read.split("(FOP)_")[1].split(".")[0].strip()

		if file_name in region_Prop_agar["agar_File_Name"]:
			df = pd.read_excel(temp_read)
			# Store the value in the temporary dictionary
			fopfiles_data[file_name] = df.iloc[0, df.columns.get_loc('mean')]
		else:
			# If the file_name is not found in stats_Actin["File_Name_Actin"], set it to 'None'
			fopfiles_data[file_name] = 'None'

	# After processing all files, append the values to the 'stats_Actin' DataFrame
	region_Prop_agar["agar_intens_FOP"] = [fopfiles_data.get(file_name, 'None') for file_name in region_Prop_agar["agar_File_Name"]]        

	# Chromocenter 
	chromo_data = {}

	for temp_read in ChromoFiles:
		file_name = temp_read.split("(Chromo)_")[1].split(".")[0].strip()

		if file_name in region_Prop_agar["agar_File_Name"]:
			df = pd.read_excel(temp_read)
			noValues = df['mean'].count()

			# Check if there are more than 1 values in 'Chromocenter Area' column
			if noValues > 1:
				# Store the mean intensity and Chromocenter Area values for the current file in the dictionary as a tuple
				chromo_data[file_name] = (df['mean'].mean(), df['median'].mean(), df.loc[0, 'Chromocenter Area'])

			else:
				# If there is only one value or no value, set both values in the tuple to the single value in the DataFrame
				mean_intensity_value = df['mean'].iloc[0] if noValues == 1 else 'None'
				median_intensity_value = df['median'].iloc[0] if noValues == 1 else 'None'
				chromocenter_area_value = df['Chromocenter Area'].iloc[0] if noValues == 1 else 'None'
				chromo_data[file_name] = (mean_intensity_value, median_intensity_value, chromocenter_area_value)

	region_Prop_agar["agar_intens_mean_Chromo"] = [chromo_data.get(file_name, ('None', 'None', 'None'))[0] for file_name in region_Prop_agar["agar_File_Name"]]
	region_Prop_agar["agar_intens_median_Chromo"] = [chromo_data.get(file_name, ('None', 'None', 'None'))[1] for file_name in region_Prop_agar["agar_File_Name"]]
	region_Prop_agar["agar_volume_Chromo"] = [chromo_data.get(file_name, ('None', 'None', 'None'))[2] for file_name in region_Prop_agar["agar_File_Name"]]

	df_agar = pd.DataFrame.from_dict(region_Prop_agar, orient='index')
	df_agar = df_agar.transpose()

	# Convert the 'agar_intens_median_Chromo' and 'agar_intens_median_Nucleus' columns to numeric type
	df_agar['agar_intens_median_Chromo'] = pd.to_numeric(df_agar['agar_intens_median_Chromo'], errors='coerce')
	df_agar['agar_intens_median_Nucleus'] = pd.to_numeric(df_agar['agar_intens_median_Nucleus'], errors='coerce')

	# Create a new column 'agar_intens_mean_Chromo' and calculate the mean of 'agar_intens_median_Chromo' column
	df_agar['agar_intens_mean_Chromo'] = df_agar['agar_intens_median_Chromo'].mean()

	# Perform element-wise division and create new columns for 'agar_intens_median_Chromo' and 'agar_volume_Chromo'
	df_agar['R'] = df_agar['agar_intens_median_Chromo'] / df_agar['agar_intens_median_Nucleus']

	pd.DataFrame(df_agar).to_excel(folder_Path + label + '_Agarose_Excel.xlsx') 

	region_Prop_agar.clear()

# Scan Cell Type Folders
rgc_Scan = folder_scan(rgc_Folder,".tif")
exportStatsAgarose(folder=rgc_Scan, label = '/RGC')
###########################
#  Halo_Folder    Reading #
###########################
# Scan Cell Type Folders
halo_Scan = folder_scan(halo_Folder,".tif")

exportStatsAgarose(folder=halo_Scan, label = '/Halo')
###########################
#  Flower_Folder  Reading #
###########################
# Scan Cell Type Folders
flower_Scan = folder_scan(flower_Folder,".tif")

exportStatsAgarose(folder=flower_Scan, label = '/Flower')
###########################
#  Multi_Folder  Reading #
###########################
# Scan Cell Type Folders
multi_Scan = folder_scan(multi_Folder,".tif")

exportStatsAgarose(folder=multi_Scan, label = '/Multi')

# Main Folder
folder_Path = "E:/Quantified/Agarose Compression/20221103 IF OF1 D1 H3K9me3 FOP FOXJ1/Result/Overnight/Ctrl"
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

###########################
#           Ctrl       #
###########################
# Scan Folder for all relevant files
chromo_Scan = folder_scan(chromo_Folder, ".xlsx", marker="(Chromo)_")
nucleus_Scan = folder_scan(nucleus_Folder, ".xlsx", marker="(Nucleus)_")
h3k27me3_Scan = folder_scan(h3k27me3_Folder, ".xlsx", marker="(H3K27me3)_")
foxj_Scan = folder_scan(foxj_Folder, ".xlsx", marker="(FOXJ)_")
fop_Scan = folder_scan(fop_Folder, ".xlsx", marker="(FOP)_")

def exportStatsCtrl(folder=None, label = None) :
    
# Find matching files
	NucleusFiles =  matchingFiles(folder, nucleus_Scan,sortedMarker=None, sourceMarker="(Nucleus)_")

	ChromoFiles =  matchingFiles(NucleusFiles, chromo_Scan, sortedMarker="(Nucleus)_", sourceMarker="(Chromo)_")
	H3k27me3Files = matchingFiles(NucleusFiles, h3k27me3_Scan, sortedMarker="(Nucleus)_", sourceMarker="(H3K27me3)_")
	FoxjFiles = matchingFiles(NucleusFiles, foxj_Scan, sortedMarker="(Nucleus)_", sourceMarker="(FOXJ)_")
	FopFiles = matchingFiles(NucleusFiles, fop_Scan, sortedMarker="(Nucleus)_", sourceMarker="(FOP)_")

	# Dictionary to store the properties
	region_Prop_ctrl = {"ctrl_File_Name":[],"ctrl_intens_H3k27me3":[],"ctrl_intens_FOP":[],
						"ctrl_volume_Chromo":[],"ctrl_intens_mean_Chromo":[], "ctrl_intens_median_Chromo":[],
						"ctrl_intens_mean_Nucleus":[], "ctrl_intens_median_Nucleus":[] ,"R":[]}

	###########################
	#  Multi_Folder  Reading #
	###########################

	for temp_read in NucleusFiles:
		df = pd.read_excel(temp_read)
		file_name = temp_read.split("(Nucleus)_")[1].split(".")[0].strip()
		region_Prop_ctrl["ctrl_File_Name"].append(file_name)    
		region_Prop_ctrl["ctrl_intens_mean_Nucleus"].append(df.loc[0,'mean'])
		region_Prop_ctrl["ctrl_intens_median_Nucleus"].append(df.loc[0,'median'])

	h3k27me3_data= {}

	for temp_read in H3k27me3Files:
		file_name = temp_read.split("(H3K27me3)_")[1].split(".")[0].strip()

		if file_name in region_Prop_ctrl["ctrl_File_Name"]:
			df = pd.read_excel(temp_read)
			# Store the value in the temporary dictionary
			h3k27me3_data[file_name] = df.iloc[0, df.columns.get_loc('mean')]
		else:
			# If the file_name is not found in stats_Actin["File_Name_Actin"], set it to 'None'
			h3k27me3_data[file_name] = 'None'

	# After processing all files, append the values to the 'stats_Actin' DataFrame
	region_Prop_ctrl["ctrl_intens_H3k27me3"] = [h3k27me3_data.get(file_name, 'None') for file_name in region_Prop_ctrl["ctrl_File_Name"]]        
	 
	fopfiles_data= {}

	for temp_read in FopFiles:
		file_name = temp_read.split("(FOP)_")[1].split(".")[0].strip()

		if file_name in region_Prop_ctrl["ctrl_File_Name"]:
			df = pd.read_excel(temp_read)
			# Store the value in the temporary dictionary
			fopfiles_data[file_name] = df.iloc[0, df.columns.get_loc('mean')]
		else:
			# If the file_name is not found in stats_Actin["File_Name_Actin"], set it to 'None'
			fopfiles_data[file_name] = 'None'

	# After processing all files, append the values to the 'stats_Actin' DataFrame
	region_Prop_ctrl["ctrl_intens_FOP"] = [fopfiles_data.get(file_name, 'None') for file_name in region_Prop_ctrl["ctrl_File_Name"]]        

	# Chromocenter 
	# Chromocenter 
	chromo_data = {}

	for temp_read in ChromoFiles:
		file_name = temp_read.split("(Chromo)_")[1].split(".")[0].strip()

		if file_name in region_Prop_ctrl["ctrl_File_Name"]:
			df = pd.read_excel(temp_read)
			noValues = df['mean'].count()

			# Check if there are more than 1 values in 'Chromocenter Area' column
			if noValues > 1:
				# Store the mean intensity and Chromocenter Area values for the current file in the dictionary as a tuple
				chromo_data[file_name] = (df['mean'].mean(), df['median'].mean(), df.loc[0, 'Chromocenter Area'])

			else:
				# If there is only one value or no value, set both values in the tuple to the single value in the DataFrame
				mean_intensity_value = df['mean'].iloc[0] if noValues == 1 else 'None'
				median_intensity_value = df['median'].iloc[0] if noValues == 1 else 'None'
				chromocenter_area_value = df['Chromocenter Area'].iloc[0] if noValues == 1 else 'None'
				chromo_data[file_name] = (mean_intensity_value, median_intensity_value, chromocenter_area_value)

	region_Prop_ctrl["ctrl_intens_mean_Chromo"] = [chromo_data.get(file_name, ('None', 'None', 'None'))[0] for file_name in region_Prop_ctrl["ctrl_File_Name"]]
	region_Prop_ctrl["ctrl_intens_median_Chromo"] = [chromo_data.get(file_name, ('None', 'None', 'None'))[1] for file_name in region_Prop_ctrl["ctrl_File_Name"]]
	region_Prop_ctrl["ctrl_volume_Chromo"] = [chromo_data.get(file_name, ('None', 'None', 'None'))[2] for file_name in region_Prop_ctrl["ctrl_File_Name"]]
		   

	df_ctrl = pd.DataFrame.from_dict(region_Prop_ctrl, orient='index')
	df_ctrl = df_ctrl.transpose()

	# Convert the 'ctrl_intens_median_Chromo' and 'ctrl_intens_median_Nucleus' columns to numeric type
	df_ctrl['ctrl_intens_median_Chromo'] = pd.to_numeric(df_ctrl['ctrl_intens_median_Chromo'], errors='coerce')
	df_ctrl['ctrl_intens_median_Nucleus'] = pd.to_numeric(df_ctrl['ctrl_intens_median_Nucleus'], errors='coerce')

	# Create a new column 'ctrl_intens_mean_Chromo' and calculate the mean of 'ctrl_intens_median_Chromo' column
	df_ctrl['ctrl_intens_mean_Chromo'] = df_ctrl['ctrl_intens_median_Chromo'].mean()

	# Perform element-wise division and create new columns for 'ctrl_intens_median_Chromo' and 'ctrl_volume_Chromo'
	df_ctrl['R'] = df_ctrl['ctrl_intens_median_Chromo'] / df_ctrl['ctrl_intens_median_Nucleus']

	pd.DataFrame(df_ctrl).to_excel(folder_Path + label + '_Ctrl_Excel.xlsx') 

	region_Prop_ctrl.clear()

###########################
#       RGC    Reading    #
###########################
# Scan Cell Type Folders
rgc_Scan = folder_scan(rgc_Folder,".tif")

exportStatsCtrl(folder=rgc_Scan, label = '/RGC')

###########################
#  Halo_Folder    Reading #
###########################
# Scan Cell Type Folders
halo_Scan = folder_scan(halo_Folder,".tif")

exportStatsCtrl(folder=halo_Scan, label = '/Halo')

###########################
#  Flower_Folder  Reading #
###########################
# Scan Cell Type Folders
flower_Scan = folder_scan(flower_Folder,".tif")

exportStatsCtrl(folder=flower_Scan, label = '/Flower')
###########################
#  Multi_Folder  Reading #
###########################
# Scan Cell Type Folders
multi_Scan = folder_scan(multi_Folder,".tif")

exportStatsCtrl(folder=multi_Scan, label = '/Multi')
