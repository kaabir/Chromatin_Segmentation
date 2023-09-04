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
import numpy as np
import re

def folder_Scan(folder_path= None, Type = None):
    # Actin after  30min  - value assigned is 3
    
    # Sorting Files
    actin_files =[]
    h3K9me2_files =[]
    h3K27me3_files =[]
    nuclei_files = [] 
    chromo_files = [] 
    
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
                    elif f_name.startswith('(NUC04)'):
                        nuclei_files.append(f_name)
                    elif f_name.startswith('(Actin)'):
                        actin_files.append(f_name)
                    elif f_name.startswith('(Chromo)'):
                        chromo_files.append(f_name)                        
                    else:
                        print('Some random file or folders in the directory',f_name)
                        
    return actin_files, h3K9me2_files, h3K27me3_files, nuclei_files, chromo_files
                                         
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
    # Find the last occurrence of "(Chromo)"or "(Actin)" or "(NUC04)"
    # 
    directory = os.path.splitext(input_string)[0]
    path = os.path.basename(directory )  
    identifier = path.split(marker)[1]
    
    return identifier

def ProcessData(FolderType= None, folder_path = None):
	Imaris_Type = {}
	nuclei_channel = 1
	
	actin_files, h3K9me2_files, h3K27me3_files, nuclei_files, chromo_files = folder_Scan(folder_path, FolderType)
	
	for nucleiread in nuclei_files:
		file_name = process_string(nucleiread, marker="(NUC04) ")
		if file_name not in Imaris_Type:
			Imaris_Type[file_name] = {}
		CT_nuclei_TO = read_CSV(nucleiread)
		
		# Use regular expression to find the pattern of "number_number" in the whole string
		duplicate_string_match = re.search(r'(\d+)_(\d+)', file_name)
		
		if duplicate_string_match:
			# Extract both parts of the pattern
			file_number_first_part = int(duplicate_string_match.group(1))
			file_number_second_part = int(duplicate_string_match.group(2))
			#print('Nuclei Actin Split', file_name)
			if file_number_first_part >= 1:
				row_index = file_number_second_part #+ nuclei_channel  # Adjust the row index in each iteration
				Imaris_Type[file_name]['Nuclei Volume'] = CT_nuclei_TO.loc['Volume', 'Sum']
				Imaris_Type[file_name]['Nuclei Area'] = CT_nuclei_TO.loc['Area', 'Sum']
				Imaris_Type[file_name]['Sphericity'] = CT_nuclei_TO.loc['Sphericity', 'Mean']
				Imaris_Type[file_name]['Nuclei Intensity'] = CT_nuclei_TO.loc['Intensity Mean'].iloc[row_index]['Mean']
			# Handle cases like "14_1", "13_1", etc.
		else:
			#print('Nuclei Actin', file_name)
			# Handle normal files
			Imaris_Type[file_name]['Nuclei Volume'] = CT_nuclei_TO.loc['Volume', 'Sum']
			Imaris_Type[file_name]['Nuclei Area'] = CT_nuclei_TO.loc['Area', 'Sum']
			Imaris_Type[file_name]['Sphericity'] = CT_nuclei_TO.loc['Sphericity', 'Mean']
			Imaris_Type[file_name]['Nuclei Intensity'] = CT_nuclei_TO.loc['Intensity Mean'].iloc[nuclei_channel]['Mean']
		
	# Process chromo Imaris_Type
	for chromoread in chromo_files:
		file_name = process_string(chromoread, marker="(Chromo) ")
		if file_name not in Imaris_Type:
			Imaris_Type[file_name] = {}
		CT_chromo_TO = read_CSV(chromoread)
		
		# Use regular expression to find the pattern of "number_number" in the whole string
		duplicate_string_match = re.search(r'(\d+)_(\d+)', file_name)
		
		if duplicate_string_match:
			# Extract both parts of the pattern
			file_number_first_part = int(duplicate_string_match.group(1))
			file_number_second_part = int(duplicate_string_match.group(2))
			#print('Chromo Actin Split', file_name)
			
			if file_number_first_part >= 1:
				row_index = file_number_second_part #+ nuclei_channel  # Adjust the row index in each iteration
				Imaris_Type[file_name]['Chromo Intensity'] = CT_chromo_TO.loc['Intensity Mean'].iloc[row_index]['Mean']
				Imaris_Type[file_name]['Chromo Volume'] = CT_chromo_TO.loc['Volume', 'Sum']
				Imaris_Type[file_name]['Chromo Count'] = CT_chromo_TO.loc['Volume', 'Count']
		else:
			#print('Chromo Actin', file_name)
			# Handle normal files
			Imaris_Type[file_name]['Chromo Intensity'] = CT_chromo_TO.loc['Intensity Mean'].iloc[nuclei_channel]['Mean']
			Imaris_Type[file_name]['Chromo Volume'] = CT_chromo_TO.loc['Volume', 'Sum']
			Imaris_Type[file_name]['Chromo Count'] = CT_chromo_TO.loc['Volume', 'Count']

	# Process H3K27me3 Imaris_Type
	for h3K27me3read in h3K27me3_files:
		file_name = process_string(h3K27me3read, marker="(H3K27me3) ")
		if file_name not in Imaris_Type:
			Imaris_Type[file_name] = {}
		CT_h3K27me3_TO = read_CSV(h3K27me3read)
		
		# Use regular expression to find the pattern of "number_number" in the whole string
		duplicate_string_match = re.search(r'(\d+)_(\d+)', file_name)
		
		if duplicate_string_match:
			# Extract both parts of the pattern
			file_number_first_part = int(duplicate_string_match.group(1))
			file_number_second_part = int(duplicate_string_match.group(2))
			
			if file_number_second_part == 1:
				print('H3K27me3 Split 1', file_name)
				row_index = 4  # Adjust the row index for file_number_second_part == 1
				Imaris_Type[file_name]['H3K27me3 Intensity'] = CT_h3K27me3_TO.loc['Intensity Mean'].iloc[row_index]['Mean']
			elif file_number_second_part == 2:
				print('H3K27me3 Split 2', file_name)
				row_index = 7  # Adjust the row index for file_number_second_part == 2
				Imaris_Type[file_name]['H3K27me3 Intensity'] = CT_h3K27me3_TO.loc['Intensity Mean'].iloc[row_index]['Mean']
		else:
			print('H3K27me3', file_name)
			# Handle normal files
			Imaris_Type[file_name]['H3K27me3 Intensity'] = CT_h3K27me3_TO.loc['Intensity Mean'].iloc[5]['Mean']

	# Process H3K9me2 Imaris_Type
	for h3K9me2read in h3K9me2_files:
		file_name = process_string(h3K9me2read, marker="(H3K9me2) ")
		if file_name not in Imaris_Type:
			Imaris_Type[file_name] = {}
		CT_h3K9me2_TO = read_CSV(h3K9me2read)
		
		# Use regular expression to find the pattern of "number_number" in the whole string
		duplicate_string_match = re.search(r'(\d+)_(\d+)', file_name)
		
		if duplicate_string_match:
			# Extract both parts of the pattern
			file_number_first_part = int(duplicate_string_match.group(1))
			file_number_second_part = int(duplicate_string_match.group(2))
			
			if file_number_second_part == 1:
				print('H3K9me2 Split 1', file_name)
				row_index = 5  # Adjust the row index for file_number_second_part == 1
				Imaris_Type[file_name]['H3K9me2 Intensity'] = CT_h3K9me2_TO.loc['Intensity Mean'].iloc[row_index]['Mean']
			elif file_number_second_part == 2:
				print('H3K9me2 Split 2', file_name)
				row_index = 8  # Adjust the row index for file_number_second_part == 2
				Imaris_Type[file_name]['H3K9me2 Intensity'] = CT_h3K9me2_TO.loc['Intensity Mean'].iloc[row_index]['Mean']
		else:
			print('H3K9me2', file_name)
			# Handle normal files
			Imaris_Type[file_name]['H3K9me2 Intensity'] = CT_h3K9me2_TO.loc['Intensity Mean'].iloc[6]['Mean']

	# Process Actin Imaris_Type
	for actinread in actin_files:
		file_name = process_string(actinread, marker="(Actin) ")
		if file_name not in Imaris_Type:
			Imaris_Type[file_name] = {}
		CT_actin_TO = read_CSV(actinread)
		Imaris_Type[file_name]['Actin Volume'] = CT_actin_TO.loc['Volume','Sum']

	for file_name in Imaris_Type:
		#print(file_name)
		nuc_area = Imaris_Type[file_name]['Nuclei Area']
		
		# Check if 'Actin Area' exists for the current Imaris_Type point
		if 'Actin Volume' in Imaris_Type[file_name]:
			actin_volume = Imaris_Type[file_name]['Actin Volume']
			
			if actin_volume == "None":
				Imaris_Type[file_name]['Actin Coverage'] = 0
			else:
				coverage = (float(actin_volume) / 2) / float(nuc_area) * 100
				Imaris_Type[file_name]['Actin Coverage'] = coverage
		else:
			# Handle the case where 'Actin Area' does not exist
			Imaris_Type[file_name]['Actin Coverage'] = 0
            
	df = pd.DataFrame.from_dict(Imaris_Type, orient='index')
    
	return df

folder_path = 'C:/Users/kaabi/Documents/Nuceli_Data/Enucleation/230207 OF1 D0 Enuc 5_Ctrl_woArp3 _Result/T0'
# C:/Users/kaabi/Documents/Nuceli_Data/Enucleation/230207 OF1 D0 Enuc 5_Ctrl_woArp3 _Result/T0
# 'C:/Users/kaabi/Documents/Nuceli_Data/Enucleation/230207 OF1 D0 Enuc 5_Ctrl_woArp3 _Result/T30/'
# C:/Users/kaabi/Documents/Nuceli_Data/Enucleation/230208 OF1 D1 Enuc/220208 OF1 D1 Enuc/

df_actin = ProcessData(FolderType= 'Actin', folder_path = folder_path)
df_ctrl = ProcessData(FolderType= 'Ctrl', folder_path = folder_path)

# # Export Ctrl Data 
Result_folder = os.path.join(folder_path, 'Result')
if not os.path.exists(Result_folder):
    os.makedirs(Result_folder)

# Export the DataFrame to an Excel file
df_ctrl.to_excel(Result_folder + '/Export_Ctrl_Excel.xlsx')  

df_actin.to_excel(Result_folder + '/Export_Actin_Excel.xlsx')
# # Sort Based on Coverage Percentage
Actin_L_30 = df_actin[df_actin['Actin Coverage'] < 30]

Actin_L_30.to_excel(Result_folder + '/Actin_Less_30.xlsx')
Actin_O_30 = df_actin[df_actin['Actin Coverage'] > 30]
Actin_O_30.to_excel(Result_folder + '/Actin_Over_30.xlsx')    

# with open (Result_folder + '/Export_Ctrl.csv','w', newline='') as csv_file:
#     csv.writer(csv_file).writerows([k, *v] for k,v in Imaris_ctrl.items())
# colsCtrl = ['nuclei_ctrl_0', 'nuclei_ctrl_1', 'nuclei_ctrl_2','nuclei_ctrl_3']
# df_ctrl[colsCtrl] = df_ctrl[colsCtrl].apply(pd.to_numeric, errors='coerce', axis=1)
# df_ctrl['Ratio1'] = df_ctrl['nuclei_ctrl_0']/df_ctrl['nuclei_ctrl_1']
# df_ctrl['Ratio2'] = df_ctrl['nuclei_ctrl_0']/df_ctrl['nuclei_ctrl_2']

# # Sort based on Nuclei shrinkage >1.2  or <1.2 if inflates
# Ctrl_O_R2 = df_ctrl[df_ctrl['Ratio2'] > 1.2]
# pd.DataFrame(Ctrl_O_R2).to_excel(Result_folder + '/Export_Ctrl_O_R2.xlsx')
# Ctrl_L_R2 = df_ctrl[df_ctrl['Ratio2'] < 1.2]
# pd.DataFrame(Ctrl_L_R2).to_excel(Result_folder + '/Export_Ctrl_L_R2.xlsx')

# pd.DataFrame(df_ctrl).to_excel(Result_folder + '/Export_Ctrl_Excel.xlsx')

# # Export Actin Data
# with open (Result_folder + '/Export_Actin.csv','w', newline='') as csv_file:
#     csv.writer(csv_file).writerows([k, *v] for k,v in Imaris_vca.items())

# df_actin = pd.DataFrame.from_dict(Imaris_vca, orient='index')
# df_actin = df_actin.transpose()

# # colsActin = ['nuclei_vca_0', 'nuclei_vca_1', 'nuclei_vca_2','nuclei_vca_3']
# # df_actin[colsActin] = df_actin[colsActin].apply(pd.to_numeric, errors='coerce', axis=1)
# # df_actin['Ratio1'] = df_actin['nuclei_vca_0']/df_actin['nuclei_vca_1']
# # df_actin['Ratio2'] = df_actin['nuclei_vca_0']/df_actin['nuclei_vca_2']

# #Actin_L_30.to_csv(index=False, path_or_buf='C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Actin_Low_30.csv')
# # Actin_L_30_R2 = Actin_L_30[Actin_L_30['Ratio2'] > 1.2]
# # pd.DataFrame(Actin_L_30_R2).to_excel(Result_folder + '/Export_Actin_Low_30_Over_SHR_1p2_R2.xlsx')
# # Actin_L_30_R1 = Actin_L_30[Actin_L_30['Ratio2'] < 1.2]
# # pd.DataFrame(Actin_L_30_R1).to_excel(Result_folder + '/Export_Actin_Low_30_Less_noShrink_1p2_R2.xlsx')

# Actin_B_30_to_70 = df_actin[df_actin['actin_area_60'].between(30, 70, inclusive="both")]
# #Actin_B_30_to_70.to_csv(index=False, path_or_buf='C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Actin_Bt_30_to_70.csv')
# Actin_O_70 = df_actin[df_actin['Actin Coverage'] > 70]
# #Actin_O_70.to_csv(index=False, path_or_buf='C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Actin_Ovr_70.csv')

# Actin_O_30.to_csv(index=False, path_or_buf=Result_folder + '/Actin_Ovr_30.csv')

# # # Sort based on Nuclei shrinkage >1.2  or <1.2 if inflates
# # Actin_O_30_O_R2 = Actin_O_30[Actin_O_30['Ratio2'] > 1.2]
# # pd.DataFrame(Actin_O_30_O_R2).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Export_Actin_Ovr_30_SHR_1p2_R2.xlsx')
# # Actin_O_30_L_R2 = Actin_O_30[Actin_O_30['Ratio2'] < 1.2]
# # pd.DataFrame(Actin_O_30_L_R2).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Export_Actin_Ovr_30_Less_noShrink_1p2_R2.xlsx')
