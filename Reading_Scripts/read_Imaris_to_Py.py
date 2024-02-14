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
    extension = 'xls'
    for g_name in global_path:
        # Scan Actin/Ctrl first
        if g_name.startswith(Type):
                os.chdir(g_name)
                path_local = glob.glob('*.{}'.format(extension))
                for f_name in path_local:
                    #print(f_name)
                    if f_name.startswith('(NUC04)'):
                        nuclei_files.append(f_name)
                    elif f_name.startswith('(Actin)'):
                        actin_files.append(f_name)
                    elif f_name.startswith('(Chromo)'):
                        chromo_files.append(f_name)
                    elif f_name.startswith('(H3k27me3)'):
                        h3K27me3_files.append(f_name)
                    elif f_name.startswith('(H3k9me2)'):
                        h3K27me3_files.append(f_name)                          
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

def combined_metric(Ch_intensity=None,Ch_inten_weight=None, Ch_vol=None,Ch_vol_weight=None, nuc_volume=None):

    #den = intens/volume
    #density = (Ch_inten_weight * Ch_intensity) /(Ch_vol_weight * Ch_vol)
    combined_metric = Ch_vol/ nuc_volume
    return combined_metric


def calculate_frequency_distr(position_input, volume_input, num_bins=None):
    """
    Calculate frequency distribution and mean volume values for given position and volume data.
    
    Parameters:
        position_values (list): List of position values.
        volume_values (list): List of volume values.
        num_bins (int): Number of bins for position values. Default is 8.
    
    Returns:
        combined_df (DataFrame): DataFrame containing frequency distribution and mean volume values.
    """
    dfs = []
    for volume_vals, position_vals in zip(volume_input, position_input):
        df = pd.DataFrame({'Position': position_vals, 'Volume': volume_vals})
        dfs.append(df)

    # Merge DataFrames based on a common identifier or index if available
    merged_df = pd.concat(dfs, ignore_index=True)

    # Define bins for position values
    position_values = merged_df['Position']
    volume_values = merged_df['Volume']

    num_bins = num_bins  # Adjust the number of bins as needed
    bins = pd.cut(position_values, bins=num_bins, include_lowest=True)

    # Create a DataFrame to store position, volume, and their corresponding bins
    data = pd.DataFrame({'Position': position_values, 'Volume': volume_values, 'Position_Bin': bins})

    # Calculate frequency distribution for each bin
    frequency_distribution = data['Position_Bin'].value_counts().sort_index()

    # Calculate relative frequency percentage for each bin
    relative_frequency_percentage = frequency_distribution / len(position_values) * 100

    # Group data by position bins and calculate mean volume value for each bin
    grouped_data = data.groupby('Position_Bin')['Volume'].mean()

    # Convert frequency distribution and mean volume values to DataFrames
    combined_df = pd.DataFrame({'Position_Bin': frequency_distribution.index, 'Frequency': frequency_distribution.values})
    
    # Convert Interval objects to strings and adjust labels
    combined_df['Position_Bin'] = combined_df['Position_Bin'].apply(lambda x: round((x.left + x.right) / 2, 3))
    combined_df['Position_Bin'] = combined_df['Position_Bin'].astype(str) 
    
    combined_df['Mean_Volume'] = grouped_data.values

    # Output the combined DataFrame
    print(combined_df)
    #return combined_df['Position_Bin'], combined_df['Frequency'], combined_df['Mean_Volume']
    return combined_df

def process_string(input_string, marker = None):
    # Find the last occurrence of "(Chromo)"or "(Actin)" or "(NUC04)"
    # 
    directory = os.path.splitext(input_string)[0]
    #print(directory)
    path = os.path.basename(directory)
    #print(path)
    identifier = path.split(marker)[1]
    #print(identifier)
    
    return identifier

def ProcessData(FolderType=None, folder_path=None):
    print(FolderType)
    Imaris_Type = {}
    # nuclei_channel = 1

    actin_files, h3K9me2_files, h3K27me3_files, nuclei_files, chromo_files = folder_Scan(folder_path, FolderType)

    # Nucleus Read
    for nucleiread in nuclei_files:
        file_name = process_string(nucleiread, marker="(NUC04)")
        
        print(nucleiread)
        if file_name not in Imaris_Type:
            Imaris_Type[file_name] = {}
        
        # Read the entire Excel file with sheets
        ## Nucleus Volume
        nuc_excel = pd.ExcelFile(nucleiread)
        df = pd.read_excel(nuc_excel, sheet_name='Volume', header=1)      
        Imaris_Type[file_name]['Nuclei Volume'] = df['Volume'].sum()
        df = pd.read_excel(nuc_excel, sheet_name='Area', header=1)      
        Imaris_Type[file_name]['Nuclei Area'] = df['Area'].sum()
        df = pd.read_excel(nuc_excel, sheet_name='Sphericity', header=1)
        Imaris_Type[file_name]['Nuclei Sphericity'] = df['Sphericity'].sum()
        
        # Nucleus Intensity Read
        sheet_names = nuc_excel.sheet_names
        #print(sheet_names)
        # Last Channel is Nuclei == Maximum channel
        max_channel = max(int(name.split('Ch=')[1].split()[0]) for name in sheet_names if 'Intensity Mean' in name)
        # Construct the target sheet name
        target_sheet_name = f'Intensity Mean Ch={max_channel} Img=1'
        # Read the specified sheet
        nuc_inten = pd.read_excel(nuc_excel, sheet_name=target_sheet_name, header=1)
        Imaris_Type[file_name]['Nuclei Intensity'] = nuc_inten['Intensity Mean'].sum()
    
    # Chromocenter Read           
    for chromoread in chromo_files:
        file_name = process_string(chromoread, marker="(Chromo)")
        #print(file_name)
        if file_name not in Imaris_Type:
            Imaris_Type[file_name] = {}
        chromo_excel = pd.ExcelFile(chromoread)    
        df = pd.read_excel(chromo_excel, sheet_name='Volume', header=1)
        Imaris_Type[file_name]['Chromocenter Volume'] = df['Volume'].sum()
        
        chroVol_values = df['Volume'].values#.flatten()#.reshape(-1, 1)
        Imaris_Type[file_name]['Chromocenter Volume Values'] = chroVol_values
        # Coordinates of chromocenters to the nucleus surface 
        # Shortest Distance to Surfac-88
        
        sheet_names = chromo_excel.sheet_names

        # Last Channel is Chromocenter == Maximum channel
        max_channel = max(int(name.split('Ch=')[1].split()[0]) for name in sheet_names if 'Intensity Mean' in name)
        # Construct the target sheet name
        target_sheet_name = f'Intensity Mean Ch={max_channel} Img=1'
        # Read the specified sheet
        chromo_inten = pd.read_excel(chromo_excel, sheet_name=target_sheet_name, header=1)
        Imaris_Type[file_name]['Chromocenter Intensity Values'] = chromo_inten['Intensity Mean'].values.flatten()#.sum()
        Imaris_Type[file_name]['Chromocenter Mean'] = chromo_inten['Intensity Mean'].sum()
        # Find the maximum channel number
        max_channel = max(int(name.split('Shortest Distance to Surfac-')[1].split()[0]) for name in sheet_names if 'Shortest Distance to Surfac-' in name)

        for target_channel in range(max_channel, 0, -1):
    # Construct the target sheet name
            target_sheet_name = f'Shortest Distance to Surfac-{target_channel}'

            try:
        # Read the specified sheet
                df_coordinates = pd.read_excel(chromo_excel, sheet_name=target_sheet_name)
        
        # Check if the desired column exists
                if 'Shortest Distance to Surfaces Surfaces=nuc' in df_coordinates.columns:
                    df_coordinates = pd.read_excel(chromo_excel, sheet_name=target_sheet_name, skiprows=1)
                    Imaris_Type[file_name]['Chromocenter Distance to Borders'] = df_coordinates['Shortest Distance to Surfaces'].values.flatten()
    
                    collected_positions = df_coordinates['Shortest Distance to Surfaces'].tolist()  # Convert to list
                    #print(collected_positions)
                    # This divides the data in three bins of range from -1 to 0
                    bins = [-1, -0.7, -0.4, 0]
                    categories = []
                    bins_aggregated = []

                    for i in range(len(bins) - 1):
                        lower_bound = bins[i]
                        upper_bound = bins[i + 1]
        
                        count = len([val for val in collected_positions if lower_bound <= val < upper_bound])
                        categories.append(f'{lower_bound}-{upper_bound}')
                        bins_aggregated.append(count)
                    # Center = (bins_aggregated[0]/ sum(bins_aggregated)) *100 Frequency percentage of each region
                    Imaris_Type[file_name]['Center'] = (bins_aggregated[0]/ sum(bins_aggregated)) *100#bins_aggregated[0]
                    Imaris_Type[file_name]['Inside'] = (bins_aggregated[1]/ sum(bins_aggregated)) *100#bins_aggregated[1]
                    Imaris_Type[file_name]['Edge'] = (bins_aggregated[2]/ sum(bins_aggregated)) *100#bins_aggregated[2]              

                    #print(bins_aggregated)
                    break

            except ValueError:
                
        # Handle the case when the sheet cannot be read
                print(f"Error reading sheet: {target_sheet_name}")

    df = pd.DataFrame.from_dict(Imaris_Type, orient='index')
    
    df['Chromo Foci Density'] = combined_metric(Ch_intensity=df['Chromocenter Mean'],Ch_inten_weight=0.6, Ch_vol=df['Chromocenter Volume'],Ch_vol_weight=0.4, nuc_volume=df['Nuclei Volume'])   
    
    return df

folder_path = 'E:/Quantified/Agarose Compression/[Imaris 3D] 240207 e230926 OF1 D6 Tf shNeg Gfop Rdsred FRfoxj1/Result/'

#df_actin = ProcessData(FolderType= 'spVCA_but_no_Actin', folder_path = folder_path)
df_Diff = ProcessData(FolderType= 'Diff', folder_path = folder_path)
freq_diff_df = calculate_frequency_distr(df_Diff['Chromocenter Distance to Borders'], df_Diff['Chromocenter Volume Values'], num_bins=22)
df_Diff = pd.concat([df_Diff, freq_diff_df], axis=1)

# Multiciliated Cells
df_Multi_EC = ProcessData(FolderType= 'Multi_EC', folder_path = folder_path)
freq_Multi_df = calculate_frequency_distr(df_Multi_EC['Chromocenter Distance to Borders'], df_Multi_EC['Chromocenter Volume Values'], num_bins=22)
df_Multi_EC = pd.concat([df_Multi_EC, freq_Multi_df], axis=1)
# Radial Glial Cells
df_RGC = ProcessData(FolderType= 'RGC', folder_path = folder_path)
freq_RGC_df = calculate_frequency_distr(df_RGC['Chromocenter Distance to Borders'], df_RGC['Chromocenter Volume Values'], num_bins=22)
df_RGC = pd.concat([df_RGC,freq_RGC_df], axis=1)

# # Export Ctrl Data 
Result_folder = os.path.join(folder_path)#, 'Result')
if not os.path.exists(Result_folder):
    os.makedirs(Result_folder)

# Export the DataFrame to an Excel file
df_Diff.to_excel(Result_folder + '/Export_Diff_Stage.xlsx')  
df_Multi_EC.to_excel(Result_folder + '/Export_Multi_Stage.xlsx')
df_RGC.to_excel(Result_folder + '/Export_RGC_Stage.xlsx')
#df_14percent.to_excel(Result_folder + '/Export_14percent_Excel.xlsx')  

# df_actin.to_excel(Result_folder + '/Export_Actin_Excel.xlsx')
# # # Sort Based on Coverage Percentage
# Actin_L_30 = df_actin[df_actin['Actin Coverage'] < 30]

# Actin_L_30.to_excel(Result_folder + '/Actin_Less_30.xlsx')
# Actin_O_30 = df_actin[df_actin['Actin Coverage'] > 30]
# Actin_O_30.to_excel(Result_folder + '/Actin_Over_30.xlsx')    
