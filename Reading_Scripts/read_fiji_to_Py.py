# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 18:23:42 2023

@author: kaabi
"""
import pandas as pd
import glob
import os
import csv

actin_Channel = {"H3k9me2":[],"H3k27me3":[],"Hoescht":[]}
ctrl_Channel = {"H3k9me2":[],"H3k27me3":[],"Hoescht":[]}

# Reading CSV
def read_CSV(fil_Nam):
     with open(fil_Nam) as temp_fil:
         col_count = [ len(l.split(",")) for l in temp_fil.readlines() ]
         column_names = [j for j in range(0, max(col_count))]
         df = pd.DataFrame()
         df = pd.read_csv(fil_Nam, delimiter=",", names=column_names)
         df.columns = df.iloc[0]
         df = df.set_index(df.columns[0])
         df = df[1:]   
     return df
     df.clear()

getFiles = [] 
H3k27me3_files = []   #1 in csv # no actin intensity taken
H3k9me2_files =[]   #2 in csv 
Hoescht = []    # 3 in csv eg. Channel 4 221222 RGC E-Nuc Hoescht B H3K9me2-A568 H3k27me3-A647 Ctrl 01.czi

# Actin channel
directory = 'C:/Users/kaabi/Documents/Nuceli_Data/Fixed Nulcei Staining/221222 RGC E-Nuc Hoescht B H3K9me2-A568 H3k27me3-A647/Output_fiji_Hoescht_H3k9me2_H3k7me3/Actin'
os.chdir(directory)
global_path = glob.glob('*')
extension = '.csv'

for g_name in global_path:
     if g_name.endswith(extension):
            get_files.append(g_name)
     else:
            print('Some random file or folders in the directory') 
        
for f_read in getFiles:
    act_Markers = read_CSV(f_read)
    actin_Channel["H3k27me3"].append(act_Markers.loc['1','Mean'])
    actin_Channel["H3k9me2"].append(act_Markers.loc['2','Mean'])
    actin_Channel["Hoescht"].append(act_Markers.loc['3','Mean'])
getFiles.clear()          

# Ctrl channel
directory = 'C:/Users/kaabi/Documents/Nuceli_Data/Fixed Nulcei Staining/221222 RGC E-Nuc Hoescht B H3K9me2-A568 H3k27me3-A647/Output_fiji_Hoescht_H3k9me2_H3k7me3/Ctrl'
os.chdir(directory)
global_path = glob.glob('*')
extension = 'csv'
for g_name in global_path:
    getFiles.append(g_name)
        
for f_read in getFiles:
    ctrl_Markers = read_CSV(f_read)
    ctrl_Channel["H3k27me3"].append(ctrl_Markers.loc['1','Mean'])
    ctrl_Channel["H3k9me2"].append(ctrl_Markers.loc['2','Mean'])
    ctrl_Channel["Hoescht"].append(ctrl_Markers.loc['3','Mean'])
getFiles.clear()     

df_actin = pd.DataFrame.from_dict(actin_Channel, orient='index')
df_actin = df_actin.transpose()
pd.DataFrame(df_actin).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/Fixed Nulcei Staining/221222 RGC E-Nuc Hoescht B H3K9me2-A568 H3k27me3-A647/Output_fiji_Hoescht_H3k9me2_H3k7me3/Export_actin.xlsx')

df_ctrl = pd.DataFrame.from_dict(ctrl_Channel, orient='index')
df_ctrl = df_ctrl.transpose()
pd.DataFrame(df_ctrl).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/Fixed Nulcei Staining/221222 RGC E-Nuc Hoescht B H3K9me2-A568 H3k27me3-A647/Output_fiji_Hoescht_H3k9me2_H3k7me3/Export_ctrl.xlsx')


