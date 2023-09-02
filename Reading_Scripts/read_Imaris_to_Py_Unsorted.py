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

Imaris_vca = {"File_Name_VCA":[],"VCA_chromo_volume_0":[],"VCA_chromo_volume_1":[],"VCA_chromo_volume_2":[],"VCA_chromo_volume_3":[],
           "VCA_chromo_count_0":[],"VCA_chromo_count_1":[],"VCA_chromo_count_2":[],"VCA_chromo_count_3":[],
          "nuclei_vca_0":[],"nuclei_vca_1":[],"nuclei_vca_2":[],"nuclei_vca_3":[],
           "actin_area_30":[],"actin_area_60":[],"sphericity_vca_0":[],"sphericity_vca_1":[],"sphericity_vca_2":[],"sphericity_vca_3":[]}


Imaris_ctrl = {"File_Name_Ctrl":[],"Ctrl_chromo_volume_0":[],"Ctrl_chromo_volume_1":[],"Ctrl_chromo_volume_2":[],"Ctrl_chromo_volume_3":[],
           "Ctrl_chromo_count_0":[],"Ctrl_chromo_count_1":[],"Ctrl_chromo_count_2":[],"Ctrl_chromo_count_3":[],
           "nuclei_ctrl_0":[],"nuclei_ctrl_1":[],"nuclei_ctrl_2":[],"nuclei_ctrl_3":[],
          "sphericity_ctrl_0":[],"sphericity_ctrl_1":[],"sphericity_ctrl_2":[],"sphericity_ctrl_3":[]}
            
# Type- Actin, Ctrl  && Time - T0,T15,T30,T60
# ===========================================

# Sorting Files
actin_files =[]
chromo_files =[]
nuclei_files = [] 
def folder_Scan(Type,Time):
    # Actin after  30min  - value assigned is 3
    directory = 'C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/'
    os.chdir(directory)
    global_path = glob.glob('*/*')
    extension = 'csv'
    for g_name in global_path:
        # Scan Actin/Ctrl first
        if g_name.startswith(Type) and g_name.endswith(Time):
                os.chdir(g_name)
                path_local = glob.glob('*/*.{}'.format(extension))
                for f_name in path_local:
                    #print(f_name)
                    if f_name.startswith('(Chromo)'):
                        chromo_files.append(f_name)
                       # return chromo_files
                    elif f_name.startswith('(NUC04)'):
                        nuclei_files.append(f_name)
                        #return nuclei_files
                    elif f_name.startswith('(Actin)'):
                        actin_files.append(f_name)
                        #return actin_files
                    else:
                        print('Some random file or folders in the directory')
                        
                        
def actin_Coverage(Nuc_Area,Actin_Area):
           # Pd to np et divide
    Nuclei_area_den = list(map(float, Nuc_Area))
    Act_area_to_Flt = list(map(float, Actin_Area))
    value = 2
    Area_Div = [x / value for x in Act_area_to_Flt] 
    # Actin coverage in perentage       
    Actin_area_num = []
    for i in Area_Div :
        Actin_area_num.append(i * 100)
    #Actin_coverage_per = []
    # get last index for the lists for iteration
    end_index = len(Actin_area_num)
    for i in range(end_index):
        Actin_coverage_per = (Actin_area_num[i]/Nuclei_area_den[i])
    return Actin_coverage_per
    #Actin_coverage_per.clear()
# ==============
# Reading CSV
def read_CSV(fil_Nam):
     filename_delimiter = ','
     largest_column_count = 0
     with open(fil_Nam) as temp_fil:
         #col_count = list(map(len(l.split(",")), temp_f.readlines())) # where do I get index
         col_count = [ len(l.split(",")) for l in temp_fil.readlines() ]
         column_names = [j for j in range(0, max(col_count))]
         df = pd.DataFrame()
         df = pd.read_csv(fil_Nam, delimiter=",", names=column_names)
         df.columns = df.iloc[2]
         df = df[3:]
         df = df.set_index(df.columns[0])
     return df
     df.clear()

def get_Filename(fil_Str):
    #filename = Path(file_path).stem
    File_Name=fil_Str[50:]
    index = File_Name.find('Average') # Find endswith key to locate the image number
    index = index - 3 
    File_Name = File_Name[index:index+2]
    File_Name = File_Name.split(',') 
    File_Name = File_Name[0] #extracts the first field
    return File_Name

# ===============
# Actin Channel
# ===============
#### T0
##
actin_files0 = folder_Scan('Actin','T0') # Scan
for chromoread in chromo_files:
    AC_Chromoread_TO = read_CSV(chromoread)
    #exec(f'AC_Chromoread_TO = read_CSV(chromoread)') 
    Imaris_vca["File_Name_VCA"].append(get_Filename(chromoread))
    Imaris_vca["VCA_chromo_volume_0"].append(AC_Chromoread_TO.loc['Volume','Sum'])
    Imaris_vca["VCA_chromo_count_0"].append(AC_Chromoread_TO.loc['Volume','Count'])
chromo_files.clear()     
# Nuclei Read
Nuclei_area_T0 = []
for nucleiread in nuclei_files:
    AC_nuclei_TO = read_CSV(nucleiread)
    Imaris_vca["nuclei_vca_0"].append(AC_nuclei_TO.loc['Volume','Sum'])
    Nuclei_area_T0.append(AC_nuclei_TO.loc['Area','Sum'])
    Imaris_vca["sphericity_vca_0"].append(AC_nuclei_TO.loc['Sphericity','Mean']) 
nuclei_files.clear()  

# ACtin Read
Actin_area_T0 = []  #No area at T0
for actinread in actin_files:
    AC_actin_TO = read_CSV(actinread)
    Actin_area_T0.append(AC_actin_TO.loc['Area','Sum'])
actin_files.clear()    

#### T30
##
actin_files30 = folder_Scan('Actin','T30')  # Scan
for chromoread in chromo_files:
    AC_Chromoread_T3O = read_CSV(chromoread)
    #Imaris_vca["File_Name_VCA"].append(File_Name)
    Imaris_vca["VCA_chromo_volume_1"].append(AC_Chromoread_T3O.loc['Volume','Sum'])
    Imaris_vca["VCA_chromo_count_1"].append(AC_Chromoread_T3O.loc['Volume','Count'])
chromo_files.clear() 

Nuclei_area_T30 = []
for nucleiread in nuclei_files:
    AC_nuclei_T3O = read_CSV(nucleiread)
    Imaris_vca["nuclei_vca_1"].append(AC_nuclei_T3O.loc['Volume','Sum'])
    Nuclei_area_T30.append(AC_nuclei_T3O.loc['Area','Sum'])
    Imaris_vca["sphericity_vca_1"].append(AC_nuclei_T3O.loc['Sphericity','Mean']) 
nuclei_files.clear()  

# Actin Read
Actin_area_T30 = []  
for actinread in actin_files:
    AC_actin_T3O = read_CSV(actinread)
    Actin_area_T30.append(AC_actin_T3O.loc['Area','Sum'])
    Imaris_vca["actin_area_30"].append(actin_Coverage(Nuclei_area_T30,Actin_area_T30))
actin_files.clear() 


#### T60  
actin_files60 = folder_Scan('Actin','T60')  # Scan
for chromoread in chromo_files:
    AC_Chromoread_T6O = read_CSV(chromoread)
    #Imaris_vca["File_Name_VCA"].append(File_Name)
    Imaris_vca["VCA_chromo_volume_2"].append(AC_Chromoread_T6O.loc['Volume','Sum'])
    Imaris_vca["VCA_chromo_count_2"].append(AC_Chromoread_T6O.loc['Volume','Count'])    
chromo_files.clear() 

Nuclei_area_T60 = []
for nucleiread in nuclei_files:
    AC_nuclei_T6O = read_CSV(nucleiread)
    Imaris_vca["nuclei_vca_2"].append(AC_nuclei_T6O.loc['Volume','Sum'])
    Nuclei_area_T60.append(AC_nuclei_T6O.loc['Area','Sum'])
    Imaris_vca["sphericity_vca_2"].append(AC_nuclei_T6O.loc['Sphericity','Mean'])    
nuclei_files.clear()  
# ACtin Read
Actin_area_T60 = []  #No area at T0
for actinread in actin_files:
    AC_actin_T6O = read_CSV(actinread)
    Actin_area_T60.append(AC_actin_T6O.loc['Area','Sum'])
    Imaris_vca["actin_area_60"].append(actin_Coverage(Nuclei_area_T60,Actin_area_T60))

actin_files.clear()  
# ================
# Control Channel
# ================
#### T0
##
ctrl_files0 = folder_Scan('Ctrl','T0') # Scan
for chromoread in chromo_files:
    CT_Chromoread_TO = read_CSV(chromoread)
    Imaris_ctrl["File_Name_Ctrl"].append(get_Filename(chromoread))
    Imaris_ctrl["Ctrl_chromo_volume_0"].append(CT_Chromoread_TO.loc['Volume','Sum'])
    Imaris_ctrl["Ctrl_chromo_count_0"].append(CT_Chromoread_TO.loc['Volume','Count']) 
chromo_files.clear()     
# Nuclei Read
for nucleiread in nuclei_files:
    CT_nuclei_TO = read_CSV(nucleiread)
    Imaris_ctrl["nuclei_ctrl_0"].append(CT_nuclei_TO.loc['Volume','Sum'])
    Nuclei_area_T0.append(CT_nuclei_TO.loc['Area','Sum'])
    Imaris_ctrl["sphericity_ctrl_0"].append(CT_nuclei_TO.loc['Sphericity','Mean']) 
nuclei_files.clear()  

#### T30
##
ctrl_files30 = folder_Scan('Ctrl','T30') # Scan
for chromoread in chromo_files:
    CT_Chromoread_T3O = read_CSV(chromoread)
    Imaris_ctrl["Ctrl_chromo_volume_1"].append(CT_Chromoread_T3O.loc['Volume','Sum'])
    Imaris_ctrl["Ctrl_chromo_count_1"].append(CT_Chromoread_T3O.loc['Volume','Count'])
chromo_files.clear()     
# Nuclei Read
for nucleiread in nuclei_files:
    CT_nuclei_T3O = read_CSV(nucleiread)
    Imaris_ctrl["nuclei_ctrl_1"].append(CT_nuclei_T3O.loc['Volume','Sum'])
    Nuclei_area_T30.append(CT_nuclei_T3O.loc['Area','Sum'])
    Imaris_ctrl["sphericity_ctrl_1"].append(CT_nuclei_T3O.loc['Sphericity','Mean']) 
nuclei_files.clear()  

#### T60
##
ctrl_files60 = folder_Scan('Ctrl','T60') # Scan
for chromoread in chromo_files:
    CT_Chromoread_T6O = read_CSV(chromoread)
    Imaris_ctrl["Ctrl_chromo_volume_2"].append(CT_Chromoread_T6O.loc['Volume','Sum'])
    Imaris_ctrl["Ctrl_chromo_count_2"].append(CT_Chromoread_T6O.loc['Volume','Count'])
chromo_files.clear()     
# Nuclei Read
for nucleiread in nuclei_files:
    CT_nuclei_T6O = read_CSV(nucleiread)
    Imaris_ctrl["nuclei_ctrl_2"].append(CT_nuclei_T6O.loc['Volume','Sum'])
    Nuclei_area_T0.append(CT_nuclei_T6O.loc['Area','Sum'])
    Imaris_ctrl["sphericity_ctrl_2"].append(CT_nuclei_T6O.loc['Sphericity','Mean']) 
nuclei_files.clear()  

# Export Ctrl Data 
with open ('C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Export_Ctrl.csv','w', newline='') as csv_file:
    csv.writer(csv_file).writerows([k, *v] for k,v in Imaris_ctrl.items())
    
df_ctrl = pd.DataFrame.from_dict(Imaris_ctrl, orient='index')
df_ctrl = df_ctrl.transpose()

colsCtrl = ['nuclei_ctrl_0', 'nuclei_ctrl_1', 'nuclei_ctrl_2','nuclei_ctrl_3']
df_ctrl[colsCtrl] = df_ctrl[colsCtrl].apply(pd.to_numeric, errors='coerce', axis=1)
df_ctrl['Ratio1'] = df_ctrl['nuclei_ctrl_0']/df_ctrl['nuclei_ctrl_1']
df_ctrl['Ratio2'] = df_ctrl['nuclei_ctrl_0']/df_ctrl['nuclei_ctrl_2']

# Sort based on Nuclei shrinkage >1.2  or <1.2 if inflates
Ctrl_O_R2 = df_ctrl[df_ctrl['Ratio2'] > 1.2]
pd.DataFrame(Ctrl_O_R2).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Export_Ctrl_O_R2.xlsx')
Ctrl_L_R2 = df_ctrl[df_ctrl['Ratio2'] < 1.2]
pd.DataFrame(Ctrl_L_R2).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Export_Ctrl_L_R2.xlsx')

pd.DataFrame(df_ctrl).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Export_Ctrl_Excel.xlsx')

# Export Actin Data
with open ('C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Export_VCA.csv','w', newline='') as csv_file:
    csv.writer(csv_file).writerows([k, *v] for k,v in Imaris_vca.items())

df_actin = pd.DataFrame.from_dict(Imaris_vca, orient='index')
df_actin = df_actin.transpose()

colsActin = ['nuclei_vca_0', 'nuclei_vca_1', 'nuclei_vca_2','nuclei_vca_3']
df_actin[colsActin] = df_actin[colsActin].apply(pd.to_numeric, errors='coerce', axis=1)
df_actin['Ratio1'] = df_actin['nuclei_vca_0']/df_actin['nuclei_vca_1']
df_actin['Ratio2'] = df_actin['nuclei_vca_0']/df_actin['nuclei_vca_2']

pd.DataFrame(df_actin).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Export_VCA_Excel.xlsx')
# Sort Based on Coverage Percentage
Actin_L_30 = df_actin[df_actin['actin_area_30'] < 30]
#Actin_L_30.to_csv(index=False, path_or_buf='C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Actin_Low_30.csv')
Actin_L_30_R2 = Actin_L_30[Actin_L_30['Ratio2'] > 1.2]
pd.DataFrame(Actin_L_30_R2).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Export_Actin_Low_30_Over_SHR_1p2_R2.xlsx')
Actin_L_30_R1 = Actin_L_30[Actin_L_30['Ratio2'] < 1.2]
pd.DataFrame(Actin_L_30_R1).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Export_Actin_Low_30_Less_noShrink_1p2_R2.xlsx')

pd.DataFrame(Actin_L_30).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Actin_Low_30.xlsx')
Actin_B_30_to_70 = df_actin[df_actin['actin_area_60'].between(30, 70, inclusive="both")]
#Actin_B_30_to_70.to_csv(index=False, path_or_buf='C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Actin_Bt_30_to_70.csv')
Actin_O_70 = df_actin[df_actin['actin_area_60'] > 70]
#Actin_O_70.to_csv(index=False, path_or_buf='C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Actin_Ovr_70.csv')
Actin_O_30 = df_actin[df_actin['actin_area_30'] > 30]
Actin_O_30.to_csv(index=False, path_or_buf='C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Actin_Ovr_30.csv')

# Sort based on Nuclei shrinkage >1.2  or <1.2 if inflates
Actin_O_30_O_R2 = Actin_O_30[Actin_O_30['Ratio2'] > 1.2]
pd.DataFrame(Actin_O_30_O_R2).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Export_Actin_Ovr_30_SHR_1p2_R2.xlsx')
Actin_O_30_L_R2 = Actin_O_30[Actin_O_30['Ratio2'] < 1.2]
pd.DataFrame(Actin_O_30_L_R2).to_excel('C:/Users/kaabi/Documents/Nuceli_Data/Imaris_to_Py/Output/Export_Actin_Ovr_30_Less_noShrink_1p2_R2.xlsx')
