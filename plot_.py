# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 15:57:13 2022

@author: kaabir
"""
import pandas as pd
import glob
import os
import re
import csv
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cbook as cbook
from pathlib import Path
plt.style.use('seaborn-notebook')


Ar_csv = []
Ar_csv1 = []
Chromo_volume_ctrl = []
Chromo_count_ctrl = []


Chromo_volume_VCA = []
Chromo_volume_VCA_Actin = []
Chromo_count_VCA = []
Chromo_count_VCA_Actin = []
Nuclei_area = []
Actin_area = []
Percentage_coverage = []

# Search in directory for a single csv within all folders  
path = 'C:/Users/adity/Documents/20221006 Nuclei vs Chromocenter/Output/'
extension = 'csv'
os.chdir(path)
path_ctrl_chromo_1 = glob.glob('Ctrl/Output_single_csv/*(Chromo)*/*.{}'.format(extension))
path_ctrl_chromo_2 = glob.glob('VCA-Actin 15min/Output_Single_CSV/*(Chromo)*/*.{}'.format(extension))
path_ctrl_chromo_3 = glob.glob('VCA-Actin 30min/Output_Single_CSV/*(Chromo)*/*.{}'.format(extension))
path_ctrl_chromo_4 = glob.glob('VCA-Actin 1H/Output_Single_CSV/*(Chromo)*/*.{}'.format(extension))
path_ctrl_nuclei_1 = glob.glob('Ctrl/Output_single_csv/*(NUC*/*.{}'.format(extension))
path_ctrl_nuclei_2 = glob.glob('VCA-Actin 15min/Output_Single_CSV/*(NUC*/*.{}'.format(extension))
path_ctrl_nuclei_3 = glob.glob('VCA-Actin 30min/Output_Single_CSV/*(NUC*/*.{}'.format(extension))
path_ctrl_nuclei_4 = glob.glob('VCA-Actin 1H/Output_Single_CSV/*(NUC*/*.{}'.format(extension))
path_ctrl_actin_3 = glob.glob('VCA-Actin 30min/Output_Single_CSV/*(Actin)*/*.{}'.format(extension))
path_ctrl_actin_4 = glob.glob('VCA-Actin 1H/Output_Single_CSV/*(Actin)*/*.{}'.format(extension))

# Ctrl without VCA - initial imaging
for filename in path_ctrl_chromo_1:
    # Delimiter
    filename_delimiter = ','
    # The max column count a line in the file could have
    largest_column_count = 0
    # Loop the data lines
    with open(filename) as temp_f:
    # get No of columns in each line
        col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
    # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
        column_names = [i for i in range(0, max(col_count))]
    # Read csv
        df = pd.read_csv(filename, delimiter=",", names=column_names)
        df.columns = df.iloc[2]
        df = df[3:]
        df = df.set_index(df.columns[0])
        Ar_csv.append(df)
        Chromo_volume_ctrl.append(df.loc['Volume','Sum'])
        Chromo_count_ctrl.append(df.loc['Volume','Count'])

# For VCA Chromo without Actin - 15mins
for filename in path_ctrl_chromo_2:
    # Delimiter
    filename_delimiter = ','
    # The max column count a line in the file could have
    largest_column_count = 0
    # Loop the data lines
    with open(filename) as temp_f:
    # get No of columns in each line
        col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
    # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
        column_names = [i for i in range(0, max(col_count))]
    # Read csv
        df1 = pd.read_csv(filename, delimiter=",", names=column_names)
        df1.columns = df1.iloc[2]
        df1 = df1[3:]
        df1 = df1.set_index(df1.columns[0])
        Ar_csv.append(df1)
        Chromo_volume_VCA.append(df1.loc['Volume','Sum'])
        Chromo_count_VCA.append(df1.loc['Volume','Count'])

# For Actin Channel

for filename in path_ctrl_actin_3:
    # Delimiter
    filename_delimiter = ','
    # The max column count a line in the file could have
    largest_column_count = 0
    # Loop the data lines
    with open(filename) as temp_f:
    # get No of columns in each line
        col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
    # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
        column_names = [i for i in range(0, max(col_count))]
    # Read csv
        df2 = pd.read_csv(filename, delimiter=",", names=column_names)
        df2.columns = df2.iloc[2]
        df2 = df2[3:]
        df2 = df2.set_index(df2.columns[0])

        #Ar_csv1.append(df2)
        Actin_area.append(df2.loc['Area','Sum'])

#Actin_area_divided = [float(x) for x in Actin_area]
Actin_area_divided = list(map(float, Actin_area))
value = 2
Actin_area_divided = [x / value for x in Actin_area_divided] 
    
# Nuclei channel for actin coverage

for filename in path_ctrl_nuclei_3:
    filename_delimiter = ','
    largest_column_count = 0
    with open(filename) as temp_f:
        col_count = [ len(l.split(",")) for l in temp_f.readlines() ]    # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
        column_names = [i for i in range(0, max(col_count))]
        df3 = pd.read_csv(filename, header=None, delimiter=",", names=column_names)
        df3.columns = df3.iloc[2]
        df3 = df3[3:]
        df3 = df3.set_index(df3.columns[0])
        #Ar_csv.append(df)
        Nuclei_area.append(df3.loc['Area','Sum'])
   
# Actin coverage in perentage       
numerator = []
for number in Actin_area_divided:
    numerator.append(number * 100)
     
Nuclei_area_conv = list(map(float,Nuclei_area))

Actin_coverage_per = []

# get last index for the lists for iteration
end_index = len(numerator)

for i in range(end_index):
    Actin_coverage_per.append(numerator[i]/Nuclei_area_conv[i])
Actin_coverage_per

    
for filename in path_ctrl_chromo_3:
    filename_delimiter = ','
    largest_column_count = 0
    with open(filename) as temp_f:
        col_count = [ len(l.split(",")) for l in temp_f.readlines() ]    # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
        column_names = [i for i in range(0, max(col_count))]
        df4 = pd.read_csv(filename, header=None, delimiter=",", names=column_names)
        df4.columns = df4.iloc[2]
        df4 = df4[3:]
        df4 = df4.set_index(df4.columns[0])
        Ar_csv1.append(df4)
        Chromo_volume_VCA_Actin.append(df4.loc['Volume','Sum'])
        Chromo_count_VCA_Actin.append(df4.loc['Volume','Count'])    
# =============================================================================
#df= [pd.read_csv(filename, sep='\t') for filename in path_ctrl_chromo]
#chromo_Ctrl = df.drop(df.columns[[0, 1, 2]])
#print(chromo_Ctrl)
# =============================================================================

# =============================================================================
# for filename in path_ctrl_chromo:
#     ##### plot U1
# 
#  	fig0= plt.figure(0)
#  	plt.plot(df["Area"],df["Sum"] ,'-')
#  	plt.xlabel("Nuclei")
#  	plt.ylabel("Chromocenter")
#  	plt.legend(loc="center left")
#  	plt.title("Nuclei vs Chromocenter")
#  	#fig0.savefig("comparison.png", dpi=400)
#plt.show()
# =============================================================================
np.random.seed(8520)
Per = range(len(Actin_coverage_per))
N = 13
colors = np.random.rand(N)
fig, ax = plt.subplots(figsize=[5,4])
plt.scatter(Per,Actin_coverage_per,c=colors, alpha=1)
#ax.plot(Actin_coverage_per,'.', ls='-.', linewidth=1, label='Actin Coverage Percentage')
#ax.plot(Chromo_number_ctrl,'.', ls='-.', linewidth=1,label='Hoescht Dots Volume')



plt.legend()
plt.xlabel(r'$Nuclei$ $Count$')
plt.ylabel(r'$Actin$ $Coverage$ $(percent)$')

#plt.savefig('1.png', dpi = 300, bbox_inches = 'tight')
plt.show()
