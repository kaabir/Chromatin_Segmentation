# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 15:57:13 2022

@author: kaabir
"""

import pandas as pd
import glob
import re
import csv
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cbook as cbook
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go
plt.style.use('seaborn-notebook')



Chromo_volume_ctrl = []
Chromo_number_ctrl = []



Nuclei_area = []
Actin_area = []
Percentage_coverage = []

# Search in directory for a single csv within all folders  
path_ctrl = "G:/Kaabir/221006_Enculeation_RGC_D0_LCbiot_RedVCA/Ctrl/Output_Stats/"
path_ctrl_chromo = glob.glob(path_ctrl + "*(Chromo)*/*.csv")

# Search in directory for a single csv within all folders for chromocenter and actin
path_VCA = "G:/Kaabir/221006_Enculeation_RGC_D0_LCbiot_RedVCA/Actin/Output-Stats-Actin/"
path_VCA_chromo = glob.glob(path_VCA + "*(Chromo)*/*.csv")
path_VCA_actin = glob.glob(path_VCA + "*(Actin)*/*.csv")
path_VCA_nuclei = glob.glob(path_VCA + "*(NUC04)*/*.csv")
# For Control
for filename in path_VCA_chromo:
    # Input
    #data_file = filename

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
        df = pd.read_csv(filename, header=None, delimiter=",", names=column_names)

        #Ar_csv.append(df)
        Chromo_volume_ctrl.append(df.loc[90,8])
        Chromo_number_ctrl.append(df.loc[90,9])

# For Actin Channel
for filename in path_VCA_actin:
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
        df1 = pd.read_csv(filename, header=None, delimiter=",", names=column_names)

        #Ar_csv.append(df)
        Actin_area.append(df1.loc[3,8])

#Actin_area_divided = [float(x) for x in Actin_area]
Actin_area_divided = list(map(float, Actin_area))
value = 2
Actin_area_divided = [x / value for x in Actin_area_divided] 
    # String search
    
# Nuclei channel for actin coverage
for filename in path_VCA_nuclei:
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
        df2 = pd.read_csv(filename, header=None, delimiter=",", names=column_names)

        #Ar_csv.append(df)
        Nuclei_area.append(df2.loc[3,8])
   
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

# =============================================================================
# for filename in path_ctrl_chromo:
 #   df = pd.read_csv(filename)
  #  df = pd.read_csv(filename, sep=';')
#     regex = re.compile(r'Area')
#     Area, Sum = regex.match('Area','Sum').groups()
#     df.rate = df.rate.str.strip('')
#     li.append(df)
# =============================================================================
    
    
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
Per = np.arange(0,9,1)
N = 9
colors = np.random.rand(N)
fig, ax = plt.subplots(figsize=[5,4])
plt.scatter(Per,Actin_coverage_per,c=colors, alpha=1)
#ax.plot(Actin_coverage_per,'.', ls='-.', linewidth=1, label='Actin Coverage Percentage')
#ax.plot(Chromo_number_ctrl,'.', ls='-.', linewidth=1,label='Hoescht Dots Volume')



plt.legend()
plt.xlabel(r'$Label$ (Count)')
plt.ylabel(r'$Number$ $Count$')

#plt.savefig('1.png', dpi = 300, bbox_inches = 'tight')
plt.show()
