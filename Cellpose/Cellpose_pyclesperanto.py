# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 18:22:51 2023

@author: kaabir
"""

import time, os, sys
import skimage.io 
from skimage.io import imread, imshow
from cellpose import models, utils, plot, io 
import matplotlib as mpl
import napari 
import pyclesperanto_prototype as cle
import numpy as np
import pandas as pd
from aicsimageio import AICSImage

from pathlib import Path
from astropy.visualization import simple_norm
from skimage import img_as_float32
from skimage.util import img_as_ubyte
from tifffile import imsave

# Read the Tif File
img_path = 'C:/Users/kaabi/Documents/Nuceli_Data/Enucleation/Cellpose/Test_Dataset/Files/3D-Q-230120 OF1 D0 Enuc Ctrl5_actin5 T30 #01-1.tif'
img = imread(img_path)
get_Zstack = img.shape[0]
get_Xresl = img.shape[1]
get_yresl = img.shape[2]

# # Read the CZI File
# img_path = AICSImage("C:/Users/kaabi/Documents/Nuceli_Data/Enucleation/Cellpose/Test_Dataset/Files/02.czi")
# # Hoescht Channel
# img = img_path.get_image_data("ZYX", C=2, S=0, T=0)
# get_Zstack = img.shape[0]

# For saving results
filename = Path(img_path).parts[-1]
Result_folder = "C:/Users/kaabi/Documents/Nuceli_Data/Enucleation/Cellpose/Test_Dataset/Result" #@param {type:"string"}

# Upscale the Z dimension to be isotropic
anisotropy = 2*get_Zstack

# channel to segment and nuclear channel numbering starts at 1 
# for your single channel - [0, 0] 
# for the multi channel - [3, 0]
channels = [0, 0] 
diameter = 120 #169.708
use_GPU = True
stitch_threshold=1
#do_3D=True
#resample=True
cellprob_threshold = 0

pretrained_model = "C:/Users/kaabi/Documents/Nuceli_Data/Enucleation/Cellpose/models_trained/CP_20230306_150448_200"

model_match_threshold = 30
flow_threshold = (31.0 - model_match_threshold) / 10.

logger = io.logger_setup()

model = models.CellposeModel(gpu=use_GPU, model_type=pretrained_model)#,residual_on=False, style_on=False, concatenation=False)#, resample=True, do_3D=True,stitch_threshold=1)
mask, flows,styles = model.eval(img,
                                channels =channels,
                                #anisotropy=anisotropy,
                                diameter=diameter, 
                                #pretrained_model=pretrained_model,
                                do_3D=False,
                                resample=True,
                                min_size = -1,
                                flow_threshold =flow_threshold,
                                cellprob_threshold = cellprob_threshold,
                                net_avg=True,
                                stitch_threshold=stitch_threshold,
                                #model_match_threshold = model_match_threshold
                                    )

# Save Labels
prediction_stack_32 = img_as_float32(mask, force_copy=False)     
os.chdir(Result_folder)
imsave(str(filename)+".tif", prediction_stack_32)
 # save results to load in GUI
io.masks_flows_to_seg(img, mask, flows, diameter, filename, channels)

outlines = utils.masks_to_outlines(mask)
viewer = napari.Viewer() 
viewer.add_image(img) 
viewer.add_labels(mask) 
viewer.add_labels(outlines)

# Merge Segmented Mask
merged_Labels = cle.connected_components_labeling_box(mask)

# Getting Original intensities
#intensity_vector = cle.read_intensities_from_map(mask, img)
statistics = cle.statistics_of_labelled_pixels(img, merged_Labels)

# Conversion to um
px_to_um_X = 0.0351435 
px_to_um_Y = 0.0351435
nuclei_Area = statistics['area']*px_to_um_Y*px_to_um_Y

intensity_vector = cle.read_intensities_from_map(img,merged_Labels)
intensity_vector_1_voxel = cle.read_intensities_from_map(img,mask)

#statistics_df = pd.DataFrame(statistics)
#statistics_df.head()

# Chromocenter Segmentation
# Get voxel of pixel intensity values inside the mask 
def intensity_vector(mask, img):
    if not (mask.shape == img.shape):
        return False
    
    mat_intensity = np.where(np.logical_and(mask,img),img,0)
    return mat_intensity

intensity_map= intensity_vector(mask, img)

mean_Intensity = statistics['mean_intensity']
stD_Intensity = statistics['standard_deviation_intensity']
factor = 2.5
threshold_Chromo = mean_Intensity + factor * stD_Intensity

# Regions of intensity greater than mean/threshold value
def Chromo_vector(intensity_vector, threshold_value):
    mat_chromo = np.where(intensity_vector>threshold_value,intensity_vector,0)
    return mat_chromo

mat_intensity = Chromo_vector(intensity_map,threshold_Chromo)
segmented_chromo = cle.gaussian_blur(mat_intensity, sigma_x=2, sigma_y=2)
segmented_chromo = cle.threshold_otsu(segmented_chromo)

# Merge Segmented Mask
chromo_Labels = cle.connected_components_labeling_box(segmented_chromo) #Chromo are seperate
# Getting stats
statistics_chromo = cle.statistics_of_labelled_pixels(img, chromo_Labels)
chromo_Area = np.sum(statistics_chromo['area'], axis=0)
chromo_Area = chromo_Area * px_to_um_Y*px_to_um_X

# Max Entropy Chromocenter Segmentation -> Orginial Voxel -> Blur (rx=ry=2) -> Max(Entropy Threshold) 
intensity_map_blur = cle.gaussian_blur(intensity_map, sigma_x=2, sigma_y=2)
entropy_Threshold_Chromo = nsitk.threshold_maximum_entropy(intensity_map_blur)

# Sort Nuclei Label Mask
#print(len(np.unique(merged_Labels)))
number_of_Nuclei = len(np.unique(merged_Labels))-1

im_obj1 = np.zeros(img.shape) 
im_obj2 = np.zeros(img.shape) 
im_obj3 = np.zeros(img.shape)
nuc_lbl_lst1 = []
nuc_lbl_lst2 = []
nuc_lbl_lst3 = []

for i in np.unique(merged_Labels_np)[1:]: # Initial (1) is mask
    if (len(np.unique(merged_Labels)-1) == 1):
        break
    elif i ==1:
        im_obj1[merged_Labels_np == i] = 1
        label1 = im_obj1
        nuc_lbl_lst1.append(label1)
        #viewer = napari.Viewer()
        #viewer.add_image(im_obj1)
    elif i ==2: 
        im_obj2[merged_Labels_np == i] = 1
        label2 = im_obj2
        nuc_lbl_lst2.append(label2)
        #viewer = napari.Viewer()
        #viewer.add_image(im_obj2)
       
    elif i ==3: 
        im_obj3[merged_Labels_np == i] = 1
        label3 = im_obj3
        nuc_lbl_lst3.append(label3)        
        #viewer = napari.Viewer()
        #viewer.add_image(im_obj3)

    else:
        print(' More than four nuclei...')    

# Not replacing back the original mask values associted with labels
# Place back the intensity values to the individual label
def chromocenter_vector(mask, img):
    if not (mask.shape == img.shape):
        return False
      
    chromocenter_intensity = np.where(np.logical_and(mask,img),img,0) # Overlap Mask on original image (as reference)
    return chromocenter_intensity

def empty_nuclei(nuclei_label):
    if nuclei_label:
        #print('Label is Not an empty list')
        return True
    #print('Label is an empty list')
    return False

for i in np.unique(merged_Labels_np)[1:]:
    # Nuclei one to be segmented for chromocenter 1
    if not empty_nuclei(nuc_lbl_lst1):
        break
    nuc_lbl_lst1 = np.vstack(nuc_lbl_lst1)
    intensity_map1= chromocenter_vector(nuc_lbl_lst1, img)
    # Check the nuclei
    viewer = napari.Viewer()
    viewer.add_image(intensity_map1)
    # Nuclei 2 to be segmented for chromocenter 2
    if not empty_nuclei(nuc_lbl_lst2):
        break
    nuc_lbl_lst2 = np.vstack(nuc_lbl_lst2)
    intensity_map2= chromocenter_vector(nuc_lbl_lst2, img)        
    # Nuclei 3 to be segmented for chromocenter 3
    if not empty_nuclei(nuc_lbl_lst3):
        break
    nuc_lbl_lst3 = np.vstack(nuc_lbl_lst3)
    intensity_map3= chromocenter_vector(nuc_lbl_lst3, img)

 # Method running chromocenters segmentation.
 # Simple algorithms description :
 # 1- Gaussian blur
 # 2- Image gradient : for each voxels computing sum of difference with
 # X (_neigh parameter) neighborhood.
 # Thresholding value : keep voxels having value higher mean nucleus
 # intensity plus X (_factor parameter) standard deviation value
 # 4- Binarization of threshold image
 # 5- Connected component computation from binarized image
