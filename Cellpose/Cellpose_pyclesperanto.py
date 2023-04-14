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
from napari_segment_blobs_and_things_with_membranes import split_touching_objects
import napari_simpleitk_image_processing as nsitk

from pathlib import Path
from astropy.visualization import simple_norm
from skimage import img_as_float32
from skimage.util import img_as_ubyte
from tifffile import imsave

# Read the Tif File
folder_path = 'C:/Users/kaabir/Documents/'
file_path = '212155 t30 #09.tif'
img_path = folder_path + file_path
img = imread(img_path)

# In case two channels
# img_actin = img[:,0,:,:] # A Channel
# img_nuclei = img[:,1,:,:] # N Channel

# In case three channels
img_actin = img[:,:,:,0] # A Channel
img_dextran = img[:,:,:,1] # Dextran Channel
img_nuclei = img[:,:,:,2] # N Channel

get_Zstack = img_nuclei.shape[0]
get_Xresl = img_nuclei.shape[1]
get_yresl = img_nuclei.shape[2]

# # Read the CZI File
# img_path = AICSImage("C:/Users/kaabi/Documents/1234.czi")
# # Hoescht Channel
# img = img_path.get_image_data("ZYX", C=2, S=0, T=0)
# get_Zstack = img.shape[0]

# For saving results
filename = Path(file_path).stem
Result_folder = folder_path + '/Result/' 
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

 # save results to load in GUI
#io.masks_flows_to_seg(img, mask, flows, diameter, filename, channels)
# Merge Segmented Mask
merged_Labels = cle.connected_components_labeling_box(mask)
merged_Labels_np = np.array(merged_Labels)

# Getting Original intensities
#intensity_vector = cle.read_intensities_from_map(mask, img_nuclei)
#statistics = cle.statistics_of_labelled_pixels(img_nuclei, merged_Labels)
#
# Conversion to um
px_to_um_X = 0.0351435
px_to_um_Y = 0.0351435
px_to_um_Z = 1
#nuclei_Area = statistics['area']*px_to_um_Y*px_to_um_Y*px_to_um_Z

# Chromocenter Segmentation
# Get voxel of pixel intensity values inside the mask 
def intensity_vector(mask, img_nuclei):
    if not (mask.shape == img_nuclei.shape):
        return False
    
    mat_intensity = np.where(np.logical_and(mask,img_nuclei),img_nuclei,0) # Overlap Mask on original image (as reference)
    return mat_intensity

intensity_map= intensity_vector(mask, img_nuclei)

def normalize_intensity(k):
    k_min = k.min(axis=(1, 2), keepdims=True)
    k_max = k.max(axis=(1, 2), keepdims=True)

    k = (k - k_min)/(k_max-k_min)
    k[np.isnan(k)] = 0
    return k

# structure for actin segmentation
from scipy import ndimage as ndi
diamond = np.zeros((3, 3, 3), dtype=bool)
diamond[1:2, 0:3, 0:3] = True
  
#viewer = napari.Viewer()
#viewer.add_image(mask)
#viewer.add_image(intensity_map)

from scipy import ndimage
from scipy.ndimage import gaussian_filter

#print(len(np.unique(merged_Labels)))
number_of_Nuclei = len(np.unique(merged_Labels))-1

# Create a matrice to place the single nuleus with intensity voxel
im_obj1 = np.zeros(img_nuclei.shape) 
im_obj2 = np.zeros(img_nuclei.shape) 
im_obj3 = np.zeros(img_nuclei.shape)
nuc_lbl_lst1 = []
nuc_lbl_lst2 = []
nuc_lbl_lst3 = []
# Create actin mask
actin_lbl_lst1 = []
actin_lbl_lst2 = []
actin_lbl_lst3 = []
# Extracting mask for each nucleus intensity and storing them
for i in np.unique(merged_Labels_np)[1:]: # Initial (1) is mask
    if (len(np.unique(merged_Labels)-1) == 1):
        break
    elif i ==1:
        im_obj1[merged_Labels_np == i] = 1
        label1 = im_obj1
        nuc_lbl_lst1.append(label1)
        # Save Labels
        prediction_stack_32 = img_as_float32(label1, force_copy=False)     
        os.chdir(Result_folder)
        imsave(filename+"_1.tif", prediction_stack_32)
        # Use nuclei mask as seed for actin segmentation
        # Pad some pixels for actin with constant value which matches unique label value
        # place back actin and segment with classical thresholding
        #lbl1_actin = np.pad(label1, (25), 'constant', constant_values=(1)) # constant_values to fill mask with 
        dilated1 = ndi.binary_dilation(label1, diamond, iterations=10).astype(label1.dtype)
        actin_img = intensity_vector(dilated1,img_actin)
        # Apply filter  
        actin_filter = nsitk.median_filter(actin_img, radius_x:=2, radius_y=2)
        # Apply Otsu Threshold
        actin_threshold_1 = nsitk.threshold_otsu(actin_filter)

        #actin_lbl_lst1.append(lbl1_actin)#
        # viewer = napari.Viewer()
        # viewer.add_image(im_obj1)
    elif i ==2: 
        im_obj2[merged_Labels_np == i] = 1
        label2 = im_obj2
        nuc_lbl_lst2.append(label2)
        # Save Labels
        prediction_stack_32 = img_as_float32(label2, force_copy=False)     
        os.chdir(Result_folder)
        imsave(filename+"_2.tif", prediction_stack_32)
        # viewer = napari.Viewer()
        # viewer.add_image(im_obj2)
       
    elif i ==3: 
        im_obj3[merged_Labels_np == i] = 1
        label3 = im_obj3
        nuc_lbl_lst3.append(label3)    
        # Save Labels
        prediction_stack_32 = img_as_float32(label3, force_copy=False)     
        os.chdir(Result_folder)
        imsave(filename+"_3.tif", prediction_stack_32)
        # viewer = napari.Viewer()
        # viewer.add_image(im_obj3)

    else:
        raise ValueError('More than three nucleus')   
    # viewer = napari.Viewer()
    # viewer.add_image(y)

# Not replacing back the original mask values associted with labels
# Place back the intensity values to the individual label
def chromocenter_vector(mask, img_nuclei):
    if not (mask.shape == img_nuclei.shape):
        return False

    chromocenter_intensity = np.where(np.logical_and(mask,img_nuclei),img_nuclei,0) # Overlap Mask on original image (as reference)
    return chromocenter_intensity

def empty_nuclei(nuclei_label):
    if nuclei_label:
        #print('Label is Not an empty list')
        return True
    #print('Label is an empty list')
    return False

for i in np.unique(merged_Labels_np)[1:]: # Background mask at 1
    # Nuclei one to be segmented for chromocenter 1
    if not empty_nuclei(nuc_lbl_lst1):
        break
    nuc_lbl_lst1 = np.vstack(nuc_lbl_lst1)
    intensity_map1= chromocenter_vector(nuc_lbl_lst1, img_nuclei)
    intensity_map1_norm = normalize_intensity(intensity_map1)
    
    # Get Nucleus Statistics
    statistics_nucleus1 = cle.statistics_of_labelled_pixels(img_nuclei, nuc_lbl_lst1)
    pd.DataFrame(statistics_nucleus1).to_excel(Result_folder+filename+'(Nucleus)Nucleus_statistics_1.xlsx')
 
    # Blur to enhance chromocenters
    #intensity_map_blur1 = nsitk.gaussian_blur(intensity_map1, variance_x=2, variance_y=2)
    intensity_map_blur1 = nsitk.median_filter(intensity_map1_norm, radius_x:=2, radius_y=2)
    # Apply Intermodes Thresholding for Chromocenters
    intermodes_Threshold_Chromo1 = nsitk.threshold_intermodes(intensity_map_blur1)
    # Caculate intermodes Area
    chromo_intermodes_Labels1 = cle.connected_components_labeling_box(intermodes_Threshold_Chromo1)
    statistics_intermodes_chromo1 = cle.statistics_of_labelled_pixels(img_nuclei, chromo_intermodes_Labels1)
    
    chromointermodes_Area1 = np.sum(statistics_intermodes_chromo1['area'], axis=0)
    chromointermodes_Area1 = chromointermodes_Area1 * px_to_um_Y*px_to_um_X*px_to_um_Z
    # Export Chromo 1 Stats 
    pd.DataFrame(statistics_intermodes_chromo1).to_excel(Result_folder+filename+'(Chromo)Chromo_statistics_1.xlsx')
    
    # Max Entropy Chromocenter Segmentation -> Orginial Voxel -> Blur (rx=ry=2) -> Max(Entropy Threshold) 
    entropy_Threshold_Chromo1 = nsitk.threshold_maximum_entropy(intensity_map_blur1)
    chromo_entropy_Labels1 = cle.connected_components_labeling_box(entropy_Threshold_Chromo1)
    # Caculate Entropy Area
    statistics_entropy_chromo1 = cle.statistics_of_labelled_pixels(img_nuclei, chromo_entropy_Labels1)
    
    chromoentropy_Area1 = np.sum(statistics_entropy_chromo1['area'], axis=0)
    chromoentropy_Area1 = chromoentropy_Area1 * px_to_um_Y*px_to_um_X*px_to_um_Z
    
    #binary_fill_holes = binary_fill_holes(intermodes_Threshold_Chromo1)
    #exclude_tiny_label = cle.exclude_small_labels(binary_fill_holes, maximum_size = 500)
    #exclude_tiny_label = cle.exclude_large_labels(binary_fill_holes, minimum_size = 40)

    # Check the nuclei
    viewer = napari.Viewer()
    #viewer.add_image(binary_fill_holes)
    viewer.add_image(intensity_map1_norm)
    viewer.add_image(chromo_entropy_Labels1)
    viewer.add_image(chromo_intermodes_Labels1)
    #viewer.add_image(chromoshanbhag_Area1)
    
    # Nuclei 2 to be segmented for chromocenter 2
    if not empty_nuclei(nuc_lbl_lst2):
        break
    nuc_lbl_lst2 = np.vstack(nuc_lbl_lst2)
    intensity_map2= chromocenter_vector(nuc_lbl_lst2, img_nuclei)
    intensity_map2_norm = normalize_intensity(intensity_map2)

    # Get Nucleus Statistics
    statistics_nucleus2 = cle.statistics_of_labelled_pixels(img_nuclei, nuc_lbl_lst2)
    pd.DataFrame(statistics_nucleus2).to_excel(Result_folder+filename+'(Nucleus)Nucleus_statistics_2.xlsx')

    # Blur to enhance chromocenters
    #intensity_map_blur2 = nsitk.gaussian_blur(intensity_map2, variance_x=2, variance_y=2)
    intensity_map_blur2 = nsitk.median_filter(intensity_map2_norm, radius_x:=2, radius_y=2)
    # Apply Intermodes Thresholding for Chromocenters
    intermodes_Threshold_chromo2 = nsitk.threshold_intermodes(intensity_map_blur2)
    # Caculate intermodes Area
    chromo_intermodes_Labels2 = cle.connected_components_labeling_box(intermodes_Threshold_chromo2)
    statistics_intermodes_chromo2 = cle.statistics_of_labelled_pixels(img_nuclei, chromo_intermodes_Labels2)

    chromointermodes_Area2 = np.sum(statistics_intermodes_chromo2['area'], axis=0)
    chromointermodes_Area2 = chromointermodes_Area2 * px_to_um_Y*px_to_um_X*px_to_um_Z
    # Export Chromo 1 Stats 
    pd.DataFrame(statistics_intermodes_chromo2).to_excel(Result_folder+filename+'(Chromo)Chromo_statistics_2.xlsx')

    # Max Entropy Chromocenter Segmentation -> Orginial Voxel -> Blur (rx=ry=2) -> Max(Entropy Threshold) 
    entropy_Threshold_chromo2 = nsitk.threshold_maximum_entropy(intensity_map_blur2)
    chromo_entropy_Labels2 = cle.connected_components_labeling_box(entropy_Threshold_chromo2)
    # Caculate Entropy Area
    statistics_entropy_chromo2 = cle.statistics_of_labelled_pixels(img_nuclei, chromo_entropy_Labels2)

    chromoentropy_Area2 = np.sum(statistics_entropy_chromo2['area'], axis=0)
    chromoentropy_Area2 = chromoentropy_Area2 * px_to_um_Y*px_to_um_X*px_to_um_Z

    #binary_fill_holes = binary_fill_holes(intermodes_Threshold_chromo2)
    #exclude_tiny_label = cle.exclude_small_labels(binary_fill_holes, maximum_size = 500)
    #exclude_tiny_label = cle.exclude_large_labels(binary_fill_holes, minimum_size = 40)

    # Check the nuclei
    viewer = napari.Viewer()
    #viewer.add_image(binary_fill_holes)
    viewer.add_image(intensity_map2_norm)
    viewer.add_image(chromo_entropy_Labels2)
    viewer.add_image(chromo_intermodes_Labels2)
    
    # Nuclei 3 to be segmented for chromocenter 3
    if not empty_nuclei(nuc_lbl_lst3):
        break
    nuc_lbl_lst3 = np.vstack(nuc_lbl_lst3)
    intensity_map3= chromocenter_vector(nuc_lbl_lst3, img_nuclei)
    intensity_map3_norm = normalize_intensity(intensity_map3)

    # Get Nucleus Statistics
    statistics_nucleus3 = cle.statistics_of_labelled_pixels(img_nuclei, nuc_lbl_lst3)
    pd.DataFrame(statistics_nucleus3).to_excel(Result_folder+filename+'(Nucleus)Nucleus_statistics_3.xlsx')

    # Blur to enhance chromocenters
    #intensity_map_blur3 = nsitk.gaussian_blur(intensity_map3, variance_x=2, variance_y=2)
    intensity_map_blur3 = nsitk.median_filter(intensity_map3_norm, radius_x:=2, radius_y=2)
    # Apply Intermodes Thresholding for Chromocenters
    intermodes_Threshold_chromo3 = nsitk.threshold_intermodes(intensity_map_blur3)
    # Caculate intermodes Area
    chromo_intermodes_Labels3 = cle.connected_components_labeling_box(intermodes_Threshold_chromo3)
    statistics_intermodes_chromo3 = cle.statistics_of_labelled_pixels(img_nuclei, chromo_intermodes_Labels3)

    chromointermodes_Area3 = np.sum(statistics_intermodes_chromo3['area'], axis=0)
    chromointermodes_Area3 = chromointermodes_Area3 * px_to_um_Y*px_to_um_X*px_to_um_Z
    # Export Chromo 1 Stats 
    pd.DataFrame(statistics_intermodes_chromo3).to_excel(Result_folder+filename+'(Chromo)Chromo_statistics_3.xlsx')

    # Max Entropy Chromocenter Segmentation -> Orginial Voxel -> Blur (rx=ry=2) -> Max(Entropy Threshold) 
    entropy_Threshold_chromo3 = nsitk.threshold_maximum_entropy(intensity_map_blur3)
    chromo_entropy_Labels3 = cle.connected_components_labeling_box(entropy_Threshold_chromo3)
    # Caculate Entropy Area
    statistics_entropy_chromo3 = cle.statistics_of_labelled_pixels(img_nuclei, chromo_entropy_Labels3)

    chromoentropy_Area3 = np.sum(statistics_entropy_chromo3['area'], axis=0)
    chromoentropy_Area3 = chromoentropy_Area2 * px_to_um_Y*px_to_um_X*px_to_um_Z

    #binary_fill_holes = binary_fill_holes(intermodes_Threshold_chromo3)
    #exclude_tiny_label = cle.exclude_small_labels(binary_fill_holes, maximum_size = 500)
    #exclude_tiny_label = cle.exclude_large_labels(binary_fill_holes, minimum_size = 40)

    # Check the nuclei
    viewer = napari.Viewer()
    #viewer.add_image(binary_fill_holes)
    viewer.add_image(intensity_map3_norm)
    viewer.add_image(chromo_entropy_Labels3)
    viewer.add_image(chromo_intermodes_Labels3)
