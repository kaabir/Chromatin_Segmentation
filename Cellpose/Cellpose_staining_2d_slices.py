# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 18:22:51 2023

@author: kaabir
"""

import time, os, sys
import skimage.io 
from skimage.io import imread, imshow
from cellpose import models, utils, plot, io, core
import matplotlib as mpl
import napari 
import pyclesperanto_prototype as cle
import numpy as np
import pandas as pd
from aicsimageio import AICSImage
import napari_simpleitk_image_processing as nsitk

from pathlib import Path
from astropy.visualization import simple_norm
from skimage import img_as_float32
from skimage.util import img_as_ubyte
from tifffile import imsave
from napari_segment_blobs_and_things_with_membranes import split_touching_objects
# model_type='cyto' or model_type='nuclei'
#model = models.Cellpose(model_type='cyto', gpu = True) #nuclei model for brdu
# Read the Tif File
#folder_path = 'C:/Users/NBMS1/Documents/Nuclei Data/Cellpose/Test_Dataset/15-12-22/T30/'
folder_path = 'C:/Users/ChimieENS/Documents/Cellpose/H3K_Staining/'

file_path = '230413 OF1 D1 Ctrl 10uMVCA Gactin RH3K9me2 FRH3K27me2 01-1.tif'
img_path = folder_path + file_path

aics_image = AICSImage(img_path)
img = aics_image.get_image_data("YX", T=0)
#img = imread(img_path)

folder_path_red = '230413 OF1 D1 Ctrl 10uMVCA FRH3K27me2 03.tif'
img_path_red = folder_path + folder_path_red 
#img_red = imread(img_path)
aics_image_red = AICSImage(img_path_red)
img_red = aics_image_red.get_image_data("YX", T=0)

# In case two channels
# img_actin = img[:,0,:,:] # A Channel
# img_nuclei = img[:,1,:,:] # N Channel

# In case three channels
# img_actin = img[:,:,:,0] # A Channel
# img_dextran = img[:,:,:,1] # Dextran Channel
# img_nuclei = img[:,:,:,2] # N Channel
#img = img[:,:,:,-1] # # Hoescht Channel

# get_Zstack = img.shape[0]
# get_Xresl = img.shape[1]
# get_yresl = img.shape[2]

# In case four channels

# # Read the CZI File
# img_path = AICSImage("C:/Users/kaabi/Documents/Nuceli_Data/Enucleation/Cellpose/Test_Dataset/Files/02.czi")
# # Hoescht Channel
# img = img_path.get_image_data("ZYX", C=2, S=0, T=0)
# get_Zstack = img.shape[0]

# For saving results
filename = Path(file_path).parts[-1]
#.stem
filename = Path(file_path).stem
Result_folder = folder_path + '/Result/' 
# Upscale the Z dimension to be isotropic
#anisotropy = 2*get_Zstack

# channel to segment and nuclear channel 
# numbering starts at 1 
# for your single channel image use [0, 0] 
# for the multi channel image it's [3, 0]
channels = [0, 0] 
diameter = 30 #169.708
use_GPU = True
cellprob_threshold = 0

model_match_threshold = 27#30
flow_threshold = (31.0 - model_match_threshold) / 10.0

logger = io.logger_setup()

model = models.CellposeModel(gpu=use_GPU, model_type='cyto')#,residual_on=False, style_on=False, concatenation=False)#, resample=True, do_3D=True,stitch_threshold=1)
mask, flows, styles = model.eval(img, 
                                            channels=channels, 
                                            diameter=diameter, 
                                            #anisotropy=anisotropy, 
                                            do_3D=False, 
                                            net_avg=True, 
                                            resample=True,
                                            #augment=False, 
                                            cellprob_threshold=0)

merged_Labels_np = np.array(mask)
#merged_Labels_np = split_touching_objects(merged_Labels_np)

    # Check the nuclei
viewer = napari.Viewer()
viewer.add_image(img)
viewer.add_image(merged_Labels_np)

# Save Labels
prediction_stack_32 = img_as_float32(mask, force_copy=False)     
os.chdir(Result_folder)
#imsave(str(filename)+".tif", prediction_stack_32)
# Getting Original intensities
#intensity_vector = cle.read_intensities_from_map(mask, img)
##statistics_hoescht = cle.statistics_of_labelled_pixels(img, merged_Labels_np)
statistics_hoescht = nsitk.label_statistics(img, merged_Labels_np)

statistics_h39me2 = cle.statistics_of_labelled_pixels(img_red, merged_Labels_np)
#
# Conversion to um
px_to_um_X = 0.0351435
px_to_um_Y = 0.0351435
px_to_um_Z = 1
#nuclei_Area = statistics['area']*px_to_um_Y*px_to_um_Y*px_to_um_Z

# Chromocenter Segmentation
# Get voxel of pixel intensity values inside the mask 
def intensity_vector(mask, img):
    if not (mask.shape == img.shape):
        return False
    
    mat_intensity = np.where(np.logical_and(mask,img),img,0) # Overlap Mask on original image (as reference)
    return mat_intensity

intensity_map= intensity_vector(mask, img)

def normalize_intensity(k):
    k_min = k.min(axis=(1, 2), keepdims=True)
    k_max = k.max(axis=(1, 2), keepdims=True)

    k = (k - k_min)/(k_max-k_min)
    k[np.isnan(k)] = 0
    return k

from scipy import ndimage
from scipy.ndimage import gaussian_filter


# Create a matrice to place the single nuleus with intensity voxel
im_obj1 = np.zeros(img.shape) 
im_obj2 = np.zeros(img.shape) 
im_obj3 = np.zeros(img.shape)
nuc_lbl_lst1 = []
nuc_lbl_lst2 = []
nuc_lbl_lst3 = []

