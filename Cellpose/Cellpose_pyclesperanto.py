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

model_match_threshold = 27 #30
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

 # Method running chromocenters segmentation.
 # Simple algorithms description :
 # 1- Gaussian blur
 # 2- Image gradient : for each voxels computing sum of difference with
 # X (_neigh parameter) neighborhood.
 # Thresholding value : keep voxels having value higher mean nucleus
 # intensity plus X (_factor parameter) standard deviation value
 # 4- Binarization of threshold image
 # 5- Connected component computation from binarized image
