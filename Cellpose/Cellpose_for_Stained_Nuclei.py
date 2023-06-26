# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 18:22:51 2023

@author: kaabi
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
from aicsimageio.writers import OmeTiffWriter
import napari_simpleitk_image_processing as nsitk

from pathlib import Path
from astropy.visualization import simple_norm
from skimage import img_as_float32
from skimage.util import img_as_ubyte
from tifffile import imsave
from napari_segment_blobs_and_things_with_membranes import split_touching_objects

device = cle.select_device("AMD Radeon Pro W6600")

from scipy import ndimage as ndi
# Structure  for Actin Dilation
diamond = np.zeros((20, 20), dtype=bool)
diamond[1:19, 1:19] = True

# model_type='cyto' or model_type='nuclei'
#model = models.Cellpose(model_type='cyto', gpu = True) #nuclei model for brdu
# # Read the Tif File

# folder_path = 'E:/Quantified/Agarose Compression/20221103 IF OF1 D1 H3K9me3 FOP FOXJ1/'

# file_path = 'Binned_z_3d-20221103_IF_OF1 D1 H3K27me3 FOP FOXJ1_1H_Ctrl_S2 _01.tif'
# #3d-20221103 IF OF1 D1 H3K27me3 FOP FOXJ1 -1H Ctrl -01-1
# img_path = folder_path + file_path


def trim_array(arr):
    """Returns a trimmed view of an n-D array excluding any outer
    regions which contain only zeros.
    https://stackoverflow.com/questions/55917328/numpy-trim-zeros-in-2d-or-3d
    """
    slices = tuple(slice(idx.min(), idx.max() + 1) for idx in np.nonzero(arr))
    return arr[slices]

def mid_Slice(img_stack):
    # Find the halfway point of the array along the first axis
    mid_slice = img_stack.shape[0] // 2

    # Calculate the start and end indices for slicing
    start_index = mid_slice - 2
    end_index = mid_slice + 2

    # Slicing the array
    img_stack = img_stack[start_index:end_index+1, :, :]
    
    return img_stack

# Read the Tif/CZI File
def folder_scan(directory):
    # Actin after 30min - value assigned is 3
    get_files = []
    extension = ".czi"  # ".czi"
    for f_name in os.listdir(directory):
        if f_name.endswith(extension):
            get_files.append(os.path.join(directory, f_name))
    return get_files

def normalize_intensity(k):
    k_min = k.min()
    k_max = k.max()

    k = (k - k_min)/(k_max-k_min)
    k[np.isnan(k)] = 0
    return k

def pad_with(vector, pad_width, iaxis, kwargs):
    pad_value = kwargs.get('padder', 0)
    vector[:pad_width[0]] = pad_value
    vector[-pad_width[1]:] = pad_value
    
# Get voxel of pixel intensity values inside the mask 
def replace_intensity(mask, img):
    if not (mask.shape == img.shape):
        return False
    
    mat_intensity = np.where(np.logical_and(mask,img),img,0) # Overlap Mask on original image (as reference)
    return mat_intensity

folder_path = 'E:/Quantified/Agarose Compression/20221103 IF OF1 D1 H3K9me3 FOP FOXJ1/'

get_files = folder_scan(folder_path)

# Conversion to um
px_to_um_X = 0.1098235
px_to_um_Y = 0.1098235
px_to_um_Z = 0.32 #1.0/0.0351435 # making pixels isotropic


for image in get_files:
    aics_image = AICSImage(image) # Read indivi
    filename = Path(image).stem
    # <Dimensions [T: 1, C: 1, Z: 256, Y: 256, X: 256]> # ("ZYX", T=0)
    no_chnls = aics_image.dims.C
    if no_chnls == 1:
        img_nuclei = aics_image.get_image_data("ZYX", T=0, C=0) 
    elif no_chnls == 2: 
        img_actin = aics_image.get_image_data("ZYX", T=0, C=0) # A Channel
        img_nuclei = aics_image.get_image_data("ZYX", T=0, C=1)
    elif no_chnls == 3:
        img_actin = aics_image.get_image_data("ZYX", T=0, C=0) # A Channel
        img_dextran_read = aics_image.get_image_data("ZYX", T=0, C=1) # Dextran Channel
        img_nuclei = aics_image.get_image_data("ZYX", T=0, C=2)
    else:
        img_fop = aics_image.get_image_data("ZYX", T=0, C=0) # Actin Channel 
        img_h3k27me3 = aics_image.get_image_data("ZYX", T=0, C=1) # yH2AX
        img_foxj = aics_image.get_image_data("ZYX", T=0, C=2) # LaminB1
        img_nuclei = aics_image.get_image_data("ZYX", T=0, C=3) # Nucleus

        img_nuclei = np.max(img_nuclei, axis=0)

    if 'img_fop' in globals():
        img_fop = np.max(img_fop, axis=0)

    if 'img_h3k27me3' in globals():
        img_h3k27me3 = np.max(img_h3k27me3, axis=0)
    
    if 'img_foxj' in globals():
        img_foxj = np.max(img_foxj, axis=0)   

    if not os.path.exists(folder_path + '/Result'):
        os.makedirs(folder_path + '/Result')
        
    Result_folder = folder_path + '/Result/' 


# In case four channels

# # Read the CZI File
# img_path = AICSImage("C:/Users/kaabi/Documents/Nuceli_Data/Enucleation/Cellpose/Test_Dataset/Files/02.czi")
# # Hoescht Channel
# img = img_path.get_image_data("ZYX", C=2, S=0, T=0)
# get_Zstack = img.shape[0]


# For saving results
#.stem
# Upscale the Z dimension to be isotropic
#anisotropy = 2*get_Zstack

# channel to segment and nuclear channel 
# numbering starts at 1 
# for your single channel image use [0, 0] 
# for the multi channel image it's [3, 0]
    channels = [0, 0] 
    diameter = 90 #169.708
    use_GPU = True
    cellprob_threshold = 0

    stitch_threshold = 0

    model_match_threshold = 27#30
    flow_threshold = (31.0 - model_match_threshold) / 10.0

    print('Segmentation begins')

    logger = io.logger_setup()

    pretrained_model = "C:/Users/kaabi/Documents/Nuceli_Data/Enucleation/Cellpose/Models_fixed_eNUC/bulk/CP_LC3_Bulk_nuclei_Main"

    model = models.CellposeModel(gpu=use_GPU, model_type=pretrained_model)#,residual_on=False, style_on=False, concatenation=False)#, resample=True, do_3D=True,stitch_threshold=1)
    mask, flows, styles = model.eval(img_nuclei, 
                                channels =channels,
                                #anisotropy=anisotropy,
                                diameter=diameter, 
                                #pretrained_model=pretrained_model,
                                do_3D=False,
                                resample=True,
                                min_size = 10, #-1
                                flow_threshold =flow_threshold,
                                cellprob_threshold = cellprob_threshold,
                                net_avg=True,
                                stitch_threshold=stitch_threshold,)

#mask_np1 = nsitk.morphological_watershed(mask)
#mask_np = np.array(mask)

    clean_mask = cle.exclude_labels_on_edges(mask)

    from skimage.segmentation import clear_border
    from scipy.ndimage import label

#merged_Labels_np = split_touching_objects(merged_Labels_np)

#composite = np.dstack((img_fop, img_h3k27me3, img_foxj, img_nuclei))


# # Perform Z-projection on each channel
# projection_method = "max"  # Choose the desired projection method ("max", "min", "mean", "sum", etc.)
# composite = []
# for channel_index in range(merged_Labels_np):
#     z_projected_channel = np.max(channel_index, axis=0)  # Replace with the desired projection method
#     composite.append(z_projected_channel)

# Save Labels
#prediction_stack_32 = img_as_float32(mask, force_copy=False)     
#os.chdir(Result_folder)
#imsave(str(filename)+".tif", prediction_stack_32)
# Getting Original intensities
#intensity_vector = cle.read_intensities_from_map(mask, img)
##statistics_hoescht = cle.statistics_of_labelled_pixels(img, merged_Labels_np)
#
# Conversion to um
    px_to_um_X = 0.1098235
    px_to_um_Y = 0.1098235
#px_to_um_Z = 1 # check for bulk
#nuclei_Area = statistics['area']*px_to_um_Y*px_to_um_Y*px_to_um_Z

# Chromocenter Segmentation


    from skimage import measure
    labels = measure.label(clean_mask)

    for i, label in enumerate(np.unique(labels)[1:]):
        print('i ---', i)
        #print('label ---', label)
        if label in labels: # if label in labels:
            print('label ---', label)
            # create a mask for the current label
            #if (len(np.unique(merged_Labels)-1) == 1):# and no_chnls>1:
                # break    
            mask_unique = labels == label
            #nuc_lbl = np.array(label)
        #label_crop = trim_array(mask)

        
        # Dilation of unqiue label
            dilated_mask = ndi.binary_dilation(mask_unique, structure = diamond, iterations=2).astype(mask_unique.dtype)
        # FoxJ Channel
            intensity_foxj = replace_intensity(dilated_mask, img_foxj)
            #foxj_arr = trim_array(intensity_foxj)
        
            statistics_Foxj = nsitk.label_statistics(intensity_image=intensity_foxj, label_image=dilated_mask,
                                            size = True, intensity = True, perimeter= True,
                                            shape= True, position= True)#, moments= True)
            pd.DataFrame(statistics_Foxj).to_excel(Result_folder+'(FoxJ)_' + filename +'_' + str(i)+'.xlsx')


        #print('nucelus_arr ---', nucelus_arr.shape)
        # Nucleus Stats
            statistics_nucleus = nsitk.label_statistics(intensity_image=img_nuclei, label_image=mask_unique,
                                            size = True, intensity = True, perimeter= True,
                                            shape= True, position= True)#, moments= True)
            pd.DataFrame(statistics_nucleus).to_excel(Result_folder+'(Nucleus)_' + filename +'_' + str(i)+'.xlsx')

        # Chromocenter Dots
            normal_nucelus= replace_intensity(mask_unique, img_nuclei)
        
            #padding_nucleus_arr = trim_array(normal_nucelus)
        # Get the shape difference between the two arrays
            #shape_diff = np.subtract(foxj_arr.shape,padding_nucleus_arr.shape)
            #shape_diff = int(shape_diff[0]/2)
        
            intensity_map_norm = normalize_intensity(normal_nucelus)
            intensity_map_blur = nsitk.median_filter(intensity_map_norm, radius_x=2, radius_y=2, radius_z=0)
            intensity_map_thb = cle.top_hat_box(intensity_map_blur, None, 10.0, 10.0, 0.0)
            
            try:
                intermodes_Threshold_Chromo = nsitk.threshold_intermodes(intensity_map_thb)
            except RuntimeError:
                print("No Chromocenter Found")
            
            
            # chromointermodes_Area = intermodes_chromo_Area *px_to_um_X* px_to_um_Y*px_to_um_Z
            #
        
            statistics_Chromo = nsitk.label_statistics(intensity_image=intermodes_Threshold_Chromo, label_image=mask_unique,
                                            size = True, intensity = True, perimeter= True,
                                            shape= True, position= True)#, moments= True)
        
                        
            intermodes_chromo_Area = statistics_Chromo.loc[0, 'number_of_pixels']
            chromointermodes_Area = intermodes_chromo_Area *px_to_um_X* px_to_um_Y
            statistics_Chromo['Chromocenter Area'] = chromointermodes_Area 
            print('Chromocenter_Area ---', chromointermodes_Area)
            pd.DataFrame(statistics_Chromo).to_excel(Result_folder+'(Chromo)_' + filename + '_' + str(i)+'.xlsx')

        
        # Threshold isodata
        # image2_T = nsitk.threshold_isodata(image1_gb)
        # image3_T = nsitk.threshold_maximum_entropy(image1_gb)
        
        # Calculate the padding values for each dimension
        #padding = [(0, diff) for diff in shape_diff]  # Pad with zeros

        # Pad array1 to match the shape of array2
        
        # Normal nucleus
        
            #nucelus_arr = np.pad(padding_nucleus_arr, shape_diff, pad_with)
        
        # H3k27me3 Channel
            intensity_h3k27me3 = replace_intensity(mask_unique, img_h3k27me3)
            #h3k27me3_trim = trim_array(intensity_h3k27me3)
            #h3k27me3_arr = np.pad(h3k27me3_trim, shape_diff, pad_with)
        
            statistics_h3k27me3 = nsitk.label_statistics(intensity_image=intensity_h3k27me3, label_image=mask_unique,
                                            size = True, intensity = True, perimeter= True,
                                            shape= True, position= True)#, moments= True)
            pd.DataFrame(statistics_h3k27me3).to_excel(Result_folder+'(H3K27me3)_' + filename+'_'+ str(i)+'.xlsx')
        
        #h3k27me3_arr = np.pad(intensity_h3k27me3, padding,  constant_values=0, mode='constant')
        #intensity_h3k27me3 = replace_intensity(dilated_mask, img_h3k27me3)
        #h3k27me3_trim = trim_array(intensity_h3k27me3)
        #h3k27me3_arr = ndi.binary_dilation(h3k27me3_trim, structure = diamond, iterations=2).astype(mask_unique.dtype)
        
        # FOP Channel
            intensity_fop = replace_intensity(mask_unique, img_fop)
            #fop_trim = trim_array(intensity_fop)
            #fop_arr = np.pad(fop_trim, shape_diff, pad_with)
        
            statistics_Fop = nsitk.label_statistics(intensity_image=intensity_fop, label_image=mask_unique,
                                            size = True, intensity = True, perimeter= True,
                                            shape= True, position= True)#, moments= True)
            pd.DataFrame(statistics_Fop).to_excel(Result_folder+'(FOP)_' + filename+'_'+ str(i)+'.xlsx')

        
        #fop_trim = trim_array(intensity_fop)
        #fop_arr = ndi.binary_dilation(fop_trim, structure = diamond, iterations=2).astype(mask_unique.dtype)
        
        # Green with Far Red and then Red and Blue channel
            merged_image = np.stack([intensity_fop, intensity_h3k27me3, intensity_foxj, normal_nucelus], axis=0)
            os.chdir(Result_folder)
            OmeTiffWriter.save(merged_image,filename+'_'+ str(i) + ".tif")
        
        #print(image3_T.shape)
        #viewer = napari.Viewer()
        # ####viewer.add_image(merged_image)
        # ####viewer.add_image(image2_T)
        # ####viewer.add_image(image3_T)
        #viewer.add_image(intermodes_Threshold_Chromo)
        # viewer.add_image(img_nuclei)
        # viewer.add_image(mask)
        # viewer.add_image(clean_mask)
        # viewer.add_image(mask_unique)
        # viewer.add_image(intensity_h3k27me3)
        # viewer.add_image(intensity_foxj)

