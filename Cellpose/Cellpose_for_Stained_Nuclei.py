# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 11:54:07 2023
@author: kaabir
"""
import time, os, sys
import glob
import skimage.io 
from skimage.io import imread, imshow
from cellpose import models, utils, plot, io 
import matplotlib as mpl
import napari 
import pyclesperanto_prototype as cle
import numpy as np
import pandas as pd
from aicsimageio import AICSImage
import napari_simpleitk_image_processing as nsitk
from skimage.morphology import skeletonize, thin, skeletonize_3d
from scipy.ndimage import label

from pathlib import Path
#from astropy.visualization import simple_norm
from skimage import img_as_float32
from skimage.util import img_as_ubyte
from tifffile import imwrite
#import SimpleITK as sitk

# structure for actin segmentation
from skimage import measure
from scipy import ndimage
from scipy.ndimage import gaussian_filter
from skimage import morphology
from skimage.morphology import binary_dilation
from scipy import ndimage as ndi

device = cle.select_device("gfx1032")
#device = cle.select_device("AMD Radeon Pro W6600")

def replace_intensity(mask, img):
    if not (mask.shape == img.shape):
        return False

    replace_intensity = np.where(np.logical_and(mask,img),img,0) # Overlap Mask on original image (as reference)
    #chromocenter_intensity = np.array(cle.replace_intensities(mask,img))
    return replace_intensity

def normalize_intensity(k):
    #k_min = k.min(axis=(1, 2), keepdims=True)
    #k_max = k.max(axis=(1, 2), keepdims=True)
    k_min = k.min()
    k_max = k.max()

    k = (k - k_min)/(k_max-k_min)
    k[np.isnan(k)] = 0
    return k

def trim_array(arr, mask):
    bounding_box = tuple(
    slice(np.min(indexes), np.max(indexes) + 1)
    for indexes in np.where(mask))
    return arr[bounding_box]

def analyze_actin(mask, img_actin, filename, i):
    act_obj = np.zeros(img_nuclei.shape)
    dilated = ndi.binary_dilation(mask, diamond, iterations=10).astype(mask.dtype)
    actin_img = replace_intensity(dilated, img_actin)
    actin_filter = nsitk.median_filter(actin_img, radius_x=2, radius_y=2, radius_z=0)
    actin_binary = nsitk.threshold_otsu(actin_filter)

    act_obj = np.zeros(img_nuclei.shape)

    # Apply thinning function to each slice of the actin binary image
    for i in range(actin_binary.shape[0]):
        thinned_slice = thin(actin_binary[i])
        act_obj[i] = thinned_slice

    statistics_surf_actin = cle.statistics_of_labelled_pixels(actin_img, act_obj)

    actin_surf_Area = np.sum(statistics_surf_actin['area'], axis=0)
    actin_surf_Area = actin_surf_Area * px_to_um_X
    print('Actin_surf_Area ---', actin_surf_Area)
    
    # Write statistics to Excel file
    statistics_df = pd.DataFrame(statistics_surf_actin)
    statistics_df.to_excel(Result_folder + '(Actin)_' + filename + '_' + str(i) + '.xlsx')
    

def calculate_surface_area(mask, threshold=None):
    # generate surface mesh using marching cubes algorithm
    verts, faces, _, _ = measure.marching_cubes(mask)
        # calculate surface area using mesh surface area function
    surface_area = measure.mesh_surface_area(verts, faces)
    #surface_area = surface_area*px_to_um_Y*px_to_um_Y*px_to_um_Z

    return surface_area

# Read the Tif/CZI File
def folder_scan(directory):
    # Actin after 30min - value assigned is 3
    get_files = []
    extension = ".czi"  # ".czi"
    for f_name in os.listdir(directory):
        if f_name.endswith(extension):
            get_files.append(os.path.join(directory, f_name))
    return get_files


folder_path = 'C:/Users/kaabi/Documents/Nuceli_Data/Enucleation/230406 OF1 D1 Actin 10uMVCA Gactin RH3K9me2 FRH3K27me2/'

get_files = folder_scan(folder_path)

# Conversion to um
px_to_um_X = 0.0351435
px_to_um_Y = 0.0351435
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
        img_actin = aics_image.get_image_data("ZYX", T=0, C=0) # Actin Channel
        img_H3K27me3 = aics_image.get_image_data("ZYX", T=0, C=1) # yH2AX
        img_H3K9me2 = aics_image.get_image_data("ZYX", T=0, C=2) # LaminB1
        img_nuclei = aics_image.get_image_data("ZYX", T=0, C=3) # Nucleus
    
    if not os.path.exists(folder_path + '/Result'):
        os.makedirs(folder_path + '/Result')
        
    Result_folder = folder_path + '/Result/' 
    # For saving results #.stem

    # Upscale the Z dimension to be isotropic

    # channel to segment and nuclear channel 
    # numbering starts at 1 
    # for your single channel image use [0, 0] 
    # for the multi channel image it's [3, 0]
    channels = [0, 0] 
    diameter = 224.397 #169.708
    use_GPU = True
    stitch_threshold=1
    #do_3D=True
    #resample=True
    cellprob_threshold = 0

    pretrained_model = "C:/Users/kaabi/Documents/Nuceli_Data/Enucleation/Cellpose/Models_fixed_eNUC/CP_20230619_150738"

    model_match_threshold = 27 #30
    flow_threshold = (31.0 - model_match_threshold) / 10.

    logger = io.logger_setup()

    model = models.CellposeModel(gpu=use_GPU, model_type=pretrained_model)#,residual_on=False, style_on=False, concatenation=False)#, resample=True, do_3D=True,stitch_threshold=1)
    mask, flows,styles = model.eval(img_nuclei,
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

    # Merge Segmented Mask

    merged_Labels = cle.connected_components_labeling_box(mask)
    merged_Labels_np = np.array(merged_Labels).astype('int32')
    print('The Filename ---', filename)
    ### Save the Mask
    prediction_stack_32 = img_as_float32(merged_Labels_np, force_copy=False)     
    os.chdir(Result_folder)
    imwrite(filename+".tif", prediction_stack_32)
    #merged_Labels_np = split_touching_objects(merged_Labels_np)

    # Structure  for Actin Dilation
    diamond = np.zeros((3, 3, 3), dtype=bool)
    diamond[1:2, 0:3, 0:3] = True

    #print('Number_of_Nuclei',len(np.unique(merged_Labels))-1) first is mask

    from skimage import measure
    labels = measure.label(merged_Labels_np)


    for i, label in enumerate(np.unique(labels)[1:]):
        print('i ---', i)
        #print('label ---', label)
        if label in labels:
            print('label ---', label)
            # create a mask for the current label
            #if (len(np.unique(merged_Labels)-1) == 1):# and no_chnls>1:
                # break    
            mask = labels == label
            ###nuc_lbl = np.array(label)
            intensity_nucelus= replace_intensity(mask, img_nuclei)
            image1_gb = nsitk.gaussian_blur(intensity_nucelus, 2.0, 2.0, 0.0)
            # threshold isodata
            image2_T = nsitk.threshold_isodata(image1_gb)
            # connected component labeling
            image3_C = nsitk.connected_component_labeling(image2_T)
        
            statistics_nucleus = nsitk.label_statistics(intensity_image=intensity_nucelus, label_image=mask,
                                                size = True, intensity = True, perimeter= True,
                                                shape= True, position= True)#, moments= True)
        
            nuclei_Area = statistics_nucleus.iloc[0, 26]
            nuclei_Area = nuclei_Area*px_to_um_X*px_to_um_Y*px_to_um_Z
            statistics_nucleus['Nucleus Area'] = nuclei_Area 
            print('nuclei_Area ---', nuclei_Area)

            pd.DataFrame(statistics_nucleus).to_excel(Result_folder+'(Nucleus)_' + filename+'_'+ str(i)+'.xlsx')
        
           
            # # Chromocenter Segmentation

            # intensity_map_norm = normalize_intensity(intensity_nucelus)
    
            # intensity_map_blur = nsitk.median_filter(intensity_map_norm, radius_x=2, radius_y=2, radius_z=0)
            # intensity_map_thb = cle.top_hat_box(intensity_map_blur, None, 10.0, 10.0, 0.0)
            # intermodes_Threshold_Chromo = nsitk.threshold_intermodes(intensity_map_thb)

            # statistics_intermodes_chromo = nsitk.label_statistics(intensity_image=img_nuclei, label_image=intermodes_Threshold_Chromo,
            #                                     size = True, intensity = True, perimeter= True,
            #                                     shape= True, position= True)#, moments= True)
            
            # pd.DataFrame(statistics_intermodes_chromo).to_excel(Result_folder+'(Chromo)Chromo_statistics_'+filename+'_'+ str(i) +'.xlsx') 
            
            # intermodes_chromo_Area = statistics_intermodes_chromo.iloc[0, 26]
            
            # chromointermodes_Area = intermodes_chromo_Area *px_to_um_X* px_to_um_Y*px_to_um_Z
            # print('Chromointermodes_Area ---', chromointermodes_Area)
            # Actin Segmentation
            if 'img_actin' in globals(): #or 'img_actin' in globals():
                # Separate Ctrl
                pass
           
            act_obj = np.zeros(img_nuclei.shape)
            dilated = ndi.binary_dilation(mask, diamond, iterations=10).astype(mask.dtype)
            actin_img = replace_intensity(dilated, img_actin)
            actin_filter = nsitk.median_filter(actin_img, radius_x=2, radius_y=2, radius_z=0)
            actin_binary = nsitk.threshold_otsu(actin_filter)

            act_obj = np.zeros(img_nuclei.shape)

            # Apply thinning function to each slice of the actin binary image
            for i in range(actin_binary.shape[0]):
                    thinned_slice = thin(actin_binary[i])
                    act_obj[i] = thinned_slice
            
            statistics_actin = nsitk.label_statistics(intensity_image=actin_img, label_image=act_obj,
                                                size = True, intensity = True, perimeter= True,
                                                shape= True, position= True)#, moments= True)

            #actin_surf_Area = np.sum(statistics_surf_actin['area'], axis=0)
            actin_surf_Area = statistics_actin.iloc[0, 26]
            actin_surf_Area = actin_surf_Area * px_to_um_X
            statistics_actin['Actin Surface Area'] = actin_surf_Area
            print('Actin_surf_Area ---', actin_surf_Area) 
            # Write statistics to Excel file
            
            statistics_actin.to_excel(Result_folder + '(Actin)_' + filename + '_' + str(i) + '.xlsx')   
            
            # H3K27me3 Segmentation

            statistics_H3K27me3 = nsitk.label_statistics(intensity_image=img_H3K27me3, label_image=dilated,
                                                size = True, intensity = True, perimeter= True,
                                                shape= True, position= True)#, moments= True)
            
            statistics_H3K27me3.to_excel(Result_folder + '(H3K27me3)_' + filename + '_' + str(i) + '.xlsx')   
            H3K27me3_Intensity = statistics_H3K27me3.iloc[0, 2]
            print('H3K27me3_Intensity ---', H3K27me3_Intensity)
            
            # H3K9me2 Segmentation
            
            statistics_H3K9me2 = nsitk.label_statistics(intensity_image=img_H3K9me2, label_image=dilated,
                                    size = True, intensity = True, perimeter= True,
                                    shape= True, position= True)#, moments= True)
            statistics_H3K9me2.to_excel(Result_folder + '(H3K9me2)_' + filename + '_' + str(i) + '.xlsx')   

            H3K9me2_Intensity = statistics_H3K9me2.iloc[0, 2]
            print('H3K9me2_Intensity ---', H3K9me2_Intensity)
