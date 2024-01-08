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
from aicsimageio.writers import OmeTiffWriter
import napari_simpleitk_image_processing as nsitk

from pathlib import Path
from astropy.visualization import simple_norm
from skimage import img_as_float32
from skimage.util import img_as_ubyte
from tifffile import imsave
from napari_segment_blobs_and_things_with_membranes import split_touching_objects
import napari_segment_blobs_and_things_with_membranes as nsbatwm
from skimage.measure import regionprops

device = cle.select_device("AMD Radeon Pro W6600")

from scipy import ndimage as ndi
# Structure  for Actin Dilation
diamond = np.zeros((20, 20), dtype=bool)
diamond[1:19, 1:19] = True


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

from scipy.ndimage import label

def connected_components_labeling_box(binary_input):
    labeled_array, num_labels = label(binary_input)
    return labeled_array

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

print(" Folder List Stage")

folder_list = [
    'E:/Quantified/Agarose Compression/231219 Compression OF1 Volume Fraction/03-Ctrl_20ul/Agarose/',
    'E:/Quantified/Agarose Compression/231219 Compression OF1 Volume Fraction/03-Ctrl_20ul/Ctrl/',
]

for folder_name in folder_list:
    print(folder_name)
    get_files = folder_scan(folder_name)
    print("Scanning Files in the folder >>>>>>>", folder_name)

    #folder_path = 'E:/Quantified/Agarose Compression/231109 e230926 OF1 D6 Tf shNeg shArp3 Pad2percent 3R Gfop Rdsred FRfoxj1/'

    get_files = folder_scan(folder_name)

    # Conversion to um
    px_to_um_X = 0.1098235
    px_to_um_Y = 0.1098235
    px_to_um_Z = 0.3054316 #1.0/0.0351435 # making pixels isotropic


    for image in get_files:
        # # Read the CZI File
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
            img_fop = aics_image.get_image_data("ZYX", T=0, C=0) # Actin Channel   # Green with Far Red and then Red and Blue channel
            img_foxj = aics_image.get_image_data("ZYX", T=0, C=1) # yH2AX
            img_phosRB = aics_image.get_image_data("ZYX", T=0, C=2) # LaminB1
            img_nuclei = aics_image.get_image_data("ZYX", T=0, C=3) # Nucleus  # Hoescht Channel

            img_nuclei = np.max(img_nuclei, axis=0)

        if 'img_fop' in globals():
            img_fop = np.max(img_fop, axis=0)

        if 'img_phosRB' in globals():
            img_phosRB = np.max(img_phosRB, axis=0)
    
        if 'img_foxj' in globals():
            img_foxj = np.max(img_foxj, axis=0)   

        if not os.path.exists(folder_name + '/Result'):
            os.makedirs(folder_name + '/Result')
        
        Result_folder = folder_name + '/Result/' 

        channels = [0, 0] 
        diameter = 90 #169.708
        use_GPU = True
        cellprob_threshold = 0

        stitch_threshold = 0

        model_match_threshold = 26#30
        flow_threshold = (31.0 - model_match_threshold) / 10.0

        print('Segmentation Begins')

        logger = io.logger_setup()
        print('Segmenting Image --->', filename)
    
        pretrained_model = "C:/Users/kaabi/Documents/Nuceli_Data/Enucleation_Quantification/Cellpose/Models_fixed_eNUC/bulk/CP_CPX_Bulk_nuclei_Main_20231130"

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

        clean_mask = cle.exclude_labels_on_edges(mask)
        
        from skimage.segmentation import clear_border
        from scipy.ndimage import label

        # Chromocenter Segmentation

        from skimage import measure
        from skimage.filters import threshold_multiotsu
        labels = measure.label(clean_mask)

        for i, mask_label in enumerate(np.unique(labels)[1:]):
            print('i ---', i)
            #print('label ---', label)
            if mask_label in labels:  # This check is not necessary
                print('mask_label ---', mask_label)
            # create a mask for the current label
                mask_unique = labels == mask_label
            # Continue with your processing using the mask_label            
            # Dilation of unqiue label
                dilated_mask = ndi.binary_dilation(mask_unique, structure = diamond, iterations=3).astype(mask_unique.dtype)
                # FoxJ Channel
                intensity_Fop = replace_intensity(dilated_mask, img_fop)
        
                statistics_Fop = nsitk.label_statistics(intensity_image=intensity_Fop, label_image=dilated_mask,
                                            size = True, intensity = True, perimeter= True,
                                            shape= True, position= True)#, moments= True)
                pd.DataFrame(statistics_Fop).to_excel(Result_folder+'(FOP)_' + filename +'_' + str(i)+'.xlsx')


            #print('nucelus_arr ---', nucelus_arr.shape)
            # Nucleus Stats
                statistics_nucleus = nsitk.label_statistics(intensity_image=img_nuclei, label_image=mask_unique,
                                            size = True, intensity = True, perimeter= True,
                                            shape= True, position= True)#, moments= True)
            
                nucleus_Area = statistics_nucleus.loc[0, 'number_of_pixels']
               
                nucleus_Area = nucleus_Area *px_to_um_X* px_to_um_Y
                statistics_nucleus['Nucleus Area'] = nucleus_Area 
            
                pd.DataFrame(statistics_nucleus).to_excel(Result_folder+'(Nucleus)_' + filename +'_' + str(i)+'.xlsx')

        
                normal_nucelus= replace_intensity(mask_unique, img_nuclei)
        
            #padding_nucleus_arr = trim_array(normal_nucelus)
            # Get the shape difference between the two arrays
            #shape_diff = np.subtract(foxj_arr.shape,padding_nucleus_arr.shape)
            #shape_diff = int(shape_diff[0]/2)
            # Napari Viewer
            # viewer = napari.Viewer()
            # viewer.add_image(normal_nucelus)
            
            # Nucleus_inten = normalize_intensity(normal_nucelus)
                nucleus_thf = cle.top_hat_sphere(normal_nucelus, None, 10.0, 10.0, 0.0)
            
                nucleus_median = nsitk.median_filter(nucleus_thf, radius_x=2, radius_y=2, radius_z=0)
            
        
            # H3k27me3 Channel
                intensity_phosRB = replace_intensity(mask_unique, img_phosRB)
            #h3k27me3_trim = trim_array(intensity_h3k27me3)
            #h3k27me3_arr = np.pad(h3k27me3_trim, shape_diff, pad_with)
        
                statistics_phosRB = nsitk.label_statistics(intensity_image=intensity_phosRB, label_image=mask_unique,
                                            size = True, intensity = True, perimeter= True,
                                            shape= True, position= True)#, moments= True)
                pd.DataFrame(statistics_phosRB).to_excel(Result_folder+'(PhosRB)_' + filename+'_'+ str(i)+'.xlsx')
        
        #h3k27me3_arr = np.pad(intensity_h3k27me3, padding,  constant_values=0, mode='constant')
        #intensity_h3k27me3 = replace_intensity(dilated_mask, img_h3k27me3)
        #h3k27me3_trim = trim_array(intensity_h3k27me3)
        #h3k27me3_arr = ndi.binary_dilation(h3k27me3_trim, structure = diamond, iterations=2).astype(mask_unique.dtype)
        
        # FOP Channel
                intensity_FoxJ = replace_intensity(mask_unique, img_foxj)
            #fop_trim = trim_array(intensity_fop)
            #fop_arr = np.pad(fop_trim, shape_diff, pad_with)
        
                statistics_FoxJ = nsitk.label_statistics(intensity_image=intensity_FoxJ, label_image=mask_unique,
                                            size = True, intensity = True, perimeter= True,
                                            shape= True, position= True)#, moments= True)
                pd.DataFrame(statistics_FoxJ).to_excel(Result_folder+'(FOXJ)_' + filename+'_'+ str(i)+'.xlsx')

        
        #fop_trim = trim_array(intensity_fop)
        #fop_arr = ndi.binary_dilation(fop_trim, structure = diamond, iterations=2).astype(mask_unique.dtype)
        
       
                # Chromocenter Dots
                try:
                    chrom_MultiOTSU = threshold_multiotsu(nucleus_median)  
                    ChromoSeg = nucleus_median > chrom_MultiOTSU[1]
                    intensity_ChromoSeg = replace_intensity(mask_unique, ChromoSeg)
                #os.chdir(Result_folder)
                #OmeTiffWriter.save(intensity_ChromoSeg,'(ChormoTif)_' + filename+'_'+ str(i) + ".tif")
                
                # viewer.add_image(ChromoSeg)
                    labeled_slice = measure.label(ChromoSeg)
                    props = regionprops(labeled_slice)
                    num_bright_spots = len(props)
                    print('Chromocenter Dots Count ---', num_bright_spots)
                    statistics_Chromo = nsitk.label_statistics(intensity_image=img_nuclei, label_image=ChromoSeg,
                                                    size = True, intensity = True, perimeter= True,
                                                    shape= True, position= True)#, moments= True)
                                               
                    intermodes_chromo_Area = statistics_Chromo.loc[0, 'number_of_pixels']
                    chromointermodes_Area = intermodes_chromo_Area *px_to_um_X* px_to_um_Y
                    statistics_Chromo['Chromocenter Area'] = chromointermodes_Area 
                    statistics_Chromo['Chromocenter Spots'] = num_bright_spots
                
                    print('Chromocenter_Area ---', chromointermodes_Area)
                    pd.DataFrame(statistics_Chromo).to_excel(Result_folder+'(Chromo)_' + filename + '_' + str(i)+'.xlsx')

                    merged_image = np.stack([intensity_phosRB, intensity_FoxJ,intensity_Fop, normal_nucelus, intensity_ChromoSeg], axis=0)
                    os.chdir(Result_folder)
                    OmeTiffWriter.save(merged_image,filename+'_'+ str(i) + ".tif")                    
                
                except ValueError:               
                    print("No Chromocenter Found")
                
                    merged_image = np.stack([intensity_phosRB, intensity_FoxJ,intensity_Fop, normal_nucelus], axis=0)
                    os.chdir(Result_folder)
                    OmeTiffWriter.save(merged_image,filename+'_'+ str(i) + ".tif")
        
