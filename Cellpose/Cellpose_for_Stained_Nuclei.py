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
        img_h3k27me3 = aics_image.get_image_data("ZYX", T=0, C=0) # Actin Channel 
        img_foxj = aics_image.get_image_data("ZYX", T=0, C=1) # yH2AX
        img_fop = aics_image.get_image_data("ZYX", T=0, C=2) # LaminB1
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
                
                statistics_Chromo = nsitk.label_statistics(intensity_image=img_nuclei, label_image=intermodes_Threshold_Chromo,
                                                size = True, intensity = True, perimeter= True,
                                                shape= True, position= True)#, moments= True)
            
                            
                intermodes_chromo_Area = statistics_Chromo.loc[0, 'number_of_pixels']
                chromointermodes_Area = intermodes_chromo_Area *px_to_um_X* px_to_um_Y
                statistics_Chromo['Chromocenter Area'] = chromointermodes_Area 
                print('Chromocenter_Area ---', chromointermodes_Area)
                pd.DataFrame(statistics_Chromo).to_excel(Result_folder+'(Chromo)_' + filename + '_' + str(i)+'.xlsx')

            except RuntimeError:
                print("No Chromocenter Found")
            
            
            # chromointermodes_Area = intermodes_chromo_Area *px_to_um_X* px_to_um_Y*px_to_um_Z
            #
        
            
        
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
            merged_image = np.stack([intensity_h3k27me3, intensity_foxj,intensity_fop, normal_nucelus], axis=0)
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

#################################
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
#from scipy.ndimage import label

from pathlib import Path
#from astropy.visualization import simple_norm
from skimage import img_as_float32
from skimage.util import img_as_ubyte
from tifffile import imwrite
#import SimpleITK as sitk

# structure for actin segmentation
from skimage import measure
from scipy.ndimage import gaussian_filter
from skimage import morphology
from skimage.morphology import binary_dilation
from scipy import ndimage
from skimage import exposure
from scipy.ndimage import binary_fill_holes
from skimage.morphology import closing, square
from scipy.ndimage import binary_closing
from skimage import measure
from scipy import ndimage as ndi
import napari_segment_blobs_and_things_with_membranes as nsbatwm
from scipy.ndimage import label, generate_binary_structure
#from topostats.tracing.tracingfuncs import getSkeleton

# Open a file for writing the log
#log_file = open('output_log.txt', 'w')

# Redirect the standard output to the log file
#sys.stdout = log_file

#device = cle.select_device("gfx1032")
device = cle.select_device("AMD Radeon Pro W6600")

def pruneSkeleton(inputImage, opImgSize):
    """Function to remove the hanging branches from the skeletons - these
    are a persistent problem in the overall tracing process."""
    
    number_of_branches = 0
    thinned_slice = thin(inputImage)
    coordinates = np.argwhere(thinned_slice == 1).tolist()
    #print('coordinates >>>>', coordinates)

    # The branches are typically short so if a branch is longer than a quarter
    # of the total points, it's assumed to be part of the real data
    length_of_trace = len(coordinates)
    max_branch_length = int(length_of_trace * 0.05)

    # first check to find all the end coordinates in the trace
    potential_branch_ends = []

    # Most of the branch ends are just points with one neighbor
    for x, y in coordinates:
        if countNeighbours(x, y, coordinates) == 1:
            potential_branch_ends.append([x, y])

    # Now check if it's a branch and if it is, delete it
    for x_b, y_b in potential_branch_ends:
        branch_coordinates = [[x_b, y_b]]
        branch_continues = True
        temp_coordinates = coordinates[:]
        temp_coordinates.pop(temp_coordinates.index([x_b, y_b]))

        is_branch = False

        while branch_continues:
            no_of_neighbours, neighbours = countandGetNeighbours(x_b, y_b, temp_coordinates)

            # If branch continues
            if no_of_neighbours == 1:
                x_b, y_b = neighbours[0]
                branch_coordinates.append([x_b, y_b])
                temp_coordinates.pop(temp_coordinates.index([x_b, y_b]))

            # If the branch reaches the edge of the main trace
            elif no_of_neighbours > 1:
                branch_coordinates.pop(branch_coordinates.index([x_b, y_b]))
                branch_continues = False
                is_branch = True
            # Weird case that happens sometimes
            elif no_of_neighbours == 0:
                is_branch = True
                branch_continues = False

            if len(branch_coordinates) > max_branch_length:
                branch_continues = False
                is_branch = False

        if is_branch:
            number_of_branches += 1
            for x, y in branch_coordinates:
                thinned_slice[x, y] = 0

    #remaining_coordinates = np.argwhere(thinned_slice)
    #purged = cle.exclude_small_labels(thinned_slice, maximum_size=25)
    # Create an image of zeros with the same shape as the input image
    result_image = np.zeros_like(opImgSize)
    result_image[np.where(thinned_slice)] = 1

    # Put the remaining coordinates back into the result image
    # for x, y in remaining_coordinates:
    #     result_image[x, y] = 1
        
    #print('resulting Image size >>>>', result_image.size)
    #print('remaining_coordinates >>>>', remaining_coordinates)
    return result_image



def countandGetNeighbours(x, y, trace_coordinates):
    """Returns the number of neighboring points for a coordinate and an
    array containing those points"""

    neighbour_array = []
    number_of_neighbours = 0
    if [x, y + 1] in trace_coordinates:
        neighbour_array.append([x, y + 1])
        number_of_neighbours += 1
    if [x + 1, y + 1] in trace_coordinates:
        neighbour_array.append([x + 1, y + 1])
        number_of_neighbours += 1
    if [x + 1, y] in trace_coordinates:
        neighbour_array.append([x + 1, y])
        number_of_neighbours += 1
    if [x + 1, y - 1] in trace_coordinates:
        neighbour_array.append([x + 1, y - 1])
        number_of_neighbours += 1
    if [x, y - 1] in trace_coordinates:
        neighbour_array.append([x, y - 1])
        number_of_neighbours += 1
    if [x - 1, y - 1] in trace_coordinates:
        neighbour_array.append([x - 1, y - 1])
        number_of_neighbours += 1
    if [x - 1, y] in trace_coordinates:
        neighbour_array.append([x - 1, y])
        number_of_neighbours += 1
    if [x - 1, y + 1] in trace_coordinates:
        neighbour_array.append([x - 1, y + 1])
        number_of_neighbours += 1
    return number_of_neighbours, neighbour_array


def countNeighbours(x, y, trace_coordinates):
    """Counts the number of neighboring points for a given coordinate in
    a list of points"""

    number_of_neighbours = 0
    if [x, y + 1] in trace_coordinates:
        number_of_neighbours += 1
    if [x + 1, y + 1] in trace_coordinates:
        number_of_neighbours += 1
    if [x + 1, y] in trace_coordinates:
        number_of_neighbours += 1
    if [x + 1, y - 1] in trace_coordinates:
        number_of_neighbours += 1
    if [x, y - 1] in trace_coordinates:
        number_of_neighbours += 1
    if [x - 1, y - 1] in trace_coordinates:
        number_of_neighbours += 1
    if [x - 1, y] in trace_coordinates:
        number_of_neighbours += 1
    if [x - 1, y + 1] in trace_coordinates:
        number_of_neighbours += 1
    return number_of_neighbours

def fill_large_hole(binary_image):
    # Fills large holes and smoothes edges
    # Iterate over each 2D slice of the 3D image
    filled_slices = []
    smoothing_radius=5
    for slice in binary_image:
        # Fill holes in the slice
        filled_slice = binary_fill_holes(slice)
        
        # Smooth the filled slice using binary closing
        struct = morphology.disk(smoothing_radius)
        smoothed_slice = morphology.binary_closing(filled_slice, footprint=struct)
        
        filled_slices.append(smoothed_slice)

    # Stack the filled slices back into a 3D image
    filled_image = np.stack(filled_slices)

    return filled_image

def replace_intensity(mask, img):
    if not (mask.shape == img.shape):
        print("Mask shape differs from image shape")
        return None
    
    # Convert the mask to the same data type as the image
    mask = mask.astype(img.dtype)
    
    # Overlap Mask on original image (as reference)
    mat_intensity = np.where(np.logical_and(mask, img), img, 0)
    
    return mat_intensity

#np.seterr(divide='ignore', invalid='ignore')
def normalize_intensity(k):
    #return np.clip(k, 0., 1., out=k)
    # ##k_min = k.min(axis=(1, 2), keepdims=True)
    # ##k_max = k.max(axis=(1, 2), keepdims=True)
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

# Read the Tif/CZI File
def folder_scan(directory):
    valid_extensions = (".czi", ".tif")  # Valid extensions: CZI and TIF
    get_files = []
    for f_name in os.listdir(directory):
        if f_name.lower().endswith(valid_extensions):
            get_files.append(os.path.join(directory, f_name))
    return get_files

# Contrast stretching
def contrast_stretching(img):
    p2, p98 = np.percentile(img, (2, 98))
    img_rescale = exposure.rescale_intensity(img, in_range=(p2, p98))
    
    return img_rescale

def connected_components_labeling_box(binary_input):
    binary_input = binary_input.astype('int32')
    labeled_array, num_labels = ndimage.label(binary_input)
    return labeled_array

#folder_path = 'E:/Quantified/eNUC Live/230628 enuc vca vs novca/Test_Center/toGenerateMask_Cellpose/'
#get_files = folder_scan(folder_path)

print(" Folder List Stage")

folder_list = [
    'C:/Users/ChimieENS/Documents/Cellpose/230406 OF1 D1 Actin 10uMVCA Gactin RH3K9me2 FRH3K27me2/'
]

for folder_name in folder_list:
    print(folder_name)
    get_files = folder_scan(folder_name)
    print("Scanning Files in the folder >>>>>>>", folder_name)

    # Conversion to um
    px_to_um_X = 0.0263576
    px_to_um_Y = 0.0263576
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
            img_dextran = aics_image.get_image_data("ZYX", T=0, C=1) # Dextran Channel
            img_nuclei = aics_image.get_image_data("ZYX", T=0, C=2)
        else:
            img_actin = aics_image.get_image_data("ZYX", T=0, C=0) # Actin Channel
            img_H3K27me3 = aics_image.get_image_data("ZYX", T=0, C=1) # yH2AX
            img_H3K9me2 = aics_image.get_image_data("ZYX", T=0, C=2) # LaminB1
            img_nuclei = aics_image.get_image_data("ZYX", T=0, C=3) # Nucleus
    
        Result_folder = os.path.join(folder_name, 'Result')
        if not os.path.exists(Result_folder):
            os.makedirs(Result_folder)
            
        # For saving results #.stem

        # Upscale the Z dimension to be isotropic

        # channel to segment and nuclear channel 
        # numbering starts at 1 
        # for your single channel image use [0, 0] 
        # for the multi channel image it's [3, 0]
        channels = [0, 0] 
        diameter = 224.397#model diam_labels  165.257 (mean diameter of training ROIs)
        use_GPU = True
        stitch_threshold=1
        #do_3D=True
        #resample=True
        cellprob_threshold = 0

        pretrained_model = "C:/Users/ChimieENS/Documents/Cellpose/Trained_model/Fixed_eNUC/CP_20230620_182329"

        model_match_threshold = 27 #30
        flow_threshold = (31.0 - model_match_threshold) / 10.

        logger = io.logger_setup()
        print('Segmenting Image ---', filename)
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
        print('The Filename ---', filename)
        merged_Labels = connected_components_labeling_box(mask)
        merged_Labels_np = np.array(merged_Labels)#.astype('int32')
    
    
        ### Save the Mask
        prediction_stack_32 = img_as_float32(merged_Labels_np, force_copy=False)     
        os.chdir(Result_folder)
        imwrite("(Mask)_"+filename+".tif", prediction_stack_32)

        # Structure  for Actin Dilation
        diamond = np.zeros((3, 3, 3), dtype=bool)
        diamond[1:2, 0:3, 0:3] = True

        label_OrgCnt = measure.label(merged_Labels_np)

        # Calculate the number of stacks for each label
        label_counts = np.bincount(label_OrgCnt.ravel())
    
        print("Number of Unique mask Found", np.unique(label_OrgCnt)[1:])
        for lbl_count in np.unique(label_OrgCnt)[1:]:
            
            if label_counts[lbl_count] >= 100000:# 100000from 120um3 min nuclei I expect / 0.0351435*2

                print('lbl_count >>>>', lbl_count)
    
                maskLBL = label_OrgCnt == lbl_count
                nucleus_Inten = replace_intensity(maskLBL,img_nuclei)
                print("Nucleus Intensity")
            
                nucelus_IntenNorm = normalize_intensity(nucleus_Inten)
                print("normalize_intensity")
                nucelus_AdaptEq = exposure.equalize_adapthist(nucelus_IntenNorm, clip_limit=0.03)
                print("equalize_adapthist")
                nucelus_GB = nsitk.gaussian_blur(nucelus_AdaptEq, 2.0, 2.0, 0.0)
                nucelus_OT = nsitk.threshold_otsu(nucelus_GB)
                print("Threshold Otsu")
                # Chromo_SGB = cle.subtract_gaussian_background(nucelus_GB, None, 10.0, 10.0, 0.0)            
                # Chromo_Thres_Inter = nsitk.threshold_intermodes(Chromo_SGB)
                
                structuring_element = generate_binary_structure(nucelus_OT.ndim, 3)
        
                nucelus_FLH1 = fill_large_hole(nucelus_OT)
                print("Nucleus Fill Large Holes")
                nucelus_FLH2 = fill_large_hole(nucelus_FLH1)#.astype(int)

                print("Nucleus Fill Large Holes 2")    
                # # Nucleus Segmentation
                
                statistics_nucleus = nsitk.label_statistics(intensity_image=nucleus_Inten, label_image=maskLBL,
                                                size = True, intensity = True, perimeter= True,
                                                shape= True, position= True)#, moments= True)
        
                nuclei_Area = statistics_nucleus.loc[0, 'number_of_pixels']
                nuclei_Area = nuclei_Area*px_to_um_X*px_to_um_Y*px_to_um_Z
                statistics_nucleus['Nucleus Area'] = nuclei_Area 
                print('nuclei_Area ---', nuclei_Area)
        
                # statistics_nucleus = cle.statistics_of_labelled_pixels(img_nuclei, nucelus_FLH2)    

                # print("nucleus Statisitcs")
                # nuclei_Area = np.sum(statistics_nucleus['area'], axis=0)
                # #nuclei_Area = statistics_nucleus.loc[0, 'number_of_pixels']
                # print("Number of Pixels", nuclei_Area)
                # nuclei_Area = nuclei_Area*px_to_um_X*px_to_um_Y*px_to_um_Z
                # statistics_nucleus['Nucleus Area'] = nuclei_Area
                # print('nuclei_Area >>>>', nuclei_Area)
            
                ### Save the Nucleus Mask
                prediction_stack_32 = img_as_float32(nucelus_FLH2, force_copy=False)     
                os.chdir(Result_folder)
                imwrite("(Nucleus)_"+filename+".tif", prediction_stack_32)              
                pd.DataFrame(statistics_nucleus).to_excel('(Nucleus)_' + filename+'_'+ str(lbl_count)+'.xlsx')
            
                # Chromocenter Segmentation # Approach 1
                # Normmalized nucleus, Top Hat box, Intermodes Threshold
                # Chromo_GB = cle.gaussian_blur(nucelus_IntenNorm, None, 2.0, 2.0, 0.0)
                # #image1_gb = cle.gaussian_blur(intensity_map_norm, None, 2.0, 2.0, 0.0)
                # Chromo_THB = cle.top_hat_box(Chromo_GB, None, 15.0, 15.0, 0.0)
                # # subtract gaussian background
                # Chromo_SGB = cle.subtract_gaussian_background(Chromo_THB, None, 10.0, 10.0, 0.0)
                # Chromo_Thres_Inter1 = nsitk.threshold_intermodes(Chromo_SGB)
                
                # Chromocenter Segmentation # Approach 2
                # Normalized Nucleus, Max Filter, Max entropy
                Chromo_GB = cle.gaussian_blur(nucelus_IntenNorm, None, 2.0, 2.0, 0.0)
                # maximum filter
                Chromo_THB = cle.subtract_gaussian_background(Chromo_GB, None, 10.0, 10.0, 0.0)
                # threshold maximum entropy
                Chromo_Thres_Inter2 = nsitk.threshold_maximum_entropy(Chromo_THB)
                
                # # Chromocenter Segmentation # Approach 3
                # # Adaptive Thresholding to Intermodes
                # Chromo_GB = cle.gaussian_blur(nucelus_AdaptEq, None, 2.0, 2.0, 0.0) # 'nucelus_AdaptEq'
                # # top hat box
                # Chromo_THB = cle.top_hat_box(Chromo_GB, None, 15.0, 15.0, 0.0)
                # # threshold intermodes
                # Chromo_Thres_Inter3 = nsitk.threshold_intermodes(Chromo_THB)
                
                Chromo_Thres_Inter_Comb = cle.connected_components_labeling_box(Chromo_Thres_Inter2)
                statistics_Chromo = cle.statistics_of_labelled_pixels(img_nuclei, Chromo_Thres_Inter_Comb)                 
                intermodes_chromo_Area = np.sum(statistics_Chromo['area'], axis=0)

                
                #intermodes_chromo_Area = statistics_Chromo.loc[0, 'number_of_pixels']
                chromointermodes_Area = intermodes_chromo_Area *px_to_um_X* px_to_um_Y*px_to_um_Z
                statistics_Chromo['Chromocenter Area'] = chromointermodes_Area 
                print('Chromocenter_Area >>>>', chromointermodes_Area)

                ### Save the Chromocenter Mask
                prediction_stack_32 = img_as_float32(Chromo_Thres_Inter2, force_copy=False)     
                os.chdir(Result_folder)
                imwrite("(Chromocenter)_"+filename+".tif", prediction_stack_32)                            
                pd.DataFrame(statistics_Chromo).to_excel('(Chromo)_' + filename + '_' + str(lbl_count)+'.xlsx')   
            
                nuc_EdgThin = np.zeros(img_nuclei.shape)
                
                for i in range(nucelus_FLH2.shape[0]):            
                    thinned_slice = cle.reduce_labels_to_label_edges(nucelus_FLH2[i]) #cle.reduce_labels_to_label_edges(image2_FH[i])            
                    nuc_EdgThin[i] = thinned_slice
                    
                print('nucleus Edge Thinning >>>>')
            
                nuc_EdgThin_Ar = np.array(nuc_EdgThin)            

                # Actin Segmentation

                act_obj = np.zeros(img_actin.shape)
                
                dilated = ndimage.binary_dilation(maskLBL, diamond, iterations=15).astype(maskLBL.dtype)
                actin_img = replace_intensity(dilated, img_actin)
                actin_img_N = normalize_intensity(actin_img)
                print(">> Actin Intensity Normalized <<")
                
                # H3K27me3 Segmentation                
                statistics_H3K27me3 = cle.statistics_of_labelled_pixels(img_H3K27me3, dilated)
                os.chdir(Result_folder)
                pd.DataFrame(statistics_H3K27me3).to_excel('(H3K27me3)_' + filename + '_' + str(lbl_count) + '.xlsx')                   
                print('H3K27me3_Intensity --- ok')
                
                # H3K9me2 Segmentation
                statistics_H3K9me2 = cle.statistics_of_labelled_pixels(img_H3K9me2, dilated)
                os.chdir(Result_folder)
                pd.DataFrame(statistics_H3K9me2).to_excel('(H3K9me2)_' + filename + '_' + str(lbl_count) + '.xlsx')                   
                print('H3K9me2_Intensity --- ok')
                
                # Adaptive Thresholding
                min_threshold_mean = 20
                min_threshold_std = 800
                
                actin_adapteq = exposure.equalize_adapthist(actin_img_N, clip_limit=0.03)
            
                img_mean = np.mean(actin_img)
                print('Mean Actin Intensity >>>>', img_mean)
                img_std = np.std(actin_img)
                print('Std Actin Intensity >>>>', img_std)
            
                #######
                # View Segmentation
                #######
                # viewer = napari.Viewer()
                # viewer.add_image(nucleus_Inten)  
                # viewer.add_image(nucelus_AdaptEq)
                # #viewer.add_image(nuc_EdgThin_Ar)
                # viewer.add_image(nucelus_FLH2)
                # #viewer.add_labels(Chromo_Thres_Inter1)
                # viewer.add_labels(Chromo_Thres_Inter2)
                # #viewer.add_labels(Chromo_Thres_Inter3)
                # viewer.add_image(actin_img_N)
                # viewer.add_image(actin_adapteq)
        
                image2_Blur = cle.gaussian_blur(actin_adapteq, None, 2.0, 2.0, 0.0)
                image2_Gaus = nsitk.white_top_hat(image2_Blur, 10, 10, 0)
                print('>> Pre Filtering Actin <<')
                if img_mean >= min_threshold_mean and img_std <= min_threshold_std:
                    actin_binary = nsitk.threshold_otsu(image2_Gaus)
                    print("Ostu")

                    act_obj = np.zeros(img_nuclei.shape)

                    # Apply Otsu thresholding
                    ####################
                    # Condition to remove background Actin segmentation
                    
                    opImBase = img_nuclei[0,:,:]        
            
                    for i in range(actin_binary.shape[0]):
                        prune = pruneSkeleton(actin_binary[i], opImBase)
                        act_obj[i] = prune
                         
                    #######
                    # View Segmentation
                    #######
                    #viewer.add_image(act_obj)
                
                    # Compute the standard deviation for each slice in the binary mask
                    std_per_slice = np.std(act_obj, axis=(1, 2))

                    # Calculate the ratio of the first and second value of the standard deviation
                    ratio = std_per_slice[0] / std_per_slice[1]

                    # Find the indices of slices that have a ratio higher than a threshold
                    threshold_ratio = 5.0  # Set your desired threshold ratio here
                    indices_to_remove = np.where(std_per_slice[1:] / std_per_slice[:-1] > threshold_ratio)[0]  + 1

                    # Set the values of the slices to zero instead of deleting them
                    act_obj[indices_to_remove] = 0
        
                    # Write statistics to Excel file
                    statistics_Actin =  cle.statistics_of_labelled_pixels(img_actin, actin_binary)    
                    
                    print("Actin Found")
                

                    Coverage = []
                    for act_Z, nucl_Z in zip(act_obj, nuc_EdgThin_Ar):
                        ratio = np.sum(act_Z) / np.sum(nucl_Z)
                        if np.sum(nucl_Z) != 0 and np.sum(act_Z) != 0:
                            if ratio > 1 and ratio < 1.6:#ratio <= 100:
                                Coverage.append(100)
                            elif ratio < 1:
                                Coverage.append(ratio * 100)                                
                            else:
                                Coverage.append(0)
            
                            coverage_values = [value for value in Coverage if value is not None and not np.isnan(value) and value != 0]

                            if coverage_values:
                                actin_coverage = sum(coverage_values) / len(coverage_values)
                            else:
                                actin_coverage = 0.0
                        
                            print("Actin Coverage:", actin_coverage)
                            statistics_Actin['Actin Coverage'] = actin_coverage 
                        
                    ### Save the Actin Mask
                    prediction_stack_32 = img_as_float32(actin_binary, force_copy=False)     
                    os.chdir(Result_folder)
                    imwrite("(Actin)_"+filename+".tif", prediction_stack_32)                    
                    pd.DataFrame(statistics_Actin).to_excel('(Actin)_' + filename + '_' + str(lbl_count) + '.xlsx')             
            
                elif img_std >= min_threshold_std:
                    actin_binary = nsitk.threshold_maximum_entropy(image2_Gaus)
                #viewer.add_image(actin_img)
                    print("Max entropy")
            
                    opImBase = img_nuclei[0,:,:]        
            
                    for i in range(actin_binary.shape[0]):
                        prune = pruneSkeleton(actin_binary[i], opImBase)
                        act_obj[i] = prune
                    
                    #######
                    # View Segmentation
                    #######
                    #viewer.add_image(act_obj)
                    # Compute the standard deviation for each slice in the binary mask
                    std_per_slice = np.std(act_obj, axis=(1, 2))

                    # Calculate the ratio of the first and second value of the standard deviation
                    ratio = std_per_slice[0] / std_per_slice[1]

                    # Find the indices of slices that have a ratio higher than a threshold
                    threshold_ratio = 5.0  # Set your desired threshold ratio here
                    indices_to_remove = np.where(std_per_slice[1:] / std_per_slice[:-1] > threshold_ratio)[0]  + 1

                    # Set the values of the slices to zero instead of deleting them
                    act_obj[indices_to_remove] = 0

                    statistics_Actin =  cle.statistics_of_labelled_pixels(img_actin, actin_binary)    

                    print("Actin Found")

                    Coverage = []
                    for act_Z, nucl_Z in zip(act_obj, nuc_EdgThin_Ar):
                        ratio = np.sum(act_Z) / np.sum(nucl_Z)
                        if np.sum(nucl_Z) != 0 and np.sum(act_Z) != 0:
                            if ratio > 1 and ratio < 1.6:#ratio <= 100:
                                Coverage.append(100)
                            elif ratio < 1:
                                Coverage.append(ratio * 100)                                
                            else:
                                Coverage.append(0)
            
                            coverage_values = [value for value in Coverage if value is not None and not np.isnan(value) and value != 0]

                            if coverage_values:
                                actin_coverage = sum(coverage_values) / len(coverage_values)
                            else:
                                actin_coverage = 0.0

                            print("Actin Coverage:", actin_coverage)
                            statistics_Actin['Actin Coverage'] = actin_coverage   
                    
                    ### Save the Actin Mask
                    prediction_stack_32 = img_as_float32(actin_binary, force_copy=False)     
                    os.chdir(Result_folder)
                    imwrite("(Actin)_"+filename+".tif", prediction_stack_32)            
                    pd.DataFrame(statistics_Actin).to_excel('(Actin)_' + filename + '_' + str(lbl_count) + '.xlsx')
                                
                else:
                    #actin_binary = nsitk.threshold_maximum_entropy(image2_Gaus)
                    print("No actin present")

# Close the log file
#log_file.close()
