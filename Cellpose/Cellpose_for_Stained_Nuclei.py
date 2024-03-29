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
from skimage.filters import threshold_multiotsu
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


def threshold_binary_image(image, min_threshold_mean, min_threshold_std):
    """Function to process Actin Coverage. In case 1. when normal actin is present
    then I can apply OTSU simply. If a big actin blob is present then I apply maximum
    entropy to segment actin region."""
    
    if img_mean >= min_threshold_mean and img_std <= min_threshold_std:
        binary_image = nsitk.threshold_otsu(image)
        return binary_image, "Otsu"
    elif img_std >= min_threshold_std:
        binary_image = nsitk.threshold_maximum_entropy(image)
        return binary_image, "Max entropy"
    else:
        return None, "No actin present"

def get_actin_coverage(act_obj, nuc_EdgThin_Ar):
    """Function to get Actin coverage. Here If the ratio of actin to nucleus is in the rang of
    1-1.6 in cases where due to nucleus membrane curling I have extra regions after thinning
    A true value is taken if the coverge is below 1 ."""
    
    Coverage = []
    actin_coverage = 0.0
    for act_Z, nucl_Z in zip(act_obj, nuc_EdgThin_Ar):
        ratio = np.sum(act_Z) / np.sum(nucl_Z)
        if np.sum(nucl_Z) != 0 and np.sum(act_Z) != 0:
            if 1 < ratio < 1.6:
                Coverage.append(100)
            elif ratio < 1:
                Coverage.append(ratio * 100)
            else:
                Coverage.append(0)

    coverage_values = [value for value in Coverage if value is not None and not np.isnan(value) and value != 0]

    if coverage_values:
        actin_coverage = sum(coverage_values) / len(coverage_values)

    return actin_coverage

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
    min_NUC_Size = px_to_um_X * px_to_um_Y
    
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
        # Anisotropic pixel scaling
        Z_Scaling = 1.0 #/ 0.0351435
        Y_Scaling = 1.0  # Assuming Y and X pixels are already isotropic
        X_Scaling = 1.0
        # Create an anisotropic dilation structuring element
        diamond = np.zeros((int(3 / Z_Scaling), int(3 / Y_Scaling), int(3 / X_Scaling)), dtype=bool)
        diamond[1:2, 0:int(3 / Y_Scaling), 0:int(3 / X_Scaling)] = True

        label_OrgCnt = measure.label(merged_Labels_np)

        # Calculate the number of stacks for each label
        label_counts = np.bincount(label_OrgCnt.ravel())
    
        print("Number of Unique mask Found", np.unique(label_OrgCnt)[1:])
        for lbl_count in np.unique(label_OrgCnt)[1:]:
            
            if label_counts[lbl_count] >= 110/min_NUC_Size:# 100000from 120um3 min nuclei I expect / 0.0351435*2

                print('lbl_count >>>>', lbl_count)
    
                maskLBL = label_OrgCnt == lbl_count
                nucleus_Inten = replace_intensity(maskLBL,img_nuclei)            
                print("Nucleus Intensity")     
                
                if 'img_dextran' in globals():# and img_dextran not in None:
                    checkPorous = replace_intensity(maskLBL,img_dextran)
                    z_slice = checkPorous.shape[0]
                    avrMaskSize = int(z_slice * 0.5)
                    dexIntenThresh = np.mean(checkPorous[avrMaskSize])
                    print('Dextran Intensity >>>>', dexIntenThresh)
                    if dexIntenThresh >= 2800:# Here I want to scan The mean Z
                    # Additional check to make sure no porous nuclei is segmented
                        continue
                else:                    
                        pass  
                #img_dextran = None                     

                #######
                # View Segmentation
                #######
                #viewer = napari.Viewer()         
                               
                nucelus_OrgNorm = normalize_intensity(nucleus_Inten)
                print("Nucleus ->  normalize_intensity")                         
                                
                nucelus_GB = nsitk.gaussian_blur(nucelus_OrgNorm, 2.0, 2.0, 0.0)
                nucelus_GC = cle.gamma_correction(nucelus_GB, None, 1.0)
                nucelus_OT = cle.threshold_otsu(nucelus_GC)
                print("Nucleus ->  Threshold Otsu")
                nucelus_CL = cle.closing_labels(nucelus_OT, None, 1.0)
                nucelus_BF = nsitk.binary_fill_holes(nucelus_CL)
                nucelus_FLH = fill_large_hole(nucelus_BF)
                print("Nucleus -> Fill Binary Mask")
                nucleus_Trim = replace_intensity(nucelus_FLH,img_nuclei)
                
                #viewer.add_image(nucleus_Trim)
                
                nucleus_TrimNorm = normalize_intensity(nucleus_Trim)   
                #######
                # View Segmentation
                #######
                #viewer.add_image(nucleus_TrimNorm)
                            
                print("Nucleus -> Normalized Intensity")                
                # # Nucleus Segmentation       
                statistics_nucleus = cle.statistics_of_labelled_pixels(img_nuclei, nucelus_FLH)    

                print("Nucleus ->  Statisitcs")
                
                nuclei_Area = np.sum(statistics_nucleus['area'], axis=0)
                #nuclei_Area = statistics_nucleus.loc[0, 'number_of_pixels']
                print("Nucleus -> Number of Pixels", nuclei_Area)
                nuclei_Area = nuclei_Area*px_to_um_X*px_to_um_Y*px_to_um_Z
                statistics_nucleus['Nucleus Area'] = nuclei_Area
                
                print("Nucleus  Area -> ", nuclei_Area)
                               
                #######
                # View Segmentation
                #######
                # viewer.add_image(multiOTSU_GB) 
                ### Save the Nucleus Mask
                prediction_stack_32 = img_as_float32(nucelus_FLH, force_copy=False)     
                os.chdir(Result_folder)
                imwrite("(Nucleus)_"+filename + str(lbl_count)+".tif", prediction_stack_32)              
                pd.DataFrame(statistics_nucleus).to_excel('(Nucleus)_' + filename+'_'+ str(lbl_count)+'.xlsx')

                
                # Chromocenter Segmentation
                
                ChromoBlur = cle.gaussian_blur(nucleus_Trim, None, 2.0, 2.0, 0.0)                
                ChromoTHB = cle.top_hat_box(ChromoBlur, None, 20.0, 20.0, 0.0)
                ChromoTHBAr = np.array(ChromoTHB) 
                Chrom_multiOTSU = threshold_multiotsu(nucleus_Trim)                    
                #ChromoSeg = ChromoTHBAr > Chrom_multiOTSU[1]
                
                try:       
                    ChromoSeg = ChromoTHBAr > Chrom_multiOTSU[1]
                    print("<- Chromocenter Segmentation -> ")
                    #######
                    # View Segmentation
                    #######
                    #viewer.add_image(Chromo_Thres_Inter2)
                    
                    Chromo_Thres_Inter_Comb = cle.connected_components_labeling_box(ChromoSeg)
                    statistics_Chromo = cle.statistics_of_labelled_pixels(img_nuclei, Chromo_Thres_Inter_Comb)                 
                    intermodes_chromo_Area = np.sum(statistics_Chromo['area'], axis=0)
                
                    #intermodes_chromo_Area = statistics_Chromo.loc[0, 'number_of_pixels']
                    chromointermodes_Area = intermodes_chromo_Area *px_to_um_X* px_to_um_Y
                    statistics_Chromo['Chromocenter Area'] = chromointermodes_Area 
                    print('Chromocenter_Area ->', chromointermodes_Area)        
                    
                    ### Save the Chromocenter Mask
                    prediction_stack_32 = img_as_float32(ChromoSeg, force_copy=False)     
                    os.chdir(Result_folder)
                    imwrite("(Chromocenter)_"+filename + str(lbl_count) +".tif", prediction_stack_32)                            
                    pd.DataFrame(statistics_Chromo).to_excel('(Chromo)_' + filename + '_' + str(lbl_count)+'.xlsx')   
                    
                except RuntimeError:
                    
                    print("No Chromocenter Found")

                nuc_EdgThin = np.zeros(img_nuclei.shape)
                
                for i in range(nucelus_FLH.shape[0]):            
                    thinned_slice = cle.reduce_labels_to_label_edges(nucelus_FLH[i]) #cle.reduce_labels_to_label_edges(image2_FH[i])            
                    nuc_EdgThin[i] = thinned_slice
                    
                print('Nucleus Edge Thinning >>>>')
            
            
                nuc_EdgThin_Ar = np.array(nuc_EdgThin)

                # Actin Segmentation
                #if 'img_actin' in globals():
                if 'img_actin' in globals() and img_actin is not None:       
                    
                    act_obj = np.zeros(img_actin.shape)
                    
                    # # To verify in some case it is not Dextran instead
                    z_slice = img_actin.shape[0]
                    avrMaskSize = int(z_slice * 0.5)
                    checkDextran = np.mean(img_actin[avrMaskSize])
                    print('checkDextran >>>>', checkDextran)
                    if checkDextran > 2800:# Here I want to scan The mean Z
                    # Additional check to dextran not in Actin channel
                        print(">> It is Dextran not Actin <<")
                        print(">> Skipping <<")
                        continue
                    else:                    
                        pass                    
                    
                    dilated = ndimage.binary_dilation(maskLBL, diamond, iterations=10).astype(maskLBL.dtype)
                    actin_img = replace_intensity(dilated, img_actin)
                    actin_img_N = normalize_intensity(actin_img)
                    print(">> Actin Intensity Normalized <<")                                    
                
                    # Adaptive Thresholding
                    min_threshold_mean = 20
                    min_threshold_std = 800
                
                    actin_adapteq = exposure.equalize_adapthist(actin_img_N, clip_limit=0.03)
            
                    img_mean = np.mean(actin_img)
                    print('Actin  -> Mean Actin Intensity >>>>', img_mean)
                    img_std = np.std(actin_img)
                    print('Actin  -> Std Actin Intensity >>>>', img_std)
                    
                    image2_Blur = cle.gaussian_blur(actin_adapteq, None, 2.0, 2.0, 0.0)
                    image2_Gaus = nsitk.white_top_hat(image2_Blur, 10, 10, 0)
                    print('>> Pre Filtering Actin <<')
                    
                    actin_binary, threshold_method = threshold_binary_image(image2_Gaus, min_threshold_mean, min_threshold_std)
                    print(threshold_method)
                    
                    opImBase = img_nuclei[0,:,:]        
                    if actin_binary is not None:
                        
                        for i in range(actin_binary.shape[0]):
                            prune = pruneSkeleton(actin_binary[i], opImBase)
                            act_obj[i] = prune
                            # All the thinning operations are not uniform in top and bottom slices
                            # I decided to ignore those cases but I consider them in the slices
                            # Considering normal slice with actin and top/bottom slices with actin will have more threads
                            # then here is                         
                        std_per_slice = np.std(act_obj, axis=(1, 2))
                        ratio = std_per_slice[0] / std_per_slice[1]
                        threshold_ratio = 5.0  #  threshold ratio 
                        indices_to_remove = np.where(std_per_slice[1:] / std_per_slice[:-1] > threshold_ratio)[0]  + 1
                        act_obj[indices_to_remove] = 0
                        
                        numSliceRem = len(indices_to_remove)

                        # Write Actin statistics to Excel file
                        statistics_Actin =  cle.statistics_of_labelled_pixels(img_actin, actin_binary)    

                        print("Actin  -> Actin Found")
                        actin_coverage = get_actin_coverage(act_obj, nuc_EdgThin_Ar)                        
                    
                        print("Actin  -> Actin Coverage:", actin_coverage)
                        statistics_Actin['Actin Coverage'] = actin_coverage 
                        
                        ### Save the Actin Mask
                        prediction_stack_32 = img_as_float32(actin_binary, force_copy=False)     
                        os.chdir(Result_folder)
                        imwrite("(Actin)_" + str(lbl_count) +filename+".tif", prediction_stack_32)                    
                        pd.DataFrame(statistics_Actin).to_excel('(Actin)_' + filename + '_' + str(lbl_count) + '.xlsx')             
                         
                    else:
                    
                        #actin_binary = nsitk.threshold_maximum_entropy(image2_Gaus)
                        print("Actin  -> No actin present")

                       
        img_actin = None
