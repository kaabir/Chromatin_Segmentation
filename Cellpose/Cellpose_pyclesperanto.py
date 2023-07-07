# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 18:22:51 2023

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

from pathlib import Path
#from astropy.visualization import simple_norm
from skimage import img_as_float32
from skimage.util import img_as_ubyte
from tifffile import imsave

device = cle.select_device("gfx1032")
#device = cle.select_device("AMD Radeon Pro W6600")

# Read the Tif/CZI File
get_files = []

def folder_Scan(directory):
    # Actin after  30min  - value assigned is 3
    os.chdir(directory)
    global_path = glob.glob('*')
    
    extension = '.tif' # '.czi'
    for g_name in global_path:
        if g_name.endswith(extension):
            get_files.append(g_name)
        else:
            print('Some random file or folders in the directory')    

folder_path = 'C:/Users/kaabi/Documents/Nuceli_Data/Enucleation/Cellpose/Test_Dataset_Nuclei/15-12-2022/'

folder_Scan(folder_path)

# Conversion to um
px_to_um_X = 0.0351435
px_to_um_Y = 0.0351435
px_to_um_Z = 1 #1.0/0.0351435 # making pixels isotropic


for image in get_files:
    aics_image = AICSImage(image)
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
        img_chnl1_read = aics_image.get_image_data("ZYX", T=0, C=0) # A Channel
        img_chnl2_read = aics_image.get_image_data("ZYX", T=0, C=1) # Dextran Channel
        img_chnl3_read = aics_image.get_image_data("ZYX", T=0, C=2)
        img_nuclei = aics_image.get_image_data("ZYX", T=0, C=3)

# For saving results #.stem
Result_folder = folder_path + '/Result/' 
# Upscale the Z dimension to be isotropic

# channel to segment and nuclear channel 
# numbering starts at 1 
# for your single channel image use [0, 0] 
# for the multi channel image it's [3, 0]
channels = [2, 0] 
diameter = 165.257 #169.708
use_GPU = True
stitch_threshold=1
#do_3D=True
#resample=True
cellprob_threshold = 0

pretrained_model = "C:/Users/ChimieENS/Documents/Cellpose/Trained_model/CP_20230306_150448_Main"

model_match_threshold = 29 #30
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
#merged_Labels_np = split_touching_objects(merged_Labels_np)

# Getting Original intensities
#intensity_vector = cle.read_intensities_from_map(mask, img)
statistics = cle.statistics_of_labelled_pixels(img_nuclei, merged_Labels)
#

px_to_um_scaled = 4.340452764123788e-5
nuclei_Area = statistics['area']*px_to_um_Y*px_to_um_Y*px_to_um_Z

# Chromocenter Segmentation
# Get voxel of pixel intensity values inside the mask 
def intensity_vector(mask, img):
    if not (mask.shape == img.shape):
        return False
    
    mat_intensity = np.where(np.logical_and(mask,img),img,0) # Overlap Mask on original image (as reference)
    ##mat_intensity = np.array(cle.replace_intensities(mask,img))
    
    return mat_intensity

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

from skimage import measure

def calculate_surface_area(mask, threshold=None):
    # generate surface mesh using marching cubes algorithm
    verts, faces, _, _ = measure.marching_cubes(mask)
    # calculate surface area using mesh surface area function
    surface_area = measure.mesh_surface_area(verts, faces)
    #surface_area = surface_area*px_to_um_Y*px_to_um_Y*px_to_um_Z

    return surface_area

# structure for actin segmentation
from scipy import ndimage
from scipy.ndimage import gaussian_filter
from skimage import morphology
from skimage.morphology import binary_dilation
from scipy import ndimage as ndi

diamond = np.zeros((3, 3, 3), dtype=bool)
diamond[1:2, 0:3, 0:3] = True

#print(len(np.unique(merged_Labels)))
number_of_Nuclei = len(np.unique(merged_Labels))-1 # first is mask

# Create an image to place a single nuleus with intensity voxel
im_obj1 = np.zeros(img_nuclei.shape) 
im_obj2 = np.zeros(img_nuclei.shape) 
im_obj3 = np.zeros(img_nuclei.shape)
nuc_lbl_lst1 = []
nuc_lbl_lst2 = []
nuc_lbl_lst3 = []
# Create actin mask
act_obj1 = np.zeros(img_nuclei.shape)
act_obj2 = np.zeros(img_nuclei.shape)
act_obj3 = np.zeros(img_nuclei.shape)
# actin_lbl_lst1 = []
# actin_lbl_lst2 = []
# actin_lbl_lst3 = []


# Extracting mask for each nucleus intensity and storing them
for i in np.unique(merged_Labels_np)[1:]: # Initial (1) is mask
    if (len(np.unique(merged_Labels)-1) == 1):# and no_chnls>1:
        break
    elif i ==1:
        #nuc_img = intensity_vector(merged_Labels_np,img_nuclei)
        nuc_lbl_lst1.append(merged_Labels_np)
        # Save Labels
        ##prediction_stack_32 = img_as_float32(label1, force_copy=False)     
        ##os.chdir(Result_folder)
        #imsave(filename+"_1.tif", prediction_stack_32)
        # Use nuclei mask as seed for actin segmentation
        # Pad some pixels for actin with constant value which matches unique label value
        # place back actin and segment with classical thresholding
    
        dilated1 = ndi.binary_dilation(merged_Labels_np, diamond, iterations=10).astype(merged_Labels_np.dtype)
        actin_img1 = intensity_vector(dilated1,img_actin)
        actin_filter1 = nsitk.median_filter(actin_img1, radius_x=2, radius_y=2, radius_z=0)
        actin_binary1 = nsitk.threshold_otsu(actin_filter1)
        
        #skeleton = skeletonize(actin_binary,method="lee")
        
        for i in range(actin_binary1.shape[0]):
            thinned_slice = thin(actin_binary1[i])
            act_obj1[i] = thinned_slice
        
        statistics_surf_actin1 = cle.statistics_of_labelled_pixels(actin_img1, act_obj1)
        actin_surf_Area1 = np.sum(statistics_surf_actin1['area'], axis=0)
        actin_surf_Area1 = actin_surf_Area1*px_to_um_Y
        #pd.DataFrame(statistics_surf_actin1).to_excel(Result_folder+filename+'(Actin)Actin_statistics_1.xlsx')

        #print(merged_Labels_np.shape)
        #print(dilated1.shape)

    elif i ==2: 
        im_obj2[merged_Labels_np == i] = 1
        label2 = im_obj2
        nuc_lbl_lst2.append(label2)
        # Save Labels
        prediction_stack_32 = img_as_float32(label2, force_copy=False)     
        os.chdir(Result_folder)
        #imsave(filename+"_2.tif", prediction_stack_32)
        
        dilated2 = ndi.binary_dilation(merged_Labels_np, diamond, iterations=10).astype(merged_Labels_np.dtype)
        actin_img2 = intensity_vector(dilated2,img_actin)
        actin_filter2 = nsitk.median_filter(actin_img2, radius_x=2, radius_y=2, radius_z=0)
        actin_binary2 = nsitk.threshold_otsu(actin_filter2)
        
        #skeleton = skeletonize(actin_binary,method="lee")
        
        for i in range(actin_binary2.shape[0]):
            thinned_slice = thin(actin_binary2[i])
            act_obj2[i] = thinned_slice
        
        statistics_surf_actin2 = cle.statistics_of_labelled_pixels(actin_img2, act_obj2)
        actin_surf_Area2 = np.sum(statistics_surf_actin2['area'], axis=0)
        actin_surf_Area2 = actin_surf_Area2*px_to_um_Y
        #pd.DataFrame(statistics_surf_actin2).to_excel(Result_folder+filename+'(Actin)Actin_statistics_2.xlsx')

        # viewer = napari.Viewer()
        # viewer.add_image(im_obj2)

       
    elif i ==3:
        im_obj3[merged_Labels_np == i] = 1
        label3 = im_obj3
        nuc_lbl_lst3.append(label3)    
        # Save Labels
        prediction_stack_32 = img_as_float32(label3, force_copy=False)     
        #os.chdir(Result_folder)
        #imsave(filename+"_3.tif", prediction_stack_32)
        
        dilated3 = ndi.binary_dilation(merged_Labels_np, diamond, iterations=10).astype(merged_Labels_np.dtype)
        actin_img3 = intensity_vector(dilated3,img_actin)
        actin_filter3 = nsitk.median_filter(actin_img3, radius_x=2, radius_y=2, radius_z=0)
        actin_binary3 = nsitk.threshold_otsu(actin_filter3)
        
        #skeleton = skeletonize(actin_binary,method="lee")
        
        for i in range(actin_binary3.shape[0]):
            thinned_slice = thin(actin_binary3[i])
            act_obj3[i] = thinned_slice
        
        statistics_surf_actin3 = cle.statistics_of_labelled_pixels(actin_img3, act_obj3)
        actin_surf_Area3 = np.sum(statistics_surf_actin3['area'], axis=0)
        actin_surf_Area3 = actin_surf_Area3*px_to_um_Y
        #pd.DataFrame(statistics_surf_actin3).to_excel(Result_folder+filename+'(Actin)Actin_statistics_3.xlsx')

        
        # viewer = napari.Viewer()
        # viewer.add_image(im_obj3)

    else:
        raise ValueError('More than three nucleus')   
    # viewer = napari.Viewer()
    # viewer.add_image(y)

# Not replacing back the original mask values associted with labels
# Place back the intensity values to the individual label
def chromocenter_vector(mask, img):
    if not (mask.shape == img.shape):
        return False

    chromocenter_intensity = np.where(np.logical_and(mask,img),img,0) # Overlap Mask on original image (as reference)
    #chromocenter_intensity = np.array(cle.replace_intensities(mask,img))

    return chromocenter_intensity

for i in np.unique(merged_Labels_np)[1:]: # Background mask at 1
    # Nuclei one to be segmented for chromocenter 1
    if not nuc_lbl_lst1:
        break
    nuc_lbl_lst1 = np.vstack(nuc_lbl_lst1)

    intensity_map1= chromocenter_vector(nuc_lbl_lst1, img_nuclei)
    intensity_map1_norm = normalize_intensity(intensity_map1)
    
    # viewer = napari.Viewer()
    # viewer.add_image(nuc_lbl_lst1)
    # viewer.add_image(intensity_map1_norm)
    # Get Nucleus Statistics
    statistics_nucleus1 = cle.statistics_of_labelled_pixels(img_nuclei, nuc_lbl_lst1)
    #pd.DataFrame(statistics_nucleus1).to_excel(Result_folder+filename+'(Nucleus)Nucleus_statistics_1.xlsx')

    # Blur to enhance chromocenters
    
    intensity_map_blur1 = nsitk.median_filter(intensity_map1_norm, radius_x=2, radius_y=2, radius_z=0)
    # Apply Intermodes Thresholding for Chromocenters
    intensity_map_thb1 = cle.top_hat_box(intensity_map_blur1, None, 10.0, 10.0, 0.0)
    #intensity_map_thb1 = nsitk.black_top_hat(image2_M, 10, 10, 0)
    
    intermodes_Threshold_Chromo1 = nsitk.threshold_intermodes(intensity_map_thb1)
    # Caculate intermodes Area
    chromo_intermodes_Labels1 = cle.connected_components_labeling_box(intermodes_Threshold_Chromo1)
    statistics_intermodes_chromo1 = cle.statistics_of_labelled_pixels(img_nuclei, chromo_intermodes_Labels1)
     
    chromointermodes_Area1 = np.sum(statistics_intermodes_chromo1['area'], axis=0)
    chromointermodes_Area1 = chromointermodes_Area1 * px_to_um_Y*px_to_um_X*px_to_um_Z
    # Export Chromo 1 Stats 
    #pd.DataFrame(statistics_intermodes_chromo1).to_excel(Result_folder+filename+'(Chromo)Chromo_statistics_1.xlsx')
    
    # Max Entropy Chromocenter Segmentation -> Orginial Voxel -> Blur (rx=ry=2) -> Max(Entropy Threshold) 
    entropy_Threshold_Chromo1 = nsitk.threshold_maximum_entropy(intensity_map_thb1)
    chromo_entropy_Labels1 = cle.connected_components_labeling_box(entropy_Threshold_Chromo1)
    # Caculate Entropy Area
    statistics_entropy_chromo1 = cle.statistics_of_labelled_pixels(img_nuclei, chromo_entropy_Labels1)
    
    chromoentropy_Area1 = np.sum(statistics_entropy_chromo1['area'], axis=0)
    chromoentropy_Area1 = chromoentropy_Area1 * px_to_um_Y*px_to_um_X*px_to_um_Z
    
    #binary_fill_holes = binary_fill_holes(intermodes_Threshold_Chromo1)
    #exclude_tiny_label = cle.exclude_small_labels(binary_fill_holes, maximum_size = 500)
    #exclude_tiny_label = cle.exclude_large_labels(binary_fill_holes, minimum_size = 40)

    #Check the nuclei
    #viewer = napari.Viewer()
    
    # viewer.add_image(intensity_map1)
    #viewer.add_image(intensity_map1_norm)
    # viewer.add_image(nuc_lbl_lst1)
    # viewer.add_image(chromo_intermodes_Labels1)
    
    # Nuclei 2 to be segmented for chromocenter 2
    if not nuc_lbl_lst2:
        break
    nuc_lbl_lst2 = np.vstack(nuc_lbl_lst2)
    intensity_map2= chromocenter_vector(nuc_lbl_lst2, img_nuclei)
    intensity_map2_norm = normalize_intensity(intensity_map2)

    # Get Nucleus Statistics
    statistics_nucleus2 = cle.statistics_of_labelled_pixels(img_nuclei, nuc_lbl_lst2)
    #pd.DataFrame(statistics_nucleus2).to_excel(Result_folder+filename+'(Nucleus)Nucleus_statistics_2.xlsx')

    # Blur to enhance chromocenters
    #intensity_map_blur2 = nsitk.gaussian_blur(intensity_map2, variance_x=2, variance_y=2)
    intensity_map_blur2 = nsitk.median_filter(intensity_map2_norm, radius_x=2, radius_y=2)
    # Apply Intermodes Thresholding for Chromocenters
    intensity_map_thb2 = cle.top_hat_box(intensity_map_blur2, None, 10.0, 10.0, 0.0)
    
    intermodes_Threshold_chromo2 = nsitk.threshold_intermodes(intensity_map_thb2)
    # Caculate intermodes Area
    chromo_intermodes_Labels2 = cle.connected_components_labeling_box(intermodes_Threshold_chromo2)
    statistics_intermodes_chromo2 = cle.statistics_of_labelled_pixels(img_nuclei, chromo_intermodes_Labels2)

    chromointermodes_Area2 = np.sum(statistics_intermodes_chromo2['area'], axis=0)
    chromointermodes_Area2 = chromointermodes_Area2 * px_to_um_Y*px_to_um_X*px_to_um_Z
    # Export Chromo 1 Stats 
    #pd.DataFrame(statistics_intermodes_chromo2).to_excel(Result_folder+filename+'(Chromo)Chromo_statistics_2.xlsx')

    # Max Entropy Chromocenter Segmentation -> Orginial Voxel -> Blur (rx=ry=2) -> Max(Entropy Threshold) 
    entropy_Threshold_chromo2 = nsitk.threshold_maximum_entropy(intensity_map_thb2)
    chromo_entropy_Labels2 = cle.connected_components_labeling_box(entropy_Threshold_chromo2)
    # Caculate Entropy Area
    statistics_entropy_chromo2 = cle.statistics_of_labelled_pixels(img_nuclei, chromo_entropy_Labels2)

    chromoentropy_Area2 = np.sum(statistics_entropy_chromo2['area'], axis=0)
    chromoentropy_Area2 = chromoentropy_Area2 * px_to_um_Y*px_to_um_X*px_to_um_Z

    #binary_fill_holes = binary_fill_holes(intermodes_Threshold_chromo2)
    #exclude_tiny_label = cle.exclude_small_labels(binary_fill_holes, maximum_size = 500)
    #exclude_tiny_label = cle.exclude_large_labels(binary_fill_holes, minimum_size = 40)

    # Check the nuclei
    # viewer = napari.Viewer()
    # #viewer.add_image(binary_fill_holes)
    # viewer.add_image(intensity_map2_norm)
    # viewer.add_image(chromo_entropy_Labels2)
    # viewer.add_image(chromo_intermodes_Labels2)
    
    # Nuclei 3 to be segmented for chromocenter 3
    if not nuc_lbl_lst3:
        break
    nuc_lbl_lst3 = np.vstack(nuc_lbl_lst3)
    intensity_map3= chromocenter_vector(nuc_lbl_lst3, img_nuclei)
    intensity_map3_norm = normalize_intensity(intensity_map3)

    # Get Nucleus Statistics
    statistics_nucleus3 = cle.statistics_of_labelled_pixels(img_nuclei, nuc_lbl_lst3)
    #pd.DataFrame(statistics_nucleus3).to_excel(Result_folder+filename+'(Nucleus)Nucleus_statistics_3.xlsx')

    # Blur to enhance chromocenters
    #intensity_map_blur3 = nsitk.gaussian_blur(intensity_map3, variance_x=2, variance_y=2)
    intensity_map_blur3 = nsitk.median_filter(intensity_map3_norm, radius_x=2, radius_y=2)
    # Apply Intermodes Thresholding for Chromocenters
    intensity_map_thb3 = cle.top_hat_box(intensity_map_blur3, None, 10.0, 10.0, 0.0)
    
    intermodes_Threshold_chromo3 = nsitk.threshold_intermodes(intensity_map_thb3)
    # Caculate intermodes Area
    chromo_intermodes_Labels3 = cle.connected_components_labeling_box(intermodes_Threshold_chromo3)
    statistics_intermodes_chromo3 = cle.statistics_of_labelled_pixels(img_nuclei, chromo_intermodes_Labels3)

    chromointermodes_Area3 = np.sum(statistics_intermodes_chromo3['area'], axis=0)
    chromointermodes_Area3 = chromointermodes_Area3 * px_to_um_Y*px_to_um_X*px_to_um_Z
    # Export Chromo 1 Stats 
    #pd.DataFrame(statistics_intermodes_chromo3).to_excel(Result_folder+filename+'(Chromo)Chromo_statistics_3.xlsx')

    # Max Entropy Chromocenter Segmentation -> Orginial Voxel -> Blur (rx=ry=2) -> Max(Entropy Threshold) 
    entropy_Threshold_chromo3 = nsitk.threshold_maximum_entropy(intensity_map_thb3)
    chromo_entropy_Labels3 = cle.connected_components_labeling_box(entropy_Threshold_chromo3)
    # Caculate Entropy Area
    statistics_entropy_chromo3 = cle.statistics_of_labelled_pixels(img_nuclei, chromo_entropy_Labels3)

    chromoentropy_Area3 = np.sum(statistics_entropy_chromo3['area'], axis=0)
    chromoentropy_Area3 = chromoentropy_Area2 * px_to_um_Y*px_to_um_X*px_to_um_Z

    #binary_fill_holes = binary_fill_holes(intermodes_Threshold_chromo3)
    #exclude_tiny_label = cle.exclude_small_labels(binary_fill_holes, maximum_size = 500)
    #exclude_tiny_label = cle.exclude_large_labels(binary_fill_holes, minimum_size = 40)

    # Check the nuclei
    # viewer = napari.Viewer()
    # #viewer.add_image(binary_fill_holes)
    # viewer.add_image(intensity_map3_norm)
    # viewer.add_image(chromo_entropy_Labels3)

# Next Optimization
# from skimage import measure
# labels = measure.label(image)
# # loop over labels
# for i, label in enumerate(np.unique(labels)[1:]):
#     # create a mask for the current label
#     mask = labels == label
#     # do some processing on the mask
#     # save the output to a file with a unique name based on the counter i
#     output_filename = f"output_{i}.tif"
#     # save the output to file

# Making Isotropic 
# voxel_size = [px_to_um_Z, 1, 1] # ZXY
# scale nuclei
# img_nuclei = cle.scale(img_nuclei_read, 
#                    factor_x=voxel_size[2],
#                    factor_y=voxel_size[1],
#                    factor_z=voxel_size[0],
#                    auto_size=True,
#                    linear_interpolation=True
#                   )

# img_nuclei = np.array(img_nuclei)
# from skimage import measure
# from vedo import Mesh, show
# # Generate the voxel using marching cubes algorithm
# verts, faces, _, _ = measure.marching_cubes(img_nuclei, 0.5)

# # Build the polygonal Mesh object from the vertices and faces
# mesh = Mesh([verts, faces])

# # Set the backcolor of the mesh to violet
# # and show edges with a linewidth of 2
# mesh.backcolor('violet').linecolor('tomato').linewidth(2)

# # # Show the mesh, vertex labels, and docstring
# show(mesh, __doc__, viewup='z', axes=1).close()

# scale nuclei
#img_actin = img_actin_read
# img_actin = cle.scale(img_actin_read, 
#                    factor_x=voxel_size[2],
#                    factor_y=voxel_size[1],
#                    factor_z=voxel_size[0],
#                    auto_size=True,
#                    linear_interpolation=True
#                   )

# img_actin = np.array(img_actin)

# Actin Surface Visualization with Napari Vedo
# viewer = napari.Viewer(ndisplay=3)
# actin_surface = nppas.label_to_surface(actin_merged_np)
# viewer.add_surface(actin_surface)








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
from scipy import ndimage as ndi

#device = cle.select_device("gfx1032")
#device = cle.select_device("AMD Radeon Pro W6600")


from topostats.tracing.tracingfuncs import getSkeleton
from skimage.morphology import skeletonize, thin

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

        count = 0
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
    
    statistics_surf_actin = nsitk.label_statistics(intensity_image=actin_img, label_image=act_obj,
                                            size = True, intensity = True, perimeter= True,
                                            shape= True, position= True)#, moments= True)
    #statistics_surf_actin = cle.statistics_of_labelled_pixels(actin_img, act_obj)

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

from scipy import ndimage
# def connected_components_labeling_box(binary_input):
#     labeled_array, nb_labels = ndimage.label(binary_input)
#     sizes = ndimage.sum(binary_input, labeled_array, range(nb_labels + 1))
#     mask = sizes == max(sizes)
#     print(mask)
#     binary_input = mask[labeled_array]
#     print(binary_input)
#     return binary_input

def connected_components_labeling_box(binary_input):
    binary_input = binary_input.astype('int32')
    labeled_array, num_labels = ndimage.label(binary_input)
    return labeled_array


import numpy as np
from scipy.ndimage import binary_erosion



def filter_images(images, threshold_mean, threshold_std):
    filtered_images = []
    for img in images:
        # Apply Otsu thresholding
        threshold_value = nsitk.threshold_otsu(img)
        binary_img = img > threshold_value

        img_mean = np.mean(img)
        img_std = np.std(img)

        if img_mean >= threshold_mean and img_std >= threshold_std:
            filtered_images.append(binary_img)
    return filtered_images


def compute_threshold_mean(binary_mask, original_image):
    masked_image = np.where(binary_mask, original_image, 0)  # Apply binary mask to the original image
    threshold_mean = np.mean(masked_image)  # Compute the mean intensity within the masked region
    return threshold_mean

def compute_threshold_std(binary_mask, original_image):
    masked_image = np.where(binary_mask, original_image, 0)  # Apply binary mask to the original image
    threshold_mean = np.std(masked_image)  # Compute the mean intensity within the masked region
    return threshold_mean


from skimage.morphology import binary_erosion

# Now, eroded_labels will have labels reduced to single pixels

def erosion_to_single_pixels(input_image):
    """
    Reduce a segmented binary image to single pixels using erosion.

    Args:
        binary_image (ndarray): Binary image where the object of interest is represented by non-zero values.

    Returns:
        ndarray: Image with only single pixels representing the object.
    """
    labeled_image = input_image
    
    # Convert the labeled image to a binary mask
    binary_mask = labeled_image > 0
    
    # Perform erosion operation on the binary image
    eroded_mask = binary_erosion(binary_mask)

    # Convert the eroded mask back to labeled image
    eroded_labels = np.zeros_like(labeled_image)
    eroded_labels[eroded_mask] = labeled_image[eroded_mask]
    
    return eroded_labels


folder_path = 'E:/Quantified/eNUC Live/230628 enuc vca vs novca/Test_Center/'

get_files = folder_scan(folder_path)

# Conversion to um
px_to_um_X = 0.0351435
px_to_um_Y = 0.0351435
px_to_um_Z = 1 #1.0/0.0351435 # making pixels isotropic



for image in get_files:

    aics_image = AICSImage(image) # Read indivi
    filename = Path(image).stem
# <Dimensions [T: 1, C: 1, Z: 256, Y: 256, X: 256]> # ("ZYX", T=0)
    no_chnls = aics_image.dims.C
    if no_chnls == 1:
        img_nuclei = aics_image.get_image_data("ZYX", T=0, C=0) 
    elif no_chnls == 2: 
        img_dextran = aics_image.get_image_data("ZYX", T=0, C=0) # A Channel
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
    
    if not os.path.exists(folder_path + '/Result'):
        os.makedirs(folder_path + '/Result')
        
    Result_folder = folder_path + '/Result/' 
    # For saving results #.stem

    # Upscale the Z dimension to be isotropic

    # channel to segment and nuclear channel 
    # numbering starts at 1 
    # for your single channel image use [0, 0] 
    # for the multi channel image it's [3, 0]
    ##############################################
    # channels = [0, 0] 
    # diameter = 170 #model diam_labels  165.257 (mean diameter of training ROIs)
    # use_GPU = True
    # stitch_threshold=1
    # #do_3D=True
    # #resample=True
    # cellprob_threshold = 0

    # pretrained_model = "C:/Users/kaabi/Documents/Nuceli_Data/Enucleation/Cellpose/Models_live_eNUC/CP_20230306_150448_200"

    # model_match_threshold = 25 #30
    # flow_threshold = (31.0 - model_match_threshold) / 10.

    # logger = io.logger_setup()
    # print('Segmenting Image ---', filename)
    # model = models.CellposeModel(gpu=use_GPU, model_type=pretrained_model)#,residual_on=False, style_on=False, concatenation=False)#, resample=True, do_3D=True,stitch_threshold=1)
    # mask, flows,styles = model.eval(img_nuclei,
    #                             channels =channels,
    #                             #anisotropy=anisotropy,
    #                             diameter=diameter, 
    #                             #pretrained_model=pretrained_model,
    #                             do_3D=False,
    #                             resample=True,
    #                             min_size = -1,
    #                             flow_threshold =flow_threshold,
    #                             cellprob_threshold = cellprob_threshold,
    #                             net_avg=True,
    #                             stitch_threshold=stitch_threshold,
    #                             #model_match_threshold = model_match_threshold
    #                                 )
    ##################################################
    # Merge Segmented Mask
    print('The Filename ---', filename)
    
    mask_path = folder_path + '230628 D0enuc dex vca 0xArp3 actin 1.tif'
    mask = imread(mask_path)

    merged_Labels = connected_components_labeling_box(mask)
    merged_Labels_np = np.array(merged_Labels)#.astype('int32')
    
    # ###############################
    # ### Save the Mask
    # prediction_stack_32 = img_as_float32(merged_Labels_np, force_copy=False)     
    # os.chdir(Result_folder)
    # imwrite(filename+".tif", prediction_stack_32)
    # #merged_Labels_np = split_touching_objects(merged_Labels_np)
    # ###############################
    
    # Structure  for Actin Dilation
    diamond = np.zeros((3, 3, 3), dtype=bool)
    diamond[1:2, 0:3, 0:3] = True

    #print('Number_of_Nuclei',len(np.unique(merged_Labels))-1) first is mask

    from skimage import measure
    #labels = measure.label(merged_Labels_np)

    viewer = napari.Viewer()
    
    
    #viewer.add_labels(labels)
    
    for i, label in enumerate(np.unique(mask)[1:]):
        print('i >>>>', i)
        #print('label ---', label)
        if label in mask:
            print('label >>>>', label)
            # create a mask for the current label
            #if (len(np.unique(merged_Labels)-1) == 1):# and no_chnls>1:
                # break    
            mask = mask == label
            ###nuc_lbl = np.array(label)
            intensity_nucelus= replace_intensity(mask, img_nuclei)
            
            viewer.add_image(intensity_nucelus)
            
            image1_gb = nsitk.gaussian_blur(intensity_nucelus, 2.0, 2.0, 0.0)
            # threshold isodata
            image2_T = nsitk.threshold_isodata(image1_gb)

            # Reducing Labels to Edges
            
            
            #from scipy.ndimage import binary_fill_holes
            #image2_FH = binary_fill_holes(image2_T)
            image2_Li = nsitk.threshold_li(image1_gb)
            image2_BF = nsitk.binary_fill_holes(image2_Li)
            
            nuc_border = np.zeros(img_nuclei.shape)
            for i in range(image2_BF.shape[0]):                
                thinned_slice = cle.reduce_labels_to_label_edges(image2_BF[i]) #cle.reduce_labels_to_label_edges(image2_FH[i])
                nuc_border[i] = thinned_slice
                
            nuc_border_EC = cle.exclude_small_labels(nuc_border, maximum_size=250)
            nuc_border_T = np.array(nuc_border_EC)
            viewer.add_image(nuc_border_T) 
            # connected component labeling
            image3_C = nsitk.connected_component_labeling(image2_T)
        
            #statistics_nucleus = nsitk.label_statistics(intensity_image=intensity_nucelus, label_image=mask,
            #                                     size = True, intensity = True, perimeter= True,
            #                                     shape= True, position= True)#, moments= True)
        
            # nuclei_Area = statistics_nucleus.loc[0, 'number_of_pixels']
            # nuclei_Area = nuclei_Area*px_to_um_X*px_to_um_Y*px_to_um_Z
            # statistics_nucleus['Nucleus Area'] = nuclei_Area 
            
            
            
            # viewer.add_image(mask)
            
                # # Chromocenter Segmentation

            intensity_map_norm = normalize_intensity(intensity_nucelus)
                #viewer.add_image(intensity_map_norm)
                # Gaussian blur
            image1_gb = cle.gaussian_blur(intensity_map_norm, None, 2.0, 2.0, 0.0)
            image2_thb = cle.top_hat_box(image1_gb, None, 15.0, 15.0, 0.0)
                # subtract gaussian background
            image3_sgb = cle.subtract_gaussian_background(image2_thb, None, 10.0, 10.0, 0.0)
            

            intermodes_Threshold_Chromo = nsitk.threshold_intermodes(image3_sgb)
            viewer.add_labels(intermodes_Threshold_Chromo)
                    
            # Chromocenter Segmentation
            #viewer.add_image(intermodes_Threshold_Chromo)
            ###############################
            # statistics_Chromo = nsitk.label_statistics(intensity_image=img_nuclei, label_image=intermodes_Threshold_Chromo,
            #                             size = True, intensity = True, perimeter= True,
            #                             shape= True, position= True)#, moments= True)
            
                            
            # intermodes_chromo_Area = statistics_Chromo.loc[0, 'number_of_pixels']
            # chromointermodes_Area = intermodes_chromo_Area *px_to_um_X* px_to_um_Y
            # statistics_Chromo['Chromocenter Area'] = chromointermodes_Area 
            # print('Chromocenter_Area >>>>', chromointermodes_Area)
            # pd.DataFrame(statistics_Chromo).to_excel(Result_folder+'(Chromo)_' + filename + '_' + str(i)+'.xlsx')
                    
            # # Save Chromocenter Mask
            # ### Save the Mask
            # chromo_stack_32 = img_as_float32(intermodes_Threshold_Chromo, force_copy=False)     
            # os.chdir(Result_folder)
            # imwrite("(Intermodes_chromocenter)" +filename+".tif", chromo_stack_32)

                    ###############################
                #  # # Chromocenter Segmentation   with Max Entropy
                # intensity_gauss = nsitk.gaussian_blur(intensity_map_norm, variance_x=2, variance_y=2, variance_z=0)
                # intensity_maxEntrp = nsitk.threshold_maximum_entropy(intensity_gauss)
                # ### Save the Mask
                # chromo_stack_32 = img_as_float32(intensity_maxEntrp , force_copy=False)     
                # os.chdir(Result_folder)
                # imwrite("(Intermodes_chromocenter)" +filename+".tif", chromo_stack_32) 
                    
                # Actin Segmentation

            act_obj = np.zeros(img_nuclei.shape)
            act_erode = np.zeros(img_nuclei.shape)
            dilated = ndi.binary_dilation(mask, diamond, iterations=10).astype(mask.dtype)
            actin_img = replace_intensity(dilated, img_actin)
            actin_img_N = normalize_intensity(actin_img)
            
            viewer.add_image(actin_img_N)
            
            # Thresholding Actin 
            # Earlier tried cle.gaussian_blur(actin_img_N, None, 2.0, 2.0, 0.0), cle.top_hat_box(actin_filter, None, 15.0, 15.0, 0.0)
            # and nsitk.threshold_otsu(actin_filter_thb)
            
            actin_filter = nsitk.median_filter(actin_img_N, radius_x=2, radius_y=2, radius_z=0)
            actin_binary_GS = cle.subtract_gaussian_background(actin_filter, None, 15.0, 15.0, 0.0)
            ####actin_filter = cle.gaussian_blur(actin_img_N, None, 2.0, 2.0, 0.0)
            ####actin_filter_thb = cle.top_hat_box(actin_filter, None, 15.0, 15.0, 0.0)
            ##actin_binary_GS = cle.subtract_gaussian_background(actin_filter_thb, None, 10.0, 10.0, 0.0)
            image2_Otsu = nsitk.threshold_otsu(actin_binary_GS)
            actin_binary =nsitk.binary_fill_holes(image2_Otsu)
            
            
            act_obj = np.zeros(img_nuclei.shape)
            #act_obj_2 = np.zeros(img_nuclei.shape)
            
            #actin_binary = nsitk.threshold_otsu(actin_filter)
            # Apply Otsu thresholding
            ####################
            # Condition to remove background Actin segmentation
            threshold_mean = compute_threshold_mean(actin_binary,actin_img)
            threshold_std = compute_threshold_std(actin_binary,actin_img)
            
            img_mean = np.mean(actin_img)
            img_std = np.std(actin_img)

            # Apply thinning function to each slice of the actin binary image
            # for i in range(actin_binary.shape[0]):
            #     if img_mean >= threshold_mean and img_std >= threshold_std:
            #         thinned_slice = thin(actin_binary[i])
            #         act_obj[i] = thinned_slice
            # from topostats.tracing.tracingfuncs import getSkeleton
            
            # smoothed_grain = ndimage.binary_dilation(actin_binary , iterations=1).astype(actin_binary.dtype)
            
            # dna_skeleton = getSkeleton(actin_binary_GS, smoothed_grain, actin_binary.shape[1],
            #                            actin_binary.shape[2],px_to_um_X)
            opImBase = img_nuclei[0,:,:]
            
            
            for i in range(actin_binary.shape[0]):
                if img_mean >= threshold_mean and img_std >= threshold_std:
                    prune = pruneSkeleton(actin_binary[i], opImBase)
                    #prune_labels = cle.exclude_small_labels(prune, maximum_size=50)
                    #thinned_slice = thin(actin_binary[i])
                    act_obj[i] = prune

            # Compute the standard deviation for each slice in the binary mask
            std_per_slice = np.std(act_obj, axis=(1, 2))

            # Calculate the ratio of the first and second value of the standard deviation
            ratio = std_per_slice[0] / std_per_slice[1]

            # Find the indices of slices that have a ratio higher than a threshold
            threshold_ratio = 5.0  # Set your desired threshold ratio here
            indices_to_remove = np.where(std_per_slice[1:] / std_per_slice[:-1] > threshold_ratio)[0]  + 1

            # Set the values of the slices to zero instead of deleting them
            act_obj[indices_to_remove] = 0

            # Filter out the zero slices from the array
            #act_obj_filtered = act_obj[np.any(act_obj, axis=(1, 2))]

                
            statistics_actin = nsitk.label_statistics(intensity_image=img_actin, label_image=act_obj,
                                                size = True, intensity = True, perimeter= True,
                                                shape= True, position= True)#, moments= True)
            
            actin_no_of_pixels = statistics_actin.loc[0, 'number_of_pixels']
            actin_Perimeter = actin_no_of_pixels*px_to_um_X
            
            # for i in range(actin_binary.shape[0]):                
            #     thinned_slice = reduce_to_single_pixels(actin_binary[i])
            #     act_erode[i] = thinned_slice
            print("Actin Found")
            viewer.add_labels(actin_binary)    
            viewer.add_image(act_obj)
            #viewer.add_image(act_obj_2)
            
            for i, (act_Z, nucl_Z) in enumerate(zip(act_obj, nuc_border_T)):
                actin = np.sum(act_Z)
                nucleus = np.sum(nucl_Z)
                
                print('Actin pixels', actin)
                print('Nucleus pixels',nucleus)
                Coverage = []
                indiv_Cov = (actin/nucleus)*100
                Coverage.append(indiv_Cov)
                #print(indiv_Cov)
                
