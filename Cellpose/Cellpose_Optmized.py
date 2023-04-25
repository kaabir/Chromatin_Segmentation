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
channels = [0, 0] 
diameter = 165.257 #169.708
use_GPU = True
stitch_threshold=1
#do_3D=True
#resample=True
cellprob_threshold = 0

pretrained_model = "C:/Users/kaabi/Documents/Nuceli_Data/Enucleation/Cellpose/models_trained_Nuclei/CP_20230306_150448_200"

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
#merged_Labels = cle.connected_components_labeling_box(mask)
merged_Labels_np = np.array(mask)
#merged_Labels_np = np.array(merged_Labels).astype('int32')
#merged_Labels_np = split_touching_objects(merged_Labels_np)

# Getting Original intensities
#intensity_vector = cle.read_intensities_from_map(mask, img)
#statistics = cle.statistics_of_labelled_pixels(img_nuclei, merged_Labels)

#nuclei_Area = statistics['area']*px_to_um_Y*px_to_um_Y*px_to_um_Z

# Chromocenter Segmentation
# Get voxel of pixel intensity values inside the mask 

def nuceleus_intensity(mask, img):
    if not (mask.shape == img.shape):
        return False

    nuceleus_intensity = np.where(np.logical_and(mask,img),img,0) # Overlap Mask on original image (as reference)
    #chromocenter_intensity = np.array(cle.replace_intensities(mask,img))
    return nuceleus_intensity

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


def calculate_surface_area(mask, threshold=None):
    # generate surface mesh using marching cubes algorithm
    verts, faces, _, _ = measure.marching_cubes(mask)
    # calculate surface area using mesh surface area function
    surface_area = measure.mesh_surface_area(verts, faces)
    #surface_area = surface_area*px_to_um_Y*px_to_um_Y*px_to_um_Z

    return surface_area

# structure for actin segmentation
from skimage import measure
from scipy import ndimage
from scipy.ndimage import gaussian_filter
from skimage import morphology
from skimage.morphology import binary_dilation
from scipy import ndimage as ndi

diamond = np.zeros((3, 3, 3), dtype=bool)
diamond[1:2, 0:3, 0:3] = True

#print(len(np.unique(merged_Labels)))
number_of_Nuclei = len(np.unique(merged_Labels_np))-1 # first is mask

# Create an image to place a single nuleus with intensity voxel
im_obj = np.zeros(img_nuclei.shape) 

nuc_lbl_lst1 = []

from skimage import measure
labels = measure.label(merged_Labels_np)

for i, label in enumerate(np.unique(labels)[1:]):
    
    if label in labels:
    # create a mask for the current label
        mask = labels == label
        #nuc_lbl = np.array(label)
        intensity_map= nuceleus_intensity(mask, img_nuclei)
        statistics_nucleus = cle.statistics_of_labelled_pixels(intensity_map, mask)
        pd.DataFrame(statistics_nucleus).to_excel(Result_folder+filename+'_'+ str(i)+'_(Nucleus)Nucleus_statistics.xlsx')
    
        nuclei_Area = statistics_nucleus['area']*px_to_um_Y*px_to_um_Y*px_to_um_Z
        print('nuclei_Area ---', nuclei_Area)
        # Chromocenter Segmentation

        intensity_map_norm = normalize_intensity(intensity_map)
    
        intensity_map_blur = nsitk.median_filter(intensity_map_norm, radius_x=2, radius_y=2, radius_z=0)
        intensity_map_thb = cle.top_hat_box(intensity_map_blur, None, 10.0, 10.0, 0.0)
        intermodes_Threshold_Chromo = nsitk.threshold_intermodes(intensity_map_thb)
        statistics_intermodes_chromo = cle.statistics_of_labelled_pixels(img_nuclei, intermodes_Threshold_Chromo)
        pd.DataFrame(statistics_intermodes_chromo).to_excel(Result_folder+filename+'_'+ str(i)+'_(Chromo)Chromo_statistics.xlsx')
    
        chromointermodes_Area = np.sum(statistics_intermodes_chromo['area'], axis=0)
        chromointermodes_Area = chromointermodes_Area *px_to_um_X* px_to_um_Y*px_to_um_Z
        print('chromointermodes_Area ---', chromointermodes_Area)
        # Actin Segmentation
        act_obj = np.zeros(img_nuclei.shape)
        if 'img_actin' in globals(): #or 'img_actin' in globals():
            # Separate Ctrl
            pass
        
            dilated = ndi.binary_dilation(mask, diamond, iterations=10).astype(mask.dtype)
            actin_img = nuceleus_intensity(dilated,img_actin)
            actin_filter = nsitk.median_filter(actin_img, radius_x=2, radius_y=2, radius_z=0)
            actin_binary = nsitk.threshold_otsu(actin_filter)
    
            for i in range(actin_binary.shape[0]):
                thinned_slice = thin(actin_binary[i])
                act_obj[i] = thinned_slice

            statistics_surf_actin = cle.statistics_of_labelled_pixels(actin_img, act_obj)
    
            actin_surf_Area = np.sum(statistics_surf_actin['area'], axis=0)
            actin_surf_Area = actin_surf_Area*px_to_um_X
            print(actin_surf_Area)
            pd.DataFrame(statistics_surf_actin).to_excel(Result_folder+filename+'_'+ str(i)+'_(Actin)Actin_statistics.xlsx')
                 
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
