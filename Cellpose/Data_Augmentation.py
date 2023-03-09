import random
import glob
import os
import cv2
import uuid
from matplotlib import pyplot as plt
import tifffile as tiff
import skimage.io

import albumentations as A
from albumentations import OneOf

def visualize(image):
    # Divide all values by 65535 so we can display the image using matplotlib
    image = image / 65535
    plt.figure(figsize=(10, 10))
    plt.axis('off')
    plt.imshow(image)


transform = A.Compose([
    A.ToFloat(max_value=65535.0),
    A.HorizontalFlip(p=0.5),
    A.RandomRotate90(),
    A.FromFloat(max_value=65535.0),
])

transform1 = A.Compose([
    A.ToFloat(max_value=65535.0),
    A.Flip(),
    A.RandomRotate90(),
    A.OneOf([
            A.MotionBlur(p=0.2),
            A.MedianBlur(blur_limit=5, p=0.1),
            A.Blur(blur_limit=5, p=0.1),
        ], p=0.2),
    A.HueSaturationValue(hue_shift_limit=40, sat_shift_limit=0.1, val_shift_limit=0.1, p=0.3),

    A.FromFloat(max_value=65535.0),
])
    
transform2 = A.Compose([
    A.ToFloat(max_value=65535.0),

    A.RandomRotate90(),
    A.Flip(),
    A.OneOf([
        #A.MotionBlur(p=0.2),
        #A.MedianBlur(blur_limit=3, p=0.1),
        A.Blur(blur_limit=1, p=0.1),
    ], p=0.2),
    A.ShiftScaleRotate(shift_limit=0.0625, scale_limit=0.2, rotate_limit=45, p=0.2),
    A.OneOf([
        A.OpticalDistortion(p=0.3),
        A.GridDistortion(p=0.1),
    ], p=0.2),        
    A.HueSaturationValue(hue_shift_limit=20, sat_shift_limit=0.1, val_shift_limit=0.1, p=0.3),

    A.FromFloat(max_value=65535.0),
])

images = []
def folder_Scan(directory):
    os.chdir(directory)
    extension = 'tif'
    path_local = glob.glob('*.{}'.format(extension))
    for g_name in path_local:
        images.append(g_name)


folder_Scan('C:/Users/kaabi/Documents/Nuceli_Data/Enucleation/Cellpose/Sorted/')
folder ='C:/Users/kaabi/Documents/Nuceli_Data/Enucleation/Cellpose/Sorted/'

for g_name in images:
    img = cv2.imread(os.path.join(folder,g_name), cv2.IMREAD_UNCHANGED)
    img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    augmented = transform(image=img)  
    #visualize(augmented['image'])
    filename = str(uuid.uuid4())
    #####skimage.io.imsave(folder+filename+'.tif', augmented['image'], plugin='tifffile')
    tiff.imwrite(folder+filename+'.tif', augmented['image'])
    augmented1 = transform1(image=img)
    filename = str(uuid.uuid4())
    tiff.imwrite(folder+filename+'.tif', augmented1['image'])
    augmented2 = transform2(image=img)
    filename = str(uuid.uuid4())
    tiff.imwrite(folder+filename+'.tif', augmented2['image'])
    #cv2.imwrite(folder+filename.tif, img)

#for g_name in range(10):
