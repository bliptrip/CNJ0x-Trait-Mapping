#!/usr/bin/env python3
# Author: Andrew Maule
# Objective: Read in berry shape templates (binary grayscale image files), determine all berry areas, and normalize all images
#            such that all berries are the same area and have their center of masses in the center of the image.  This is
#            done to allow distance transform 'averaging' of berry shapes based on berry shape categories input by the user.
#
# Output: Derived binary grayscale image files of all input berry template images, with identical areas, the same center of mass,
#         and equivalent resolutions.
# import the necessary packages
import argparse
import cv2
import math
import numpy as np
from pathlib import Path
from shapely.geometry import Polygon
from shapely.affinity import scale, translate
import sys

def convertToPolygon(mask):
    contours, hierarchy = cv2.findContours(mask.astype('uint8'), cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)
    contours_stacked    = np.vstack(np.vstack(contours[0]))
    polygon             = Polygon(contours_stacked)
    return(polygon)

def convertToContours(polygon):
    (x,y)     = polygon.exterior.coords.xy
    contours    = np.zeros((len(x),1,2), dtype=int)
    for i in range(0,len(x)):
        contours[i,0,:] = np.array([x[i],y[i]], dtype=int)
    return [contours]

# construct the argument parser and parse the arguments
def parse_args():
    parser = argparse.ArgumentParser(prog=sys.argv[0], description="A command-line utility for normalizing areas and centering template berry shapes.")
    parser.add_argument('-i', '--input', action='append', required=True, help="Input binary grayscale image files containing berry shape template.")
    parser.add_argument('-o', '--output', action='store', required=True, help="Output folder for normalized/centered/resized images.  Images will have the same name in output folder as in input folder, so they need to be different paths.")
    parsed = parser.parse_args(sys.argv[1:])
    return(parsed)


if __name__ == '__main__':
    parsed      = parse_args()
    images      = []
    polys       = []
    for i,image_file in enumerate(parsed.input):
        #Snap polygons with each other
        image       = cv2.imread(image_file, cv2.IMREAD_GRAYSCALE)
        images.append(image)
        polys.append(convertToPolygon(image))
    image_dims          = np.array([i.shape for i in images])
    image_dims_padded   = 2 * np.max(image_dims, axis=0)
    areas               = np.array([p.area for p in polys])
    areas_mean          = np.mean(areas)
    scale_factors       = np.sqrt(areas_mean/areas)
    polys_scaled        = []
    for i,p in enumerate(polys):
        polys_scaled.append(scale(p, xfact=scale_factors[i], yfact=scale_factors[i], origin='center'))
    centroids       = np.array([[p.centroid.y,p.centroid.x] for p in polys_scaled])
    shift = (image_dims_padded/2) - centroids #Take halfway point of new padded images
    polys_translated = []
    for i,p in enumerate(polys_scaled):
        polys_translated.append(translate(p, xoff=shift[i,1], yoff=shift[i,0]))
    images_out = []
    for i,p in enumerate(polys_translated):
        contours    = convertToContours(p)
        iout        = cv2.drawContours(np.zeros(image_dims_padded, dtype='uint8'), contours, -1, (255), cv2.FILLED)
        basename    = Path(parsed.input[i]).name
        newname     = Path(parsed.output) / basename
        cv2.imwrite(str(newname), iout); 
        
        

    
#    images_padded       = []
#    for i,image in enumerate(images):
#        x_left          = int((image_dims_padded[0]-image_dims[i][0])/2)
#        x_right         = image_dims_padded[0]-image_dims[i][0]-x_left
#        y_left          = int((image_dims_padded[1]-image_dims[i][1])/2)
#        y_right         = image_dims_padded[1]-image_dims[i][1]-y_left
#        image_padded    = np.pad(images[i], ((x_left, x_right), (y_left, y_right)))
#        images_padded.append(image_padded)
