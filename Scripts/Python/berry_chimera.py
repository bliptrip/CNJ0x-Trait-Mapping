#!/usr/bin/env python
# Author: Andrew Maule
# Objective: Using normalized berry templates, merge a set of categorical berry shape descriptors to create a chimera shape.  Uses 
# distance transforms to generate the chimera.
#
# Output: A chain code version of chimera.
#         
# import the necessary packages
import argparse
import cv2
from glob import *
import matplotlib.pyplot as plt
import numpy as np
from skimage import measure
import sys

shape_choices=['round', 'oblong', 'oval', 'pyriform', 'spindle']
stdout_default='-'

# construct the argument parser and parse the arguments
def parse_args():
    parser = argparse.ArgumentParser(prog=sys.argv[0], description="A command-line utility for generating a chain code representation of berry shapes.")
    parser.add_argument('-p', '--path', action='store', type=str, required=True, help="Input folder path containing normalized berry template representations of different berry shapes.")
    parser.add_argument('-m', '--map', action='store', default="{'round': '*round*.small.png', 'oblong': '*oblong*.small.png', 'oval': '*oval*.small.png', 'pyriform': '*pyriform*.small.png', 'spindle': '*spindle*.small.png'}", help="Dictionary mapping berry shapes to berry template filename (glob-patterns allowed).")
    parser.add_argument('-o', '--output', action='store', default=stdout_default, help="Output file to write chaincode representation to ('{}' for stdout)".format(stdout_default))
    parser.add_argument('shapes', metavar='shapes', type=str, nargs='+', choices=shape_choices, help="Sequence of input shapes to merge into chimera.  Valid values: {}".format(shape_choices))
    parsed = parser.parse_args(sys.argv[1:])
    return(parsed)


if __name__ == '__main__':
    templates   = {}
    dtransforms = {}
    parsed = parse_args()
    shapes = parsed.shapes
    parsed.map = eval(parsed.map)
    for k in parsed.map:
        names           = glob("{}/{}".format(parsed.path,parsed.map[k]))
        templates[k]    = cv2.imread(names[0], cv2.IMREAD_GRAYSCALE)
        dtransforms[k]  = cv2.distanceTransform(templates[k], distanceType=cv2.DIST_L2, maskSize=cv2.DIST_MASK_PRECISE) - cv2.distanceTransform(np.bitwise_not(templates[k]), distanceType=cv2.DIST_L2, maskSize=cv2.DIST_MASK_PRECISE)
    dsums = dtransforms[shapes[0]]
    for s in shapes[1:]:
        dsums += dtransforms[s]
    dsums = ((dsums / len(shapes)) > 0).astype(dtype='uint8')
    dsums[dsums > 0] = 255
    if( parsed.output != stdout_default ):
        cv2.imwrite(parsed.output, dsums)
    contours = measure.find_contours(dsums)[0]
    print("Number of contours: {}".format(len(contours)))
    plt.imshow(dsums)
    plt.plot(contours[:, 1], contours[:, 0], '-g', linewidth=2)
    plt.show()
