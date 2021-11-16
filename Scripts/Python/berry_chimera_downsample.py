#!/usr/bin/env python3
# Author: Andrew Maule
# Objective: Downsamples a chimera berry image in order to reduce the number of contours needed for describing the outline.
#
# Output: A downsampled version of the chimera berry image.
#         
# import the necessary packages
import argparse
import cv2
import sys

# construct the argument parser and parse the arguments
def parse_args():
    parser = argparse.ArgumentParser(prog=sys.argv[0], description="A command-line utility for downsampling a chimera berry image in order to reduce number of contours needed for describing the berry contour.")
    parser.add_argument('-i', '--input', action='store', type=str, required=True, help="Path to input berry chimera image file.")
    parser.add_argument('-o', '--output', action='store', type=str, required=True, help="Path to output file containing downsampled berry chimera image.")
    parser.add_argument('-d', '--downsample', action='store', type=int, default=10, help="Factor to downsample each axis by.")
    parsed = parser.parse_args(sys.argv[1:])
    return(parsed)


if __name__ == '__main__':
    parsed = parse_args()
    iimage = cv2.imread(parsed.input, cv2.IMREAD_GRAYSCALE)
    oimage = cv2.resize(iimage, dsize=(0,0), fx=1/parsed.downsample, fy=1/parsed.downsample)
    cv2.imwrite(parsed.output, oimage)
