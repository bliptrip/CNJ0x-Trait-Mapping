#!/usr/bin/env python3
# Author: Andrew Maule
# Objective: Given a set of input unsigned manhattan chain code (UMCC) numpy array files, generate the numerical values for
# UMCCx and UMCCy for each shape type, and plot differences on 2D cartesian plot.
#
# Output: A chain code version of chimera.
#         
# import the necessary packages
import argparse
import logging
import math
import matplotlib.pyplot as plt
import numpy as np
import sys

from berry_umcc import UMCCEncoder

def parse_args():
    parser = argparse.ArgumentParser(prog=sys.argv[0], description="A command-line utility for generating a chain code representation of berry shapes.")
    parser.add_argument('-i', '--input', action='store', type=str, required=True, help="Input folder prefix path containing UMCC numpy array files to open.")
    parser.add_argument('--names', '-n', type=str, nargs='+', help="Labels to apply to each input file that is plotted on cartesian plane -- same order as input files.")
    parser.add_argument('--files', '-f', type=str, nargs='+', help="Sequence of UMCC numpy array files to open, one per shape type.")
    parsed = parser.parse_args(sys.argv[1:])
    return(parsed)

if __name__ == '__main__':
    coords_x = []
    coords_y = []
    coords_both = []
    coords_div = []
    logger = logging.getLogger()
    parsed = parse_args()
    for f in parsed.files:
        umcc    = np.load(parsed.input + '/' + f, allow_pickle=True)
        umccO   = UMCCEncoder(logger,U=umcc)
        umccOn  = umccO.numeric()
        coords_x.append(math.log(umccOn['Ux']))
        coords_y.append(math.log(umccOn['Uy']))
        coords_both.append(math.log(umccOn['Ux'] * umccOn['Uy']))
        coords_div.append(math.log(umccOn['Uy']) - math.log(umccOn['Ux']))
    #plt.scatter(coords_both, [1]*len(coords_both))
    #plt.scatter(coords_x, coords_y)
    plt.scatter(coords_div, [1]*len(coords_div))
    for i,n in enumerate(parsed.names):
        #plt.text(coords_x[i]+0.2,coords_y[i]+0.2,n)
        plt.text(coords_div[i],1.002,n)
    plt.show()
