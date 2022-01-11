#!/usr/bin/env python3
# Author: Andrew Maule
# Objective: Generate precursor images for demonstrating the berry chimera transformation process.

import argparse
import cv2
import logging
from matplotlib.colors import Normalize
import numpy as np
import seaborn as sns
import sys
from ztools.dip.chimera import Chimera

# construct the argument parser and parse the arguments
def parse_args():
    parser = argparse.ArgumentParser(prog=sys.argv[0], description="A command-line utility for generating GiNA-like berry parameters from a set of individual berry images and merging into collated dataset.")
    parser.add_argument('-p', '--path', action='store', type=str, required=True, help="Input folder path containing normalized berry template representations of different berry shapes.")
    parser.add_argument('-o', '--output', action='store', type=str, required=True, help="Output folder to write images to.")
    parser.add_argument('-m', '--map', action='store', default="{'round': '*fruit_template_round_binary.png', 'oblong': '*fruit_template_oblong_binary.png', 'oval': '*fruit_template_oval_binary.png', 'pyriform': '*fruit_template_pyriform_binary.png', 'spindle': '*fruit_template_spindle_binary.png'}", help="Dictionary mapping shape categories to template binary image file (glob-patterns allowed).")
    parser.add_argument('--omap', action='store', default="{'round': 'fruit_template_round_diffdistance.png', 'oblong': 'fruit_template_oblong_diffdistance.png', 'oval': 'fruit_template_oval_diffdistance.png', 'pyriform': 'fruit_template_pyriform_diffdistance.png', 'spindle': 'fruit_template_spindle_diffdistance.png'}", help="Dictionary mapping shape categories to template binary image file (glob-patterns allowed).")
    parser.add_argument('-l', '--level', type=str, default="WARNING", choices=["DEBUG","INFO","WARNING","ERROR","CRITICAL"], help="Set the debug log level.")
    parser.add_argument('shapes', metavar='shapes', type=str, nargs='+', help="Sequence of input shape categories to merge into chimera.  Valid values: {}")
    parsed = parser.parse_args(sys.argv[1:])
    return(parsed)

def diffDistanceColorize(diffDistance):
    c = np.delete(palette(nmlz(diffDistance)), 3, axis=2) #Converts to sequence of rgb values for create a color-gradient image - remove alpha channel
    ci = np.array(c * 255.0, dtype="uint8")
    return(ci)

if __name__ == '__main__':
    logger      = logging.getLogger()
    parsed      = parse_args()
    decoded_level = eval("logging.{}".format(parsed.level))
    parsed = parse_args()
    palette = sns.diverging_palette(h_neg=50,h_pos=240,s=100,l=50,center="dark",as_cmap=True)
    nmlz = Normalize()
    omap = eval(parsed.omap)
    #First generate unitary diff distance transform images
    for shape in np.unique(parsed.shapes):
        berryc = Chimera()
        berryc.loadMap(parsed.path,parsed.map)
        berryc.compose([shape])
        composite_raw_normalized = diffDistanceColorize(berryc.composite_raw)
        cv2.imwrite(parsed.output + "/" + omap[shape], composite_raw_normalized)
    #Now generate composite diff distance transform image and composite binary image
    berryc = Chimera()
    berryc.loadMap(parsed.path,parsed.map)
    berryc.compose(parsed.shapes)
    composite_raw_normalized = diffDistanceColorize(berryc.composite_raw)
    cv2.imwrite(parsed.output + "/fruit_chimera_binary.png", berryc.composite)
    cv2.imwrite(parsed.output + "/fruit_chimera_diffdistance.png", composite_raw_normalized)

