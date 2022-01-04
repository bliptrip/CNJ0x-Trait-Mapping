#!/usr/bin/env python3
# Author: Andrew Maule
# Objective: Given a set of berry shape classifications and their templates, generates a chimera of these
#               shapes and calculates a subset of GiNA parameters on them, in addition to the slope chain
#               code parameter tortuosity, along with unsigned manhattan chain code x/y outputs as binary
#               strings.

import argparse
import cv2
import logging
import math
import numpy as np
import pandas as pd
import re
import sys
from ztools.dip.binary_segment import BinaryThresholdSegment
from ztools.dip.chaincodes.umcc import UMCCEncoder
from ztools.dip.chaincodes.scc import SCCEncoder
from ztools.dip.chimera import Chimera

# construct the argument parser and parse the arguments
def parse_args():
    parser = argparse.ArgumentParser(prog=sys.argv[0], description="A command-line utility for generating GiNA-like berry parameters from a set of individual berry images and merging into collated dataset.")
    parser.add_argument('-d', '--data', action='store', type=str, required=True, help="Input collated upright trait dataset.")
    parser.add_argument('-o', '--odata', action='store', type=str, required=True, help="Output collated upright trait dataset with chimeric shape parameters injected.")
    parser.add_argument('-l', '--level', type=str, default="WARNING", choices=["DEBUG","INFO","WARNING","ERROR","CRITICAL"], help="Set the debug log level.")
    parser.add_argument('--audit-path', action='store', type=str, help="If --audit-path is defined, this is where the chimeric binary images and their contours will be stored for subsequent auditing.")
    parser.add_argument('--min_area', '--mina', '--minArea', dest='mina', type=int, default='100', help="Remove foreground blobs (fruit) less than this size.")
    parser.add_argument('--max_area', '--maxa', '--maxArea', dest='maxa', type=int, default='3000000', help="Remove foreground blobs (fruit) greather than this size.")
    parsed = parser.parse_args(sys.argv[1:])
    return(parsed)

def cleanupShapes(s):
    try:
        if np.nan == s:
            return(None)
        elif re.search(r'\*|-',s):
            return(None)
        elif re.search('^\s*ob',s,re.IGNORECASE):
            return('oblong')
        elif re.search('^\s*o',s,re.IGNORECASE):
            return('oval')
        elif re.search('^\s*p',s,re.IGNORECASE):
            return('pyriform')
        elif re.search('^\s*r',s,re.IGNORECASE):
            return('round')
        elif re.search('^\s*s',s,re.IGNORECASE):
            return('spindle')
        else:
            return(None)
    except Exception as e:
        print(e)
        pass
    return(None)

def mixupShapes(df,chimera,segmenter,logger):
    shapes = df['berry shape'].dropna()
    chimera.compose(shapes)
    segment_df = segmenter.segment(chimera.composite)[0]
    umcc    = UMCCEncoder(logger=logger, contours=chimera.contours)
    tolerance = 0.5
    highest_tolerance = tolerance
    scc        = SCCEncoder(logger=logger,contours=chimera.contours, tolerance=tolerance)
    highest_tortuosity = scc.tortuosity()
    while( True ):
        tolerance = tolerance/2
        scc.encode(contours=chimera.contours, tolerance=tolerance)
        new_tortuosity = scc.tortuosity()
        if( scc.tortuosity() > 1.025 * highest_tortuosity ): #Give ourselves some buffer room on this
            highest_tolerance = tolerance
            highest_tortuosity = new_tortuosity
        else:
            break
    scc.encode(contours=chimera.contours, tolerance=highest_tolerance)
    df.insert(df.shape[1], 'chimera_LvsW', segment_df['LvsW'][0])
    df.insert(df.shape[1], 'chimera_eccentricity', segment_df['blobEccentricity'][0])
    df.insert(df.shape[1], 'chimera_solidity', segment_df['blobSolidity'][0])
    df.insert(df.shape[1], 'chimera_tortuosity', scc.tortuosity())
    df.insert(df.shape[1], 'chimera_umccX', str(umcc.x))
    df.insert(df.shape[1], 'chimera_umccLogX', math.log(umcc.x))
    df.insert(df.shape[1], 'chimera_umccY', str(umcc.y))
    df.insert(df.shape[1], 'chimera_umccLogY', math.log(umcc.y))
    return(df)


if __name__ == '__main__':
    logger      = logging.getLogger()
    parsed      = parse_args()
    decoded_level = eval("logging.{}".format(parsed.level))
    logger.setLevel(decoded_level)

    collated_df = pd.read_excel(parsed.data, header=1)
    #Preprocess berry shape categories, fixing shortened abbreviations and removing anything that is invalid (mark as None, or NA).
    collated_shapes = collated_df['berry shape']
    collated_df['berry shape'] = collated_shapes.apply(cleanupShapes)

    segmenter = BinaryThresholdSegment( channel=0, resize=1.0, minArea=parsed.mina, maxArea=parsed.maxa, numRefs=0 )

    chimera = Chimera()
    chimera.loadMap(parsed.mapp, parsed.map)
    collated_shapes_df = collated_df[['population','year','accession name','up-right no.','berry shape']]
    collated_shapes_condensed_df = collated_shapes_df.groupby(['population','year','accession name'],sort=False,dropna=True)
    collated_shapes_new_df = collated_shapes_condensed_df.apply(mixupShapes, chimera, segmenter, logger)
    collated_df = collated_df.merge(collated_shapes_new_df,on=['population','year','accession name','up-right no.'], how='left')
    collated_df.to_csv(parsed.odata, index=False)
