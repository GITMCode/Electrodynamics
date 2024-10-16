#!/usr/bin/env python3

import numpy as np
from datetime import datetime
from amie_routines import *
from rim_routines import *
import argparse
from glob import glob

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(description = 'Convert SWMF/RIM to AMIE')

    parser.add_argument('-dir', default = './', \
                        help = 'directory of files to process')

    args = parser.parse_args()

    return args

#------------------------------------------------------------------------------
# main code:
#------------------------------------------------------------------------------

# Get the input arguments
args = parse_args()

filelist = sorted(glob(args.dir + '/it*.idl'))

if (len(filelist) <= 1):
    print(' --> Checking for gz files:')
    filelist = sorted(glob(args.dir + '/it*.gz'))

dataToWriteN, dataToWriteS = read_all_rim_files(filelist)

ymd = dataToWriteN['times'][0].strftime('%Y%m%d')
outfile = 'swmf' + ymd + dataToWriteN['hem'] + '.bin'
print(' -> Writing AMIE-style file : ', outfile)
amie_write_binary(outfile, dataToWriteN)

outfile = 'swmf' + ymd + dataToWriteS['hem'] + '.bin'
print(' -> Writing AMIE-style file : ', outfile)
amie_write_binary(outfile, dataToWriteS)
