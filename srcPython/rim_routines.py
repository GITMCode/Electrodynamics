#!/usr/bin/env python

import re
import datetime as dt
from numpy import zeros, reshape
import numpy as np

#------------------------------------------------------------------------------
# SWMF / RIM files have 0 / 360 as noon and the 0th element in the array.
# To fix this, we need to:
#  - take off the 360
#  - shift the data so that midnight is at 0th element
#  - copy the 0th element on the end to give the wrap around 
#------------------------------------------------------------------------------

def fix_mlts(data):
    nMlts = len(data[0,:])
    nLats = len(data[:,0])
    nShift = int((nMlts-1)/2)
    newData = np.zeros((nLats, nMlts))
    newData[:, 0:nShift] = data[:, nShift:-1]
    newData[:, nShift:-1] = data[:, 0:nShift]
    newData[:, -1] = newData[:, 0]
    return newData


#------------------------------------------------------------------------------
# Read a single SWMF / RIM file
#------------------------------------------------------------------------------

def read_rim_ascii_file(filename):

    '''
    Read an ascii ".idl" output file and load the contents into the object.
    '''

    import gzip

    # slurp entire file.
    print(' --> Reading File : ', filename)

    if (filename[-3:] == '.gz'):
        infile = gzip.open(filename, 'rt')
    else:
        infile = open(filename, 'r')
    raw = infile.readlines()
    infile.close()

    dataOut = {}
    
    # Parse header
    title = raw[raw.index('TITLE\n')+1]
    dataOut['title'] = title[title.index('"')+1:title.rindex('"')]

    i = raw.index('NUMERICAL VALUES\n')
    dataOut['nvars'] = int(raw[i+1].split()[0])
    dataOut['ntheta']= int(raw[i+2].split()[0])
    dataOut['nphi']  = int(raw[i+3].split()[0])

    # Convenience:
    nphi, ntheta = dataOut['nphi'], dataOut['ntheta']

    i = raw.index('TIME\n')
    dataOut['time'] = dt.datetime(
        int(raw[i+1].split()[0]),      #year
        int(raw[i+2].split()[0]),      #month
        int(raw[i+3].split()[0]),      #day
        int(raw[i+4].split()[0]),      #hour
        int(raw[i+5].split()[0]),      #min
        int(raw[i+6].split()[0]),      #sec
        int(raw[i+7].split()[0])*1000  #microsec
    )

    i = raw.index('SIMULATION\n')
    dataOut['iter']    =   int(raw[i+1].split()[0])
    dataOut['simtime'] = float(raw[i+2].split()[0])

    i = raw.index('DIPOLE TILT\n')
    dataOut['tilt'] = zeros(2)
    dataOut['tilt'] = float(raw[i+1].split()[0])
    dataOut['tilt'] = float(raw[i+2].split()[0])

    i = raw.index('VARIABLE LIST\n')
    namevar = []
    units   = {}
    for j in range(i + 1, i + dataOut['nvars'] + 1):
        match = re.match(r'\s*\d+\s+([\w\s\W]+)\[([\w\s\W]+)\]',raw[j])
        if match:
            name = (match.group(1).strip()).lower()
            namevar.append(name)
            units[name] = match.group(2).strip()
        else:
            raise ValueError('Could not parse %s' % raw[j])

    ### Read all data ###

    # Create data arrays
    nPts = dataOut['ntheta'] * dataOut['nphi']
    for key in namevar:
        dataOut['n_'+key] = zeros(nPts)
        dataOut['s_'+key] = zeros(nPts)
    i = raw.index('BEGIN NORTHERN HEMISPHERE\n')+1

    # Some compilers insert line breaks automatically when fortran format
    # string is not adequately specified.  Let's see if that's the
    # case here: how many lines does it take to cover all variables?
    nvars, nvarline, nwrap = len(namevar), 0, 0
    while nvarline<nvars:
        nvarline += len(raw[i+nwrap].split())
        nwrap += 1

    # Fill data arrays:
    for j in range(nPts):
        # Build list of variables; accounting for line-wrapping:
        parts = []
        iLine = i + j*nwrap
        for iwrap in range(nwrap):
            parts += raw[iLine+iwrap].split()
        for k in range(dataOut['nvars']):
            dataOut['n_'+namevar[k]][j] = parts[k]
    i = raw.index('BEGIN SOUTHERN HEMISPHERE\n')+1
    for j in range(nPts):
        parts = []
        iLine = i + j*nwrap
        for iwrap in range(nwrap):
            parts += raw[iLine+iwrap].split()
        for k in range(dataOut['nvars']):
            dataOut['s_'+namevar[k]][j] = parts[k]

    # Create 2-D arrays.
    for key in namevar:
        nkey, skey = 'n_'+key, 's_'+key
        dataNorth = reshape(dataOut[nkey], (ntheta, nphi), order='F')
        dataSouth = reshape(dataOut[skey], (ntheta, nphi), order='F')
        dataOut[nkey] = fix_mlts(dataNorth)
        dataOut[skey] = fix_mlts(dataSouth)

    # Some extra grid info:
    dataOut['dlon'] = dataOut['n_psi'  ][0,3] - dataOut['n_psi'  ][0,2]
    dataOut['dlat'] = dataOut['n_theta'][3,0] - dataOut['n_theta'][2,0]
    
    return dataOut


#------------------------------------------------------------------------------
# Read and extract needed information from SWMF / RIM files
#------------------------------------------------------------------------------

def read_all_rim_files(filelist):

    # Let's get some information before starting this:
    nTimes = len(filelist)
    dataOneFile = read_rim_ascii_file(filelist[0])

    nLats = dataOneFile['ntheta']
    nMlts = dataOneFile['nphi']

    dataToWriteN = {}
    dataToWriteS = {}

    dataToWriteN["nLats"] = nLats
    dataToWriteN["nMlts"] = nMlts
    dataToWriteN["nTimes"] = nTimes

    dataToWriteS["nLats"] = nLats
    dataToWriteS["nMlts"] = nMlts
    dataToWriteS["nTimes"] = nTimes

    LatsNorth = 90.0 - dataOneFile['n_theta']
    LatsSouth = np.flip(dataOneFile['s_theta'] - 90.0, 0)

    MltsNorth = (dataOneFile['n_psi'] + 180.0) % 360.0
    MltsSouth = (np.flip(dataOneFile['s_psi'], 0) + 180.0) % 360.0
    MltsNorth[:, -1] = MltsNorth[:, -1] + 360.0
    MltsSouth[:, -1] = MltsSouth[:, -1] + 360.0

    dataToWriteN["lats"] = LatsNorth[:, 0]
    dataToWriteN["mlts"] = MltsNorth[0, :] / 15.0

    dataToWriteS["lats"] = LatsSouth[:, 0]
    dataToWriteS["mlts"] = MltsSouth[0, :] / 15.0

    dLat = dataToWriteN["lats"][1] - dataToWriteN["lats"][2]
    dLon = (dataToWriteN["mlts"][2] - dataToWriteN["mlts"][1]) * 15.0

    dtom = 2 * np.pi * (6371.0 + 120.0) * 1000.0 / 360.0
    area = dtom * dtom * dLat * dLon * np.cos(LatsNorth * np.pi / 180.0)

    # These are the RIM / SWMF variable names:
    rimVars = ['phi', 'e-flux', 'ave-e']

    # These are the AMIE variable names to write that IE understands:
    dataToWriteN["Vars"] = ['Potential (V)', \
                            'Electron Energy Flux (ergs/cm2/s)', \
                            'Electron Mean Energy (keV)']
    dataToWriteN["nVars"] = len(dataToWriteN["Vars"])
    dataToWriteN["version"] = 0.9
    dataToWriteN["imf"] = []
    dataToWriteN["ae"] = []
    dataToWriteN["dst"] = []
    dataToWriteN["hp"] = []
    dataToWriteN["cpcp"] = []
    dataToWriteN["hem"] = 'N'
    dataToWriteN["times"] = []
    for var in dataToWriteN["Vars"]:
        dataToWriteN[var] = []

    dataToWriteS["Vars"] = ['Potential (V)', \
                            'Electron Energy Flux (ergs/cm2/s)', \
                            'Electron Mean Energy (keV)']
    iPot_ = 0
    iEflux_ = 1
    dataToWriteS["nVars"] = len(dataToWriteS["Vars"])
    dataToWriteS["version"] = 0.9
    dataToWriteS["imf"] = []
    dataToWriteS["ae"] = []
    dataToWriteS["dst"] = []
    dataToWriteS["hp"] = []
    dataToWriteS["cpcp"] = []
    dataToWriteS["hem"] = 'S'
    dataToWriteS["times"] = []
    for var in dataToWriteS["Vars"]:
        dataToWriteS[var] = []

    cPot = dataToWriteN["Vars"][iPot_]
    cEflux = dataToWriteN["Vars"][iEflux_]

    for file in filelist:
        dataOneFile = read_rim_ascii_file(file)

        # Move variables over to export dictionary:
        for iVar, rimVar in enumerate(rimVars):
            amieVar = dataToWriteN["Vars"][iVar]
            dataToWriteN[amieVar].append(dataOneFile['n_' + rimVar])
            dataToWriteS[amieVar].append(np.flip(dataOneFile['s_' + rimVar], 0))
        dataToWriteN['times'].append(dataOneFile['time'])
        dataToWriteS['times'].append(dataOneFile['time'])
        
        cpcpN = np.max(dataToWriteN[cPot][-1]) - np.min(dataToWriteN[cPot][-1])
        cpcpS = np.max(dataToWriteS[cPot][-1]) - np.min(dataToWriteS[cPot][-1])
        hpN = np.sum(area * dataToWriteN[cEflux][-1]) / 1e9
        hpS = np.sum(area * dataToWriteS[cEflux][-1]) / 1e9
        
        dataToWriteN["imf"].append([-1e32, -1e32, -1e32, -1e32])
        dataToWriteN["ae"].append([-1e32, -1e32, -1e32, 0.0])
        dataToWriteN["dst"].append([0.0, 0.0])
        dataToWriteN["hp"].append([hpN, 0.0])
        dataToWriteN["cpcp"].append(cpcpN)

        dataToWriteS["imf"].append([-1e32, -1e32, -1e32, -1e32])
        dataToWriteS["ae"].append([-1e32, -1e32, -1e32, 0.0])
        dataToWriteS["dst"].append([0.0, 0.0])
        dataToWriteS["hp"].append([hpS, 0.0])
        dataToWriteS["cpcp"].append(cpcpS)


        # Convert from W/m2 to mW/m2 or ergs/cm2/s:
        dataToWriteN[cEflux][-1] = dataToWriteN[cEflux][-1] * 1000.0
        dataToWriteS[cEflux][-1] = dataToWriteS[cEflux][-1] * 1000.0
        # Convert from kV to V:
        dataToWriteN[cPot][-1] = dataToWriteN[cPot][-1] * 1000.0
        dataToWriteS[cPot][-1] = dataToWriteS[cPot][-1] * 1000.0
        
    dataToWriteN['nTimes'] = len(dataToWriteN['times'])
    dataToWriteS['nTimes'] = len(dataToWriteS['times'])
        
    return dataToWriteN, dataToWriteS


