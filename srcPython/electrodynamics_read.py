#!/usr/bin/env python

import netCDF4 as nc
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from pylab import cm
from datetime import datetime
from datetime import timedelta
import sys

args = sys.argv

file = args[1]
ncfile = nc.Dataset(file, 'r')

lats = np.array(ncfile.variables['MagneticLatitude'])
mlts = np.array(ncfile.variables['MagneticLocalTime'])
potential = np.array(ncfile.variables['Potential'])/1000.0
eflux = np.array(ncfile.variables['EFlux'])
times_in_s = np.array(ncfile.variables['Time'])

hasFac = 'FAC' in ncfile.variables
if hasFac:
    fac = np.array(ncfile.variables['FAC'])

nTimes = len(times_in_s)
base = datetime(1965,1,1,0,0,0)
times = []
for t in times_in_s:
    times.append(base + timedelta(seconds = t))

theta, r = np.meshgrid(mlts * np.pi/12.0 - np.pi/2.0, 90.0 - lats)

nLevels = 21
maxp = np.max(np.abs(potential))
dr = np.ceil(maxp/(nLevels-1)*2)
minp = -dr * (nLevels-1)/2
maxp = np.abs(minp)
levelp = np.arange(minp,maxp+dr,dr)

maxi = np.max(np.abs(eflux))
mini = 0.0

norm = cm.colors.Normalize(vmax=mini, vmin=maxi)
if (mini >= 0):
    cmap = cm.plasma
else:
    cmap = cm.bwr

di = int(np.ceil(nTimes/9))

fig = plt.figure(figsize = (10,10))

iPlot = 1
ax = []
cax = []
for iT in np.arange(0,nTimes,di):

    subplot = 330 + iPlot
    iPlot += 1
    print(subplot, iT)
    ax.append(fig.add_subplot(subplot, projection = 'polar'))

    pot2d = potential[iT]
    eflux2d = eflux[iT]

    eflux2d[eflux2d < 0.1] = np.nan

    cax.append(ax[-1].pcolor(theta, r, eflux2d, \
                             vmin = mini, vmax = maxi, cmap = cmap))
    ax[-1].contour(theta, r, pot2d, levelp, colors = 'k')
    fig.colorbar(cax[-1])

fig.savefig('netcdf_test.png')
plt.close()
mlt_ticks = np.arange(0, 24, 3)
if hasFac:

    maxFac = np.max(np.abs(fac))
    print('maxFac = ', maxFac)
    normFac = cm.colors.Normalize(vmax=maxFac, vmin=-maxFac)

    fig = plt.figure(figsize = (10,10))

    iPlot = 1
    ax = []
    cax = []
    for iT in np.arange(0,nTimes,di):

        subplot = 330 + iPlot
        iPlot += 1
        ax.append(fig.add_subplot(subplot, projection = 'polar'))

        pot2d = potential[iT]
        fac2d = fac[iT]

        cax.append(ax[-1].pcolor(theta, r, fac2d, \
                                 vmin = -maxFac, vmax = maxFac, cmap = cm.bwr))
        # ax[-1].contour(theta, r, pot2d, levelp, colors = 'k')
        ax[-1].set_xticklabels([f"{val:.1f}" for val in mlt_ticks])
        fig.colorbar(cax[-1])

    fig.savefig('netcdf_test_fac.png')
    plt.close()
