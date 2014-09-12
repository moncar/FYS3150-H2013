#!/usr/bin/env python
"""
Created on Mon 9 Dec 2013 01:06

Created energy plots and saves.

@author Benedicte Emilie Braekken
"""
import os
from dataReader import ClusterReader,ClusterVisualizer
import matplotlib.pyplot as plt

"""CONSTANTS"""
_DATADIR = 'data'
_OUTDIR = 'plots'

if __name__ == '__main__':
    if not os.path.exists(_OUTDIR):
        os.mkdir(_OUTDIR)

    for datafile in os.listdir(_DATADIR + '/'):
        data = ClusterReader(_DATADIR + '/' + datafile)
        viz = ClusterVisualizer(data)

        viz.plotBoundObjects()
        filename = _OUTDIR + '/' + 'bound_'
        for part in datafile.split('.')[0:-1]:
            filename += part
        filename += '.eps'
        plt.savefig(filename, format='eps')

        viz.plotTotalAndBoundEng()
        filename = _OUTDIR + '/' + 'boundAndTotEng_'
        for part in datafile.split('.')[0:-1]:
            filename += part
        filename += '.eps'
        plt.savefig(filename, format='eps')

    print 'Made all plots.'
