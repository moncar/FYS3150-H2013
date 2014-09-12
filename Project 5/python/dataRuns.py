#!/usr/bin/env python
"""
Created on Sun 8 Dec 2013

Script is used to run the correct sequences for producing data for energy
simulations.

@author Benedicte Emilie Braekken
"""
import os

"""CONSTANTS"""
_EXECUTABLE = 'build/ClusterModel.x'
_METHODS = ['lf', 'rk4']
_DATADIR = 'data'
_ULTIMATE_DT = 0.01

"""FUNCTIONS"""
def run(dt, T, eps, N, method):
    """
    Run the cluster simulation with the correct parameters.
    """
    if not method in _METHODS:
      print 'Not a method.'
      return

    # Get optimal save each
    counter = dt
    SE = 1
    while counter < _ULTIMATE_DT:
      counter *= 10
      SE *= 10

    filename = '%s_N%g_dt%g_T%g_epsilon%g_saveEach%g.dat' % \
                          (method,N,dt,T,eps,SE)

    runstring = './' + _EXECUTABLE + ' '
    runstring += '--dt %g ' % dt
    runstring += '--T %g ' % T
    runstring += '--epsilon %g ' % eps
    runstring += '--N %g ' % N
    runstring += '--saveEach %g ' % SE
    runstring += '--savefile %s/%s' % (_DATADIR,filename)

    os.system('clear')

    print '$', runstring
    os.system(runstring)

if __name__ == '__main__':
    # Check that datadir is present, if not create
    if not os.path.exists(_DATADIR):
        os.mkdir(_DATADIR)

    Ns = [100, 300, 500, 1000]
    dt = 0.001
    T = 5
    epsilon = 0.15

    for method in _METHODS:
        for N in Ns:
            run(dt,T,epsilon,N,method)

    print 'Ran all simulations.'
