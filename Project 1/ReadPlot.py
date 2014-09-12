#!/usr/bin/env python
"""
Created on Fri 30 Aug 2013

Script for reading and plotting the data produced by the ThomasAlgorithm.

author: Benedicte Emilie Braekken
"""
from matplotlib import pyplot as plt
import sys

try:
    global _infile
    _infile = sys.argv[1];
except IndexError:
    print 'Please specify datafile to read.'
    print 'Usage: python %s datafile.dat' % sys.argv[0]
    sys.exit(1)

def read(infile=False):
    '''
    Reads the data from given datafile.
    '''
    if infile: infile = open(infile,'r')
    else: infile = open(_infile,'r')

    global _x; _x = []
    global _y; _y = []

    for line in infile:
        columns = line.split()
        _x.append(float(columns[0]))
        _y.append(float(columns[1]))

    infile.close()

    return _x,_y

def plot():
    '''
    Plots data.
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(_x,_y)

if __name__ == '__main__':
    print 'Sucessfully started.'
    read()
    print 'File read.'
    plot()

    plt.show()
