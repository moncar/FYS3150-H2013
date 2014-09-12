#!/usr/bin/env python
"""
Created on Wed 4 September 2013

FYS3150: Project 1
Contains the closed form of the differential equation for checking relative
error and plotting function.

@author: Benedicte Emilie Braekken
"""
import numpy as np
import ReadPlot as rp
from matplotlib import pyplot as plt

def cf(x):
    """
    Solution of general one-dimension Poisson equation.
    """
    return 1 - (1 - np.exp(-10))*x - np.exp(-10*x)

def source_term(x):
    """
    Source term used in solution of differential equation. One-dimension
    Poisson equation.
    """
    return 100*np.exp(-10*x)

def plot_source_term(a=0,b=1,n=1e3):
    """
    Plots the source term.
    """
    x = np.linspace(a,b,n)

    fig = plt.figure()
    fig.suptitle('The source term')
    ax = fig.add_subplot(111)
    ax.plot(x,source_term(x))
    ax.set_xlabel('$x$')
    ax.set_ylabel('$f(x)$')

def plot_closed_form(a=0,b=1,n=1e3):
    """
    Plots the closed form.
    """
    x = np.linspace(a,b,n)

    fig = plt.figure()
    fig.suptitle('Closed form')
    ax = fig.add_subplot(111)
    ax.plot(x,cf(x))
    ax.set_xlabel('$x$')
    ax.set_ylabel('$f(x)$')

def plot_data_against_cf(datafile,a=0,b=1,n=1e3):
    """
    Plots the given datafile against the closed form of the solution.
    """
    x = np.linspace(a,b,n)
    d_x,d_y = rp.read(datafile)

    fig = plt.figure()
    fig.suptitle('Analytical vs. Numerical')
    ax = fig.add_subplot(111)
    ax.grid('on')
    ax.plot(x,cf(x),label='Analytical')
    ax.plot(d_x,d_y,label='Numerical')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$f(x)$')
    ax.legend(loc='best')

if __name__ == '__main__':
    import sys
    plot_data_against_cf(sys.argv[1])
    plt.show()
