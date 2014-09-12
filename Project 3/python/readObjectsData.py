#!/usr/bin/env python
"""
Created on Tir 22 Oct 2013

Script for reading and plotting data from celestial objects outputted from C++
simulations of the solar system.

@author Benedicte Emilie Braekken
"""
from matplotlib import pyplot as plt
import numpy as np
import os,sys

"""Constants"""

_OBJECTS_PATH = '../data/objects/'
_PLOT_EVERY = 100
_XLIM = 50
_YLIM = 50

"""Classes"""

class ObjectData(object):
    """
    Contains the data for one single object.
    """
    def __init__(self, datafile):
        """
        Name for object is found in datafile.

        @param datafile The file containing positional data.
        """
        self._file = datafile
        self.name = self._findName()
        self.currentLine = 0
        self.allDataRead = False

    def getNext(self):
        """
        @return The next x,y pair from datafile.
        """
        inData = open(self._file,'r')

        # Read til end of header
        while True:
            if inData.readline().startswith('#syntax'):
                break

        # Read til the next line and save x,y pair
        counter = 0
        while counter < self.currentLine:
            inData.readline()
            counter += 1

        # Save correct x,y pair
        x,y = np.asarray(inData.readline().strip().split(),float)

        # Close file
        inData.close()

        # Update what is next line
        self.currentLine += 1

        return x,y

    def readAllData(self):
        """
        Reads datafile and puts into arrays.
        """
        if self.allDataRead:
            print 'All data already read. Returning.'
            return

        inData = open(self._file, 'r')
        x = []
        y = []

        # First read till end of header
        for line in inData:
            if line.startswith('#syntax'):
                break

        for line in inData:
            columns = line.strip().split()
            x.append(columns[0])
            y.append(columns[1])

        self.x,self.y = np.asarray(x), np.asarray(y)
        inData.close()
        self.allDataRead = True

    def hasNextPoint(self):
        """
        Says whether there are more points left in datafile. Tries to fetch
        more data from file. If that fails because of something, returns False.
        Else returns True.
        """
        try:
            self.getNext()
        except:
            return False

        return True

    def _findName(self):
        """
        Finds the name of the object from the datafile.

        @param datafile The file containing data.
        """
        inData = open(self._file, 'r')

        for line in inData:
            if line.startswith('Positions'):
                inData.close()
                return line.strip().split(':')[1].strip()
        
        # Rapid err msg, no time dude
        print 'If you see this the datafile has wrong syntax.'
        sys.exit(1)

    def __len__(self):
        """
        Length means number of x,y points.
        """
        if not self.allDataRead:
            print "Cannot find length when all data isn't read. Returning."

        return len(self.x)

"""Methods"""

if __name__ == '__main__':
    """
    Running script will loop through all files in objects folder and add them
    to the plot along with name as label.
    """
    animate = False
    try:
        animate = eval(sys.argv[1])
    except:
        print 'Want animation?\nUsage: ./%s True/False' % sys.argv[0]
        sys.exit(1)
    if animate != False and animate != True:
        print 'Command line argument must be True or False.'
        sys.exit(1)

    files = os.listdir(_OBJECTS_PATH)

    plt.ion()
    fig = plt.figure()
    fig.suptitle('The solar system')
    ax = fig.add_subplot(111)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.grid('on')
    ax.axis('equal')
    ax.set_ylim([-_YLIM,_YLIM])
    ax.set_xlim([-_XLIM,_XLIM])

    if animate:
        datafiles = {}
        lines = {}
        for datafile in files:
            newDataObject = ObjectData(os.path.join(_OBJECTS_PATH,datafile))
            datafiles[newDataObject.name] = newDataObject
            x,y = newDataObject.getNext()
            lines[newDataObject.name] = ax.plot(x,y,'o',label=newDataObject.name)[0]

        plt.draw()
        ax.legend(loc='best')

        while datafiles[datafiles.keys()[0]].hasNextPoint():
            for name in lines.keys():
                x,y = datafiles[name].getNext()
                lines[name].set_data(x,y)

            plt.draw()

    else:
        for datafile in files:
            data = ObjectData(os.path.join(_OBJECTS_PATH,datafile))
            ax.plot(data.x,data.y,label=data.name)

    plt.ioff()
    plt.show()
