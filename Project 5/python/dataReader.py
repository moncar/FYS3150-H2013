#!/usr/bin/env python
"""
Created on Tir 26 Nov 2013

Script is used to read data produced by the open cluster simulation. Datafile
has the following syntax:

t=1 E=<total energy> (x y z bound(True/False))(x y z bound(True/False)) ...
t=2 E=<total energy> (x y z bound(True/False))(x y z bound(True/False)) ...

Where each line represents one state of the open cluster and each object has
three spatial coordinates. This can later be extended to incorporate all
possible dimensions.

@author Benedicte Emilie Braekken
"""
import re,sys,argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

"""CONSTANTS"""

_UNIVERSAL_NUMBER = r'[+\-]?\d+\.?\d*[Ee]?[+\-]?\d*'
_FONTSIZE = 26
_ENG_UNIT = '$\mathrm{ly}^2 \mathrm{M}_\odot \\tau_{\mathrm{crunch}}^{-2}$'
_TITLE = '$N=%g, \ R_0=%g, \ \mathrm{d}t=%g, \ \epsilon=%g$'

"""CLASSES"""

class ClusterReader(object):
    """
    Automatic fields:
    
    @field t 1D Array of time steps
    @field positions 3D array of positions, (#timesteps, #particles,
                     #dimensions)
    @field n Integer total number of time steps.
    @field N Integer total number of particles.
    @field dimensions Integer describing number of dimensions.
    @field E Array of total energies.
    @field R0 The radius of the cluster.
    @field epsilon The smoothing factor, if not included is false.
    @field dt Float describing time step.
    @field SE Indicating frequency of writing to file.
    @field T Length of simulation (time).
    @field kinetic Array with all kinetic energies per particle.
    @field potential Array with all potential energies per particle.
    """
    def __init__(self, datafile):
        """
        @param datafile The file containing data from cluster simulation.
        """
        self.datafile = datafile
        self._readFile()

    def _readFile(self):
        """
        Reads the file. Data is stored in instance.
        """
        datafile = self.datafile

        # Fetch universal number regex for shorter code
        un = _UNIVERSAL_NUMBER

        # Expression for finding times
        regexp_t = r't=(' + un + ')'

        # Expression for finding energies
        regexp_E = r'E=(' + un + ')'

        # Expression for finding bound / not
        regexp_bound_t = r'bound=(False|True)'

        # Expression for finding coordinates
        regexp_pos = r'x=(' + un + ')\s+y=(' + un + ')\s+z=(' + un + ')\s+'
        regexp_x = r'x=(' + un + ')'
        regexp_y = r'y=(' + un + ')'
        regexp_z = r'z=(' + un + ')'

        # Expression for finding kinetic and potential energy
        regexp_kin = r'kinetic=(' + un + ')'
        regexp_pot = r'potential=(' + un + ')'
        regexp_pot_b = r'potential_bound=(' + un + ')'

        # Read whole file
        data = open(datafile,'r')
        contents = data.read()
        data.close()

        # Fetch general values
        self.t = np.asarray(re.findall(regexp_t, contents), dtype=np.float64)
        self.T = self.t[-1]
        self.E = np.asarray(re.findall(regexp_E, contents), dtype=np.float64)
        self.R0 = float(re.search(r'R0=(' + un + ')', contents).groups(0)[0])
        self.epsilon = float(re.search(r'epsilon=(' + un + ')', contents).groups(0)[0])
        self.N = float(re.search(r'N=(' + un + ')', contents).groups(0)[0])
        self.SE = float(re.search(r'saveEach=(' + un + ')', contents).groups(0)[0])
        self.dt = ( self.t[1]-self.t[0] ) / self.SE

        # Read all positions
        data = open(datafile,'r')

        """
        Need empty lists for creating (#timesteps, #particles, #dimensions)
        sized array

        And other relevant data
        """
        kinetic = []
        potential = []
        potentialBound = []
        positions = []
        bound = []

        # Some values need temporary storage
        newline = []

        # Switch for knowing when past meta-part
        pastMeta = False

        # A line is the same as a time step
        for line in data:
            # Do not record unless past meta header
            if not pastMeta:
                if line.strip() == '********/META********':
                    pastMeta = True
            else:
                # Fetch lists of position for each particle for this timestep
                X = np.asarray(re.findall(regexp_x, line), dtype=np.float64)
                Y = np.asarray(re.findall(regexp_y, line), dtype=np.float64)
                Z = np.asarray(re.findall(regexp_z, line), dtype=np.float64)

                # Same with kinetic and potential energy
                K = np.asarray(re.findall(regexp_kin, line), dtype=np.float64)
                P = np.asarray(re.findall(regexp_pot, line), dtype=np.float64)
                PB = np.asarray(re.findall(regexp_pot_b, line), dtype=np.float64)

                # Bound particles
                B = np.asarray(re.findall(regexp_bound_t, line))

                # Each timestep is a new line in the matrices
                for x,y,z in zip(X,Y,Z):
                    newline.append((x,y,z))

                # Add lines to correct matrices
                positions.append(newline)
                kinetic.append(K)
                potential.append(P)
                potentialBound.append(PB)
                bound.append(B)

                # Ready for next line
                newline = []

        data.close()

        # Save as fields
        self.positions = np.asarray(positions, dtype=np.float64)
        self.n,self.N,self.dimensions = self.positions.shape
        self.kinetic,self.potential = np.asarray(kinetic),np.asarray(potential)
        self.PBound = potentialBound

        """
        Ugly hack for converting to boolean matrix. Does not change values but
        casts from string to boolean. Had problems with built-in numpy way ---
        it converted everything to true.
        """
        self.bound = np.asarray(bound) == 'True'

        KBound = []
        # Loop through to find kinetic energy only from bound objects
        for K,B in zip(self.kinetic,self.bound):
            kin = 0

            for k,b in zip(K,B):
                if b: kin += k

            KBound.append(kin)

        self.KBound = np.asarray(KBound,dtype=np.float64)
        self.EBound = np.sum(self.PBound,axis=1) + self.KBound

class ClusterVisualizer(object):
    """
    Can use a ClusterReader to visualize properties of the system.
    """
    def __init__(self, reader):
        """
        @param reader The reader to animate.
        """
        self.reader = reader

    def animateCluster(self):
        """
        Animation is done using matplotlib and interactive mode by redrawing.
        """
        plt.ion()

        fig = plt.figure()
        ax3 = fig.add_subplot(111, projection='3d')

        ax3.set_xlabel('$x$ [l.y.]')
        ax3.set_ylabel('$y$ [l.y.]')
        ax3.set_zlabel('$z$ [l.y.]')

        dots = []

        for j in range(self.reader.N):
            x,y,z = self.reader.positions[0,j]
            # Last index [0] because plot returns tuple with one element
            dots.append(ax3.plot([x],[y],[z],'o')[0])

        ax3.set_xlim(-1.5*self.reader.R0, 1.5*self.reader.R0)
        ax3.set_ylim(-1.5*self.reader.R0, 1.5*self.reader.R0)
        ax3.set_zlim(-1.5*self.reader.R0, 1.5*self.reader.R0)

        plt.draw()

        for i in range(1,self.reader.n):
            ax3.set_title('$t = %g$' % self.reader.t[i])
            for j in range(self.reader.N):
                x,y,z = self.reader.positions[i,j]

                """
                Necessity two have two lines because there is no set_data for
                3 dim data in matplotlib.

                http://matplotlib.org/examples/animation/simple_3danim.html
                """
                dots[j].set_data([x],[y])
                dots[j].set_3d_properties([z])

            plt.draw()

        plt.ioff()

    def plotTotalEnergy(self, bound=False):
        """
        Plots the total energy against time.

        @param bound Only plot for the bound particles.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.grid('on')

        ax.set_xlabel('$t$ [$\\tau_{\mathrm{crunch}}$]', fontsize=_FONTSIZE)
        ax.set_ylabel('$E$ [%s]' % _ENG_UNIT, fontsize=_FONTSIZE)

        if not bound:
            ax.plot(self.reader.t,self.reader.E)
        else:
            ax.plot(self.reader.t,self.reader.EBound)

        dt,N,eps,R0 = \
                self.reader.dt,self.reader.N,self.reader.epsilon,self.reader.R0
        ax.set_title(_TITLE % (N,R0,dt,eps),fontsize=_FONTSIZE)
        fig.tight_layout()

    def plotKineticEnergy(self, bound=False):
        """
        Plots the kinetic energy against time.

        @param bound Only plot for the bound particles.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.grid('on')

        ax.set_xlabel('$t$ [$\\tau_{\mathrm{crunch}}$]', fontsize=_FONTSIZE)
        ax.set_ylabel('$E$ [%s]' % _ENG_UNIT, fontsize=_FONTSIZE)

        if not bound:
            ax.plot(self.reader.t,np.sum(self.reader.kinetic,axis=1))
        else:
            ax.plot(self.reader.t,self.reader.KBound)

        dt,N,eps,R0 = \
                self.reader.dt,self.reader.N,self.reader.epsilon,self.reader.R0
        ax.set_title(_TITLE % (N,R0,dt,eps),fontsize=_FONTSIZE)
        fig.tight_layout()

    def plotPotentialEnergy(self, bound=False):
        """
        Plots the potential energy against time.

        @param bound Only plot for the bound particles.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.grid('on')

        ax.set_xlabel('$t$ [$\\tau_{\mathrm{crunch}}$]', fontsize=_FONTSIZE)
        ax.set_ylabel('$E$ [%s]' % _ENG_UNIT, fontsize=_FONTSIZE)

        if not bound:
            ax.plot(self.reader.t,np.sum(self.reader.potential,axis=1))
        else:
            ax.plot(self.reader.t,np.sum(self.reader.PBound,axis=1))

        dt,N,eps,R0 = \
                self.reader.dt,self.reader.N,self.reader.epsilon,self.reader.R0
        ax.set_title(_TITLE % (N,R0,dt,eps),fontsize=_FONTSIZE)
        fig.tight_layout()

    def plotAllEnergies(self, bound=False):
        """
        Plots all energies in same plot.

        @param bound Only plot for the bound particles.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.grid('on')

        ax.set_xlabel('$t$ [$\\tau_{\mathrm{crunch}}$]', fontsize=_FONTSIZE)
        ax.set_ylabel('$E$ [%s]' % _ENG_UNIT, fontsize=_FONTSIZE)

        if not bound:
            ax.plot(self.reader.t,self.reader.E, label='Total')
            ax.plot(self.reader.t,np.sum(self.reader.kinetic,axis=1), label='Kinetic')
            ax.plot(self.reader.t,np.sum(self.reader.potential,axis=1), label='Potential')
        else:
            ax.plot(self.reader.t,self.reader.EBound, label='Total')
            ax.plot(self.reader.t,self.reader.KBound, label='Kinetic')
            ax.plot(self.reader.t,np.sum(self.reader.PBound,axis=1), label='Potential')

        ax.legend(loc='best')

        dt,N,eps,R0 = \
                self.reader.dt,self.reader.N,self.reader.epsilon,self.reader.R0
        ax.set_title(_TITLE % (N,R0,dt,eps),fontsize=_FONTSIZE)
        fig.tight_layout()

    def plotBoundObjects(self):
        """
        Plots number of bounds objects as time goes.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.grid('on')

        ax.set_xlabel('$t$ [$\\tau_{\mathrm{crunch}}$]', fontsize=_FONTSIZE)
        ax.set_ylabel('# Particles', fontsize=_FONTSIZE)

        ax.plot(self.reader.t[1:], np.sum(self.reader.bound,axis=1)[1:])

        dt,N,eps,R0 = \
                self.reader.dt,self.reader.N,self.reader.epsilon,self.reader.R0
        ax.set_title(_TITLE % (N,R0,dt,eps),fontsize=_FONTSIZE)
        fig.tight_layout()

    def plotTotalAndBoundEng(self):
        """
        Plots the total energy for the whole system vs. total energy for bound
        particles.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.grid('on')

        ax.set_xlabel('$t$ [$\\tau_{\mathrm{crunch}}$]',fontsize=_FONTSIZE)
        ax.set_ylabel('$E$ [%s]' % _ENG_UNIT,fontsize=_FONTSIZE)

        ax.plot(self.reader.t,self.reader.E,label='Total')
        ax.plot(self.reader.t,self.reader.EBound,label='Bound')

        ax.legend(loc='best')

        dt,N,eps,R0 = \
                self.reader.dt,self.reader.N,self.reader.epsilon,self.reader.R0
        ax.set_title(_TITLE % (N,R0,dt,eps),fontsize=_FONTSIZE)
        fig.tight_layout()

    def show(self):
        """
        Shows all ordered plots.
        """
        plt.show()

if __name__ == '__main__':
    """
    Commandline arguments:
    * -a / --animate    Animates given file
    * -E / --energy     Plots total energy.
    * -s / --savefile   The savefile
    * --bound           Plots bound objects against time.
    * -K / --kinetic    Plots kinetic energy against time.
    * -P / --potential  Plots potential energy against time.
    * -b                Signal that for energy plots only bound
                        particles are to be considered.
    * --allE            Plots all energies against time.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-a','--animate',action='store_true',default=False)
    parser.add_argument('-E','--energy',action='store_true',default=False)
    parser.add_argument('-s','--savefile',type=str,default=False)
    parser.add_argument('--bound',action='store_true',default=False)
    parser.add_argument('-b',action='store_true',default=False)
    parser.add_argument('-K','--kinetic',action='store_true',default=False)
    parser.add_argument('-P','--potential',action='store_true',default=False)
    parser.add_argument('--totE',action='store_true',default=False)
    parser.add_argument('--allE',action='store_true',default=False)
    args = parser.parse_args()

    if not args.savefile:
        print 'File with data must be given.'
        sys.exit(1)

    my_reader = ClusterReader(args.savefile)
    viz = ClusterVisualizer(my_reader)

    if args.animate:
        viz.animateCluster()
    if args.energy:
        viz.plotTotalEnergy(args.b)
    if args.kinetic:
        viz.plotKineticEnergy(args.b)
    if args.potential:
        viz.plotPotentialEnergy(args.b)
    if args.allE:
        viz.plotAllEnergies(args.b)
    if args.bound:
        viz.plotBoundObjects()
    if args.totE:
        viz.plotTotalAndBoundEng()

    viz.show()
