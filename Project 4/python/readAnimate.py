#!/usr/bin/env python
"""
"""
from matplotlib import pyplot as plt
from closedForm import cfs
import numpy as np
import sys,re,os

"""Classes"""
class DiffusionProblem(object):
    """
    """
    def __init__(self, datafile):
        self.path = datafile

    def animate(self):
        inData = open(self.path,'r')

        plt.ion()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.grid('on')

        pastHeader = False
        firstFrame = True
        
        for line in inData:
            if not pastHeader:
                if line.startswith('#/HEADER#'):
                    pastHeader = True
                if line.startswith('x:'):
                    x = np.asarray([float(e) for e in line.strip().split()[1:]])
            else:
                u = np.asarray([float(e) for e in line.strip().split()])

                if firstFrame:
                    graph, = ax.plot(x, u)
                    firstFrame = False
                else:
                    graph.set_ydata(u)

                plt.draw()

        plt.ioff()
        plt.show()

        inData.close()

    def getState(self, t):
        """
        Returns state at time t.

        @param t Time of state.
        """
        tCur = 0

        dtfile = open(self.path, 'r')
        dt = float(re.search(r'dt: (.+)', dtfile.read()).groups(0)[0])
        dtfile.close()

        inData = open(self.path, 'r')

        pastHeader = False
        for line in inData:
            if tCur >= t:
                break
            if not pastHeader:
                if line.startswith('#/HEADER#'):
                    pastHeader = True
                if line.startswith('x:'):
                    x = np.asarray([float(e) for e in line.strip().split()[1:]])
            else:
                u = np.asarray([float(e) for e in line.strip().split()])
                tCur += dt

        inData.close()

        return x,u

def compareCases(t1=17*5e-3, t2=90*5e-3):
    """
    Makes plot for comparing the two test cases.
    """
    DATAPATH = '../data/'

    dt_1 = 5e-3
    dt_2 = 5e-5

    dx_1 = 0.1
    dx_2 = 0.01

    pFE_case1 = DiffusionProblem(os.path.join(DATAPATH,'fe_case1.dat'))
    pFE_case2 = DiffusionProblem(os.path.join(DATAPATH,'fe_case2.dat'))

    pBE_case1 = DiffusionProblem(os.path.join(DATAPATH,'be_case1.dat'))
    pBE_case2 = DiffusionProblem(os.path.join(DATAPATH,'be_case2.dat'))

    cn1 = open(os.path.join(DATAPATH,'cn_case1_t1.dat'))
    cn_1_t1 = np.asarray([float(e) for e in cn1.read().strip().split()])
    cn2 = open(os.path.join(DATAPATH,'cn_case1_t2.dat'))
    cn_1_t2 = np.asarray([float(e) for e in cn2.read().strip().split()])
    cn3 = open(os.path.join(DATAPATH,'cn_case2_t1.dat'))
    cn_2_t1 = np.asarray([float(e) for e in cn3.read().strip().split()])
    cn4 = open(os.path.join(DATAPATH,'cn_case2_t2.dat'))
    cn_2_t2 = np.asarray([float(e) for e in cn4.read().strip().split()])

    cn1.close()
    cn2.close()
    cn3.close()
    cn4.close()

    fig_case1_t1 = plt.figure()
    fig_case1_t2 = plt.figure()
    fig_case2_t1 = plt.figure()
    fig_case2_t2 = plt.figure()

    fig_case1_t1.suptitle('Comparing schemes, $ t = t_1 \\approx %g$, $\Delta t = %g$, $\Delta x = %g$' % (t1, dt_1, dx_1))
    fig_case1_t2.suptitle('Comparing schemes, $ t = t_2 \\approx %g$, $\Delta t = %g$, $\Delta x = %g$' % (t2, dt_1, dx_1))
    fig_case2_t1.suptitle('Comparing schemes, $ t = t_1 \\approx %g$, $\Delta t = %g$, $\Delta x = %g$' % (t1, dt_2, dx_2))
    fig_case2_t2.suptitle('Comparing schemes, $ t = t_2 \\approx %g$, $\Delta t = %g$, $\Delta x = %g$' % (t2, dt_2, dx_2))

    ax_case1_t1 = fig_case1_t1.add_subplot(111)
    ax_case1_t2 = fig_case1_t2.add_subplot(111)
    ax_case2_t1 = fig_case2_t1.add_subplot(111)
    ax_case2_t2 = fig_case2_t2.add_subplot(111)

    ax_case1_t1.set_xlabel('$x$ [dimensionless]')
    ax_case1_t1.set_ylabel('$y$ [dimensionless]')
    ax_case1_t2.set_xlabel('$x$ [dimensionless]')
    ax_case1_t2.set_ylabel('$y$ [dimensionless]')
    ax_case2_t1.set_xlabel('$x$ [dimensionless]')
    ax_case2_t1.set_xlabel('$x$ [dimensionless]')
    ax_case2_t2.set_ylabel('$y$ [dimensionless]')
    ax_case2_t2.set_ylabel('$y$ [dimensionless]')

    ax_case1_t1.grid('on')
    ax_case1_t2.grid('on')
    ax_case2_t1.grid('on')
    ax_case2_t2.grid('on')

    x,fe = pFE_case1.getState(t1)
    ax_case1_t1.plot(x,fe,label='ForwardEuler')
    x,fe = pFE_case2.getState(t1)
    ax_case2_t1.plot(x,fe,label='ForwardEuler')

    x,be = pBE_case1.getState(t1)
    ax_case1_t1.plot(x,be,label='BackwardEuler')
    x,be = pBE_case1.getState(t1)
    ax_case2_t1.plot(x,be,label='BackwardEuler')

    x = np.linspace(0,1,11)
    ax_case1_t1.plot(x,cfs(x, t1), label='$cfs$')
    x = np.linspace(0,1,101)
    ax_case2_t1.plot(x,cfs(x, t1), label='$cfs$')
    
    x = np.linspace(0,1,11)
    ax_case1_t1.plot(x,cn_1_t1,label='CrankNicolson')
    x = np.linspace(0,1,101)
    ax_case2_t1.plot(x,cn_2_t1,label='CrankNicolson')

    x,fe = pFE_case2.getState(t2)
    ax_case1_t2.plot(x,fe,label='ForwardEuler')
    x,fe = pFE_case2.getState(t2)
    ax_case2_t2.plot(x,fe,label='ForwardEuler')

    x,be = pBE_case2.getState(t2)
    ax_case1_t2.plot(x,be,label='BackwardEuler')
    x,be = pBE_case2.getState(t2)
    ax_case2_t2.plot(x,be,label='BackwardEuler')

    x = np.linspace(0,1,11)
    ax_case1_t2.plot(x,cfs(x, t2), label='$cfs$')
    x = np.linspace(0,1,101)
    ax_case2_t2.plot(x,cfs(x, t2), label='$cfs$')

    x = np.linspace(0,1,11)
    ax_case1_t2.plot(x,cn_1_t2,label='CrankNicolson')
    x = np.linspace(0,1,101)
    ax_case2_t2.plot(x,cn_2_t2,label='CrankNicolson')

    ax_case1_t1.legend(loc='best')
    ax_case2_t1.legend(loc='best')
    ax_case1_t2.legend(loc='best')
    ax_case2_t2.legend(loc='best')

    plt.show()

if __name__ == '__main__':
    # p = DiffusionProblem(sys.argv[1])
    # p.animate()
    compareCases()
