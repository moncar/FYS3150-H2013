#!/usr/bin/env python
import sys
from ReadPlot import read
from matplotlib import pyplot as plt

x,y = read(sys.argv[1])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x,y)
fig.suptitle('Relative error')
ax.set_xlabel('$log_{10}(h)$')
ax.set_ylabel('$log_{10}(\epsilon_{max})$')
ax.grid('on')

plt.show()
