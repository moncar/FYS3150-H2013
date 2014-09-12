#!/usr/bin/env python
"""
"""
import numpy as np

from numpy import pi, sin, exp

def cfs(x, t, n=1e3):
    """
    Closed form solution of diffusion equation.
    """
    if isinstance(x, (int,float)):
        _sum = 0
        for i in range(1,int(n)):
            _sum += 2. / (i*pi) * sin(i*pi*x) * exp(-i*i*pi*pi*t)
        return 1 - x - _sum
    elif isinstance(x, np.ndarray):
        res = np.zeros(len(x))
        for i in range(len(res)):
            res[i] = cfs(x[i], t)
        return res
