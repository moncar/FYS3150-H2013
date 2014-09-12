#!/usr/bin/env python
"""
Created on Thu 31 Oct 2013

Script has function for converting standard kilo mass into a scaled
mass with certain unit shown below. This can be used together with
simulation for solar system in order to get rid of the
gravitational constant.

@author Benedicte Emilie Braekken
"""
import unittest
from numpy import pi

"""Constants"""

_SUN_MASS = 2e30 # [kg]
_EARTH_KILO_MASS = 6e24 # [kg]
_EARTH_SCALED_MASS = 0.0001185323 # [G*AU^3/yr^2]
_TOLERANCE = 1e-7 # Tolerance for correct scaled mass

"""Test class"""

class TestMassFunction(unittest.TestCase):
    """
    Testcase for testing scaleMass function.
    """
    def test_scaleMass(self):
        """
        Tests the scaleMass function that it returns correctly
        scaled mass for earth.
        """
        self.assertTrue(abs(scaleMass(_EARTH_KILO_MASS)-_EARTH_SCALED_MASS)
                < _TOLERANCE)

"""Functions"""

def scaleMass(kiloMass):
    """
    Returns a scaled mass relative to sun mass. Relieves you of
    necessity to use gravitational constant.

    @param kiloMass Mass of object in kilos.
    """
    return (kiloMass / _SUN_MASS) * 4 * pi * pi

if __name__ == '__main__':
    unittest.main()
