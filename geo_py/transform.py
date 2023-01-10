# Import necessary packages here
import sys
import os
from pathlib import PurePath
from typing import Tuple
from math import radians, cos, sin, sqrt
p = PurePath(__file__).parent
sys.path.insert(1, os.path.abspath(p))
from geo_py.datum import Datum
# ================================================================================
# ================================================================================
# File:    transform.py
# Date:    January 09, 2023
# Purpose: Describe the purpose of functions of this file

# Source Code Metadata
__author__ = "Jonathan A. Webb"
__copyright__ = "Copyright 2023, Jon Webb Inc."
__version__ = "1.0"
# ================================================================================
# ================================================================================
# Insert Code here


class Transformations:
    """
    This class acts as a model container for several data geodetic to
    Cartesian, and Cartesian to geodetic coordinate frames.

    :param datum: A geodetic datum of type Datum.  This Protocol dataclass
                  can either be created by the user, or an existing datum
                  can be imported from the geo_py.datum file
    """
    def __init__(self, datum: Datum):
        self.datum = datum
# --------------------------------------------------------------------------------

    def llh_to_ecef(self, lat: float, lon: float,
                    alt: float) -> Tuple[float, float, float]:
        """
        :param lat: The latitude in decimal degrees
        :param lon: The longitude in decimal degrees
        :param alt: The altitude above sea-level in units of meters
        :return x, y, z: The ECEF coordinates in units of meters

        This method will transform a data from a LLH (Latitude, Longitude,
        Height) coordinate frame to an ECEF (Earth Centered Earth Fixed)
        coordinate frame.  The LLH frame is a polar coordinate frame with its
        origin at the center of the Earth.  The ECEF frame is a Cartesian
        coordinate system with its origin at the center of the Earth.

        This method solves the following equations for :math:`X`, :math:`Y`,
        and :math:`Z` assuming inputs of latitude (:math:`\\lambda`),
        longitude (:math:`\\phi`), and altitude (:math:`h`).  The terms
        :math:`a`, :math:`b`, and :math:`e` represents the semi-major axis,
        the semi-minor axis and the eccentricity of the spheroid, and are
        derived from the Datum dataclass

        .. math::

           \\begin{eqnarray}
           X = \\left(N\\left(\\phi\\right) + h\\right)cos\\phi \\: cos\\lambda \\\\
           Y = \\left(N\\left(\\phi\\right) + h\\right)cos\\phi \\: sin\\lambda \\\\
           Z = \\left(\\frac{b^2}{a^2}N\\left(\\phi\\right) + h \\right)sin\\:\\phi
          \\end{eqnarray}

        where

        .. math::

           N\\left(\\phi\\right) = \\frac{a}{\\sqrt[]{1-e^2sin^2\\:\\phi}}
        """
        lon = radians(lon)
        lat = radians(lat)
        N = self.datum._a / (sqrt(1.0 - \
                self.datum._e ** 2.0 * (sin(lat)) ** 2.0))
        xe = (N + alt) * cos(lat) * cos(lon)
        ye = (N + alt) * cos(lat) * sin(lon)
        ze = (((self.datum._b ** 2.0 / self.datum._a ** 2.0) * N) + \
                alt) * sin(lat)
        return xe, ye, ze
# ================================================================================
# ================================================================================
# eof
