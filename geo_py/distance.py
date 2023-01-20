# Import necessary packages here
import sys
import os
from pathlib import PurePath
from typing import Tuple
from math import radians, sqrt, cos, sin, atan2, acos, tan, atan, pi, asin
p = PurePath(__file__).parent
sys.path.insert(1, os.path.abspath(p))
from geo_py.datum import Datum, WGS84
# ================================================================================
# ================================================================================
# File:    distance.py
# Date:    January 17, 2023
# Purpose: This file contains functions and classes used to calculate distances
#          between two points

# Source Code Metadata
__author__ = "Jonathan A. Webb"
__copyright__ = "Copyright 2023, Jon Webb Inc."
__version__ = "1.0"
# ================================================================================
# ================================================================================
# Insert Code here


def haversine(lla1: Tuple[float, float, float],
              lla2: Tuple[float, float, float],
              cruise_alt: float=0.0, dat: Datum=WGS84()):
    """
    :param lla1: A tuple containing latitude, longitude and altitude
                 of a point on the surface of a body.  Latitude and longitude
                 are in units of decimal degrees, altitude in units of meters
    :param lla2: A tuple containing latitude, longitude and altitude
                 of a point on the surface of a body.  Latitude and longitude
                 are in units of decimal degrees, altitude in units of meters
    :param cruise_alt: The cruising altitude, if this method is used for an
                       aircraft.  Units in meters. Defaulted to 0
    :param dat: A Datum dataclass, defaulted to WGS84
    :return dist: The curvelinear distance in units of meters.

    This function will calculate the distance between to points using the Haversine
    method, assuming they are
    on the surface of a spherical body.  The altitude in lla1 and lla2 represent
    surface points such as hills, and mountains that would elevate the observer
    above the mean radius.  If this method is used to predict travel distance
    for an aircraft, the altitudes in lla1 and lla2 should be set to zero,
    and the cruise_alt variable should be envoked. This method solves the
    following equations where :math:`\\phi` and :math:`\\theta` represent
    latitude and longitude and :math:`r` represents the average radius
    of the Earth.

    .. math::

        \\begin{equation}
            d = 2rsin^{-1}\\sqrt[]{sin^2\\left(\\frac{\\phi_2-\\phi_1}{2}\\right) +
            cos\\left(\\phi_1\\right)cos\\left(\\phi_2\\right)sin^2\\left(\\frac{\\theta_2-\\theta_1}{2}\\right)}
        \\end{equation}

    Code Example

    .. code-block::

        from geo_py.distance import haversine

        lla1 = (40.98197, 111.9026, 0.0)
        lla2 = (45.01443, 113.9278, 0.0)
        d = haversine(lla1, lla2)
        print(d)
        >>> 477637.629
    """
    lat1, lon1, alt1 = radians(lla1[0]), radians(lla1[1]), lla1[2]
    lat2, lon2, alt2 = radians(lla2[0]), radians(lla2[1]), lla2[2]
    earth_radius = dat.R_avg + cruise_alt
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a_var = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c_var = 2 * atan2(sqrt(a_var), sqrt(1 - a_var))
    distance = earth_radius * c_var
    altitude = abs(alt2 - alt1)
    return sqrt(distance**2 + altitude**2)
# --------------------------------------------------------------------------------


def great_circle(lla1: Tuple[float, float, float],
                 lla2: Tuple[float, float, float],
                 cruise_alt: float=0.0, dat: Datum=WGS84()):
    """
    :param lla1: A tuple containing latitude, longitude and altitude
                 of a point on the surface of a body.  Latitude and longitude
                 are in units of decimal degrees, altitude in units of meters
    :param lla2: A tuple containing latitude, longitude and altitude
                 of a point on the surface of a body.  Latitude and longitude
                 are in units of decimal degrees, altitude in units of meters
    :param cruise_alt: The cruising altitude, if this method is used for an
                       aircraft.  Units in meters. Defaulted to 0
    :param dat: A Datum dataclass, defaulted to WGS84
    :return dist: The curvelinear distance in units of meters.

    This function will calculate the distance between two points using the
    great circle method.  If this function is used for calculating the distance
    of an aircraft traveled while in level flight, the altitude terms in the
    two tuples should be zero, and the cruise altitude term should be used
    instead.  The great cricle method impliments a spherical triangle method
    where :math:`\\phi`, and :math:`\\theta` represent latitude and longitude.
    The equation below describes the mathematics of the great circle method.

    .. math::

       \\begin{equation}
           d=cos^{-1}\\left[sin\\left(\\phi_1\\right)\\:sin\\left(\\phi_2\\right) +
           cos\\left(\\phi_1\\right)\\:cos\\left(\\phi_2\\right)
           \\:cos\\left(\\theta_1-\\theta_2\\right)\\right]
       \\end{equation}

    Code Example

    .. code-block::

        from geo_py.distance import great_circle
        lla1 = (40.98197, 111.9026, 0.0)
        lla2 = (45.01443, 113.9278, 0.0)
        d = great_circle(lla1, lla2)
        print(d)
        >>> 477637.629
    """
    lat1, lon1, alt1 = lla1[0], lla1[1], lla1[2]
    lat2, lon2, alt2 = lla2[0], lla2[1], lla2[2]
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])
    earth_radius = dat.R_avg + cruise_alt  # in meters
    distance = earth_radius * acos(sin(lat1) * sin(lat2) +
                                   cos(lat1) * cos(lat2) * cos(lon1 - lon2))
    distance = sqrt(distance**2 + (alt2-alt1)**2)
    return distance
# --------------------------------------------------------------------------------


def vincenty(lla1: Tuple[float, float, float],
             lla2: Tuple[float, float, float],
             dat: Datum=WGS84()):
    """
    :param lla1: A tuple containing latitude, longitude and altitude
                 of a point on the surface of a body.  Latitude and longitude
                 are in units of decimal degrees, altitude in units of meters
    :param lla2: A tuple containing latitude, longitude and altitude
                 of a point on the surface of a body.  Latitude and longitude
                 are in units of decimal degrees, altitude in units of meters
    :param cruise_alt: The cruising altitude, if this method is used for an
                       aircraft.  Units in meters. Defaulted to 0
    :param dat: A Datum dataclass, defaulted to WGS84
    :return dist: The curvelinear distance in units of meters.

    This function calculates the distance between to points in a geodetic
    frame using the Vincenty method, which is considered to be the most
    accurate method for determining distance. The Vincenty mathematical
    formula is to complex to document here; however, it can be found at
    `Vincenty <https://en.wikipedia.org/wiki/Vincenty%27s_formulae>`_

    Code Example

    .. code-block::

        from geo_py.distance import vincenty

        lla1 = (40.98197, 111.9026, 0.0)
        lla2 = (45.01443, 113.9278, 0.0)
        d = vincenty(lla1, lla2)
        print(d)
        >>> 477402.505
    """
    lat1, lon1, alt1 = lla1[0], lla1[1], lla1[2]
    lat2, lon2, alt2 = lla2[0], lla2[1], lla2[2]
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])
    tolerance = 1e-12 # to stop iteration

    phi1, phi2 = lat1, lat2
    U1 = atan((1-dat.f)*tan(phi1))
    U2 = atan((1-dat.f)*tan(phi2))
    L1, L2 = lon1, lon2
    L = L2 - L1

    lambda_old = L + 0

    while True:

        t = (cos(U2)*sin(lambda_old))**2
        t += (cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lambda_old))**2
        sin_sigma = t**0.5
        cos_sigma = sin(U1)*sin(U2) + cos(U1)*cos(U2)*cos(lambda_old)
        sigma = atan2(sin_sigma, cos_sigma)

        sin_alpha = cos(U1)*cos(U2)*sin(lambda_old) / sin_sigma
        cos_sq_alpha = 1 - sin_alpha**2
        cos_2sigma_m = cos_sigma - 2*sin(U1)*sin(U2)/cos_sq_alpha
        C = dat.f*cos_sq_alpha*(4 + dat.f*(4-3*cos_sq_alpha))/16

        t = sigma + C*sin_sigma*(cos_2sigma_m + C*cos_sigma*(-1 + 2*cos_2sigma_m**2))
        lambda_new = L + (1 - C)*dat.f*sin_alpha*t
        if abs(lambda_new - lambda_old) <= tolerance:
            break
        else:
            lambda_old = lambda_new

    u2 = cos_sq_alpha*((dat.a**2 - dat.b**2)/dat.b**2)
    A = 1 + (u2/16384)*(4096 + u2*(-768+u2*(320 - 175*u2)))
    B = (u2/1024)*(256 + u2*(-128 + u2*(74 - 47*u2)))
    t = cos_2sigma_m + 0.25*B*(cos_sigma*(-1 + 2*cos_2sigma_m**2))
    t -= (B/6)*cos_2sigma_m*(-3 + 4*sin_sigma**2)*(-3 + 4*cos_2sigma_m**2)
    delta_sigma = B * sin_sigma * t
    s = dat.b*A*(sigma - delta_sigma)

    distance = sqrt(s**2.0 + (alt2 - alt1)**2.0)
    return distance
# --------------------------------------------------------------------------------


def linear_dist(lla1: Tuple[float, float, float],
                lla2: Tuple[float, float, float],
                dat: Datum=WGS84()):
    """
    :param lla1: A tuple containing latitude, longitude and altitude
                 of a point on the surface of a body.  Latitude and longitude
                 are in units of decimal degrees, altitude in units of meters
    :param lla2: A tuple containing latitude, longitude and altitude
                 of a point on the surface of a body.  Latitude and longitude
                 are in units of decimal degrees, altitude in units of meters
    :param cruise_alt: The cruising altitude, if this method is used for an
                       aircraft.  Units in meters. Defaulted to 0
    :param dat: A Datum dataclass, defaulted to WGS84
    :return dist: The linear distance between two points in units of meters.

    This function will determine the linear distance between two points of
    latitudes, longitudes, and altitudes.  This method converts each
    latitude, longitude, altitude pair into ECEF coordinates and determines
    the distance as the magnitude of the distance between the X, Y, and Z
    points.

    Code Example

    .. code-block::

        from geo_py.distance import linear_dist

        lla1 = (40.98197, 111.9026, 0.0)
        lla2 = (45.01443, 113.9278, 0.0)
        d = linear_dist(lla1, lla2)
        print(d)
        >>> 477965.401
    """
    lat1, lon1, alt1 = lla1[0], lla1[1], lla1[2]
    lat2, lon2, alt2 = lla2[0], lla2[1], lla2[2]
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])

    N1 = dat.a / (sqrt(1.0 - \
        dat.e ** 2.0 * (sin(lat1)) ** 2.0))
    xe1 = (N1 + alt1) * cos(lat1) * cos(lon1)
    ye1 = (N1 + alt1) * cos(lat1) * sin(lon1)
    ze1 = (((dat.b ** 2.0 / dat.a ** 2.0) * N1) + \
            alt1) * sin(lat1)

    N2 = dat.a / (sqrt(1.0 - \
        dat.e ** 2.0 * (sin(lat2)) ** 2.0))
    xe2 = (N1 + alt2) * cos(lat2) * cos(lon2)
    ye2 = (N1 + alt2) * cos(lat2) * sin(lon2)
    ze2 = (((dat.b ** 2.0 / dat.a ** 2.0) * N2) + \
            alt2) * sin(lat2)
    dist = sqrt((xe1-xe2)**2 + (ye1-ye2)**2 + (ze1-ze2)**2)
    return dist
# ================================================================================
# ================================================================================


class Distance:
    """
    :param lla1: A tuple containing latitude, longitude and altitude
                 of a point on the surface of a body.  Latitude and longitude
                 are in units of decimal degrees, altitude in units of meters
    :param lla2: A tuple containing latitude, longitude and altitude
                 of a point on the surface of a body.  Latitude and longitude
                 are in units of decimal degrees, altitude in units of meters
    :param method: 'haversine', 'great circle', or 'vincenti'. Not case
                   sensitive. Defaulted to haversine.
    :param cruise_alt: The cruising altitude, if this method is used for an
                       aircraft.  Units in meters. Defaulted to 0
    :param dat: A Datum dataclass, defaulted to WGS84

    This class will calculate the distance between two points using either the
    Haversine, Great Circle, or Vincenti methods.

    Code Example

    .. code-block::

        from geo_py.distance import Haversine
        from geo_py.datums import ITRF

        lla1 = (40.98197, 111.9026, 0.0)
        lla2 = (45.01443, 113.9278, 0.0)
        dist = Haversine(lla1, lla2, method='haversine').km
        print(dist)
        >>> 477.637

        # Use a different coordinate system and a flight altitude
        dist = Distance(lla1, lla2, method='haversine',
                        cruise_alt=3400.0, dat=ITRF()).miles
        print(dist)
        >>> 296.948
    """
    def __init__(self, lla1: Tuple[float, float, float],
                 lla2: Tuple[float, float, float], method: str = "HAVERSINE",
                 cruise_alt: float = 0.0, dat: Datum = WGS84()):
        method = method.upper()
        if method == "HAVERSINE":
            self.distance = haversine(lla1, lla2, cruise_alt, dat)
        elif method == "GREAT CIRCLE":
            self.distance = great_circle(lla1, lla2, cruise_alt, dat)
        elif method == "VINCENTY":
            self.distance = vincenty(lla1, lla2, dat)
        elif method == "LINEAR":
            self.distance = linear_dist(lla1, lla2, dat)
        else:
            raise ValueError("Invalid method for Distance class")
# --------------------------------------------------------------------------------

    def __str__(self) -> str:
        """
        This method will returna  string representation of the object
        """
        return f"{self.distance/1000} km"
# --------------------------------------------------------------------------------

    def __setattr__(self, name, value):
        """
        This method prevents a developer from adding attributes to this class
        """
        if name in ["distance"]:
            object.__setattr__(self, name, value)
        else:
            raise AttributeError(f"Attribute {name} not allowed")
# --------------------------------------------------------------------------------

    @property
    def km(self) -> float:
        """
        This method will return the distance in units of kiometers
        """
        return self.distance/1000.0
# --------------------------------------------------------------------------------

    @property
    def m(self) -> float:
        """
        This method will return the distance in units of meters
        """
        return self.distance
# --------------------------------------------------------------------------------

    @property
    def miles(self) -> float:
        """
        This method will return the distance in units of miles
        """
        return self.distance * 0.000621371
# --------------------------------------------------------------------------------

    @property
    def feet(self) -> float:
        """
        This method will return the distance in units of feet
        """
        return self.distance * 3.28084
# --------------------------------------------------------------------------------

    def _calculate(self, lla1: Tuple[float, float, float],
                   lla2: Tuple[float, float, float],
                   cruise_alt: float, dat: Datum):
        """
        :param lla1: A tuple containing latitude, longitude and altitude
                     of a point on the surface of a body.  Latitude and longitude
                     are in units of decimal degrees, altitude in units of meters
        :param lla2: A tuple containing latitude, longitude and altitude
                     of a point on the surface of a body.  Latitude and longitude
                     are in units of decimal degrees, altitude in units of meters
        :param cruise_alt: The cruising altitude, if this method is used for an
                           aircraft.  Units in meters. Defaulted to 0
        :param dat: A Datum dataclass, defaulted to WGS84
        """
        lat1, lon1, alt1 = radians(lla1[0]), radians(lla1[1]), lla1[2]
        lat2, lon2, alt2 = radians(lla2[0]), radians(lla2[1]), lla2[2]
        earth_radius = dat.R_avg + cruise_alt
        dlat = lat2 - lat1
        dlon = lon2 - lon1
        a_var = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
        c_var = 2 * atan2(sqrt(a_var), sqrt(1 - a_var))
        distance = earth_radius * c_var
        altitude = abs(alt2 - alt1)
        return sqrt(distance**2 + altitude**2)
# ================================================================================
# ================================================================================
# eof
