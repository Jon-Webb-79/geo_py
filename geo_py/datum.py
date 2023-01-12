# Import necessary packages here
from typing import Protocol
from dataclasses import dataclass
# ================================================================================
# ================================================================================
# File:    datum.py
# Date:    January 09, 2023
# Purpose: This file contains data classes each representing different datums

# Source Code Metadata
__author__ = "Jonathan A. Webb"
__copyright__ = "Copyright 2023, Jon Webb Inc."
__version__ = "1.0"
# ================================================================================
# ================================================================================
# Insert Code here

@dataclass
class Datum(Protocol):
    """
    This class serves as a template for datum data classes.  Each implementation
    of this class, must at a minimum contain the same attribute names and
    data types listed in this class.  Each attribute is created as a
    semi-private attribute in an effort to make it more difficult for future
    developers to over-ride the attributes.

    Attributes
        **a:** The semi-major axis of the oblate spheroid in units of meters\n
        **b:** The semi-minor axis of the oblate spheriod in units of meters\n
        **f:** The flattening factor\n
        **e:** The eccentricity of the oblate spheroid\n
        **e2:** The eccentricity raised to the second power\n
        **GM:** The gravitational constant in units of m3/s-2\n
        **R_avg:** The average earth radius in units of meters\n
        **e24:** The square of the second eccentricity\n
        **g_avg:** The average gravity or gravity at a latitude of 45 degreesn
    """
    a: float
    b: float
    f: float
    e: float
    e2: float
    GM: float
    R_avg : float
    e24 : float
    g_avg: float
# --------------------------------------------------------------------------------

@dataclass(slots=True)
class WGS77(Datum):
    """
    A dataclass containing the core elements of the World Geodetic Survey
    of 1977. (*The Department of Defense World Geodetic Survey*, Defense
    Mapping Agency, Washington D.C., May 1974)
    """
    a: float  = 6378135.0  # semi-major axis (meters)
    b: float  = 6356750.5200160944  # semi-minor axis (meters)
    f: float  = 1/298.26 # flattening
    e: float  = 0.08181919084262157  # first eccentricity
    e2: float = 0.0066943177722 # first eccentricity squared
    GM: float = 3.986004418e14 #  Gravitational parameter
    R_avg : float = 6371000 # Average radius of the earth
    e24 : float = 0.00669438002290342 # second eccentricity squared
    g_avg: float = 9.780318 # average gravitational acceleration (m/s^2)
# --------------------------------------------------------------------------------

@dataclass(slots=True)
class WGS80(Datum):
    """
    A dataclass containing the core elements of the World Geodetic Survey
    of 1980. (H. Moritz, *Geodetic Reference System*)
    """
    a: float  = 6378137.0  # semi-major axis (meters)
    b: float  = 6356752.314140346  # semi-minor axis (meters)
    f: float  = 1/298.257223563  # flattening
    e: float  = 0.081819190842622  # first eccentricity
    e2: float = 0.00669437999014132  # first eccentricity squared
    GM: float = 3986005e8 #  Gravitational parameter
    R_avg : float = 6371008.7714150598 # Average radius of the earth
    e24 : float = 0.006694380022902712  # second eccentricity squared
    g_avg: float = 9.7803267715 # average gravitational acceleration (m/s^2)
# --------------------------------------------------------------------------------

@dataclass(slots=True)
class WGS84(Datum):
    """
    A dataclass containing the core elements of the World Geodetic Survey
    of 1984. (*Department of Defense World Geodetic System 1984, Its definition
    and its relationship with local geodetic systems*, Defense Mapping Agency,
    DMA TR 8350.2, September 1, 1991)
    """
    a: float  = 6378137.0  # semi-major axis (meters)
    b: float  = 6356752.314245  # semi-minor axis (meters)
    f: float  = 1/298.257223563  # flattening
    e: float  = 0.08181919084262157  # first eccentricity
    e2: float = 0.00669437999014132  # first eccentricity squared
    GM: float = 3986005e8 #  Gravitational parameter
    R_avg : float = 6371008.8 # Average radius of the earth
    e24 : float = 0.006694317779487594 # second eccentricity squared
    g_avg: float = 9.7803253359 # average gravitational acceleration (m/s^2)
# --------------------------------------------------------------------------------

@dataclass(slots=True)
class NAD83(Datum):
    """
    A dataclass containing core elements of the NAD83 datum, which should
    only be used for North America
    """
    a: float  = 6378137.0  # semi-major axis (meters)
    b: float  = 6356752.3141409  # semi-minor axis (meters)
    f: float  = 1/298.257222101  # flattening
    e: float  = 0.08181919084261  # first eccentricity
    e2: float = 0.006694380022900  # first eccentricity squared
    GM: float = 3986004.418e8 #  Gravitational parameter
    R_avg : float = 6371000 # Average radius of the earth
    e24 : float = 0.0066943800229034  # second eccentricity squared
    g_avg: float = 9.780318  # average gravitational acceleration (m/s^2)
# --------------------------------------------------------------------------------

@dataclass(slots=True)
class ETRS89(Datum):
    """
    A dataclass containing core elements of the ETRS89 datum, which should
    only be used for Europe
    """
    a: float  = 6378137.0  # semi-major axis (meters)
    b: float  = 6356752.314245  # semi-minor axis (meters)
    f: float  = 1/298.257223563  # flattening
    e: float  = 0.0818191910428158  # first eccentricity
    e2: float = 0.00669438002290343  # first eccentricity squared
    GM: float = 3986005e8 #  Gravitational parameter
    R_avg : float = 6371008.8  # Average radius of the earth
    e24 : float = 0.00669438002290343  # second eccentricity squared
    g_avg: float = 9.7803253359 # average gravitational acceleration (m/s^2)
# --------------------------------------------------------------------------------

@dataclass(slots=True)
class ITRF(Datum):
    a: float  = 6378136.6  # semi-major axis (meters)
    b: float  = 6356751.9  # semi-minor axis (meters)
    f: float  = 1/298.2572236 # flattening
    e: float  = 0.08181919084 # first eccentricity
    e2: float = 0.006694380023 # first eccentricity squared
    GM: float = 3.986004418e14 #  Gravitational parameter
    R_avg : float = 6371000 # Average radius of the earth
    e24 : float = 0.00669438002290342 # second eccentricity squared
    g_avg: float = 9.780318 # average gravitational acceleration (m/s^2)
# ================================================================================
# ================================================================================
# eof
