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


@dataclass(frozen=False)
class Datum(Protocol):
    """
    This class serves as a template for datum data classes.  Each implementation
    of this class, must at a minimum contain the same attribute names and
    data types listed in this class.  Each attribute is created as a
    semi-private attribute in an effort to make it more difficult for future
    developers to over-ride the attributes.

    Attributes
        **_a:** The semi-major axis of the oblate spheroid in units of meters\n
        **_b:** The semi-minor axis of the oblate spheriod in units of meters\n
        **_ro:** The average earth radius in units of meters\n
        **_gm:** The gravitational constant in units of m3/s-2\n
        **_f:** The flattening factor\n
        **_e:** The eccentricity of the oblate spheroid\n
        **_e2:** The eccentricity raised to the second power\n
        **_e24:** The square of the second eccentricity\n
        **_g:** The average gravity or gravity at a latitude of 45 degreesn
    """
    _a: float  # m
    _b: float  # m
    _ro: float  # m
    _gm: float  # m3/s-2
    _f: float # unitless
    _e: float  # unitless
    _e2: float  # unitless
    _e24: float  # unitless
    _g: float  # m/s2
# --------------------------------------------------------------------------------


@dataclass(frozen=False)
class WGS77(Datum):
    """
    A dataclass containing the core elements of the World Geodetic Survey
    of 1977. (*The Department of Defense World Geodetic Survey*, Defense
    Mapping Agency, Washington D.C., May 1974)
    """
    _a: float = 6378135.0  # m
    _b: float = 6356750.5  # m
    _ro: float = 6371008.7714  # m
    _gm: float = 3986008.0e8  # m3/s-2
    _f: float = 0.003352779 # unitless
    _e: float = 0.08181881066  # unitless
    _e2: float = 0.006694317778  # unitless
    _e24: float = 0.00673949674227  # unitless
    _g: float = 9.80665  # m/s2
# --------------------------------------------------------------------------------


@dataclass(frozen=False)
class WGS80(Datum):
    """
    A dataclass containing the core elements of the World Geodetic Survey
    of 1980. (H. Moritz, *Geodetic Reference System*)
    """
    _a: float = 6378137.0  # m
    _b: float = 6356752.314  # m
    _ro: float = 6371008.7714  # m
    _gm: float = 3986005.0e8  # m3/s-2
    _f: float = 0.00335281068118 # unitless
    _e: float = 0.0818191908426  # unitless
    _e2: float = 0.0066943800229  # unitless
    _e24: float = 0.00673949677548  # unitless
    _g: float = 9.80665  # m/s2
# --------------------------------------------------------------------------------


@dataclass(frozen=False)
class WGS84(Datum):
    """
    A dataclass containing the core elements of the World Geodetic Survey
    of 1984. (*Department of Defense World Geodetic System 1984, Its definition
    and its relationship with local geodetic systems*, Defense Mapping Agency,
    DMA TR 8350.2, September 1, 1991)

   :param _f: The flattening factor
    """
    _a: float = 6378137.0  # m
    _b: float = 6356752.3142  # m
    _ro: float = 6371008.7714  # m
    _gm: float = 3986005.0e8  # m3/s-2
    _f: float = 0.00335281066474 # unitless
    _e: float = 0.0818191908426  # unitless
    _e2: float = 0.0066943799013  # unitless
    _e24: float = 0.00673949674227  # unitless
    _g: float = 9.80665  # m/s2
# ================================================================================
# ================================================================================
# eof
