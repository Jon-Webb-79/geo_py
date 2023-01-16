import sys
import os
from pathlib import PurePath
from typing import Tuple
from math import radians, cos, sin, sqrt, atan, atan2, pi, degrees
import numpy as np
p = PurePath(__file__).parent
sys.path.insert(1, os.path.abspath(p))
from geo_py.datum import Datum, WGS84
from geo_py.rotations import rotation_matrix, aircraft_rotation
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


def llh_to_ecef(lat: float, lon: float,
                alt: float, dat: Datum = WGS84()) -> Tuple[float, float, float]:
    """
    :param lat: The latitude in decimal degrees
    :param lon: The longitude in decimal degrees
    :param alt: The altitude above sea-level in units of meters
    :param datum: A Datum dataclass, defaulted to WGS84
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

    Code example

    .. code-block::

        from geo_py.datum import ITRF
        from geo_py.frames import llh_to_ecef

        lat = 46.826
        lon = 107.321
        alt = 6096.0
        x, y, z = llh_to_ecef(lat, lon, alt)
        print(x, y, z)
        >>> -1302839.38, 4177542.21, 4632996.83

        # Use a different datum
        x, y, z = llh_to_ecef(lat, lon, alt, dat=ITRF())
        print(x, y, z)
        >>> -1303839.30, 4177541.94, 4632996.51
    """
    lon = radians(lon)
    lat = radians(lat)
    N = dat.a / (sqrt(1.0 - \
        dat.e ** 2.0 * (sin(lat)) ** 2.0))
    xe = (N + alt) * cos(lat) * cos(lon)
    ye = (N + alt) * cos(lat) * sin(lon)
    ze = (((dat.b ** 2.0 / dat.a ** 2.0) * N) + \
       alt) * sin(lat)
    return xe, ye, ze
# --------------------------------------------------------------------------------


def llh_to_enu(origin_lat: float, origin_lon: float, origin_alt: float,
               lat: float, lon: float,
               alt: float, dat: Datum = WGS84()) -> Tuple[float, float, float]:
    """
    :param origin_lat: The latitude for the coordinate origin in units of
                       decimal degrees
    :param origin_lon: The longitude for the coordinate origin in units of
                       decimal degrees
    :param origin_alt: The altitude for the coordinate of the origin in units
                       of meters
    :param lat: The latitude for the point of interest in units of decimal
                degrees
    :param lon: The latitude for the point of interest in units of decimal
                degrees
    :param alt: The altitude for the point of interest in units od meters
    :return E, N, U: Cartesian coordinates in the East, North and Up directions.

    This function will transfomr coordinate from a LLH (Latitude, Longitude,
    Height) from to an ENU (East, North, Up) frame.  This function first
    converts the LLH frame to an ECEF frame, then to an ENU frame.

    Code Example

    .. code-block::

        from geo_py.datum import ITRF
        from geo_py.frames import llh_to_enu

        radar_lat = 46.017
        radar_lon = 7.750
        radar_alt = 1673.0

        x = 45.976
        y = 7.658
        z = 4531.0
        E, N, U = llh_to_enu(radar_lat, radar_lon, radar_alt, x, y, z)
        print(E, N, U)
        >>> -7134.757, -4556.321, 2852.390

        # Try another coordinate frame
        E, N, U = llh_to_enu(radar_lat, radar_lon, radar_alt, x, y, z, dat=ITRF())
        print(E, N, U)
        >>> -7134.757, -4556.307, 2852.806
    """
    # ref is aircraft
    x, y, z = llh_to_ecef(lat, lon, alt)
    return ecef_to_enu(origin_lat, origin_lon, origin_alt, x, y, z, dat=dat)
# --------------------------------------------------------------------------------


def llh_to_ned(origin_lat: float, origin_lon: float,
               origin_alt: float, lat: float, lon: float,
               alt: float, dat: Datum = WGS84()) -> Tuple[float, float, float]:
    """
    :param origin_lat: The latitude of the origin in units of decimal degrees
    :param origin_lon: The longitude of the origin in units of decimal degrees
    :param origin_alt: The altitude of the origin in units of meters
    :param lat: The latitude of the point of interest in units of decimal
                degrees
    :param lon: The longitude of the point of interest in units of decimal
                degrees
    :param alt: The altitude of the point of interest in units of meters
    :param dat: A Datum dataclass, defaulted to WGS84
    :return N, E, D: The North, East, Down coordinates in units of meters

    This function transforms coordinates from a LLH (Latitude, Longitude,
    Height) frame to a NED (North, East, Down) Frame.

    Code Example

    .. code-block::

        from geo_py.datum import ITRF
        from geo_py.frames import llh_to_enu

        radar_lat = 46.017
        radar_lon = 7.750
        radar_alt = 1673.0

        x = 45.976
        y = 7.658
        z = 4531.0
        N, E, D = llh_to_ned(radar_lat, radar_lon, radar_alt, x, y, z)
        print(N, E, D)
        >>> -4556.321, -7134.752, -2852.39
    """
    x, y, z = llh_to_ecef(lat, lon, alt, dat=dat)
    return ecef_to_ned(origin_lat, origin_lon, origin_alt, x, y, z, dat=dat)
# ================================================================================
# ================================================================================
# ECEF FUNCTIONS


def ecef_to_llh(x: float, y: float,
                z: float, dat: Datum = WGS84()) -> Tuple[float, float, float]:
    """
    :param x: The x-coordinate in ECEF coordinate frame
    :param y: The y-coordinate in ECEF coordinate frame
    :param z: The z-coordinate in ECEF coordinate frame
    :param dat: A Datum dataclass, defaulted to WGS84
    :return lat, lon, alt: The latitude, longitude, and altitude in the
                           LLH coordinate frame. Latitude and longitude
                           are in units of decimal degrees, altitude
                           in units of meters

    This method will transform a data from an ECEF (Earth Centered Earth
    Fixed) coordinate frame to a LLH (Latitude, Longitude, Height)
    coordinate frame.  The LLH frame is a polar coordinate frame with its
    origin at the center of the Earth.  The ECEF frame is a Cartesian
    coordinate system with its origin at the center of the Earth.

    This method takes the semi-major and semi-minor variables
    (:math:`a`, :math:`b`) from the Datum dataclass and use them with the
    ECEF inputs :math:`X`, :math:`Y`, :math:`Z` to calculate
    latitude, longitude, and height :math:`\\phi`, :math:`\\lambda`,
    :math:`h`.  This method is iterative and starts by making the following
    assumptions;

    .. math::

        \\begin{eqnarray}
            \\lambda = atan2\\left(Y,X\\right) \\\\
            p = \\sqrt[]{X^2+Y^2} \\\\
            \\phi = atan2\\left(Z,P\\left(1-e^2\\right)\\right)
        \\end{eqnarray}


    The algorithm then iterates over the following variables until the difference
    between :math:`\\phi_1` and :math:`\\phi` is :math:`10^{-14}` or less.

    .. math::

        \\begin{eqnarray}
            N = \\frac{a}{\\sqrt[]{1-e^2sin^2\\left(\\phi\\right)}} \\\\
            \\phi_1 = atan2\\left(Z+Ne^2sin^2\\left(\\phi\\right),p\\right) \\\\
        \\end{eqnarray}

    Once the solution has converged, which is usually within 4 to 5 iterations
    the altitude can eb deterimined via;

    .. math::

        h=\\frac{p}{cos\\left(\\phi\\right)}-N

    Finally the latitude and longitude are converted from radians to degrees.

    Code example

    .. code-block::

        from geo_py.datum import ITRF
        from geo_py.frames import ecef_to_llh

        tran = Transformations(WGS84())
        x = -1302839.38
        y = 4177542.21
        z = 4632996.83
        lat, lon, alt = ecef_to_llh(x, y, z)
        print(lat, lon, alt)
        >>> 46.826, 107.321, 6095.999

        # Use a different datum
        lat, lon, alt = ecef_to_llh(x, y, z, dat=ITRF())
        >>> 46.826, 107.321, 6095.999
    """
    lon = atan2(y, x)
    p_var = np.sqrt(x**2 + y**2)
    lat = atan2(z, p_var * (1 - dat.e**2))
    N = 0.0  # Initialize variable
    for _ in range(10):
        N = dat.a / np.sqrt(1 - dat.e**2 * sin(lat)**2)
        lat1 = atan2(z + N * dat.e**2 * sin(lat), p_var)
        if abs(lat1 - lat) < 1e-14:
            break
        lat = lat1
    alt = p_var / cos(lat) - N
    lat = np.rad2deg(lat)
    lon = np.rad2deg(lon)
    return lat, lon, alt
# --------------------------------------------------------------------------------

def ecef_to_enu(origin_lat: float, origin_lon: float, origin_alt: float,
                X: float, Y: float,
                Z: float, dat: Datum = WGS84()) -> Tuple[float, float, float]:
    """
    :param origin_lat: The latitude of the point of interest in units of
                       decimal degrees
    :param origin_lon: The longitude of the point of interest in units of
                       decimal degrees
    :param origin_alt: The altitude of the point of interest in units of decimal
                       degrees
    :param dat: A Datum dataclasse, defaulted to WGS84
    :param X: The x location of the craft in ECEF coordinates
    :param Y: The y location of the craft in ECEF coordinates
    :param Z: The z position of the craft in ECEF coordinats
    :return E, N, U: The x, y, and z position in the East, North, Up coordinate
            frame in units of metersy

    This method transforms a three dimensional location from an ECEF (Earth
    Centered Earth Fixed) coordinate frame to an ENU (East, North, Up)
    coordinate frame.  The ENU frame is determined relative to a location,
    which is the point of interest.  An example would be an aircraft as the
    reference point, and a runway as the point of interest.  This method
    solves the following equation where :math:`X_p`, :math:`Y_p`, and :math:`Z_p`
    are the origin point in ECEF coordinates, :math:`\\lambda`,
    :math:`\\lambda`, and :math:`\\phi` are the latitude, longitude and altitude
    of the reference point in decimal degrees and meters.  The term
    :math:`x`, :math:`y`, and :math:`z` are the NED coordinates in units
    of meters.  In addition :math:`X_r`, :math:`Y_r`, and :math:`Z_r` are
    the reference point expressed in ECEF coordinates

    .. math::

        \\begin{bmatrix}
            x \\\\
            y \\\\
            z \\\\
        \\end{bmatrix}
        =
        \\begin{bmatrix}
            -sin\\lambda & cos\\lambda & 0\\\\
            -sin\\phi\\:\\cos\\lambda & -sin\\phi\\:\\sin\\lambda & cos\\phi\\\\
            cos\\phi\\:cos\\lambda & cos\\phi\\:sin\\lambda & sin\\phi
        \\end{bmatrix}
        \\cdot
        \\begin{bmatrix}
            X_p-X_r \\\\
            Y_p-Y_r \\\\
            Z_p-Z_r
        \\end{bmatrix}

    Code example

    .. code-block::

        from geo_py.datum import ITRF
        from geo_py.frames import ecef_to_enu

        radar_lat = 46.017
        radar_lon = 7.750
        radar_alt = 1673.0

        x = 440357.328
        y = 592124.546
        z = 4566651.751
        E, N, U = ecef_to_enu(radar_lat, radar_lon, radar_alt, x, y, z)
        print(E, N, U)
        >>> -7134.757, -4556.321, 2852.390

        # Try another coordinate frame
        E, N, U = ecef_to_enu(radar_lat, radar_lon, radar_alt, x, y, z, dat=ITRF())
        print(E, N, U)
        >>> -7134.757, -4556.307, 2852.806
    """
    x_r, y_r, z_r = llh_to_ecef(origin_lat, origin_lon, origin_alt, dat=dat)

    lat = np.deg2rad(origin_lat)
    lon = np.deg2rad(origin_lon)
    E = -sin(lon)*(X-x_r) + cos(lon)*(Y-y_r);
    N = -sin(lat)*cos(lon)*(X-x_r) - sin(lat)*sin(lon)*(Y-y_r) + cos(lat)*(Z-z_r);
    U = cos(lat)*cos(lon)*(X-x_r) + cos(lat)*sin(lon)*(Y-y_r) + sin(lat)*(Z-z_r);
    return E, N, U
# --------------------------------------------------------------------------------


def ecef_to_ned(origin_lat: float, origin_lon: float, origin_alt: float,
                x: float, y: float, z: float,
                dat: Datum = WGS84()) -> Tuple[float, float, float]:
    """
    :param origin_lat: The latitude of the point of interest in units of
                       decimal degrees
    :param origin_lon: The longitude of the point of interest in units of
                       decimal degrees
    :param origin_alt: The altitude of the point of interest in units of decimal
                       degrees
    :param dat: A Datum dataclasse, defaulted to WGS84
    :param X: The x location of the craft in ECEF coordinates
    :param Y: The y location of the craft in ECEF coordinates
    :param Z: The z position of the craft in ECEF coordinats
    :return N, E, D: The x, y, and z position in the North, East, Down coordinate
            frame in units of meters

    This function converts the ECEF (Earth Centered Earth Fixed) coordinate
    frame to the NED (North, East, Down) coordinate frame. This function
    solves the following equations where :math:`e` represents the body's
    eccentricity, :math:`a` represents the bodies semi-major axis length,
    :math:`X`, :math:`Y', and :math:`Z` represents the ECEF coordinates
    of the origin point.  Finally :math:`\\phi`, :math:`\\lambda`, and
    :math:`h` represent the latitude, longitude and height of the origin
    point above sea level.

    .. math::

        \\begin{bmatrix}
            N \\\\
            E \\\\
            D \\\\
        \\end{bmatrix}
        =
        \\begin{bmatrix}
            -sin\\lambda\\:cos\\phi & sin\\lambda\\:sin\\phi & cos\\lambda\\\\
            -sin\\phi & cos\\phi & 0\\\\
            cos\\lambda\\:cos\\phi & -cos\\lambda\\:sin\\phi & -sin\\lambda
        \\end{bmatrix}
        \\cdot
        \\begin{bmatrix}
            X-X_r \\\\
            Y-Y_r \\\\
            Z-Z_r
        \\end{bmatrix}

    Code Example

    .. code-block::

        from geo_py.datum import ITRF
        from geo_py.frames import ecef_to_ned

        radar_lat = 46.017
        radar_lon = 7.750
        radar_alt = 1673.0

        x = -7134.757
        y = -4556.321
        z = 2852.390

        N, E, D = ecef_to_ned(radar_lat, radar_lon, radar_alt, x, y, z)
        print(x, y, z)
        >>> 28882.282, -3552.574, 6372030.823

        # Use another datum
        N, E, D = ecef_to_ned(radar_lat, radar_lon, radar_alt, x, y, z, ITRF())
        >>> 28882.297, -3552.574, 6372030.407
    """
    # Convert latitude, longitude to radians
    lat_o = np.deg2rad(origin_lat)
    lon_o = np.deg2rad(origin_lon)

    x_o, y_o, z_o = llh_to_ecef(origin_lat, origin_lon, origin_alt, dat=dat)

    # Compute the offset of the point from the origin
    dx = x - x_o
    dy = y - y_o
    dz = z - z_o

    # Compute the rotation matrix from ECEF to NED
    R_e2n = np.array([[-sin(lat_o)*cos(lon_o), -sin(lat_o)*sin(lon_o),  cos(lat_o)],
                      [-sin(lon_o),               cos(lon_o),           0        ],
                      [-cos(lat_o)*cos(lon_o), -cos(lat_o)*sin(lon_o), -sin(lat_o)]])

    # Rotate the offset vector from ECEF to NED frame
    offset_ned = np.dot(R_e2n, np.array([dx, dy, dz]))

    return tuple(offset_ned)
# ================================================================================
# ================================================================================
# ENU FUNCTIONS


def enu_to_ecef(ref_lat: float, ref_lon: float, ref_alt: float,
                E: float, N: float,
                U: float, dat: Datum = WGS84()) -> Tuple[float, float, float]:
    """
    :param origin_lat: The latitude of the point of interest in units of
                       decimal degrees
    :param origin_lon: The longitude of the point of interest in units of
                       decimal degrees
    :param origin_alt: The altitude of the point of interest in units of decimal
                       degrees
    :param E: The x location of the craft in ENU coordinates
    :param N: The y location of the craft in ENU coordinates
    :param U: The z position of the craft in ENU coordinats
    :param dat: A Datum dataclasse, defaulted to WGS84
    :return X, Y, Z: The x, y, and z position in ECEF coordinats

    This method transforms a three dimensional location from a ENU (East,
    North, Up) coordinate frame to an ECEF (Eearth Centered Earth Fixed)
    coordinate frame.  The ENU frame is determined relative to a location,
    which is the point of interest.  An example would be an aircraft as the
    reference point, and a runway as the point of interest.  This method
    solves the following equation where :math:`X_p`, :math:`Y_p`, and :math:`Z_p`
    are the origin point in ECEF coordinates, :math:`\\lambda`,
    :math:`\\lambda`, and :math:`\\phi` are the latitude, longitude and altitude
    of the reference point in decimal degrees and meters.  The term
    :math:`x`, :math:`y`, and :math:`z` are the NED coordinates in units
    of meters.  In addition :math:`X_r`, :math:`Y_r`, and :math:`Z_r` are
    the reference point expressed in ECEF coordinates

    .. math::

        \\begin{bmatrix}
            X \\\\
            Y \\\\
            Z \\\\
        \\end{bmatrix}
        =
        \\begin{bmatrix}
            -sin\\lambda & cos\\lambda & 0\\\\
            -sin\\phi\\:\\cos\\lambda & -sin\\phi\\:\\sin\\lambda & cos\\phi\\\\
            cos\\phi\\:cos\\lambda & cos\\phi\\:sin\\lambda & sin\\phi
        \\end{bmatrix}^T
        \\cdot
        \\begin{bmatrix}
            x \\\\
            y \\\\
            z \\\\
        \\end{bmatrix}
        +
        \\begin{bmatrix}
            X_r \\\\
            Y_r \\\\
            Z_r
        \\end{bmatrix}

    Code example

    .. code-block::

        from geo_py.datum import ITRF
        from geo_py.frames import enu_to_ecef

        radar_lat = 46.017
        radar_lon = 7.750
        radar_alt = 1673.0

        E = -7134.757
        N = -4556.321
        U = 2852.390

        x, y, z = enu_to_ecef(radar_lat, radar_lon, radar_alt, E, N, U)
        print(x, y, z)
        >>> 4403757.605, 592124.579, 45666.52

        # Use another datum
        x, y, z = enu_to_ecef(radar_lat, radar_lon, radar_alt, E, N, U, ITRF())
        >>> 440357.328, 592124.546, 4566651.751
    """
    x_ref, y_ref, z_ref = llh_to_ecef(ref_lat, ref_lon, ref_alt, dat=dat)

    ref_lat = np.deg2rad(ref_lat)
    ref_lon = np.deg2rad(ref_lon)

    X = -sin(ref_lon)*E - cos(ref_lon)*sin(ref_lat)*N + \
            cos(ref_lon)*cos(ref_lat)*U + x_ref;
    Y = cos(ref_lon)*E - sin(ref_lon)*sin(ref_lat)*N + \
            cos(ref_lat)*sin(ref_lon)*U + y_ref;
    Z = cos(ref_lat)*N + sin(ref_lat)*U + z_ref;
    return X, Y, Z
# --------------------------------------------------------------------------------


def enu_to_llh(origin_lat: float, origin_lon: float, origin_alt: float,
               E: float, N: float, U: float,
               dat: Datum = WGS84()) -> Tuple[float, float, float]:
    """
    :param origin_lat: The latitude of the point of interest in units of
                       decimal degrees
    :param origin_lon: The longitude of the point of interest in units of
                       decimal degrees
    :param origin_alt: The altitude of the point of interest in units of decimal
                       degrees
    :param E: The x location of the craft in ENU coordinates
    :param N: The y location of the craft in ENU coordinates
    :param U: The z position of the craft in ENU coordinats
    :param dat: A Datum dataclasse, defaulted to WGS84
    :return lat, lon, alt: The latitude, longitude, and altitude position
                           in LLH coordinats

    Code Example

    .. code-block::

        from geo_py.datum import ITRF
        from geo_py.frames import enu_to_llh

        radar_lat = 46.017
        radar_lon = 7.750
        radar_alt = 1673.0

        E = -7134.757
        N = -4556.321
        U = 2852.390

        x, y, z = enu_to_llh(radar_lat, radar_lon, radar_alt, E, N, U)
        print(x, y, z)
        >>> 45.976, 7.658, 4531.989

        # Use another datum
        x, y, z = enu_to_ecef(radar_lat, radar_lon, radar_alt, E, N, U, ITRF())
        >>> 45.976, 7.678, 4532.004
    """
    x, y, z = enu_to_ecef(origin_lat, origin_lon, origin_alt, E, N, U, dat)
    return ecef_to_llh(x, y, z, dat)
# --------------------------------------------------------------------------------


def enu_to_ned(origin_lat: float, origin_lon: float, origin_alt: float,
               E: float, N: float, U: float,
               dat: Datum = WGS84()) -> Tuple[float, float, float]:
    """
    :param origin_lat: The latitude of the point of interest in units of
                       decimal degrees
    :param origin_lon: The longitude of the point of interest in units of
                       decimal degrees
    :param origin_alt: The altitude of the point of interest in units of decimal
                       degrees
    :param E: The x location of the craft in ENU coordinates
    :param N: The y location of the craft in ENU coordinates
    :param U: The z position of the craft in ENU coordinats
    :param dat: A Datum dataclasse, defaulted to WGS84
    :return lat, lon, alt: The latitude, longitude, and altitude position
                           in LLH coordinats

    Code Example

    .. code-block::

        from geo_py.frames import enu_to_ned

        radar_lat = 46.017
        radar_lon = 7.750
        radar_alt = 1673.0
        E = -7134.757
        N = -4556.321
        U = 2852.390
        x, y, z = enu_to_ned(radar_lat, radar_lon, radar_alt, E, N, U)
        print(x, y, z)
        >>> -4556.321, -7134.757, -2852.39
    """
    x, y, z = enu_to_ecef(origin_lat, origin_lon, origin_alt, E, N, U, dat=dat)
    return ecef_to_ned(origin_lat, origin_lon, origin_alt, x, y, z)
# ================================================================================
# ================================================================================
# NED FUNCTIONS


def ned_to_ecef(origin_lat: float, origin_lon: float, origin_alt: float,
                N: float, E: float, D: float,
                dat: Datum = WGS84()) -> Tuple[float, float, float]:
    """
    :param origin_lat: The latitude of the point of interest in units of
                       decimal degrees
    :param origin_lon: The longitude of the point of interest in units of
                       decimal degrees
    :param origin_alt: The altitude of the point of interest in units of decimal
                       degrees
    :param N: The x location of the craft in NED coordinates
    :param E: The y location of the craft in NED coordinates
    :param D: The z position of the craft in NED coordinats
    :param dat: A Datum dataclasse, defaulted to WGS84
    :return x, y, z: The x, y, and z position in the in ECEF coordinats

    This function converts NED (North, East, Down) coordinates into ECEF
    (Earth Centered Earth Fixed) coordinates by solving the following equations
    where :math:`X_r`, :math:`Y_r`, :math:`Z_r` represent the ECEF coordinates of
    the origin, :math:`\\lambda`, :math:`\\phi`, and :math:`h` represents the
    latitude, longitude and height of the origin, and :math:`N`, :math:`E`,
    and :math:`D` represents the NED coordinats.

    .. math::

        \\begin{bmatrix}
            X \\\\
            Y \\\\
            Z \\\\
        \\end{bmatrix}
        =
        \\begin{bmatrix}
            -sin\\lambda\\:cos\\phi & -sin\\lambda & -cos\\lambda\\:\\cos\\phi\\\\
            -sin\\lambda\\:\\sin\\phi & cos\\phi & -cos\\lambda\\:sin\\phi\\\\
            cos\\lambda & 0 & -sin\\lambda
        \\end{bmatrix}
        \\cdot
        \\begin{bmatrix}
            N \\\\
            E \\\\
            D \\\\
        \\end{bmatrix}
        +
        \\begin{bmatrix}
            X_r \\\\
            Y_r \\\\
            Z_r
        \\end{bmatrix}

    Code Example

    .. code-block::

        from geo_py.datum import ITRF
        from geo_py.frames import ned_to_ecef

        radar_lat = 46.017
        radar_lon = 7.750
        radar_alt = 1673.0

        N = 28882,283
        E = -3552.574
        D = 6372030.823

        x, y, z = ned_to_ecef(radar_lat, radar_lon, radar_alt, N, E, D)
        print(x, y, z)
        >>> -7135.757, -4556.321, 2852.390

        # Use another datum
        x, y, z = ned_to_ecef(radar_lat, radar_lon, radar_alt, N, E, D, ITRF())
        >>> -7135.038, -4556.358, 2852.081
    """
    # Convert latitude, longitude to radians
    lat_o = np.deg2rad(origin_lat)
    lon_o = np.deg2rad(origin_lon)

    x_o, y_o, z_o = llh_to_ecef(origin_lat, origin_lon, origin_alt, dat=dat)

    # Compute the rotation matrix from NED to ECEF
    R_n2e = np.array([[-sin(lat_o)*cos(lon_o), -sin(lon_o), -cos(lat_o)*cos(lon_o)],
                      [-sin(lat_o)*sin(lon_o),  cos(lon_o), -cos(lat_o)*sin(lon_o)],
                      [ cos(lat_o),            0,           -sin(lat_o)           ]])

    # Rotate the NED vector to ECEF frame
    offset_ecef = np.dot(R_n2e, np.array([N, E, D]))

    # Compute the ECEF coordinates of the point
    x = x_o + offset_ecef[0]
    y = y_o + offset_ecef[1]
    z = z_o + offset_ecef[2]

    return x, y, z
# --------------------------------------------------------------------------------


def ned_to_llh(origin_lat: float, origin_lon: float, origin_alt: float,
               N: float, E: float, D: float,
               dat: Datum = WGS84()) -> Tuple[float, float, float]:
    """
    :param origin_lat: The origin latitude in units of decimal degrres
    :param origin_lat: The origin longitude in units of decimal degrees
    :param origin_alt: The origin altitude in units of meters
    :param N: The northing direction of the NED frame
    :param E: The easting direction of the NED frame
    :param D: The down direction of the NED frame
    :param dat: A Datum dataclasse, defaulted to WGS84
    :return lat, lon, alt: The latitude, longitude, and altitude of the
                           point of interest

    This function translates a NED (North, East, Down) coordinate frame
    to a LLH (Latitude, Longitude, Height) frame.

    Code Example

    .. code-block::

        from geo_py.frames import ned_to_llh

        radar_lat = 46.017
        radar_lon = 7.750
        radar_alt = 1673.0

        N = 28882,283
        E = -3552.574
        D = 6372030.823

        lat, lon, alt = ned_to_llh(radar_lat, radar_lon, radar_alt, N, E, D)
        print(lat, lon, alt)
        >>> 45.976, -7.658, 4531.0
    """
    E, N, U = ned_to_ecef(origin_lat, origin_lon, origin_alt,
                          N, E, D, dat=dat)
    return ecef_to_llh(E, N, U, dat=dat)
# --------------------------------------------------------------------------------


def ned_to_enu(origin_lat: float, origin_lon: float, origin_alt,
               N: float, E: float, D: float,
               dat: Datum = WGS84()) -> Tuple[float, float, float]:
    """
    :param origin_lat: The origin latitude in units of decimal degrres
    :param origin_lat: The origin longitude in units of decimal degrees
    :param origin_alt: The origin altitude in units of meters
    :param N: The northing direction of the NED frame
    :param E: The easting direction of the NED frame
    :param D: The down direction of the NED frame
    :param dat: A Datum dataclasse, defaulted to WGS84
    :return lat, lon, alt: The latitude, longitude, and altitude of the
                           point of interest

    This function translates a NED (North, East, Down) coordinate frame
    to a LLH (Latitude, Longitude, Height) frame.  WARNING: The values
    of origin_lat, origin_lon, and origin_alt should be the same ones
    used to transform to ENU coordinates

    Code Example

    .. code-block::

        from geo_py.frames import ned_to_enu

        radar_lat = 46.017
        radar_lon = 7.750
        radar_alt = 1673.0
        N = 28882,283
        E = -3552.574
        D = 6372030.823

        lat, lon, alt = ned_to_enu(radar_lat, radar_lon, radar_alt, N, E, D)
        print(lat, lon, alt)
        >>> -7134.757, -4556.321, 2852.390
    """
    x, y, z = ned_to_ecef(origin_lat, origin_lon, origin_alt, N, E, D, dat=dat)
    return ecef_to_enu(origin_lat, origin_lon, origin_alt, x, y, z, dat=dat)
# ================================================================================
# ================================================================================


def ned_vector(lat: float, lon: float, alt: float,
               dat: Datum = WGS84()) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    :param lat: The latitude in units of decimal degrees
    :param lon: The longitude in units of decimal degrees
    :param alt: The altitude in units of meters
    :param dat: A datum data class.  Defaulted to WGS84
    :return NED: A tuple containing the North, East, and Down vectors

    This function converts a latitude, longitude and altitude to
    a North, East, and Down vectors via the following method

    .. math::

        \\begin{align}
            N = \\frac{a}{\\sqrt[]{1-e^2sin^2\\left(\\phi\\right)}} \\\\
            x = \\left(N+h\\right)cos\\left(\\phi\\right)\\:cos\\left(\\theta\\right) \\\\
            y = \\left(N+h\\right)cos\\left(\\phi\\right)\\:sin\\left(\\theta\\right) \\\\
            z = sin\\left(\\phi\\right)\\left[sin\\left(\\phi\\right)\\left(1-e^2\\right)N+h \\right]
        \\end{align} \\\\

    .. math::

        \\begin{bmatrix}
            N \\\\
            E \\\\
            D
        \\end{bmatrix}
        =
        \\begin{bmatrix}
            -x & -y & -z \\\\
            -y & x & 0 \\\\
            -z\\:cos\\left(\\phi\\right) & -z\\:cos\\left(\\phi\\right)\\:sin\\left(\\theta\\right) &
            \\left(1-e^2\\right)z\\:cos\\left(\\phi\\right)\\:cos\\left(\\theta\\right)
        \\end{bmatrix}

    Code Example

   .. code-block::

       from geo_py.frames import ned_vector

       lat = 45.976
       lon = 7.658
       alt = 4531.0

       N, E, D = ned_vector(lat, lon, alt)
       print(N)
       print(E)
       print(D)
       >>> -4403757.60452593 -592124.57913994 -4566652.06017423
       >>> -592124.57913994 4403757.6045293 0.0
       >>> -3282645.49857941 -472918.22155005 3124277.54183501
    """
    lat = radians(lat)
    lon = radians(lon)
    N = dat.a / sqrt(1 - dat.e**2 * sin(lat)**2)
    x = (N + alt) * cos(lat) * cos(lon)
    y = (N + alt) * cos(lat) * sin(lon)
    z = ((1 - dat.e**2) * N + alt) * sin(lat)
    north = np.array([-x, -y, -z])
    east = np.array([-y, x, 0])
    down = np.array([-z * sin(lat), -z * cos(lat) * sin(lon),
                     (1 - dat.e**2) * z * cos(lat) * cos(lon)])
    return north, east, down
# --------------------------------------------------------------------------------


def enu_vector(lat: float, lon: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    :param lat: The latitude in units of decimal degrees
    :param lon: The longitude in units of decimal degrees
    :param alt: The altitude in units of meters
    :return ENU: A tuple containing the East, North, and Up vectors

    WARNING: Function still in development, do noy use

    The ENU vector is calculated via the following method where
    :math:`\\phi` and :math:`\\theta` represent latitude and longitude

    .. math::

        \\begin{bmatrix}
           E \\\\
           N \\\\
           U
        \\end{bmatrix}
        =
        \\begin{bmatrix}
            -sin\\left(\\theta\\right) & cos\\left(\\theta\\right) & 0 \\\\
            -sin\\left(\\phi\\right)\\:cos\\left(\\theta\\right) &
            -sin\\left(\\phi\\right)\\:sin\\left(\\theta\\right) &
            cos\\left(\\phi\\right) \\\\
            cos\\left(\\phi\\right)\\:cos\\left(\\theta\\right)&
            cos\\left(\\phi\\right)\\:sin\\left(\\theta\\right) &
            sin\\left(\\phi\\right)
        \\end{bmatrix}
    """
    lat = radians(lat)
    lon = radians(lon)
    east = np.array([-sin(lon), cos(lon), 0])
    north = np.array([-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat)])
    up = np.array([cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)])
    return east, north, up
# --------------------------------------------------------------------------------


# def body_to_local(lat: float, lon: float, alt: float, pitch: float,
#                   roll: float, yaw: float, cos_x: float, cos_y: float,
#                   cos_z: float, ned: bool = True,
#                   dat: Datum = WGS84()) -> np.ndarray:

#     lat = radians(lat)
#     lon = radians(lon)
#     R_body_local = aircraft_rotation(yaw, pitch, roll)

#     if ned:
#         # Rotation matrix from NED to ECEF
#         R_ecef_local = np.array([[-sin(lat)*cos(lon), -sin(lon), -cos(lat)*cos(lon)],
#                                [-sin(lat)*sin(lon), cos(lon), -cos(lat)*sin(lon)],
#                                [cos(lat), 0, -sin(lat)]])
#     else:
#         # Rotation matrix from ENU to ECEF
#         R_ecef_local = np.array([[-sin(lon), cos(lon), 0],
#                                [-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat)],
#                                [cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)]])

#     R_ecef_body = np.matmul(R_ecef_local, R_body_local)

#     # Position vector in ECEF
#     pos_ecef = llh_to_ecef(degrees(lat), degrees(lon), alt)

#     # Vector in ECEF
#     vec_ecef = np.matmul(R_ecef_body, np.array([cos_x, cos_y, cos_z])) + pos_ecef
#     return vec_ecef
# ================================================================================
# ================================================================================
# eof
