# Import necessary packages here
import sys
import os
from pathlib import PurePath
from typing import Tuple
from math import radians, cos, sin, sqrt, atan, atan2, pi, degrees
import numpy as np
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


class CFTrans:
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

        Code example

        .. code-block::

            from geo_py.datum import WGS84
            from geo_py.transform import CFTrans

            tran = Transformations(WGS84())
            lat = 46.826
            lon = 107.321
            alt = 6096.0
            x, y, z = tran.ecef_to_llh(lat, lon, alt)
            print(x, y, z)
            >>> -1302839.38, 4177542.21, 4632996.83
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
# --------------------------------------------------------------------------------

    def llh_to_enu(self, lat: float, lon: float, alt: float,
                   ref_lat: float, ref_lon: float,
                   ref_alt: float) -> Tuple[float, float, float]:
        """
        """
       # ref is aircraft
        x, y, z = self.llh_to_ecef(lat, lon, alt)
        xr, yr, zr = self.llh_to_ecef(ref_lat, ref_lon, ref_alt)

        # convert reference latitude and longitude to radians
        ref_lat = np.deg2rad(ref_lat)
        ref_lon = np.deg2rad(ref_lon)

        # calculate transformation matrix
        R = np.array([[-np.sin(ref_lon), np.cos(ref_lon), 0],
                  [-np.sin(ref_lat)*np.cos(ref_lon),
                   -np.sin(ref_lat)*np.sin(ref_lon),
                   np.cos(ref_lat)],
                  [np.cos(ref_lat)*np.cos(ref_lon),
                   np.cos(ref_lat)*np.sin(ref_lon),
                   np.sin(ref_lat)]])

        # transform ECEF coordinates to NED coordinats
        ned = np.dot(R, np.array([x-xr, y-yr, z-zr]))
        return tuple(ned)
# --------------------------------------------------------------------------------

    def ecef_to_llh(self, x: float, y: float,
                    z: float) -> Tuple[float, float, float]:
        """
        :param x: The x-coordinate in ECEF coordinate frame
        :param y: The y-coordinate in ECEF coordinate frame
        :param z: The z-coordinate in ECEF coordinate frame
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
        :math:`h` in accordance with the following equations.

        .. math::

           \\begin{eqnarray}
           e^2 = \\frac{a^2-b^2}{a^2} \\\\
           e^{\\prime 2} = \\frac{a^2-b^2}{b^2} \\\\
           p = \\sqrt[]{X^2+Y^2} \\\\
           F = 54b^2\\:Z^2 \\\\
           G = p^2 + \\left(1 - e^2\\right)Z^2 - e^2 \\left(a^2-b^2\\right) \\\\
           c = \\frac{e^4Fp^2}{G^3} \\\\
           s = \\sqrt[3]{1 + c + \\sqrt[]{c^2 + 2c}} \\\\
           k = s + 1 + \\frac{1}{s} \\\\
           P = \\frac{F}{3k^2G^2} \\\\
           Q = \\sqrt[]{1+2e^4P} \\\\
           r_o = \\frac{-Pe^2p}{1+Q} \\\\
           r_o2 = \\frac{1}{2}a\\left(1+\\frac{1}{Q}\\right) \\\\
           r_o3 = -\\frac{P\\left(1-e^2\\right)Z^2}{Q\\left(1+Q\\right)} - \\frac{1}{2}Pp^2 \\\\
           r_o = r_o + \\sqrt[]{r_o2 + r_o3} \\\\
           U = \\sqrt[]{\\left(p-e^2r_o\\right)^2+Z^2} \\\\
           V = \\sqrt[]{\\left(p-e^2r_o\\right)^2+\\left(1-e^2\\right)Z^2} \\\\
           z_o = \\frac{b^2Z}{aV} \\\\
           h = U\\left(1-\\frac{b^2}{aV}\\right) \\\\
           \\phi = arctan\\left[\\frac{Z+e^{\\prime2}z_o}{p}\\right] \\\\
           \\lambda = arctan2\\left[Y, X\\right]
          \\end{eqnarray}

        Code example

        .. code-block::

            from geo_py.datum import WGS84
            from geo_py.transform import CFTrans

            tran = Transformations(WGS84())
            x = -1302839.38
            y = 4177542.21
            z = 4632996.83
            lat, lon, alt = tran.ecef_to_llh(x, y, z)
            print(lat, lon, alt)
            >>> 46.826, 107.321, 6097.1
        """
        # calculations:
        r = sqrt(x**2 + y**2)
        ep_sq  = (self.datum._a**2-self.datum._b**2)/self.datum._b**2
        ee = (self.datum._a**2-self.datum._b**2)
        f = (54*self.datum._b**2)*(z**2)
        g = r**2 + (1 - self.datum._e2)*(z**2) - self.datum._e2*ee*2
        c = (self.datum._e2**2)*f*r**2/(g**3)
        s = (1 + c + sqrt(c**2 + 2*c))**(1/3.)
        p = f/(3.*(g**2)*(s + (1./s) + 1)**2)
        q = sqrt(1 + 2*p*self.datum._e2**2)
        r_0 = -(p*self.datum._e2*r)/(1+q) + \
                sqrt(0.5*(self.datum._a**2)*(1+(1./q)) - \
                p*(z**2)*(1-self.datum._e2)/(q*(1+q)) - 0.5*p*(r**2))
        u = sqrt((r - self.datum._e2*r_0)**2 + z**2)
        v = sqrt((r - self.datum._e2*r_0)**2 + (1 - self.datum._e2)*z**2)
        z_0 = (self.datum._b**2)*z/(self.datum._a*v)
        h = u*(1 - self.datum._b**2/(self.datum._a*v))
        phi = atan((z + ep_sq*z_0)/r)
        lambd = atan2(y, x)
        return phi*180/pi, lambd*180/pi, h
# --------------------------------------------------------------------------------

    def ecef_to_enu(self, origin_lat: float, origin_lon: float, origin_alt: float,
                    X: float, Y: float,
                    Z: float) -> Tuple[float, float, float]:
        """
        :param origin_lat: The latitude of the point of interest in units of
                           decimal degrees
        :param origin_lon: The longitude of the point of interest in units of
                           decimal degrees
        :param origin_alt: The altitude of the point of interest in units of decimal
                           degrees
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

            from geo_py.datum import WGS84
            from geo_py.transform import CFTrans

            tran = Transformations(WGS84())
            radar_lat = 46.017
            radar_lon = 7.750
            radar_alt = 1673.0

            x = 45.976
            y = 7.658
            z = 4531.0
            E, N, U = tran.ecef_to_enu(radar_lat, radar_lon, radar_alt, x, y, z)
            print(E, N, U)
            >>> -7126.303, -4562.512, 2863.687
        """
        x_r, y_r, z_r = self.llh_to_ecef(origin_lat, origin_lon, origin_alt)

        lat = np.deg2rad(origin_lat)
        lon = np.deg2rad(origin_lon)
        E = -sin(lon)*(X-x_r) + cos(lon)*(Y-y_r);
        N = -sin(lat)*cos(lon)*(X-x_r) - sin(lat)*sin(lon)*(Y-y_r) + cos(lat)*(Z-z_r);
        U = cos(lat)*cos(lon)*(X-x_r) + cos(lat)*sin(lon)*(Y-y_r) + sin(lat)*(Z-z_r);
        return E, N, U
# --------------------------------------------------------------------------------

    def enu_to_ecef(self, ref_lat: float, ref_lon: float, ref_alt: float,
                    E: float, N: float,
                    U: float) -> Tuple[float, float, float]:
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

            from geo_py.datum import WGS84
            from geo_py.transform import CFTrans

            tran = Transformations(WGS84())
            radar_lat = 46.017
            radar_lon = 7.750
            radar_alt = 1673.0

            E = -7126.303
            N = -4562.512
            U = 2863.687

            x, y, z = tran.enu_to_ecef(radar_lat, radar_lon, radar_alt, x, y, z)
            print(x, y, z)
            >>> 4403757.60, 592124.58, 4566652.06
        """
        x_ref, y_ref, z_ref = self.llh_to_ecef(ref_lat, ref_lon, ref_alt)

        ref_lat = np.deg2rad(ref_lat)
        ref_lon = np.deg2rad(ref_lon)

        X = -sin(ref_lon)*E - cos(ref_lon)*sin(ref_lat)*N + cos(ref_lon)*cos(ref_lat)*U + x_ref;
        Y = cos(ref_lon)*E - sin(ref_lon)*sin(ref_lat)*N + cos(ref_lat)*sin(ref_lon)*U + y_ref;
        Z = cos(ref_lat)*N + sin(ref_lat)*U + z_ref;

        # r_val = np.array([[-np.sin(ref_lon), np.cos(ref_lon), 0],
        #                 [-np.sin(ref_lat)*np.cos(ref_lon),
        #                  -np.sin(ref_lat)*np.sin(ref_lon),
        #                  np.cos(ref_lat)],
        #                [np.cos(ref_lat)*np.cos(ref_lon),
        #                 np.cos(ref_lat)*np.sin(ref_lon),
        #                 np.sin(ref_lat)]])
        # enu = np.dot(r_val.T, np.array([E, N, U])) + np.array([x_ref, y_ref, z_ref])
        return X, Y, Z
# ================================================================================
# ================================================================================
# eof
