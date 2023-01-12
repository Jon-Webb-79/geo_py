******
frames
******

The **frames** file contains methods that can be used to transform back and forth
between geodetic and Cartesian coordinates.  The methods in this file reley on the
us eod **Datum** dataclasses that are wither created by the user or are selected
from pre-existing dataclasses in the ``geo_py.datum`` file.

Coordinate Transformations
==========================

This file contains several classes and functions that are used to transform back
and forth between different geodetic coordinate reference frames.

llh_to_ecef()
*************
This function will transform a LLH (Latitude, Longitude, Height) coordinate frame
to an ECEF (Earth Centered Earth Fixed) frame.

.. autofunction:: geo_py.frames.llh_to_ecef

llh_to_enu()
*************
This function will transform a LLH (Latitude, Longitude, Height) coordinate frame
to an ENU (East, North, Up) frame.

.. autofunction:: geo_py.frames.llh_to_enu

llh_to_ned()
*************
This function will transform a LLH (Latitude, Longitude, Height) coordinate frame
to an NED (North, East, Down) frame.

.. autofunction:: geo_py.frames.llh_to_ned

ecef_to_llh()
*************
This function will transform a ECEF (Earth Centered Earth Fixed) coordinate frame
to a LLH (Latitude, Longitude, Height) frame.

.. autofunction:: geo_py.frames.ecef_to_llh

ecef_to_enu()
*************
This function will transform an ECEF (Earth Centered Earth Fixed) coordinate frame
to an ENU (East, North, Up) frame.

.. autofunction:: geo_py.frames.ecef_to_enu

ecef_to_ned()
*************
This function will transform an ECEF (Earth Centered Earth Fixed) coordinate frame
to an NED (North, East, Down) frame.

.. autofunction:: geo_py.frames.ecef_to_ned

enu_to_ecef()
*************
This function will transform an ENU (East, North, Up) coordinate frame
to an ECEF (Earth Centered Earth Fixed) frame.

.. autofunction:: geo_py.frames.enu_to_ecef

enu_to_llh()
*************
This function will transform an ENU (East, North, Up) coordinate frame
to an LLH (Latitude, Longitude, Height) frame.

.. autofunction:: geo_py.frames.enu_to_llh

enu_to_ned()
*************
This function will transform an ENU (East, North, Up) coordinate frame
to an NED (North, East, Down) frame.

.. autofunction:: geo_py.frames.enu_to_ned

ned_to_ecef()
*************
This function will transform an NED (North, East, Down) coordinate frame
to an ECEF (Earth Centered Earth Fixed) frame.

.. autofunction:: geo_py.frames.ned_to_ecef

ned_to_llh()
*************
This function will transform an NED (North, East, Down) coordinate frame
to a LLH (Latitude, Longitude, Height) frame.

.. autofunction:: geo_py.frames.ned_to_llh

ned_to_enu()
*************
This function will transform an NED (North, East, Down) coordinate frame
to a ENU (East, North, UP) frame.

.. autofunction:: geo_py.frames.ned_to_enu
