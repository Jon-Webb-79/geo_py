******
frames
******

The **frames** file contains methods that can be used to transform back and forth
between geodetic and Cartesian coordinates.  The methods in this file reley on the
us eod **Datum** dataclasses that are wither created by the user or are selected
from pre-existing dataclasses in the ``geo_py.datum`` file.

Coordinate Transformations
==========================

The folloing functions can be used to transform back and forth between
LLH (Latitude, Longitude, Height), ECEF (Earth Centered Earth Fixed),
ENU (East, North, Up), and NED (North, East, Down) coordinate frames.

.. autofunction:: geo_py.frames.llh_to_ecef

.. autofunction:: geo_py.frames.llh_to_enu

.. autofunction:: geo_py.frames.llh_to_ned

.. autofunction:: geo_py.frames.ecef_to_llh

.. autofunction:: geo_py.frames.ecef_to_enu

.. autofunction:: geo_py.frames.ecef_to_ned

.. autofunction:: geo_py.frames.enu_to_ecef

.. autofunction:: geo_py.frames.enu_to_llh

.. autofunction:: geo_py.frames.enu_to_ned

.. autofunction:: geo_py.frames.ned_to_ecef

.. autofunction:: geo_py.frames.ned_to_llh

.. autofunction:: geo_py.frames.ned_to_enu

.. autofunction:: geo_py.frames.ned_vector

