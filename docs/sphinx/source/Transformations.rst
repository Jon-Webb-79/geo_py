*********
transform
*********

The **transform** file contains methods that can be used to transform back and forth
between geodetic and Cartesian coordinates.  The methods in this file reley on the
us eod **Datum** dataclasses that are wither created by the user or are selected
from pre-existing dataclasses in the ``geo_py.datum`` file.

Coordinate Transformations
==========================

This file contains several classes and functions that are used to transform back
and forth between different geodetic coordinate reference frames.

The ``CFTrans`` class contains several coordinate transformation methods that
are useful in the aviation and maratime professions

.. autoclass:: geo_py.transform.CFTrans
   :members:
