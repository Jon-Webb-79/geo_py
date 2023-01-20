********
distance
********

The **distance.py** file contains methods and classes that can be used to calculate
the distance between two points.  In addition, this file contains methods that will
determine the bearing between two points.  Finally this file contains classes
and methods that will determine the final location after traversing a specific
heading, a specific distance, across a known distance.  The methods in this file reley on the
us eod **Datum** dataclasses that are wither created by the user or are selected
from pre-existing dataclasses in the ``geo_py.datum`` file.

Distance Functions
==========================

The following functions can be used to determine the distance between two points.

.. autofunction:: geo_py.distance.haversine

.. autofunction:: geo_py.distance.great_circle

.. autofunction:: geo_py.distance.vincenty

.. autofunction:: geo_py.distance.linear_dist
