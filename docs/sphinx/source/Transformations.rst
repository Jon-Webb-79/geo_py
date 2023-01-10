*********
transform
*********

The **transform** file contains methods that can be used to transform back and forth
between geodetic and Cartesian coordinates.  The methods in this file reley on the
us eod **Datum** dataclasses that are wither created by the user or are selected
from pre-existing dataclasses in the ``geo_py.datum`` file.

Transformations
===============

The ``Transformations`` class contains severla methods shown below.

.. autoclass:: geo_py.transform.Transformations
   :members:
