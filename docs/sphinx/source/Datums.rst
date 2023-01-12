******
Datums
******

The datums sub-library contains several geodetic datums that were collected as part
of the World Geodetic Survey (WGS).  The Datums were drafted in the ``datum.py`` file
as a **Protocol** dataclass of the following format.

.. autoclass:: geo_py.datum.Datum
   :members:

The **Protocall** class titled ``Datum`` is a template that provides a guide describing how
other Classes of a similar type should be structured.  The protocal acts as an
implicit interface for interface driven development with the ``typing`` module.
This enables the developer of this software suite to implement pre-made classes of the correct type, but
also enables users to develop their own similar classes.  Currently there
are pre-programmed versions of this class titled ``WGS77``, ``WGS80``, and
``WGS84``.  For a refresher, the user is encouraged to 
review `dataclasses <https://docs.python.org/3/library/dataclasses.html>`_ for the
implementation of dataclasses in Python.  In addition the user can also
review `Protocal <https://www.pythontutorial.net/python-oop/python-protocol/>`_ for a
refresher on how **Protocol** classes should be implimented. For instance, if a user
wishes to set up a geodetic datum for Saturn, they would build a dataclass with
the following format. **Warning** These parameters are made up and are only for
the purpose of demonstrating how to instantiate a Datum Prototype class.  All
attributes listed in this example, must be in the user implementation, with
the same units.  While the user defined dataclass must at a minimum contain
a implementation of the attributes shown below, the user defined class
can contain extra attributes and methods.

.. code-block::
   :caption: User Defined Datum

   from geo_py.datum import Datum
   from dataclasses import dataclasses

   @dataclass(slots=True)
   class SaturnDatum(Datum):
       a: float = 60268000.0  # meters
       b: float = 54363000.0  # meters
       ro: float = 58232000.0  # meters
       gm: float = 37931000.0  # m3/s-2
       f: float = 0.09796
       e: float = 0.43169221
       e2: float = 0.18635
       e24: float = 0.18621
       g: float = 10.44  # m/s2

Pre Existing Datums
===================

While the user can define their own datums as a dataclass, this library comes with the following
pre-developed Datum dataclasses

WGS77
*****
The World Geodetic Survey 1977 datum results are implement in a dataclass titled ``WGS77``.
The datum can be retrieved with the following command. **Provide Reference**

.. code-block::
   :caption: User Defined Datum

   from geopy.datum import WGS77

WGS80
*****
The World Geodetic Survey 1980 datum results are implement in a dataclass titled ``WGS80``.
The datum can be retrieved with the following command. `WGS80 <https://geodesy.geology.ohio-state.edu/course/refpapers/00740128.pdf>`_

.. code-block::
   :caption: User Defined Datum

   from geo_py.datum import WGS80

WGS84
*****
The World Geodetic Survey 1984 datum results are implement in a dataclass titled ``WGS84``.
The datum can be retrieved with the following command. `WGS84 <https://apps.dtic.mil/sti/pdfs/ADA280358.pdf>`_

.. code-block::
   :caption: User Defined Datum

   from geo_py.datum import WGS84

NAD83
*****
The North American Datum 1983 datum results are implement in a dataclass titled ``NAD83``.
The datum can be retrieved with the following command.

.. code-block::
   :caption: User Defined Datum

   from geo_py.datum import NAD83

ETRS89
******
The European 1989 datum results are implement in a dataclass titled ``ETRS89``.
The datum can be retrieved with the following command.

.. code-block::
   :caption: User Defined Datum

   from geo_py.datum import ETRS89

ITRF
****
The International Terrestrial Reference Frame datum results are implement in a dataclass titled ``ITRF``.
The datum can be retrieved with the following command.

.. code-block::
   :caption: User Defined Datum

   from geo_py.datum import WGS84

   
