*********
rotations
*********

The **rotations.py** file contains methods that allow a programmer
to conduct coordinate frame rotations with dirction cosines as
well as quaternians.  Presently this library is only set up
to handle Tait-Bryan angles.

Euler Angles
============
The following functions will conduct passive rotations based on Direction Cosines.

.. autofunction:: geo_py.rotations.rotation_matrix

.. autofunction:: geo_py.rotations.dcm_euler_angles

.. autofunction:: geo_py.rotations.aircraft_rotation


