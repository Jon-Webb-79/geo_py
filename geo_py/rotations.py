# Import necessary packages here
from math import sin, cos, sqrt, radians, degrees, pi
import numpy as np
from typing import Tuple
# ================================================================================
# ================================================================================
# File:    rotations.py
# Date:    January 14, 2023
# Purpose: Describe the purpose of functions of this file

# Source Code Metadata
__author__ = "Jonathan A. Webb"
__copyright__ = "Copyright 2023, Jon Webb Inc."
__version__ = "1.0"
# ================================================================================
# ================================================================================
# DIRECTION COSINE FUNCTIONS


def rotation_matrix(x_rot: float, y_rot: float, z_rot: float,
                    order: str = "ZYX", deg: bool = False,
                    extrinsic: bool = False) -> np.ndarray:
    """
    :param x_rot: The rotation about the x-axis in units of degrees or radians.
                  This axis is also known as the yaw axis.
    :param y_rot: The rotation about the y-axis in units of degrees or radians.
                  This axis is also known as the pitch axis.
    :param z_rot: The rotation about the z-axis in units of degrees or radians.
                  This axis is also known as the roll axis.
    :param order: The order of rotation as a string.  Defaulted to 'ZYX', which
                  is the traditional order for aircraft rotations.
                  Currently this function only support 'XYZ', 'XZY', 'YXZ',
                  'YZX', 'ZXY', and 'ZYX'.
    :param deg: True if the angles x_rot, y_rot, and z_rot are in units of
                degrees, False for radians.  Defaulted to False
    :param extrinsic: True if the direction cosine matrix is to be calculated as
                      an extrinsic rotation matrix, False for an intrinsic
                      matrix.  Defaulted to False.

    This function will return a direction cosine matrix according to the Euler
    angles fed to the function.  This function solves and multiplies the
    following rotation matrices where :math:`\\theta`, :math:`\\phi`,
    and :math:`\\psi` corresponse to roll, pitch, and yaw rotations.

    .. math::

       \\begin{align}
           R_x =
           \\begin{bmatrix}
               1 & 0 & 0 \\\\
               0 & cos\\left(\\theta\\right) & sin\\left(\\theta\\right) \\\\
               0 & -sin\\left(\\theta\\right) & cos\\left(\\theta\\right) \\\\
           \\end{bmatrix} \\\\
           \\\\
           R_y =
           \\begin{bmatrix}
               cos\\left(\\phi\\right) & 0 & -sin\\left(\\phi\\right) \\\\
               0 & 1 & 0 \\\\
               sin\\left(\\phi\\right) & 0 & cos\\left(\\phi\\right)
           \\end{bmatrix} \\\\
           \\\\
           R_z =
           \\begin{bmatrix}
               cos\\left(\\psi\\right) & sin\\left(\\psi\\right) & 0 \\\\
               -sin\\left(\\psi\\right) & cos\\left(\\psi\\right) & 0 \\\\
               0 & 0 & 1 \\\\
           \\end{bmatrix}
       \\end{align}

    As an example, if you were conducting a rotation about 'XZY', the multiplication
    would occur as such to obtain the intrinsic matrix

    .. math::

        R_i = R_x\\left(\\theta\\right)R_z\\left(\\psi\\right)R_y\\left(\\phi\\right)

    The extrinsic matrix is determined by taking the transpose of the intrinstic
    matrix.

    .. math::

        R_t = R_i^T

    Code examples

    .. code-block::

        from geo_py.rotations import rotation_matrix

        pitch = 0.83  # radians
        roll = 0.21  # radians
        yaw = 1.18  # radians

        # rotate about the ZYX axis, wich is default to get intrinsic matrix
        dcm = rotation_matrix(yaw, pitch, roll)
        print(dcm)
        >>> [[ 0.25707693, -0.62399419, 0.73793137],
             [0.96289059, 0.23032494, -0.14068453],
             [-0.08217767, 0.74671392, 0.66004936]]

    .. code-block::

        from geo_py.rotations import rotation_matrix

        # Find the extrinsic matrix with inputs in degrees, rotated against xyz
        pitch = 22.1  # degrees
        roll = 8.75  # degrees
        yaw = 12.8  # degrees
        dcm = rotation_matrix(ya, pitch, roll, order="XYZ", deg=True, extrinsic=True)
        print(dcm)
        >>> [[ 0.90335038, 0.27478025, -0.32890227],
             [-0.20527103, 0.95112031, 0.23072485],
             [0.37622426, -0.14094667, 0.9157424]]
    """
    order = order.upper()
    if deg:
        x_rot = radians(x_rot)
        y_rot = radians(y_rot)
        z_rot = radians(z_rot)

    # Create rotation matrices for each axis
    rot_x = np.array([[1, 0, 0],
                     [0, np.cos(z_rot), -np.sin(z_rot)],
                     [0, np.sin(z_rot), np.cos(z_rot)]])
    rot_y = np.array([[np.cos(y_rot), 0, np.sin(y_rot)],
                     [0, 1, 0],
                     [-np.sin(y_rot), 0, np.cos(y_rot)]])
    rot_z = np.array([[np.cos(x_rot), -np.sin(x_rot), 0],
                     [np.sin(x_rot), np.cos(x_rot), 0],
                     [0, 0, 1]])

    # Multiply rotation matrices according to the order of rotations
    if order == "XYZ":
        rot_final = np.matmul(np.matmul(rot_x, rot_y), rot_z)
    elif order == "XZY":
        rot_final = np.matmul(np.matmul(rot_x, rot_z), rot_y)
    elif order == "YXZ":
        rot_final = np.matmul(np.matmul(rot_y, rot_x), rot_z)
    elif order == "YZX":
        rot_final = np.matmul(np.matmul(rot_y, rot_z), rot_x)
    elif order == "ZXY":
        rot_final = np.matmul(np.matmul(rot_z, rot_x), rot_y)
    elif order == "ZYX":
        rot_final = np.matmul(np.matmul(rot_z, rot_y), rot_x)
    else:
        raise ValueError("Invalid rotation order")
    if extrinsic:
        rot_final = rot_final.T
    return rot_final
# --------------------------------------------------------------------------------


def dcm_euler_angles(dcm: np.ndarray, order: str = "ZYX", deg: bool = False,
                     extrinsic: bool = False) -> Tuple[float, float, float]:
    """
    :param dcm: A direction cosine matrix of order 3x3
    :param order: The order of rotation as a string.  Defaulted to 'ZYX', which
                  is the traditional order for aircraft rotations.
                  Currently this function only support 'XYZ', 'XZY', 'YXZ',
                  'YZX', 'ZXY', and 'ZYX'.
    :param deg: True if the angles x_rot, y_rot, and z_rot are in units of
                degrees, False for radians.  Defaulted to False
    :param extrinsic: True if the direction cosine matrix is to be calculated as
                      an extrinsic rotation matrix, False for an intrinsic
                      matrix.  Defaulted to False.

    This function will determine hte Euler angles associated with a direction
    cosine matrix assuming knowledge of the rotation order.  Assuming,
    a rotation in the 'XYZ' orientations, the Euler angles can be determined
    from the following formula where :math:`R_{x,y}` reprsents an attribute
    of the rotation matrix from postion :math:`x`, and :math:`y`.

    .. math::

        \\begin{align}
            \\alpha=arctan2\\left(\\frac{R_{32}}{R_{22}}\\right) \\\\
            \\beta=arcsin\\left(-R_{12}\\right) \\\\
            \\gamma=arctan2\\left(\\frac{R_{13}}{R_{11}}\\right)
        \\end{align}

    However, the exact nature of the equations change as a function of
    the orientation angle.  A further description of the applicable equations
    can be found here `rotations <https://en.wikipedia.org/wiki/Euler_angles>`_.

    Code Example

    .. code-block::

        from geo_py.rotations import rotation_matrix, dcm_euler_angles

        # Find the extrinsic matrix with inputs in degrees, rotated against xyz
        pitch = 22.1  # degrees
        roll = 8.75  # degrees
        yaw = 12.8  # degrees
        dcm = rotation_matrix(ya, pitch, roll, order="XYZ", deg=True,
                              extrinsic=False)
        yaw, pitch, roll = dcm_euler_angles(dcm, order="XYZ", deg=True)
        print(yaw, pitch, roll)
        >>> 12.8, 22.1000,8.75
    """
    order = order.upper()
    if extrinsic:
        dcm = dcm.T
    if order == 'ZYX':
        yaw = np.arctan2(dcm[1, 0], dcm[0, 0])
        pitch = np.arcsin(-dcm[2, 0])
        roll = np.arctan2(dcm[2, 1], dcm[2, 2])
    elif order == 'XYZ':
        roll = np.arctan2(-dcm[1, 2], dcm[2, 2])
        pitch = np.arcsin(dcm[0, 2])
        yaw = np.arctan2(-dcm[0, 1], dcm[0, 0])
    elif order == 'XZY':
        roll = np.arctan2(dcm[2, 1], dcm[1, 1])
        pitch = np.arctan2(dcm[0, 2], dcm[0, 0])
        yaw = np.arcsin(-dcm[0, 1])
    elif order == 'YXZ':
        roll = np.arcsin(-dcm[1, 2])
        pitch = np.arctan2(dcm[0, 2], dcm[2, 2])
        yaw = np.arctan2(dcm[1, 0], dcm[1, 1])
    elif order == 'YZX':
        roll = np.arctan2(-dcm[1, 2], dcm[1, 1])
        pitch = np.arctan2(-dcm[2, 0], dcm[0, 0])
        yaw = np.arcsin(dcm[1, 0])
    elif order == 'ZXY':
        yaw = np.arctan2(-dcm[0, 1], dcm[1, 1])
        roll = np.arcsin(dcm[2, 1])
        pitch = np.arctan2(-dcm[2, 0], dcm[2, 2])
    else:
        raise ValueError("Invalid order of rotation.")
    if deg:
        yaw = degrees(yaw)
        pitch = degrees(pitch)
        roll = degrees(roll)
    return yaw, pitch, roll
# ================================================================================
# ================================================================================
# eof
