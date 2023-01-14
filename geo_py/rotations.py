# Import necessary packages here
from math import sin, cos, sqrt, radians
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


def extrinsic_dir_cos_mat(pitch: float, roll: float, yaw: float,
                          deg: bool = False) -> np.ndarray:
    """
    :param pitch: The angle of attack of the craft relative to the x-axis in
                  units of radians or degrees.
    :param roll: The angle of the craft relative to the y-axis in units of
                 radians or degrees.
    :param yaw: The angle of the craft relative to the z-axis, also
                the same as heading, in units of radians or degrees.
    :param deg: True if angles are in units od degrees, False is in radians.
                Defaulted to radians or degrees.
    :return dcm: The direction cosine matrix specific for rotations
                 around the Y, X, and Z axis, or the pitch, roll,
                 and yaw axis.

    This function will build an extrinsic direction cosine matrix.
    This function conducts an extrsinsic Euler rotation about all axis
    This function returns the
    matrix of the following format where :math:`\\theta`, :math:`\\phi`,
    and :math:`\\psi` represent pitch, roll, and yaw respectively.
    An extrinsic rotation matrix rotates all axis about a fixed axis,
    as might be the case when with an aircraft that is referenced to the
    earth.

    .. math::

        R_{xyz}=
        \\begin{bmatrix}
            cos\\theta\\:cos\\psi & cos\\theta\\:sin\\psi & -sin\\theta \\\\
            sin\\phi\\:sin\\theta\\:cos\\psi & sin\\phi\\:sin\\theta\\:sin\\psi & cos\\theta\\:sin\\phi \\\\
            cos\\phi\\:sin\\theta\\:cos\\psi+sin\\phi\\:sin\\psi &
            cos\\phi\\:sin\\theta\\:sin\\psi-sin\\phi\\:cos\\psi & cos\\theta\\:cos\\phi \\\\
        \\end{bmatrix}

    Code Example

    .. code-block::

        from geopy.rotations import extrinsic_dir_cos_mat
        pitch = 0.1
        roll = 0.0
        yaw = 0.7854
        dcm = extrinsic_dir_cos_mat(pitch, roll, yaw)
        print(dcm)
        >>> [[ 0.7035729 0.70357548 -0.09983342]
             [-0.70710808 0.70710548 0.]
             [0.07059276 0.07059302 0.9950417]]
    """
    if deg:
        pitch = radians(pitch)
        roll = radians(roll)
        yaw = radians(yaw)
    x_axis = [cos(pitch)*cos(yaw), cos(pitch)*sin(yaw), -sin(pitch)]
    y_axis = [sin(roll)*sin(pitch)*cos(yaw)-cos(roll)*sin(yaw),
         sin(roll)*sin(pitch)*sin(yaw)+cos(roll)*cos(yaw),
         cos(pitch)*sin(roll)]
    z_axis = [cos(roll)*sin(pitch)*cos(yaw)+sin(roll)*sin(yaw),
         cos(roll)*sin(pitch)*sin(yaw)-sin(roll)*cos(yaw),
         cos(pitch)*cos(roll)]
    return np.array([x_axis, y_axis, z_axis])
# --------------------------------------------------------------------------------


def intrinsic_dir_cos_mat(alpha: float, beta: float,
                          gamma: float, order: str = "XYZ",
                          deg: bool = False):
    """
    :param alpha: The rotation angle for the x rotation matrix in units of
                  degrees or radians
    :param beta: The rotation angle for the y rotation matrix in units of
                 degrees or radians
    :param gamma: The rotation angle for the z rotation matrix in units of
                  degrees or raidans
    :param order: The order of rotation, such as 'XYZ', 'XZY',
                  'YXZ', 'YZX', 'ZXY', or 'ZYX'. Defaulted to
                  'XYZ'
    :param deg: True if angles are in degrees, False, otherwise.  Defaulted
                to False.

    This function will produce an intrinsic direction cosine
    matrix in the order specified.  An intrinsic rotation matrix
    conducts each rotation along the axis produced by the
    previous rotation, as might be the case for a robotic
    arm.  This function returns
    the following matrices in the user specified order
    of operations.

    .. math::

       \\begin{align}
           R_x =
           \\begin{bmatrix}
               1 & 0 & 0 \\\\
               0 & cos\\alpha & sin\\alpha \\\\
               0 & -sin\\alpha & cos\\alpha \\\\
           \\end{bmatrix} \\\\
           \\\\
           R_y =
           \\begin{bmatrix}
               cos\\beta & 0 & -sin\\beta \\\\
               0 & 1 & 0 \\\\
               sin\\beta & 0 & cos\\beta
           \\end{bmatrix} \\\\
           \\\\
           R_z =
           \\begin{bmatrix}
               cos\\gamma & sin\\gamma & 0 \\\\
               -sin\\gamma & cos\\gamma & 0 \\\\
               0 & 0 & 1 \\\\
           \\end{bmatrix}
       \\end{align}

    Code Example

    .. code-block::

        from geopy.rotations import intrinsic_cos_mat
        pitch = 0.1
        roll = 0.0
        yaw = 0.7854
        dcm = intrinsic_dir_cos_mas(pitch, roll, yaw, "ZYX")
        print(dcm)
        >>> [[0.70710548, -0.70357548, 0.07059302],
             [0.70710808, 0.7035729, -0.07059276],
             [0.0, 0.09983342, 0.99500417]]
    """
    if deg:
        alpha = radians(alpha)
        beta = radians(beta)
        gamma = radians(gamma)
    # Create rotation matrices for each axis
    Rx = np.array([[1, 0, 0],
                   [0, np.cos(alpha), -np.sin(alpha)],
                   [0, np.sin(alpha), np.cos(alpha)]])
    Ry = np.array([[np.cos(beta), 0, np.sin(beta)],
                   [0, 1, 0],
                   [-np.sin(beta), 0, np.cos(beta)]])
    Rz = np.array([[np.cos(gamma), -np.sin(gamma), 0],
                   [np.sin(gamma), np.cos(gamma), 0],
                   [0, 0, 1]])

    # Multiply rotation matrices according to the order of rotations
    if order == "XYZ":
        R = np.matmul(np.matmul(Rx, Ry), Rz)
    elif order == "XZY":
        R = np.matmul(np.matmul(Rx, Rz), Ry)
    elif order == "YXZ":
        R = np.matmul(np.matmul(Ry, Rx), Rz)
    elif order == "YZX":
        R = np.matmul(np.matmul(Ry, Rz), Rx)
    elif order == "ZXY":
        R = np.matmul(np.matmul(Rz, Rx), Ry)
    elif order == "ZYX":
        R = np.matmul(np.matmul(Rz, Ry), Rx)
    else:
        raise ValueError("Invalid rotation order")
    return R
# --------------------------------------------------------------------------------


def direction_cosines(vector: np.ndarray) -> Tuple[float, float, float]:
    """
    :param vector: A three dimensional vector
    :return cos_x, cos_y, cos_z: The direction cosines for a vector

    This function returns the direction cosines of a vector
    (:math:`\\alpha`, :math:`\\beta`, :math:`\\gamma`) via the
    following equation.

    .. math::

        v = ai + bj + ck \\\\
        cos\\alpha = \\frac{a}{\\sqrt[]{a^2+b^2+c^2}} \\\\
        cos\\beta = \\frac{b}{\\sqrt[]{a^2+b^2+c^2}} \\\\
        cos\\gamma = \\frac{c}{\\sqrt[]{a^2+b^2+c^2}} \\\\

    Code Example

    .. code-block::

        from geo_py.rotations import direction_cosines
        vec = np.array([1., 2., 3.])
        cos_x, cos_y, cos_z = direction_cosines(vec)
        print(cos_x, cos_y, cos_z)
        >>> 0.26726, 0.53452, 0.801783
    """
    cos_x = vector[0] / sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
    cos_y = vector[1] / sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
    cos_z = vector[2] / sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
    return cos_x, cos_y, cos_z
# --------------------------------------------------------------------------------


def dcm_to_quaternion(dcm):
    """
    :param dcm: A Direction cosine matrix of size 3x3
    :return q: A quaternion of size 1x4

    This function transforms a direction cosine matrix to a quaternion.

    Code Example

    .. code-block::

        from geo_py.rotations import dcm_to_quaternion
        rot_mat = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        q = dcm_to_quaternion(rot_mat)
        print(q)
        >>> [0.01912624, -0.04617471, -0.38220603, 0.92272457],
    """
    tr = np.trace(dcm)
    if tr > 0:
        s = np.sqrt(tr + 1.0) * 2
        q = np.array([(dcm[2, 1] - dcm[1, 2]) / s,
                      (dcm[0, 2] - dcm[2, 0]) / s,
                      (dcm[1, 0] - dcm[0, 1]) / s,
                      0.25 * s])
    elif dcm[0, 0] > dcm[1, 1] and dcm[0, 0] > dcm[2, 2]:
        s = np.sqrt(1.0 + dcm[0, 0] - dcm[1, 1] - dcm[2, 2]) * 2
        q = np.array([0.25 * s,
                      (dcm[0, 1] + dcm[1, 0]) / s,
                      (dcm[2, 0] + dcm[0, 2]) / s,
                      (dcm[2, 1] - dcm[1, 2]) / s])
    elif dcm[1, 1] > dcm[2, 2]:
        s = np.sqrt(1.0 + dcm[1, 1] - dcm[0, 0] - dcm[2, 2]) * 2
        q = np.array([(dcm[0, 1] + dcm[1, 0]) / s,
                      0.25 * s,
                      (dcm[1, 2] + dcm[2, 1]) / s,
                      (dcm[0, 2] - dcm[2, 0]) / s])
    else:
        s = np.sqrt(1.0 + dcm[2, 2] - dcm[0, 0] - dcm[1, 1]) * 2
        q = np.array([(dcm[0, 2] + dcm[2, 0]) / s,
                      (dcm[1, 2] + dcm[2, 1]) / s,
                      0.25 * s,
                      (dcm[1, 0] - dcm[0, 1]) / s])
    return q
# ================================================================================
# ================================================================================
# QUATERNIANS


def quaternion_to_dcm(q):
    """
    :param q: A quaternion of size 1x4
    :return dcm: A direction cosine matrix of size 3x3

    This function convert a quaternian to a direction cosine matrix
    with the following method where :math:`q_x` is a component to
    the quaternion and :math:`C` is the direction cosine matrix.

    .. math::

       C =
       \\begin{bmatrix}
          q_4^2+q_1^2-q2^2-q_3^2 & 2\\left(q_1q_2+q_3q_4\\right) & 2\\left(q_1q_3-q_2q_4\\right) \\\\
          2\\left(q_1q_2-q_3q_4\\right) & q_4^2-q_1^2+q_2^2-q_3^2 & 2\\left(q_2q_3+q_1q_4\\right) \\\\
          2\\left(q_1q_3+q_2q_4\\right) & 2\\left(q_2q_3-q_1q_4\\right) & q_4^2-q_1^2-q_2^2+q_3^2 \\\\
       \\end{bmatrix}

    .. code-block::

        from geo_py.rotations import quaternion_to_dcm
        q = np.array([-0.1677489, -0.7369231, -0.3682588, 0.5414703]])
        dcm = quaternion_to_dcm(q)
        print(dcm)
        >>> [[0.1423907, 0.72441893, -0.67449393],
             [0.36109474, -0.67249153, -0.64603848],
             [-0.9215396, -0.15156633, -0.35734044]]
    """
    q0, q1, q2, q3 = q[0], q[1], q[2], q[3]

    # First row of the rotation matrix
    r00 = 2 * (q0 * q0 + q1 * q1) - 1
    r01 = 2 * (q1 * q2 - q0 * q3)
    r02 = 2 * (q1 * q3 + q0 * q2)

    # Second row of the rotation matrix
    r10 = 2 * (q1 * q2 + q0 * q3)
    r11 = 2 * (q0 * q0 + q2 * q2) - 1
    r12 = 2 * (q2 * q3 - q0 * q1)

    # Third row of the rotation matrix
    r20 = 2 * (q1 * q3 - q0 * q2)
    r21 = 2 * (q2 * q3 + q0 * q1)
    r22 = 2 * (q0 * q0 + q3 * q3) - 1
    # 3x3 rotation matrix
    rot_matrix = np.array([[r00, r01, r02],
                           [r10, r11, r12],
                           [r20, r21, r22]])
    return rot_matrix
# ================================================================================
# ================================================================================
# eof
