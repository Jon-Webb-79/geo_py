# Import necessary packages here
from math import isclose, radians
import os
import sys
from pathlib import PurePath
import numpy as np
p = PurePath(__file__).parent
sys.path.insert(1, os.path.abspath(p))
from geo_py.rotations import intrinsic_dir_cos_mat, extrinsic_dir_cos_mat
from geo_py.rotations import direction_cosines, dcm_to_quaternion
from geo_py.rotations import quaternion_to_dcm
# ================================================================================
# ================================================================================
# File:    test_rotations.py
# Date:    January 14, 2023
# Purpose: Describe the types of testing to occur in this file.
# Instruction: This code can be run in hte following ways
#              - pytest # runs all functions beginnning with the word test in the
#                         directory
#              - pytest file_name.py # Runs all functions in file_name beginning
#                                      with the word test
#              - pytest file_name.py::test_func_name # Runs only the function
#                                                      titled test_func_name in
#                                                      the file_name.py file
#              - pytest -s # Runs tests and displays when a specific file
#                            has completed testing, and what functions failed.
#                            Also displays print statments
#              - pytest -v # Displays test results on a function by function
#              - pytest -p no:warnings # Runs tests and does not display warning
#                          messages
#              - pytest -s -v -p no:warnings # Displays relevant information and
#                                supports debugging
#              - pytest -s -p no:warnings # Run for record

# Source Code Metadata
__author__ = "Jonathan A. Webb"
__copyright__ = "Copyright 2023, Jon Webb Inc."
__version__ = "1.0"
# ================================================================================
# ================================================================================
# Insert Code here


def test_pry_to_dcm():
    pitch = 0.1  # radians
    roll = 0.0  # radians
    yaw = 0.7854  # radians
    dcm = extrinsic_dir_cos_mat(pitch, roll, yaw)
    one = [0.7035729, 0.70357548, -0.09983342, -0.70710808, 0.70710548, 0.0,
           0.07059276, 0.0705930155, 0.99500417]
    new_dcm = list(dcm.flat)
    for count, value in enumerate(new_dcm):
        assert isclose(one[count], value, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_direction_cosine_matrix():
    x = 0.1
    y = 0.0
    z = 0.7854
    dcm = intrinsic_dir_cos_mat(0.1, 0.0, 0.7854, "ZYX")
    one = [0.70710548, -0.70357548, 0.07059302, 0.70710808, 0.7035729,
           -0.07059276, 0.0, 0.09983342, 0.99500417]
    new_dcm = list(dcm.flat)
    for count, value in enumerate(new_dcm):
        assert isclose(one[count], value, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_direction_cosines():
    a_vec = np.array([1., 2., 3.])
    cos_x, cos_y, cos_z = direction_cosines(a_vec)
    assert isclose(cos_x, 0.2672612, rel_tol=1.0e-3)
    assert isclose(cos_y, 0.5345224, rel_tol=1.0e-3)
    assert isclose(cos_z, 0.8017837, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_dcm_to_quaternion():
    pitch = 0.1  # radians
    roll = 0.0  # radians
    yaw = 0.7854  # radians
    dcm = extrinsic_dir_cos_mat(pitch, roll, yaw)
    q = dcm_to_quaternion(dcm)
    assert isclose(q[0], 0.01912624, rel_tol=1.0e-3)
    assert isclose(q[1], -0.04617471, rel_tol=1.0e-3)
    assert isclose(q[2], -0.38220603, rel_tol=1.0e-3)
    assert isclose(q[3], 0.92272457, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_quaternion_to_dcm():
    q = np.array([-0.1677489, -0.7369231, -0.3682588, 0.5414703])
    dcm = quaternion_to_dcm(q)
    new_dcm = list(dcm.flat)
    result = [0.142309, 0.72441893, -0.67449393, 0.36109474, -0.67249153,
              -0.64603848, -0.92159396, -0.15156633, -0.35734044]
    for count, value in enumerate(new_dcm):
        assert isclose(result[count], value, rel_tol=1.0e-3)
# ================================================================================
# ================================================================================
# eof
