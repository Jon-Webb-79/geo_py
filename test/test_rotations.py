# Import necessary packages here
from math import isclose
import os
import sys
from pathlib import PurePath
from scipy.spatial.transform import Rotation
p = PurePath(__file__).parent
sys.path.insert(1, os.path.abspath(p))
from geo_py.rotations import rotation_matrix, dcm_euler_angles, aircraft_rotation
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


def test_rot_mat_xyz():
    pitch = 22.1  # degrees
    roll = 8.75  # degrees
    yaw = 12.8  # degrees
    sequence = "XYZ"  # Z: yaw, Y: pitch, X: roll
    dcm = rotation_matrix(yaw, pitch, roll, order=sequence,
                          deg=True)
    rot = Rotation.from_euler(sequence, [roll, pitch, yaw], degrees=True)
    scipy_dcm = rot.as_matrix()
    dcm = list(dcm.flat)
    scipy_dcm = list(scipy_dcm.flat)
    for one, two in zip(dcm, scipy_dcm):
        assert isclose(one, two, rel_tol=1.0e-3)

    pitch = 22.1
    roll = 8.75
    yaw = 12.8
    dcm = rotation_matrix(yaw, pitch, roll, order="XYZ", deg=True)
    print(dcm)
# --------------------------------------------------------------------------------


def test_rot_mat_xzy():
    pitch = 22.1  # degrees
    roll = 8.75  # degrees
    yaw = 12.8  # degrees
    sequence = "XZY"  # Z: yaw, Y: pitch, X: roll
    dcm = rotation_matrix(yaw, pitch, roll, order=sequence,
                          deg=True)
    rot = Rotation.from_euler(sequence, [roll, yaw, pitch], degrees=True)
    scipy_dcm = rot.as_matrix()
    dcm = list(dcm.flat)
    scipy_dcm = list(scipy_dcm.flat)
    for one, two in zip(dcm, scipy_dcm):
        assert isclose(one, two, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_rot_mat_yxz():
    pitch = 22.1  # degrees
    roll = 8.75  # degrees
    yaw = 12.8  # degrees
    sequence = "YXZ"  # Z: yaw, Y: pitch, X: roll
    dcm = rotation_matrix(yaw, pitch, roll, order=sequence,
                          deg=True)
    rot = Rotation.from_euler(sequence, [pitch, roll, yaw], degrees=True)
    scipy_dcm = rot.as_matrix()
    dcm = list(dcm.flat)
    scipy_dcm = list(scipy_dcm.flat)
    for one, two in zip(dcm, scipy_dcm):
        assert isclose(one, two, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_rot_mat_zyx():
    pitch = 22.1  # degrees
    roll = 8.75  # degrees
    yaw = 12.8  # degrees
    sequence = "ZYX"  # Z: yaw, Y: pitch, X: roll
    dcm = rotation_matrix(yaw, pitch, roll, order=sequence,
                          deg=True)
    rot = Rotation.from_euler(sequence, [yaw, pitch, roll], degrees=True)
    scipy_dcm = rot.as_matrix()
    dcm = list(dcm.flat)
    scipy_dcm = list(scipy_dcm.flat)
    for one, two in zip(dcm, scipy_dcm):
        assert isclose(one, two, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_rot_mat_zxy():
    pitch = 22.1  # degrees
    roll = 8.75  # degrees
    yaw = 12.8  # degrees
    sequence = "ZXY"  # Z: yaw, Y: pitch, X: roll
    dcm = rotation_matrix(yaw, pitch, roll, order=sequence,
                          deg=True)
    rot = Rotation.from_euler(sequence, [yaw, roll, pitch], degrees=True)
    scipy_dcm = rot.as_matrix()
    dcm = list(dcm.flat)
    scipy_dcm = list(scipy_dcm.flat)
    for one, two in zip(dcm, scipy_dcm):
        assert isclose(one, two, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_rot_mat_extrinsic():
    pitch = 22.1  # degrees
    roll = 8.75  # degrees
    yaw = 12.8  # degrees
    sequence = "ZXY"  # Z: yaw, Y: pitch, X: roll
    dcm = rotation_matrix(yaw, pitch, roll, order=sequence,
                          deg=True, extrinsic=True)
    rot = Rotation.from_euler(sequence, [yaw, roll, pitch], degrees=True)
    scipy_dcm = rot.as_matrix().T
    dcm = list(dcm.flat)
    scipy_dcm = list(scipy_dcm.flat)
    for one, two in zip(dcm, scipy_dcm):
        assert isclose(one, two, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_dcm_to_euler_xyz():
    pitch = 22.1  # degrees
    roll = 8.75  # degrees
    yaw = 12.8  # degrees
    sequence = "XYZ"  # Z: yaw, Y: pitch, X: roll
    dcm = rotation_matrix(yaw, pitch, roll, order=sequence,
                           deg=True)
    nyaw, npitch, nroll = dcm_euler_angles(dcm, order=sequence, deg=True)
    assert isclose(yaw, nyaw, rel_tol=1.0e-3)
    assert isclose(pitch, npitch, rel_tol=1.0e-3)
    assert isclose(roll, nroll, rel_tol=1.0e-3)
# # --------------------------------------------------------------------------------


def test_dcm_to_euler_xzy():
    pitch = 22.1  # degrees
    roll = 8.75  # degrees
    yaw = 12.8  # degrees
    sequence = "XZY"  # Z: yaw, Y: pitch, X: roll
    dcm = rotation_matrix(yaw, pitch, roll, order=sequence,
                          deg=True)
    nyaw, npitch, nroll = dcm_euler_angles(dcm, order=sequence, deg=True)
    assert isclose(yaw, nyaw, rel_tol=1.0e-3)
    assert isclose(pitch, npitch, rel_tol=1.0e-3)
    assert isclose(roll, nroll, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_dcm_to_euler_yxz():
    pitch = 22.1  # degrees
    roll = 8.75  # degrees
    yaw = 12.8  # degrees
    sequence = "YXZ"  # Z: yaw, Y: pitch, X: roll
    dcm = rotation_matrix(yaw, pitch, roll, order=sequence,
                          deg=True)
    nyaw, npitch, nroll = dcm_euler_angles(dcm, order=sequence, deg=True)
    assert isclose(yaw, nyaw, rel_tol=1.0e-3)
    assert isclose(pitch, npitch, rel_tol=1.0e-3)
    assert isclose(roll, nroll, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_dcm_to_euler_yzx():
    pitch = 22.1  # degrees
    roll = 8.75  # degrees
    yaw = 12.8  # degrees
    sequence = "YZX"  # Z: yaw, Y: pitch, X: roll
    dcm = rotation_matrix(yaw, pitch, roll, order=sequence,
                          deg=True)
    nyaw, npitch, nroll = dcm_euler_angles(dcm, order=sequence, deg=True)
    assert isclose(yaw, nyaw, rel_tol=1.0e-3)
    assert isclose(pitch, npitch, rel_tol=1.0e-3)
    assert isclose(roll, nroll, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_dcm_to_euler_zyx():
    pitch = 22.1  # degrees
    roll = 8.75  # degrees
    yaw = 12.8  # degrees
    sequence = "ZYX"  # Z: yaw, Y: pitch, X: roll
    dcm = rotation_matrix(yaw, pitch, roll, order=sequence,
                          deg=True)
    nyaw, npitch, nroll = dcm_euler_angles(dcm, order=sequence, deg=True)
    assert isclose(yaw, nyaw, rel_tol=1.0e-3)
    assert isclose(pitch, npitch, rel_tol=1.0e-3)
    assert isclose(roll, nroll, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------

def test_dcm_to_euler_zxy():
    pitch = 22.1  # degrees
    roll = 8.75  # degrees
    yaw = 12.8  # degrees
    sequence = "ZXY"  # Z: yaw, Y: pitch, X: roll
    dcm = rotation_matrix(yaw, pitch, roll, order=sequence,
                          deg=True)
    nyaw, npitch, nroll = dcm_euler_angles(dcm, order=sequence, deg=True)
    assert isclose(yaw, nyaw, rel_tol=1.0e-3)
    assert isclose(pitch, npitch, rel_tol=1.0e-3)
    assert isclose(roll, nroll, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_ac_rotation():
    pitch = 22.1  # degrees
    roll = 8.75  # degrees
    yaw = 12.8  # degrees
    sequence = "ZYX"  # Z: yaw, Y: pitch, X: roll
    dcm = rotation_matrix(yaw, pitch, roll, order=sequence,
                          deg=True)
    new_dcm = aircraft_rotation(yaw, pitch, roll, deg=True)
    dcm = list(dcm.flat)
    new_dcm = list(new_dcm.flat)
    for one, two in zip(dcm, new_dcm):
        assert isclose(one, two, rel_tol=1.0e-3)
# ================================================================================
# ================================================================================
# QUATERNIANS TEST


# ================================================================================
# ================================================================================
# eof
