# Import necessary packages here
import pytest
import sys
import os
from math import isclose
import numpy as np
from pathlib import PurePath
p = PurePath(__file__).parent
sys.path.insert(1, os.path.abspath(p))
from geo_py.datum import ITRF
from geo_py.frames import llh_to_ecef, ecef_to_llh, ecef_to_enu, enu_to_ecef
from geo_py.frames import llh_to_enu, enu_to_llh, ecef_to_ned, ned_to_ecef
from geo_py.frames import llh_to_ned, ned_to_llh, ned_to_enu, enu_to_ned
from geo_py.frames import ned_vector, enu_vector
# ================================================================================
# ================================================================================
# File:    test.py
# Date:    January 09, 2023
# Purpose: This file will test methods from the transform.py file
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
# TEST TRANSOFMRATIONS CLASS


@pytest.fixture
def lla1():
    """
    This method returns data to be used in other functions
    """
    lat = 46.826
    lon = 107.321
    alt = 6096.0
    return lat, lon, alt


@pytest.fixture
def lla2():
    """
    This method returns lla data to be used in other functions
    """
    lat = 45.976
    lon = 107.750
    alt = 1673.0
    return lat, lon, alt


def cmp_float_arrays(arrone, arrtwo, tol=1.0e-3):
    for one, two in zip(arrone, arrtwo):
        assert isclose(one, two, rel_tol = tol)


def test_llh_to_ecef(lla1):
    """
    This method tests the llh_to_ecef method for correct results
    """
    x, y, z = llh_to_ecef(lla1[0], lla1[1], lla1[2])
    assert isclose(x, -1302839.38, rel_tol=1.0e-3)
    assert isclose(y, 4177542.21, rel_tol=1.0e-3)
    assert isclose(z, 4632996.83, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_llh_to_enu(lla1, lla2):
    """
    This method will test the llh_to_enu function to ensure it returns a
    correct value
    """
    new_x, new_y, new_z = llh_to_enu(lla1[0], lla1[1], lla1[2],
                                     lla2[0], lla2[1], lla2[2])
    assert isclose(new_x, 33254.515, rel_tol=1.0e-3)
    assert isclose(new_y, -94415.837, rel_tol=1.0e-3)
    assert isclose(new_z, -5209.195, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_llh_to_ned(lla1, lla2):
    """
    This function will test the llh_to_ned function to ensure it returns
    a correct set of values
    """
    north, east, down = llh_to_ned(lla1[0], lla1[1], lla1[2],
                                   lla2[0], lla2[1], lla2[2])
    assert isclose(north, -94415.837, rel_tol=1.0e-3)
    assert isclose(east, 33254.515, rel_tol=1.0e-3)
    assert isclose(down, 5209.195, rel_tol=1.0e-3)
# ================================================================================
# ================================================================================
# TEST ECEF FUNCTIONS


def test_ecef_to_llh(lla1):
    """
    This method tests the llh_to_ecef method for correct results
    """
    x_val, y_val, z_val = llh_to_ecef(lla1[0], lla1[1], lla1[2])
    new_lat, new_lon, new_alt = ecef_to_llh(x_val, y_val, z_val)
    assert isclose(new_lat, lla1[0], rel_tol=1.0e-3)
    assert isclose(new_lon, lla1[1], rel_tol=1.0e-3)
    assert isclose(new_alt, lla1[2], rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_ecef_to_enu(lla1, lla2):
    """
    This function tests the ecef_to_enu method for correct results
    """
    x_val, y_val, z_val = llh_to_ecef(lla1[0], lla1[1], lla1[2])
    new_x, new_y, new_z = ecef_to_enu(lla2[0], lla2[1], lla2[2], x_val, y_val, z_val)
    assert isclose(new_x, -32764.721, rel_tol=1.0e-3)
    assert isclose(new_y, 94660.445, rel_tol=1.0e-3)
    assert isclose(new_z, 3636.220, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_ecef_to_ned(lla1):
    """
    This function tests the ecef_to_ned function to ensure it yields the correct
    results
    """
    x_val = -7134.757
    y_val = -4556.321
    z_val = 2852.39
    new_x, new_y, new_z = ecef_to_ned(lla1[0], lla1[1], lla1[2], x_val, y_val, z_val)
    assert isclose(new_x, 24918.180, rel_tol=1.0e-3)
    assert isclose(new_y, 8167.738, rel_tol=1.0e-3)
    assert isclose(new_z, 6372311.075, rel_tol=1.0e-3)
# ================================================================================
# ================================================================================
# TEST ENU FUNCTIONS


def test_enu_to_ecef(lla1, lla2):
    """
    This function tests the enu_to_ecef function to ensure it yields the
    correct result
    """
    x_val, y_val, z_val = llh_to_ecef(lla1[0], lla1[1], lla1[2])
    new_x, new_y, new_z = ecef_to_enu(lla2[0], lla2[1], lla2[2], x_val, y_val, z_val)
    x_new, y_new, z_new = enu_to_ecef(lla1[0], lla1[1], lla1[2], new_x, new_y, new_z)
    new_lat, new_lon, new_alt = ecef_to_llh(x_new, y_new, z_new)
    assert isclose(new_lat, 47.675, rel_tol=1.0e-3)
    assert isclose(new_lon, 106.885, rel_tol=1.0e-3)
    assert isclose(new_alt, 10518.35, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_enu_to_llh(lla1):
    """
    This function tests the enu_to_llh function to ensure it yields the
    correct result
    """
    east = -7134.757
    north = -4556.321
    up = 2852.39
    lat, lon, alt = enu_to_llh(lla1[0], lla1[1], lla1[2], east, north, up)
    assert isclose(lat, 46.785, rel_tol=1.0e-3)
    assert isclose(lon, 107.227, rel_tol=1.0e-3)
    assert isclose(alt, 8953.995, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_enu_to_ned(lla1):
    """
    This function tests the enu_to_ned function to ensure it yields the
    correct result
    """
    east = -7134.757
    north = -4556.321
    up = 2852.39
    new_north, new_east, new_up = enu_to_ned(lla1[0], lla1[1],
                                               lla1[2], east, north, up)
    assert isclose(new_north, -4556.321, rel_tol=1.0e-3)
    assert isclose(new_east, -7134.757, rel_tol=1.0e-3)
    assert isclose(new_up, -2852.39, rel_tol=1.0e-3)

# ================================================================================
# ================================================================================
# TEST NED FUNCTIONS


def test_ned_to_ecef(lla1):
    """
    This function tests the ned_to_ecef function to ensure it yields the
    correct results
    """
    north = 28882.282
    east = -3552.574
    down = 6372030.823
    x_val, y_val, z_val = ned_to_ecef(lla1[0], lla1[1], lla1[2], north, east, down)
    assert isclose(x_val, 4857.67, rel_tol=1.0e-3)
    assert isclose(y_val, -3643.67, rel_tol=1.0e-3)
    assert isclose(z_val, 5769.07, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_ned_to_llh(lla1, lla2):
    """
    This function tests teh ned_to_llh function to ensure it yields the
    correct reuslt
    """
    north, east, down = llh_to_ned(lla2[0], lla2[1], lla2[2],
                                  lla1[0], lla1[1], lla1[2])
    lat, lon, alt = ned_to_llh(lla2[0], lla2[1], lla2[2],
                               north, east, down)
    assert isclose(lat, lla1[0], rel_tol=1.0e-3)
    assert isclose(lon, lla1[1], rel_tol=1.0e-3)
    assert isclose(alt, lla1[2], rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_ned_to_enu(lla1, lla2):
    """
    This function tests the ned_to_enu function to ensure it yields the correct
    value
    """
    north, east, down = llh_to_ned(lla2[0], lla2[1], lla2[2],
                         lla1[0], lla1[1], lla1[2])
    new_east, new_north, new_up = ned_to_enu(lla2[0], lla2[1], lla2[2],
                                             north, east, down)
    assert isclose(new_east, -32764.72, rel_tol=1.0e-3)
    assert isclose(new_north, 94660.44, rel_tol=1.0e-3)
    assert isclose(new_up, 3636.22, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_ned_vector():
    """
    This function tests the ned_vector function to ensure it yields the correct
    results
    """
    lat = 45.976
    lon = 7.658
    alt = 4531.0
    nvec = np.array([-4403757.60, -592124.57, -4566652.06])
    evec = np.array([-592124.58, 4403757.60, 0.0])
    dvec = np.array([-3283645.50, -422918.22, 3124277.54])
    north, east, down = ned_vector(lat, lon, alt)
    cmp_float_arrays(nvec, north)
    cmp_float_arrays(evec, east)
    cmp_float_arrays(dvec, down)
# --------------------------------------------------------------------------------


def test_enu_vector():
    """
    Function development in progress
    """
    lat = 17.4114
    lon = 78.27
    E, N, U = enu_vector(lat, lon)
# --------------------------------------------------------------------------------

# def test_body():
#     pitch = radians(8.)
#     roll = 0.
#     yaw = radians(13.0)
#     lat = 57.14
#     lon = 112.3
#     alt = 6000.0
#     cos_x = 0.13
#     cos_y = -0.28
#     cos_z = 0.0717
#     vec = body(lat, lon, alt, pitch, roll, yaw, cos_x, cos_y, cos_z)
#     print(vec)
# ================================================================================
# ================================================================================
# eof
