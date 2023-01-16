# Import necessary packages here
import sys
import os
from math import isclose, radians
import numpy as np
from pathlib import PurePath
p = PurePath(__file__).parent
sys.path.insert(1, os.path.abspath(p))
from geo_py.datum import ITRF
from geo_py.frames import llh_to_ecef, ecef_to_llh, ecef_to_enu, enu_to_ecef
from geo_py.frames import llh_to_enu, enu_to_llh, ecef_to_ned, ned_to_ecef
from geo_py.frames import llh_to_ned, ned_to_llh, ned_to_enu, enu_to_ned
from geo_py.frames import ned_vector, body, body_to_ecef
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


def test_llh_to_ecef():
    """
    This method tests the llh_to_ecef method for correct results
    """
    lat = 46.826
    lon = 107.321
    alt = 6096.0
    x, y, z = llh_to_ecef(lat, lon, alt)
    assert isclose(x, -1302839.38, rel_tol=1.0e-3)
    assert isclose(y, 4177542.21, rel_tol=1.0e-3)
    assert isclose(z, 4632996.83, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_llh_to_enu():
    radar_lat = 46.017
    radar_lon = 7.750
    radar_alt = 1673.0

    craft_lat = 45.976
    craft_lon = 7.658
    craft_alt = 4531.0

    new_x, new_y, new_z = llh_to_enu(radar_lat, radar_lon, radar_alt,
                                     craft_lat, craft_lon, craft_alt)
    assert isclose(new_x, -7134.757, rel_tol=1.0e-3)
    assert isclose(new_y, -4556.321, rel_tol=1.0e-3)
    assert isclose(new_z, 2852.39, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_llh_to_ned():
    radar_lat = 46.017
    radar_lon = 7.750
    radar_alt = 1673.0

    craft_lat = 45.976
    craft_lon = 7.658
    craft_alt = 4531.0

    N, E, D = llh_to_ned(radar_lat, radar_lon, radar_alt,
                         craft_lat, craft_lon, craft_alt)
    assert isclose(N, -4556.321, rel_tol=1.0e-3)
    assert isclose(E, -7134.752, rel_tol=1.0e-3)
    assert isclose(D, -2852.39, rel_tol=1.0e-3)
# ================================================================================
# ================================================================================
# TEST ECEF FUNCTIONS


def test_ecef_to_llh():
    """
    This method tests the llh_to_ecef method for correct results
    """
    lat = 46.826
    lon = 107.321
    alt = 6096.0
    x, y, z = llh_to_ecef(lat, lon, alt)
    new_lat, new_lon, new_alt = ecef_to_llh(x, y, z)
    assert isclose(new_lat, lat, rel_tol=1.0e-3)
    assert isclose(new_lon, lon, rel_tol=1.0e-3)
    assert isclose(new_alt, alt, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_ecef_to_enu():
    """
    This function tests the ecef_to_enu method for correct results
    """
    radar_lat = 46.017
    radar_lon = 7.750
    radar_alt = 1673.0

    craft_lat = 45.976
    craft_lon = 7.658
    craft_alt = 4531.0

    x, y, z = llh_to_ecef(craft_lat, craft_lon, craft_alt)
    new_x, new_y, new_z = ecef_to_enu(radar_lat, radar_lon, radar_alt, x, y, z)
    assert isclose(new_x, -7134.757, rel_tol=1.0e-3)
    assert isclose(new_y, -4556.321, rel_tol=1.0e-3)
    assert isclose(new_z, 2852.39, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_ecef_to_ned():
    radar_lat = 46.017
    radar_lon = 7.750
    radar_alt = 1673.0

    x = -7134.757
    y = -4556.321
    z = 2852.39

    new_x, new_y, new_z = ecef_to_ned(radar_lat, radar_lon, radar_alt, x, y, z)
    assert isclose(new_x, 28882.282, rel_tol=1.0e-3)
    assert isclose(new_y, -3552.574, rel_tol=1.0e-3)
    assert isclose(new_z, 6372030.823, rel_tol=1.0e-3)
# ================================================================================
# ================================================================================
# TEST ENU FUNCTIONS


def test_enu_to_ecef():
    radar_lat = 46.017
    radar_lon = 7.750
    radar_alt = 1673.0

    craft_lat = 45.976
    craft_lon = 7.658
    craft_alt = 4531.0

    x, y, z = llh_to_ecef(craft_lat, craft_lon, craft_alt)
    new_x, new_y, new_z = ecef_to_enu(radar_lat, radar_lon, radar_alt, x, y, z)
    xn, yn, zn = enu_to_ecef(radar_lat, radar_lon, radar_alt, new_x, new_y, new_z)
    new_lat, new_lon, new_alt = ecef_to_llh(xn, yn, zn)
    assert isclose(new_lat, 45.976, rel_tol=1.0e-3)
    assert isclose(new_lon, 7.6518, rel_tol=1.0e-3)
    assert isclose(new_alt, 4531.0, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_enu_to_llh():
    radar_lat = 46.017
    radar_lon = 7.750
    radar_alt = 1673.0

    E = -7134.757
    N = -4556.321
    U = 2852.39

    lat, lon, alt = enu_to_llh(radar_lat, radar_lon, radar_alt, E, N, U)
    assert isclose(lat, 45.976, rel_tol=1.0e-3)
    assert isclose(lon, 7.6518, rel_tol=1.0e-3)
    assert isclose(alt, 4531.0, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_enu_to_ned():
    radar_lat = 46.017
    radar_lon = 7.750
    radar_alt = 1673.0

    E = -7134.757
    N = -4556.321
    U = 2852.39

    NN, EE, DD = enu_to_ned(radar_lat, radar_lon, radar_alt, E, N, U)
    assert isclose(NN, -4556.321, rel_tol=1.0e-3)
    assert isclose(EE, -7134.757, rel_tol=1.0e-3)
    assert isclose(DD, -2852.39, rel_tol=1.0e-3)

# ================================================================================
# ================================================================================
# TEST NED FUNCTIONS


def test_ned_to_ecef():
    radar_lat = 46.017
    radar_lon = 7.750
    radar_alt = 1673.0

    N = 28882.282
    E = -3552.574
    D = 6372030.823

    x, y, z = ned_to_ecef(radar_lat, radar_lon, radar_alt, N, E, D)
    assert isclose(x, -7134.757, rel_tol=1.0e-3)
    assert isclose(y, -4556.321, rel_tol=1.0e-3)
    assert isclose(z, 2852.389, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_ned_to_llh():
    radar_lat = 46.017
    radar_lon = 7.750
    radar_alt = 1673.0

    craft_lat = 45.976
    craft_lon = 7.658
    craft_alt = 4531.0

    N, E, D = llh_to_ned(radar_lat, radar_lon, radar_alt,
                         craft_lat, craft_lon, craft_alt)
    lat, lon, alt = ned_to_llh(radar_lat, radar_lon, radar_alt,
                               N, E, D)
    assert isclose(lat, craft_lat, rel_tol=1.0e-3)
    assert isclose(lon, craft_lon, rel_tol=1.0e-3)
    assert isclose(alt, craft_alt, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_ned_to_enu():
    radar_lat = 46.017
    radar_lon = 7.750
    radar_alt = 1673.0

    craft_lat = 45.976
    craft_lon = 7.658
    craft_alt = 4531.0

    N, E, D = llh_to_ned(radar_lat, radar_lon, radar_alt,
                         craft_lat, craft_lon, craft_alt)
    EE, NN, UU = ned_to_enu(radar_lat, radar_lon, radar_alt,
                            N, E, D)
    assert isclose(EE, -7134.757, rel_tol=1.0e-3)
    assert isclose(NN, -4556.321, rel_tol=1.0e-3)
    assert isclose(UU, 2852.390, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_ned_vector():
    lat = 45.976
    lon = 7.658
    alt = 4531.0
    N, E, D = ned_vector(lat, lon, alt)
# --------------------------------------------------------------------------------


def test_body():
    pitch = radians(8.)
    roll = 0.
    yaw = radians(13.0)
    lat = 57.14
    lon = 112.3
    alt = 6000.0
    cos_x = 0.13
    cos_y = -0.28
    cos_z = 0.0717
    vec = body(lat, lon, alt, pitch, roll, yaw, cos_x, cos_y, cos_z)
    print(vec)
# ================================================================================
# ================================================================================
# eof
