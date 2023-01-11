# Import necessary packages here
import pytest
import sys
import os
from math import isclose
from pathlib import PurePath
p = PurePath(__file__).parent
sys.path.insert(1, os.path.abspath(p))
from geo_py.datum import WGS84, NAD83, ITRF
from geo_py.transform import CFTrans
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


def test_transform_instantiate():
    """
    This function tests the ability to properly instantiate the Transformations
    class
    """
    assert isinstance(CFTrans(WGS84()), CFTrans)
# --------------------------------------------------------------------------------


def test_llh_to_ecef():
    """
    This method tests the llh_to_ecef method for correct results
    """
    lat = 46.826
    lon = 107.321
    alt = 6096.0
    tran = CFTrans(WGS84())
    x, y, z = tran.llh_to_ecef(lat, lon, alt)
    assert isclose(x, -1302839.38, rel_tol=1.0e-3)
    assert isclose(y, 4177542.21, rel_tol=1.0e-3)
    assert isclose(z, 4632996.83, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_ecef_to_llh():
    """
    This method tests the llh_to_ecef method for correct results
    """
    lat = 46.826
    lon = 107.321
    alt = 6096.0
    tran = CFTrans(WGS84())
    x, y, z = tran.llh_to_ecef(lat, lon, alt)
    new_lat, new_lon, new_alt = tran.ecef_to_llh(x, y, z)
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

    tran = CFTrans(NAD83())
    x, y, z = tran.llh_to_ecef(craft_lat, craft_lon, craft_alt)
    new_x, new_y, new_z = tran.ecef_to_enu(radar_lat, radar_lon, radar_alt, x, y, z)
    assert isclose(new_x, -7134.757, rel_tol=1.0e-3)
    assert isclose(new_y, -4556.321, rel_tol=1.0e-3)
    assert isclose(new_z, 2852.39, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_enu_to_ecef():
    radar_lat = 46.017
    radar_lon = 7.750
    radar_alt = 1673.0

    craft_lat = 45.976
    craft_lon = 7.658
    craft_alt = 4531.0

    tran = CFTrans(WGS84())
    x, y, z = tran.llh_to_ecef(craft_lat, craft_lon, craft_alt)
    new_x, new_y, new_z = tran.ecef_to_enu(radar_lat, radar_lon, radar_alt, x, y, z)
    xn, yn, zn = tran.enu_to_ecef(radar_lat, radar_lon, radar_alt, new_x, new_y, new_z)
    new_lat, new_lon, new_alt = tran.ecef_to_llh(xn, yn, zn)
    assert isclose(new_lat, 45.976, rel_tol=1.0e-3)
    assert isclose(new_lon, 7.6518, rel_tol=1.0e-3)
    assert isclose(new_alt, 4531.0, rel_tol=1.0e-3)
# ================================================================================
# ================================================================================
# eof
