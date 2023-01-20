# Import necessary packages here
import pytest
import sys
import os
from math import isclose
from geopy.distance import great_circle
from pathlib import PurePath
p = PurePath(__file__).parent
sys.path.insert(1, os.path.abspath(p))
from geo_py.distance import Distance, vincenty
from geo_py.datum import ITRF
# ================================================================================
# ================================================================================
# File:    test_distance.py
# Date:    January 17, 2023
# Purpose: This file contains functions that will test the capabilities of the
#          distance.py file
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


@pytest.fixture
def dist_lla1():
    """
    This method contains data formatted for local functions
    """
    lat1 = 40.98197
    lon1 = 111.9026
    alt1 = 0.0  # 1306.0
    lat2 = 45.01443
    lon2 = 113.9278
    alt2 = 0.0  # 1292.0
    return (lat1, lon1, alt1), (lat2, lon2, alt2)
# --------------------------------------------------------------------------------


@pytest.fixture
def dist_lla2():
    """
    This method contains data formatted for GeoPy functions
    """
    lat1 = 40.98197
    lon1 = 111.9026
    lat2 = 45.01443
    lon2 = 113.9278
    return (lat1, lon1), (lat2, lon2)
# ================================================================================
# ================================================================================


def test_haversine(dist_lla1, dist_lla2):
    """
    This function tests the Haversine method for correctness against the Geopy
    great_circle function for units of feet, meters, km, and miles
    """
    assert isclose(Distance(dist_lla1[0], dist_lla1[1]).km,
                   great_circle(dist_lla2[0], dist_lla2[1]).km, rel_tol=1.0e-3)
    assert isclose(Distance(dist_lla1[0], dist_lla1[1]).m,
                   great_circle(dist_lla2[0], dist_lla2[1]).m, rel_tol=1.0e-3)
    assert isclose(Distance(dist_lla1[0], dist_lla1[1]).miles,
                   great_circle(dist_lla2[0], dist_lla2[1]).miles, rel_tol=1.0e-3)
    assert isclose(Distance(dist_lla1[0], dist_lla1[1]).feet,
                   great_circle(dist_lla2[0], dist_lla2[1]).feet, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_great_circle(dist_lla1, dist_lla2):
    """
    This function tests the great_circle method for correctness against the Geopy
    great_circle function for units of feet, meters, km, and miles
    """
    assert isclose(Distance(dist_lla1[0], dist_lla1[1], method="GREAT CIRCLE").km,
                   great_circle(dist_lla2[0], dist_lla2[1]).km, rel_tol=1.0e-3)
    assert isclose(Distance(dist_lla1[0], dist_lla1[1], method="GREAT CIRCLE").m,
                   great_circle(dist_lla2[0], dist_lla2[1]).m, rel_tol=1.0e-3)
    assert isclose(Distance(dist_lla1[0], dist_lla1[1], method="GREAT CIRCLE").miles,
                   great_circle(dist_lla2[0], dist_lla2[1]).miles, rel_tol=1.0e-3)
    assert isclose(Distance(dist_lla1[0], dist_lla1[1], method="GREAT CIRCLE").feet,
                   great_circle(dist_lla2[0], dist_lla2[1]).feet, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_vincenti(dist_lla1, dist_lla2):
    """
    This function tests the vincenty for correctness against the Geopy
    great_circle function for units of feet, meters, km, and miles
    """
    assert isclose(Distance(dist_lla1[0], dist_lla1[1], method="VINCENTY").km,
                   great_circle(dist_lla2[0], dist_lla2[1]).km, rel_tol=1.0e-3)
    assert isclose(Distance(dist_lla1[0], dist_lla1[1], method="VINCENTY").m,
                   great_circle(dist_lla2[0], dist_lla2[1]).m, rel_tol=1.0e-3)
    assert isclose(Distance(dist_lla1[0], dist_lla1[1], method="VINCENTY").miles,
                   great_circle(dist_lla2[0], dist_lla2[1]).miles, rel_tol=1.0e-3)
    assert isclose(Distance(dist_lla1[0], dist_lla1[1], method="VINCENTY").feet,
                   great_circle(dist_lla2[0], dist_lla2[1]).feet, rel_tol=1.0e-3)
# --------------------------------------------------------------------------------


def test_linear(dist_lla1, dist_lla2):
    """
    This function tests the linear_dist method for correctness against
    the Geopy great_circle function for units of feet, meters, km, and miles
    """
    assert isclose(Distance(dist_lla1[0], dist_lla1[1], method="LINEAR").km,
                   great_circle(dist_lla2[0], dist_lla2[1]).km, rel_tol=1.0e-3)
    assert isclose(Distance(dist_lla1[0], dist_lla1[1], method="LINEAR").m,
                   great_circle(dist_lla2[0], dist_lla2[1]).m, rel_tol=1.0e-3)
    assert isclose(Distance(dist_lla1[0], dist_lla1[1], method="LINEAR").miles,
                   great_circle(dist_lla2[0], dist_lla2[1]).miles, rel_tol=1.0e-3)
    assert isclose(Distance(dist_lla1[0], dist_lla1[1], method="LINEAR").feet,
                   great_circle(dist_lla2[0], dist_lla2[1]).feet, rel_tol=1.0e-3)
# ================================================================================
# ================================================================================
# eof
