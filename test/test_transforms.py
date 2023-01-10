# Import necessary packages here
import pytest
import sys
import os
from math import isclose
from pathlib import PurePath
p = PurePath(__file__).parent
sys.path.insert(1, os.path.abspath(p))
from geo_py.datum import WGS84
from geo_py.transform import Transformations
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
    assert isinstance(Transformations(WGS84()), Transformations)
# --------------------------------------------------------------------------------


def test_llh_to_ecef():
    """
    This method tests the llh_to_ecef method for correct results
    """
    lat = 46.826
    lon = 107.321
    alt = 6096.0
    tran = Transformations(WGS84())
    x, y, z = tran.llh_to_ecef(lat, lon, alt)
    assert isclose(x, -1302839.38, rel_tol=1.0e-3)
    assert isclose(y, 4177542.21, rel_tol=1.0e-3)
    assert isclose(z, 4632996.83, rel_tol=1.0e-3)
# ================================================================================
# ================================================================================
# eof
