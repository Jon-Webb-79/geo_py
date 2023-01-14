# Import necessary packages here
import numpy as np
# - If a package and a module within the package is to be imported
#   uncomment the following lines where dir is the directory containing
#   the source files.  These lines should go above the module imports
# import sys
# import os
# sys.path.insert(1, os.path.abspath(dir))
# ================================================================================
# ================================================================================
# File:    main.py
# Date:    January 09, 2023
# Purpose: Describe the purpose of functions of this file

# Source Code Metadata
__author__ = "Jonathan A. Webb"
__copyright__ = "Copyright 2023, Jon Webb Inc."
__version__ = "1.0"
# ================================================================================
# ================================================================================
# Insert Code here
# Create rotation matrices for each axis

def direction_cosine_matrix(alpha, beta, gamma, order):
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
rot = direction_cosine_matrix(0.1, 0.0, 0.7854, "AYX")
print(rot)
# ================================================================================
# ================================================================================
# eof
