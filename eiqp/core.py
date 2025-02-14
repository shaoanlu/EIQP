import numpy as np
import ctypes
from pathlib import Path
import os
import sys
import glob

class EIQP:
    def __init__(self):
        # Look for the compiled library
        pkg_dir = Path(__file__).parent
        parent_dir = pkg_dir.parent
        
        # Pattern to match the compiled library name
        if os.name == 'nt':
            pattern = 'eiqp_solver*.dll'
        else:
            pattern = 'eiqp_solver*.so'
            
        # Search in several possible locations
        possible_locations = [
            pkg_dir,  # /path/to/eiqp/
            pkg_dir / 'lib',  # /path/to/eiqp/lib/
            parent_dir,  # /path/to/site-packages/
        ]
        
        # Find all matching library files
        lib_paths = []
        for loc in possible_locations:
            lib_paths.extend(glob.glob(str(loc / pattern)))
            
        if not lib_paths:
            raise OSError(
                f"Could not find the compiled EIQP library. Searched in:\n"
                f"{[str(p) for p in possible_locations]}\n"
                "Make sure the package is installed correctly with 'pip install -e .'"
            )

        # Try loading the library
        for lib_path in lib_paths:
            try:
                self.lib = ctypes.CDLL(lib_path)
                break
            except OSError:
                continue
        else:
            raise OSError(
                f"Found libraries at {lib_paths} but none could be loaded.\n"
                "This might be due to missing dependencies."
            )
        
        # Define function signature
        self.lib.EIQP.argtypes = [
            ctypes.POINTER(ctypes.c_double),  # Q
            ctypes.POINTER(ctypes.c_double),  # c
            ctypes.POINTER(ctypes.c_double),  # A
            ctypes.POINTER(ctypes.c_double),  # b
            ctypes.c_double,                  # epsilon
            ctypes.c_int,                     # nc
            ctypes.c_int,                     # nb
            ctypes.POINTER(ctypes.c_double)   # z
        ]
        self.lib.EIQP.restype = ctypes.c_int

    def solve(self, Q, c, A, b, epsilon=1e-8):
        """
        Solve a quadratic programming problem using EIQP.
        
        Args:
            Q (np.ndarray): Quadratic term matrix (n x n)
            c (np.ndarray): Linear term vector (n x 1)
            A (np.ndarray): Constraint matrix (m x n)
            b (np.ndarray): Constraint vector (m x 1)
            epsilon (float): Optimality level
            
        Returns:
            tuple: (z, status)
                z (np.ndarray): Solution vector
                status (bool): True if optimal solution found, False if infeasible
        """
        # Input validation and preparation
        Q = np.ascontiguousarray(Q, dtype=np.float64)
        c = np.ascontiguousarray(c.flatten(), dtype=np.float64)
        A = np.ascontiguousarray(A, dtype=np.float64)
        b = np.ascontiguousarray(b.flatten(), dtype=np.float64)
        
        nc = Q.shape[0]
        nb = A.shape[0]
        
        # Prepare output array
        z = np.zeros(nc, dtype=np.float64)
        
        # Convert to C pointers
        Q_ptr = Q.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_ptr = c.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        A_ptr = A.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        b_ptr = b.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        z_ptr = z.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        
        # Call EIQP solver
        status = self.lib.EIQP(Q_ptr, c_ptr, A_ptr, b_ptr, epsilon, nc, nb, z_ptr)
        
        return z, bool(status)