# setup.py
from setuptools import setup, Extension, find_packages
import numpy as np
import os
import platform
import sys

# Detect OS
IS_WINDOWS = platform.system() == 'Windows'

if IS_WINDOWS:
    # On Windows, use Intel MKL
    try:
        import sysconfig
        python_lib = sysconfig.get_paths()['data']
        mkl_path = os.path.join(python_lib, 'Library')
        
        library_dirs = [
            mkl_path,
            os.path.join(mkl_path, 'lib'),
        ]
        include_dirs = [
            os.path.join(mkl_path, 'include'),
            np.get_include()
        ]
        libraries = ['mkl_rt']
        extra_compile_args = ['/O2']  # Windows optimization flag
        extra_link_args = []

    except ImportError:
        print("On Windows, this package requires numpy with Intel MKL.")
        print("Please install it with: pip install numpy")
        sys.exit(1)
else:
    # Unix-like systems use system LAPACK
    library_dirs = []
    include_dirs = [
        '/usr/include/',
        '/usr/include/x86_64-linux-gnu/',
        '/usr/local/include/',
        '/usr/include/lapacke/',
        np.get_include()
    ]
    libraries = ['blas', 'lapack', 'lapacke']
    extra_compile_args = ['-O3', '-fPIC']
    extra_link_args = ['-llapacke']

# Define the extension module
eiqp_module = Extension(
    'eiqp_solver',
    sources=['eiqp/lib/EIQP_p.c'],
    include_dirs=include_dirs,
    library_dirs=library_dirs,
    libraries=libraries,
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args
)

# Create an empty __init__.py in the lib directory if it doesn't exist
lib_init = os.path.join('eiqp', 'lib', '__init__.py')
if not os.path.exists(lib_init):
    with open(lib_init, 'w') as f:
        pass

setup(
    name='eiqp',
    version='0.1.0',
    packages=find_packages(),
    ext_modules=[eiqp_module],
    install_requires=[
        'numpy',
    ],
    include_package_data=True,
    package_data={
        'eiqp.lib': ['*.c'],
    },
    author='abc',
    author_email='abc@mail.com',
    description='EIQP solver for real-time convex quadratic programming',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/liangwu2019/EIQP',
)