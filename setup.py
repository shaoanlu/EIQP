from setuptools import setup, Extension, find_packages
import numpy as np
import os

# Try to find lapacke.h
possible_include_dirs = [
    '/usr/include/',
    '/usr/include/x86_64-linux-gnu/',
    '/usr/local/include/',
    '/usr/include/lapacke/',
    np.get_include()
]

eiqp_module = Extension(
    'eiqp_solver',  # Plain shared library
    sources=['eiqp/lib/EIQP_p.c'],
    include_dirs=possible_include_dirs,
    libraries=['blas', 'lapack', 'lapacke'],
    extra_compile_args=['-O3', '-fPIC'],
    extra_link_args=['-llapacke']
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
    install_requires=['numpy'],
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