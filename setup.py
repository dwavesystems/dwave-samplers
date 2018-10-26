from __future__ import absolute_import

import sys
import os

from setuptools import setup
from distutils.extension import Extension
from distutils.command.build_ext import build_ext
from distutils.errors import CCompilerError, DistutilsExecError, DistutilsPlatformError

# # add __version__, __author__, __authoremail__, __description__ to this namespace
# _PY2 = sys.version_info.major == 2
# my_loc = os.path.dirname(os.path.abspath(__file__))
# os.chdir(my_loc)
# if _PY2:
#     execfile(os.path.join(".", "orang", "package_info.py"))
# else:
#     exec(open(os.path.join(".", "orang", "package_info.py")).read())

install_requires = []

packages = []

try:
    from Cython.Build import cythonize
except ImportError:
    USE_CYTHON = False
else:
    USE_CYTHON = True

ext = '.pyx' if USE_CYTHON else '.cpp'

extensions = [Extension("orang._orang",
                        ["orang/orang"+ext,
                         "orang/src/conversions.cpp",
                         "orang/src/solve.cpp",
                         "orang/src/sample.cpp"],
                        include_dirs=['orang/src/include',
                                      # 'orang/src/include/operations'
                                      ])]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

setup(
    name='orang',
    # version=__version__,
    # author=__author__,
    # author_email=__authoremail__,
    # description=__description__,
    # long_description=open('README.rst').read(),
    # url='https://github.com/dwavesystems/dimod',
    # download_url='https://github.com/dwavesystems/dimod/releases',
    # license='Apache 2.0',
    packages=packages,
    install_requires=install_requires,
    # extras_require=extras_require,
    # include_package_data=True,
    # classifiers=classifiers,
    zip_safe=False,
    # python_requires=python_requires,
    ext_modules=extensions
)
