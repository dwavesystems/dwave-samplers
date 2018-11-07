from __future__ import absolute_import

import sys
import os

from setuptools import setup
from setuptools.command.build_ext import build_ext
from distutils.extension import Extension
from distutils.errors import CCompilerError, DistutilsExecError, DistutilsPlatformError

classifiers = [
    'License :: OSI Approved :: Apache Software License',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7'
]

python_requires = '>=2.7,!=3.0.*,!=3.1.*,!=3.2.*,!=3.3.*'

_PY2 = sys.version_info.major == 2

# add __version__, __author__, __authoremail__, __description__ to this namespace
if _PY2:
    execfile("./orang/package_info.py")
else:
    exec(open("./orang/package_info.py").read())

install_requires = ['numpy>=1.15.0,<2.0.0',
                    'dimod>=0.7.9,<0.8.0',
                    'dwave_networkx>=0.6.6,<0.7.0',
                    ]

setup_requires = ['numpy>=1.15.0,<2.0.0']

packages = ['orang',
            ]

try:
    from Cython.Build import cythonize
except ImportError:
    USE_CYTHON = False
else:
    USE_CYTHON = True

ext = '.pyx' if USE_CYTHON else '.cpp'

extensions = [Extension("orang.orang",
                        ["orang/orang"+ext,
                         "orang/src/conversions.cpp",
                         "orang/src/solve.cpp",
                         "orang/src/sample.cpp",
                         ],
                        include_dirs=['orang/src/include']),
              Extension("orang.poly",
                        ["orang/poly"+ext,
                         "orang/src/conversions.cpp",
                         "orang/src/solve.cpp",
                         "orang/src/sample.cpp",
                         ],
                        include_dirs=['orang/src/include']),
              ]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)


class build_ext_compiler_check(build_ext):
    def run(self):
        import numpy

        self.include_dirs.append(numpy.get_include())

        build_ext.run(self)


setup(
    name='orang',
    version=__version__,
    author=__author__,
    author_email=__authoremail__,
    description=__description__,
    long_description=open('README.rst').read(),
    url='https://github.com/dwavesystems/dwave-orang',
    download_url='https://github.com/dwavesystems/dwave-orang/releases',
    license='Apache 2.0',
    packages=packages,
    install_requires=install_requires,
    cmdclass={'build_ext': build_ext_compiler_check},
    # extras_require=extras_require,
    include_package_data=True,
    classifiers=classifiers,
    zip_safe=False,
    python_requires=python_requires,
    setup_requires=setup_requires,
    ext_modules=extensions
)
