# Copyright 2020 D-Wave Systems Inc.
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

from setuptools import setup
from setuptools.command.build_ext import build_ext
from distutils.extension import Extension

classifiers = [
    'License :: OSI Approved :: Apache Software License',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
]

python_requires = '>=3.6'

# add __version__, __author__, __authoremail__, __description__ to this namespace
exec(open("./orang/package_info.py").read())

install_requires = ['numpy>=1.15.0,<2.0.0',
                    'dimod>=0.8.0,<0.9.0',
                    'dwave_networkx>=0.7.0,<0.8.0',
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

extensions = [Extension("orang.conversions",
                        ["orang/conversions"+ext,
                         "orang/src/conversions.cpp",
                         "orang/src/solve.cpp",
                         "orang/src/sample.cpp",
                         ],
                        include_dirs=['orang/src/include']),
              Extension("orang.solve",
                        ["orang/solve"+ext,
                         "orang/src/conversions.cpp",
                         "orang/src/solve.cpp",
                         "orang/src/sample.cpp",
                         ],
                        include_dirs=['orang/src/include']),
              Extension("orang.sample",
                        ["orang/sample"+ext,
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
