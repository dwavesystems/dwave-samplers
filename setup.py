# Copyright 2020 D-Wave Systems Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from setuptools import setup
from setuptools.command.build_ext import build_ext as _build_ext
from distutils.extension import Extension

import dimod
import numpy as np
from Cython.Build import cythonize

classifiers = [
    'License :: OSI Approved :: Apache Software License',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
]

python_requires = '>=3.7'

# add __version__, __author__, __authoremail__, __description__ to this namespace
exec(open("./orang/package_info.py").read())

install_requires = ['numpy>=1.15.0,<2.0.0',
                    'dimod>=0.8.0,<0.11.0',
                    'dwave_networkx>=0.7.0,<0.9.0',
                    ]

setup_requires = ['numpy>=1.15.0,<2.0.0']

packages = ['orang',
            ]

extensions = [Extension("orang.solve",
                        ["orang/solve.pyx",
                         ],
                        include_dirs=['orang/src/include', np.get_include(), dimod.get_include()]),
              Extension("orang.sample",
                        ["orang/sample.pyx",
                         ],
                        include_dirs=['orang/src/include', np.get_include(), dimod.get_include()]),
              ]

class build_ext(_build_ext):
    extra_compile_args = {
        'msvc': ['/EHsc'],
        'unix': ['-std=c++11'],
    }

    extra_link_args = {
        'msvc': [],
        'unix': ['-std=c++11'],
    }

    def build_extensions(self):
        compiler = self.compiler.compiler_type

        compile_args = self.extra_compile_args[compiler]
        for ext in self.extensions:
            ext.extra_compile_args.extend(compile_args)

        link_args = self.extra_link_args[compiler]
        for ext in self.extensions:
            ext.extra_compile_args.extend(link_args)

        super().build_extensions()

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
    cmdclass={'build_ext': build_ext},
    include_package_data=True,
    classifiers=classifiers,
    zip_safe=False,
    python_requires=python_requires,
    setup_requires=setup_requires,
    ext_modules=cythonize(extensions),
)
