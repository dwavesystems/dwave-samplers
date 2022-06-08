# Copyright 2022 D-Wave Systems Inc.
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

from Cython.Build import cythonize
from setuptools import setup

import dimod
import numpy


setup(
    ext_modules=cythonize(
        ['dwave/samplers/test.pyx',
         'dwave/samplers/greedy/descent.pyx',
         'dwave/samplers/orang/sample.pyx',
         'dwave/samplers/orang/solve.pyx',
         ],
        ),
    include_dirs=[
        dimod.get_include(),
        numpy.get_include(),
        ],

    )
