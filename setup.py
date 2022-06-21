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

from setuptools import setup
from Cython.Build import cythonize
from setuptools.command.build_ext import build_ext

import dimod
import numpy


class build_ext_with_args(build_ext):
    """Add compiler-specific compile/link flags."""

    extra_compile_args = {
        'msvc': ['/std:c++14'],
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
            ext.extra_compile_args = compile_args

        link_args = self.extra_link_args[compiler]
        for ext in self.extensions:
            ext.extra_compile_args = link_args

        super().build_extensions()



setup(
    cmdclass={'build_ext': build_ext_with_args},
    ext_modules=cythonize(
        ['dwave/samplers/greedy/descent.pyx',
         'dwave/samplers/random/*.pyx',
         'dwave/samplers/sa/*.pyx',
         'dwave/samplers/tabu/tabu_search.pyx',
         'dwave/samplers/tree/*.pyx',
         ],
        ),
    include_dirs=[
        dimod.get_include(),
        numpy.get_include(),
        ],

    )
