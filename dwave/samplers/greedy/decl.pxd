# Copyright 2019 D-Wave Systems Inc.
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

from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "descent.h":

    unsigned int steepest_gradient_descent(
        char* states,
        double* energies,
        unsigned* num_steps,
        const int num_samples,
        const vector[double]& linear_biases,
        const vector[int]& coupler_starts,
        const vector[int]& coupler_ends,
        const vector[double]& coupler_weights,
        bool large_sparse_opt
    ) nogil
