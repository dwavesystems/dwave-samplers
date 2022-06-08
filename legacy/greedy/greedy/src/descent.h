// Copyright 2019 D-Wave Systems Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef _DESCENT_H
#define _DESCENT_H

#include <vector>

using std::vector;


double get_flip_energy(
    int var,
    char *state,
    const vector<double>& linear_biases,
    const vector<vector<int>>& neighbors,
    const vector<vector<double>>& neighbour_couplings
);

double get_state_energy(
    char* state,
    const vector<double>& linear_biases,
    const vector<int>& coupler_starts,
    const vector<int>& coupler_ends,
    const vector<double>& coupler_weights
);

unsigned steepest_gradient_descent_solver(
    char* state,
    const vector<double>& linear_biases,
    const vector<vector<int>>& neighbors,
    const vector<vector<double>>& neighbour_couplings
);

unsigned steepest_gradient_descent_ls_solver(
    char* state,
    const vector<double>& linear_biases,
    const vector<vector<int>>& neighbors,
    const vector<vector<double>>& neighbour_couplings
);

void steepest_gradient_descent(
    char* states,
    double* energies,
    unsigned* num_steps,
    const int num_samples,
    const vector<double>& linear_biases,
    const vector<int>& coupler_starts,
    const vector<int>& coupler_ends,
    const vector<double>& coupler_weights,
    bool large_sparse_opt=false
);

#endif
