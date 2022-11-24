// Copyright 2019 D-Wave Systems Inc.

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

#include <vector>

#include "../Catch2/single_include/catch2/catch.hpp"
#include "descent.h"


TEST_CASE("Test steepest_gradient_descent") {
    int num_samples = 1;

    int8_t states[3] = {1, 1, 1}, min_states[3] = {-1, 1, 1};
    double energies[1] = {0}, min_energies[1] = {-1};
    unsigned num_steps[1] = {0};

    // bqm ~ {(0, 1): 1, (1, 2): 1, (2, 0): 1})
    vector<double> linear_biases {0, 0, 0};
    vector<int> coupler_starts {0, 1, 2};
    vector<int> coupler_ends {1, 2, 0};
    vector<double> coupler_weights {1.0, 1.0, 1.0};

    steepest_gradient_descent(
        states, energies, num_steps, num_samples,
        linear_biases, coupler_starts, coupler_ends, coupler_weights
    );

    // assert correct solution
    for (auto i = 0; i < 3; i++) {
        CHECK(states[i] == min_states[i]);
    }
    for (auto i = 0; i < 1; i++) {
        CHECK(energies[i] == min_energies[i]);
    }
}
