// Copyright 2018 D-Wave Systems Inc.
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
//
// ===========================================================================

#ifndef _cpu_sa_h
#define _cpu_sa_h

#include <cstdint>

#ifdef _MSC_VER
    // add uint64_t definition for windows
    typedef __int64 int64_t;
    typedef unsigned __int64 uint64_t;

    // add thread_local (C++11) support for MSVC 9 and 10 (py2.7 and py3.4 on win)
    // note: thread_local support was added in MSVC 14 (Visual Studio 2015, MSC_VER 1900)
    #if _MSC_VER < 1900
    #define thread_local __declspec(thread)
    #endif
#endif

enum Proposal {
    Gibbs,       /// Gibbs acceptance critera on each variable update
    Metropolis,  /// Metropolis acceptance criteria one each variable update
};
enum VariableOrder {
    Sequential,   /// Variables updated sequentially on each sweep
    Random,       /// Variable updated uniformly at random per spin-update
};

double get_flip_energy(
    int var, std::int8_t *state, const std::vector<double> & h,
    const std::vector<int>& degrees,
    const std::vector<std::vector<int>>& neighbors,
    const std::vector<std::vector<double>>& neighbour_couplings
);

void simulated_annealing_run(
    std::int8_t *state, const std::vector<double>& h,
    const std::vector<int>& degrees,
    const std::vector<std::vector<int>>& neighbors,
    const std::vector<std::vector<double>>& neighbour_couplings,
    const int sweeps_per_beta,
    const std::vector<double>& beta_schedule,
    double * log_weight,
    const double init_energy
);

typedef bool (*const callback)(void * const function);

int general_simulated_annealing(
    std::int8_t *states,
    double *energies,
    const int num_samples,
    const std::vector<double> h,
    const std::vector<int> coupler_starts,
    const std::vector<int> coupler_ends,
    const std::vector<double> coupler_values,
    const int sweeps_per_beta,
    const std::vector<double> beta_schedule,
    const uint64_t seed,
    const VariableOrder varorder,
    const Proposal proposal_acceptance_criteria,
    callback interrupt_callback,
    void * const interrupt_function,
    double * log_weights
);

#endif
