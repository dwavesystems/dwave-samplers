/*
Copyright 2024 D-Wave

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

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
    GibbsNonErgodic,       /// Gibbs reflect on each variable update
    MetropolisNonErgodic,  /// Metropolis reflect one each variable update
    MetropolisUniform,     /// Metropolis on uniformly distributed angles
    MetropolisTF,          /// Metropolis on transverse field dependent angle.
};

double get_marginal_state_energyC(
    int var, std::uint8_t *state, const std::vector<double> & h,
    const std::vector<int>& degrees,
    const std::vector<std::vector<int>>& neighbors,
    const std::vector<std::vector<double>>& neighbour_couplings,
    const std::vector<double>& state_to_costheta
);

void simulated_annealing_run(
    std::uint8_t *state, const std::vector<double>& h,
    const std::vector<int>& degrees,
    const std::vector<std::vector<int>>& neighbors,
    const std::vector<std::vector<double>>& neighbour_couplings,
    const std::vector<double>& trans_fields,
    const int sweeps_per_beta,
    const std::vector<double>& Hp_field,
    const std::vector<double>& Hd_field
);

typedef bool (*const callback)(void * const function);

int general_simulated_annealing(
    std::uint8_t *states,
    double *energies,
    const int num_samples,
    const std::vector<double> h,
    const std::vector<int> coupler_starts,
    const std::vector<int> coupler_ends,
    const std::vector<double> coupler_values,
    const std::vector<double> trans_fields,
    const int sweeps_per_beta,
    const std::vector<double> Hp_field,
    const std::vector<double> Hd_field,
    const uint64_t seed,
    const bool randomize_order,
    const Proposal proposal,
    std::uint8_t *statistics,
    const int schedule_sample_interval,
    callback interrupt_callback,
    void * const interrupt_function
);

#endif
