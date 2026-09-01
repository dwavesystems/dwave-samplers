/*
Copyright 2025 D-Wave

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

#ifndef _cpu_dsbm_h
#define _cpu_dsbm_h

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


void discrete_simulated_bifurcation_run(
    double* state_x,
    double* state_y,
    const int num_vars,
    const std::vector<std::vector<int>>& neighbors,
    const std::vector<std::vector<double>>& neighbour_couplings,
    const std::vector<double>& a_schedule,
    const double a0,
    const double c0,
    const double _Delta_t
);

typedef bool (*const callback)(void * const function);

int general_discrete_simulated_bifurcation_machine(
    double* states_x,
    double* states_y,
    const int num_samples,
    const int num_vars,
    const std::vector<int> coupler_starts,
    const std::vector<int> coupler_ends,
    const std::vector<double> coupler_weights,
    const std::vector<double> a_schedule,
    const double a0,
    const double c0,
    const double Delta_t,
    callback interrupt_callback,
    void * const interrupt_function
);

#endif
