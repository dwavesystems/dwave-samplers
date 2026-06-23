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

#include <cstdint>
#include <math.h>
#include <atomic>
#include <thread>
#include <vector>
#include <stdexcept>
#include "fast_cpu_sa.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// xorshift128+ as defined https://en.wikipedia.org/wiki/Xorshift#xorshift.2B
#define FASTRAND(rand) do {                       \
    uint64_t x = fast_rng_state[0];                    \
    uint64_t const y = fast_rng_state[1];              \
    fast_rng_state[0] = y;                             \
    x ^= x << 23;                                 \
    fast_rng_state[1] = x ^ y ^ (x >> 17) ^ (y >> 26); \
    rand = fast_rng_state[1] + y;                      \
} while (0)

#define RANDMAX ((uint64_t)-1L)

// OPENMP_THRESHOLD_SAMPLES (32):
//   Each sample costs tens of ms even for small problems.
//   std::thread break-even is ~3-5 samples; OpenMP ~1 sample.
//   32 ensures OpenMP's static schedule wins over std::thread
//   work-stealing coordination cost at moderate sample counts.
//   For vars and couplers, OpenMP is always used unconditionally
//   since its wakeup cost (~10 us) is always worth it vs spawning
//   fresh std::threads on every single sample call.
static constexpr int OPENMP_THRESHOLD_SAMPLES = 32;

using namespace std;

// this holds the state of our thread-safe/local RNG
thread_local uint64_t fast_rng_state[2];

// Returns the energy delta from flipping variable at index `var`
double fast_get_flip_energy(
    int var,
    std::int8_t *state,
    const vector<double>& h,
    const vector<int>& degrees,
    const vector<vector<int>>& neighbors,
    const vector<vector<double>>& neighbour_couplings
) {
    double energy = h[var];
    for (int n_i = 0; n_i < degrees[var]; n_i++) {
        energy += state[neighbors[var][n_i]] * neighbour_couplings[var][n_i];
    }
    return -2 * state[var] * energy;
}

// ---------------------------------------------------------------------------
// Parallel helpers
// ---------------------------------------------------------------------------

// Parallel delta-energy initialisation.
// Uses OpenMP for large num_vars, std::thread otherwise.
// Each element is independent so no synchronisation is needed.
static void parallel_init_delta_energy(
    double* delta_energy,
    const int num_vars,
    std::int8_t* state,
    const vector<double>& h,
    const vector<int>& degrees,
    const vector<vector<int>>& neighbors,
    const vector<vector<double>>& neighbour_couplings
) {
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    for (int var = 0; var < num_vars; var++) {
        delta_energy[var] = fast_get_flip_energy(
            var, state, h, degrees, neighbors, neighbour_couplings);
    }
#else
    // OpenMP unavailable: fall back to sequential
    for (int var = 0; var < num_vars; var++) {
        delta_energy[var] = fast_get_flip_energy(
            var, state, h, degrees, neighbors, neighbour_couplings);
    }
#endif
}

// Parallel state-energy computation.
// Uses OpenMP reduction for large coupler counts, std::thread otherwise.
static double parallel_get_state_energy(
    std::int8_t* state,
    const vector<double>& h,
    const vector<int>& coupler_starts,
    const vector<int>& coupler_ends,
    const vector<double>& coupler_weights
) {
    double energy = 0.0;

    // h-field contribution: sequential, typically much cheaper than coupler sum
    for (unsigned int var = 0; var < h.size(); var++) {
        energy += state[var] * h[var];
    }

    const int num_couplers = static_cast<int>(coupler_starts.size());

#ifdef _OPENMP
    double coupler_energy = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:coupler_energy)
    for (int c = 0; c < num_couplers; c++) {
        coupler_energy += state[coupler_starts[c]] *
                          coupler_weights[c] *
                          state[coupler_ends[c]];
    }
    return energy + coupler_energy;
#else
    // OpenMP unavailable: fall back to sequential
    for (int c = 0; c < num_couplers; c++) {
        energy += state[coupler_starts[c]] *
                  coupler_weights[c] *
                  state[coupler_ends[c]];
    }
    return energy;
#endif
}

// ---------------------------------------------------------------------------
// Core SA run (unchanged except delta_energy init now uses parallel helper)
// ---------------------------------------------------------------------------

template <VariableOrder varorder, Proposal proposal_acceptance_criteria>
void fast_simulated_annealing_run(
    std::int8_t* state,
    const vector<double>& h,
    const vector<int>& degrees,
    const vector<vector<int>>& neighbors,
    const vector<vector<double>>& neighbour_couplings,
    const int sweeps_per_beta,
    const vector<double>& beta_schedule
) {
    const int num_vars = h.size();

    double *delta_energy = (double*)malloc(num_vars * sizeof(double));

    uint64_t rand;

    // --- PARALLELIZED: delta energy initialisation ---
    parallel_init_delta_energy(
        delta_energy, num_vars, state, h, degrees, neighbors, neighbour_couplings);

    bool flip_spin;
    for (int beta_idx = 0; beta_idx < (int)beta_schedule.size(); beta_idx++) {
        const double beta = beta_schedule[beta_idx];
        for (int sweep = 0; sweep < sweeps_per_beta; sweep++) {

            const double threshold = 44.36142 / beta;
            for (int varI = 0; varI < num_vars; varI++) {
                int var;
                if constexpr (varorder == Random) {
                    FASTRAND(rand);
                    var = rand % num_vars;
                } else {
                    var = varI;
                }
                if (delta_energy[var] >= threshold) continue;

                flip_spin = false;

                if constexpr (proposal_acceptance_criteria == Metropolis) {
                    if (delta_energy[var] <= 0.0) {
                        flip_spin = true;
                    } else {
                        FASTRAND(rand);
                        if (exp(-delta_energy[var]*beta) * RANDMAX > rand) {
                            flip_spin = true;
                        }
                    }
                } else {
                    FASTRAND(rand);
                    if (RANDMAX > rand * (1 + exp(delta_energy[var]*beta))) {
                        flip_spin = true;
                    }
                }

                if (flip_spin) {
                    const std::int8_t multiplier = 4 * state[var];
                    for (int n_i = 0; n_i < degrees[var]; n_i++) {
                        int neighbor = neighbors[var][n_i];
                        delta_energy[neighbor] += multiplier *
                            neighbour_couplings[var][n_i] * state[neighbor];
                    }
                    state[var] *= -1;
                    delta_energy[var] *= -1;
                }
            }
        }
    }

    free(delta_energy);
}

// ---------------------------------------------------------------------------

template <VariableOrder varorder, Proposal proposal_acceptance_criteria>
void run_one_sample(
    std::int8_t* state,
    double* energy,
    const vector<double>& h,
    const vector<int>& degrees,
    const vector<vector<int>>& neighbors,
    const vector<vector<double>>& neighbour_couplings,
    const vector<int>& coupler_starts,
    const vector<int>& coupler_ends,
    const vector<double>& coupler_weights,
    const int sweeps_per_beta,
    const vector<double>& beta_schedule
) {
    fast_simulated_annealing_run<varorder, proposal_acceptance_criteria>(
        state, h, degrees, neighbors, neighbour_couplings, sweeps_per_beta, beta_schedule);

    // --- PARALLELIZED: state energy computation ---
    *energy = parallel_get_state_energy(
        state, h, coupler_starts, coupler_ends, coupler_weights);
}

// ---------------------------------------------------------------------------
// Public entry point (unchanged)
// ---------------------------------------------------------------------------

int fast_cpu_general_simulated_annealing(
    std::int8_t* states,
    double* energies,
    const int num_samples,
    const vector<double> h,
    const vector<int> coupler_starts,
    const vector<int> coupler_ends,
    const vector<double> coupler_weights,
    const int sweeps_per_beta,
    const vector<double> beta_schedule,
    const uint64_t seed,
    const VariableOrder varorder,
    const Proposal proposal_acceptance_criteria,
    callback interrupt_callback,
    void * const interrupt_function
) {
    const int num_vars = h.size();
    if (!((coupler_starts.size() == coupler_ends.size()) &&
                (coupler_starts.size() == coupler_weights.size()))) {
        throw runtime_error("coupler vectors have mismatched lengths");
    }

    vector<int> degrees(num_vars, 0);
    vector<vector<int>> neighbors(num_vars);
    vector<vector<double>> neighbour_couplings(num_vars);

    for (unsigned int cplr = 0; cplr < coupler_starts.size(); cplr++) {
        int u = coupler_starts[cplr];
        int v = coupler_ends[cplr];

        if ((u < 0) || (v < 0) || (u >= num_vars) || (v >= num_vars)) {
            throw runtime_error("coupler indexes contain an invalid variable");
        }

        neighbors[u].push_back(v);
        neighbors[v].push_back(u);
        neighbour_couplings[u].push_back(coupler_weights[cplr]);
        neighbour_couplings[v].push_back(coupler_weights[cplr]);

        degrees[u]++;
        degrees[v]++;
    }

    if (interrupt_function != nullptr) return -2;

    auto run_sample = [&](int sample) {
        fast_rng_state[0] = seed ? (seed + static_cast<uint64_t>(sample) + 1) : RANDMAX;
        fast_rng_state[1] = static_cast<uint64_t>(sample) ^ 0x9E3779B97F4A7C15ULL;
        std::int8_t* state = states + sample * num_vars;

        if (varorder == Random) {
            if (proposal_acceptance_criteria == Metropolis) {
                run_one_sample<Random, Metropolis>(
                    state, energies + sample, h, degrees, neighbors, neighbour_couplings,
                    coupler_starts, coupler_ends, coupler_weights, sweeps_per_beta, beta_schedule);
            } else {
                run_one_sample<Random, Gibbs>(
                    state, energies + sample, h, degrees, neighbors, neighbour_couplings,
                    coupler_starts, coupler_ends, coupler_weights, sweeps_per_beta, beta_schedule);
            }
        } else {
            if (proposal_acceptance_criteria == Metropolis) {
                run_one_sample<Sequential, Metropolis>(
                    state, energies + sample, h, degrees, neighbors, neighbour_couplings,
                    coupler_starts, coupler_ends, coupler_weights, sweeps_per_beta, beta_schedule);
            } else {
                run_one_sample<Sequential, Gibbs>(
                    state, energies + sample, h, degrees, neighbors, neighbour_couplings,
                    coupler_starts, coupler_ends, coupler_weights, sweeps_per_beta, beta_schedule);
            }
        }
    };

#ifdef _OPENMP
    if (num_samples > OPENMP_THRESHOLD_SAMPLES) {
        #pragma omp parallel for schedule(static)
        for (int sample = 0; sample < num_samples; ++sample) {
            run_sample(sample);
        }
        return num_samples;
    }
#endif

    int thread_count = static_cast<int>(std::thread::hardware_concurrency());
    if (thread_count < 1) thread_count = 1;
    if (thread_count > num_samples) thread_count = num_samples;

    if (thread_count <= 1) {
        for (int sample = 0; sample < num_samples; ++sample) run_sample(sample);
    } else {
        std::atomic<int> next_sample{0};
        vector<std::thread> workers;
        workers.reserve(thread_count);
        for (int t = 0; t < thread_count; t++) {
            workers.emplace_back([&]() {
                while (true) {
                    int sample = next_sample.fetch_add(1);
                    if (sample >= num_samples) break;
                    run_sample(sample);
                }
            });
        }
        for (auto& worker : workers) worker.join();
    }

    return num_samples;
}
