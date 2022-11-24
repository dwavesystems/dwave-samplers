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

#include <cstdint>
#include <vector>
#include <set>
#include <cassert>
#include <stdexcept>
#include "descent.h"

using std::vector;
using std::set;
using std::runtime_error;


// Returns the energy delta from flipping a SPIN variable at index `var`.
//
// @param var the index of the variable to flip
// @param state the current state of all variables
// @param linear_biases vector of h or field value on each variable
// @param neighbors lists of the neighbors of each variable, such that
//     neighbors[i][j] is the jth neighbor of variable i.
// @param neighbour_couplings same as neighbors, but instead has the J value.
//     neighbour_couplings[i][j] is the J value or weight on the coupling
//     between variables i and neighbors[i][j]. 
//
// @return delta energy
double get_flip_energy(
    int var,
    std::int8_t *state,
    const vector<double>& linear_biases,
    const vector<vector<int>>& neighbors,
    const vector<vector<double>>& neighbour_couplings
) {
    // var -1 flips to +1 => delta is +2
    // var +1 flips to -1 => delta is -2
    double delta = -2 * state[var];

    // flip energy = delta * h[var]
    //             + sum_over_var_neighbors(delta * J[var][neigh] * state[neigh]))
    double contrib = linear_biases[var];
    for (int idx = 0; idx < neighbors[var].size(); idx++) {
        contrib += state[neighbors[var][idx]] * neighbour_couplings[var][idx];
    }

    return delta * contrib;
}


// Returns the energy of a given state for the input problem.
//
// @param state a int8 array containing the spin state to compute the energy of
// @param linear_biases vector of h or field value on each variable
// @param coupler_starts an int vector containing the variables of one side of
//        each coupler in the problem
// @param coupler_ends an int vector containing the variables of the other side 
//        of each coupler in the problem
// @param coupler_weights a double vector containing the weights of the 
//        couplers in the same order as coupler_starts and coupler_ends
//
// @return A double corresponding to the energy for `state` on the problem
//        defined by linear_biases and the couplers passed in
double get_state_energy(
    std::int8_t* state,
    const vector<double>& linear_biases,
    const vector<int>& coupler_starts,
    const vector<int>& coupler_ends,
    const vector<double>& coupler_weights
) {
    double energy = 0.0;

    // sum the energy due to local fields on variables
    for (unsigned int var = 0; var < linear_biases.size(); var++) {
        energy += state[var] * linear_biases[var];
    }

    // sum the energy due to coupling weights
    for (unsigned int c = 0; c < coupler_starts.size(); c++) {
        energy += state[coupler_starts[c]] * coupler_weights[c] * state[coupler_ends[c]];
    }

    return energy;
}


// One run of the steepest gradient descent on the input Ising model.
//
// Linear search for the steepest descent variable. Fastest approach for
// complete and/or dense problem graphs.
//
// @param state a int8 array where each int8 holds the state of a
//        variable. Note that this will be used as the initial state of the
//        run.
// @param linear_biases vector of h or field value on each variable
// @param neighbors lists of the neighbors of each variable, such that
//        neighbors[i][j] is the jth neighbor of variable i. Note
// @param neighbour_couplings same as neighbors, but instead has the J value.
//        neighbour_couplings[i][j] is the J value or weight on the coupling
//        between variables i and neighbors[i][j].
// @param flip_energies vector used for caching of variable flip delta
//        energies
//
// @return number of downhill steps; `state` contains the result of the run.
unsigned int steepest_gradient_descent_solver(
    std::int8_t* state,
    const vector<double>& linear_biases,
    const vector<vector<int>>& neighbors,
    const vector<vector<double>>& neighbour_couplings,
    vector<double>& flip_energies
) {
    const int num_vars = linear_biases.size();

    // short-circuit on empty models
    if (num_vars < 1) {
        return 0;
    }

    // calculate flip energies for all variables, based on the current
    // state (loop invariant) ~ O(num_vars * max_degree)
    for (int var = 0; var < num_vars; var++) {
        flip_energies[var] = get_flip_energy(
            var, state, linear_biases, neighbors, neighbour_couplings
        );
    }

    // descend ~ O(downhill_steps * num_vars)
    unsigned int steps = 0;
    while (true) {

        // calculate the gradient: on binary models this translates to finding
        // a dimension with the greatest flip energy
        int best_var = -1;
        double best_flip_energy = 0;

        // find the variable flipping of which results with the steepest
        // descent in energy landscape
        for (int var = 0; var < num_vars; var++) {
            double flip_energy = flip_energies[var];

            if (flip_energy < best_flip_energy) {
                best_flip_energy = flip_energy;
                best_var = var;
            }
        }

        // are we in a minimum already?
        if (best_var == -1) {
            break;
        }

        // otherwise, we can improve the solution by descending down the
        // `best_var` dimension
        state[best_var] *= -1;

        // but to maintain the `flip_energies` invariant (after flipping
        // `best_var`), we need to update flip energies for the flipped var and
        // all neighbors of the flipped var
        flip_energies[best_var] *= -1;

        for (int n_idx = 0; n_idx < neighbors[best_var].size(); n_idx++) {
            int n_var = neighbors[best_var][n_idx];
            double w = neighbour_couplings[best_var][n_idx];
            // flip energy for each `neighbor` includes the
            // `2 * state[neighbor] * coupling weight * state[best_var]` term.
            // the change of the flip energy due to flipping `best_var` is
            // twice that, hence the factor 4 below
            flip_energies[n_var] -= 4 * state[best_var] * w * state[n_var];
        }

        steps++;
    }

    return steps;
}


struct EnergyVar {
    double energy;
    int var;
};

struct EnergyVarCmp {
    bool operator()(const EnergyVar& lhs, const EnergyVar& rhs) const {
        return lhs.energy < rhs.energy || (lhs.energy <= rhs.energy && lhs.var < rhs.var);
    }
};


// One run of the steepest gradient descent on the input Ising model.
//
// Flip energies are kept in an ordered set (balanced binary tree) resulting in
// faster steepest descent variable lookup, but higher constant overhead of
// maintaining the order. Scaling advantage only for very *large* (and *sparse*)
// problem graphs.
//
// @param state a int8 array where each int8 holds the state of a
//        variable. Note that this will be used as the initial state of the
//        run.
// @param linear_biases vector of h or field value on each variable
// @param neighbors lists of the neighbors of each variable, such that
//        neighbors[i][j] is the jth neighbor of variable i. Note
// @param neighbour_couplings same as neighbors, but instead has the J value.
//        neighbour_couplings[i][j] is the J value or weight on the coupling
//        between variables i and neighbors[i][j].
// @param flip_energies_vector vector used for caching of variable flip delta
//        energies
//
// @return number of downhill steps; `state` contains the result of the run.
unsigned int steepest_gradient_descent_ls_solver(
    std::int8_t* state,
    const vector<double>& linear_biases,
    const vector<vector<int>>& neighbors,
    const vector<vector<double>>& neighbour_couplings,
    vector<double>& flip_energies_vector
) {
    const int num_vars = linear_biases.size();

    // short-circuit on empty models
    if (num_vars < 1) {
        return 0;
    }

    // calculate flip energies for all variables, based on the current
    // state (loop invariant); store them in:
    // (1) vector for fast var-to-energy look-up
    // (2) ordered set (rb-tree) for fast best var energy look-up
    // ~ O(N * (max_degree + 1 + logN))
    // => ~ O(N^2) for dense graphs, ~ O(N*logN) for sparse
    set<EnergyVar, EnergyVarCmp> flip_energies_set;

    for (int var = 0; var < num_vars; var++) {
        double energy = get_flip_energy(
            var, state, linear_biases, neighbors, neighbour_couplings
        );
        flip_energies_vector[var] = energy;
        flip_energies_set.insert({energy, var});
    }

    // descend ~ O(downhill_steps * max_degree * logN)
    unsigned int steps = 0;
    while (true) {
        // find the variable flipping of which results with the steepest
        // descent in energy landscape ~ O(1)
        auto best_energy_var_iter = flip_energies_set.begin();
        int best_var = best_energy_var_iter->var;
        double best_energy = best_energy_var_iter->energy;

        // are we in a minimum already?
        if (best_energy >= 0) {
            break;
        }

        // otherwise, we can improve the solution by descending down the
        // `best_var` dimension

        // but to maintain the `flip_energies` invariant (after flipping
        // `best_var`), we need to update flip energies for the flipped var and
        // all neighbors of the flipped var

        // update flip energies (and their ordered set) of all `best_var`'s
        // neighbors ~ O(max_degree * logN)
        for (int n_idx = 0; n_idx < neighbors[best_var].size(); n_idx++) {
            int n_var = neighbors[best_var][n_idx];
            double w = neighbour_couplings[best_var][n_idx];
            // flip energy for each `neighbor` includes the
            // `2 * state[neighbor] * coupling weight * state[best_var]` term.
            // the change of the flip energy due to flipping `best_var` is
            // twice that, hence the factor 4 below
            double n_energy_inc = 4 * state[best_var] * w * state[n_var];

            // to update the neighbor in our ordered set:
            // 1) remove the neighbor from the set
            double n_energy = flip_energies_vector[n_var];
            auto search = flip_energies_set.find({n_energy, n_var});
            assert(search != flip_energies_set.end());
            flip_energies_set.erase(search);

            // 2) insert new (energy, var) element to reflect changed flip energy
            n_energy += n_energy_inc;
            flip_energies_set.insert({n_energy, n_var});

            // also update energy in the vector
            flip_energies_vector[n_var] = n_energy;
        }

        // finally, descend down the `var` dim (flip it)
        state[best_var] *= -1;

        // and update flip energy for the flipped best_var, in both vector and set
        best_energy *= -1;
        flip_energies_vector[best_var] = best_energy;
        flip_energies_set.erase(best_energy_var_iter);
        flip_energies_set.insert({best_energy, best_var});

        steps++;
    }

    return steps;
}


// Perform `num_samples` runs of steepest gradient descent on a general problem.
//
// @param states int8 array of size num_samples * number of variables in the
//        problem. Will be overwritten by this function as samples are filled
//        in. The initial state of the samples are used to seed the gradient
//        descent runs.
// @param energies a double array of size num_samples. Will be overwritten by
//        this function as energies are filled in.
// @param steps an unsigned int array of size num_samples. Will be overwritten
//        by this function as number of downhill steps per sample are received.
// @param num_samples the number of samples to get
// @param linear_biases vector of linear bias or field value on each variable
// @param coupler_starts an int vector containing the variables of one side of
//        each coupler in the problem
// @param coupler_ends an int vector containing the variables of the other side 
//        of each coupler in the problem
// @param coupler_weights a double vector containing the weights of the couplers
//        in the same order as coupler_starts and coupler_ends
// @param large_sparse_opt boolean that determines
//        if linear search or balanced tree search should be used for descent.
//
// @return Nothing. Results are in `states` buffer.
void steepest_gradient_descent(
    std::int8_t* states,
    double* energies,
    unsigned* num_steps,
    const int num_samples,
    const vector<double>& linear_biases,
    const vector<int>& coupler_starts,
    const vector<int>& coupler_ends,
    const vector<double>& coupler_weights,
    bool large_sparse_opt
) {
    // the number of variables in the problem
    const int num_vars = linear_biases.size();
    if (coupler_starts.size() != coupler_ends.size() ||
        coupler_starts.size() != coupler_weights.size()
    ) {
        throw runtime_error("coupler vectors have mismatched lengths");
    }
    
    // neighbors is a vector of vectors, such that neighbors[i][j] is the jth
    // neighbor of variable i
    vector<vector<int>> neighbors(num_vars);
    // neighbour_couplings is another vector of vectors with the same structure
    // except neighbour_couplings[i][j] is the weight on the coupling between i
    // and its jth neighbor
    vector<vector<double>> neighbour_couplings(num_vars);

    // build the neighbors, and neighbour_couplings vectors by iterating over
    // the input coupler vectors
    for (unsigned coupler = 0; coupler < coupler_starts.size(); coupler++) {
        int u = coupler_starts[coupler];
        int v = coupler_ends[coupler];

        if (u < 0 || v < 0 || u >= num_vars || v >= num_vars) {
            throw runtime_error("coupler indexes contain an invalid variable");
        }

        // add v to u's neighbors list and vice versa
        neighbors[u].push_back(v);
        neighbors[v].push_back(u);
        // add the weights
        neighbour_couplings[u].push_back(coupler_weights[coupler]);
        neighbour_couplings[v].push_back(coupler_weights[coupler]);
    }

    // variable flip energies cache
    vector<double> flip_energies_vector(num_vars);

    // run the steepest descent for `num_samples` times,
    // each time seeded with the initial state from `states`
    for (int sample = 0; sample < num_samples; sample++) {
        // get initial state from states buffer; the solution overwrites the same buffer
        std::int8_t *state = states + sample * num_vars;

        if (large_sparse_opt) {
            num_steps[sample] = steepest_gradient_descent_ls_solver(
                state, linear_biases, neighbors, neighbour_couplings, flip_energies_vector
            );
        } else {
            num_steps[sample] = steepest_gradient_descent_solver(
                state, linear_biases, neighbors, neighbour_couplings, flip_energies_vector
            );
        }

        // compute the energy of the sample
        energies[sample] = get_state_energy(
            state, linear_biases, coupler_starts, coupler_ends, coupler_weights
        );
    }
}
