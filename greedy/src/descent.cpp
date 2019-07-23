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

#include <vector>
#include <stdexcept>
#include "descent.h"

using std::vector;
using std::runtime_error;


// Returns the energy delta from flipping a SPIN variable at index `var`.
//
// @param var the index of the variable to flip
// @param state the current state of all variables
// @param linear_biases vector of h or field value on each variable
// @param degrees the degree of each variable
// @param neighbors lists of the neighbors of each variable, such that
//     neighbors[i][j] is the jth neighbor of variable i.
// @param neighbour_couplings same as neighbors, but instead has the J value.
//     neighbour_couplings[i][j] is the J value or weight on the coupling
//     between variables i and neighbors[i][j]. 
//
// @return delta energy
double get_flip_energy(
    int var,
    char *state,
    const vector<double>& linear_biases,
    const vector<int>& degrees,
    const vector<vector<int>>& neighbors,
    const vector<vector<double>>& neighbour_couplings
) {
    // var -1 flips to +1 => delta is +2
    // var +1 flips to -1 => delta is -2
    double delta = -2 * state[var];

    // flip energy = delta * h[var]
    //             + sum_over_var_neighbors(delta * J[var][neigh] * state[neigh]))
    double contrib = linear_biases[var];
    for (int idx = 0; idx < degrees[var]; idx++) {
        contrib += state[neighbors[var][idx]] * neighbour_couplings[var][idx];
    }

    return delta * contrib;
}


// Returns the energy of a given state for the input problem.
//
// @param state a char array containing the spin state to compute the energy of
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
    char* state,
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
// @param state a signed char array where each char holds the state of a
//        variable. Note that this will be used as the initial state of the
//        run.
// @param linear_biases vector of h or field value on each variable
// @param degrees the degree of each variable
// @param neighbors lists of the neighbors of each variable, such that 
//        neighbors[i][j] is the jth neighbor of variable i. Note
// @param neighbour_couplings same as neighbors, but instead has the J value.
//        neighbour_couplings[i][j] is the J value or weight on the coupling
//        between variables i and neighbors[i][j]. 
//
// @return Nothing, but `state` now contains the result of the run.
void steepest_gradient_descent_solver(
    char* state,
    const vector<double>& linear_biases,
    const vector<int>& degrees,
    const vector<vector<int>>& neighbors,
    const vector<vector<double>>& neighbour_couplings
) {
    const int num_vars = linear_biases.size();

    // short-circuit on empty models
    if (num_vars < 1) {
        return;
    }

    // flip energies for all variables, based on the current state (invariant)
    vector<double> flip_energies(num_vars);
    for (int var = 0; var < num_vars; var++) {
        flip_energies[var] = get_flip_energy(
            var, state,
            linear_biases, degrees, neighbors, neighbour_couplings
        );
    }

    bool minimum_reached = false;
    while (!minimum_reached) {

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
            minimum_reached = true;
            break;
        }

        // otherwise, we can improve the solution by descending the `var` dim
        state[best_var] *= -1;

        // to maintain the `flip_energies` invariant, we need to update
        // flip energies for the flipped var and all neighbors of the flipped var
        flip_energies[best_var] *= -1;

        for (int n_idx = 0; n_idx < neighbors[best_var].size(); n_idx++) {
            int n_var = neighbors[best_var][n_idx];
            flip_energies[n_var] = get_flip_energy(
                n_var, state,
                linear_biases, degrees, neighbors, neighbour_couplings
            );
        }
    }
}


// Perform `num_samples` runs of steepest gradient descent on a general problem.
//
// @param states char array of size num_samples * number of variables in the
//        problem. Will be overwritten by this function as samples are filled
//        in. The initial state of the samples are used to seed the gradient
//        descent runs.
// @param energies a double array of size num_samples. Will be overwritten by
//        this function as energies are filled in.
// @param num_samples the number of samples to get
// @param linear_biases vector of linear bias or field value on each variable
// @param coupler_starts an int vector containing the variables of one side of
//        each coupler in the problem
// @param coupler_ends an int vector containing the variables of the other side 
//        of each coupler in the problem
// @param coupler_weights a double vector containing the weights of the couplers
//        in the same order as coupler_starts and coupler_ends
//
// @return Nothing. Results are in `states` buffer.
void steepest_gradient_descent(
    char* states,
    double* energies,
    const int num_samples,
    const vector<double>& linear_biases,
    const vector<int>& coupler_starts,
    const vector<int>& coupler_ends,
    const vector<double>& coupler_weights
) {
    // the number of variables in the problem
    const int num_vars = linear_biases.size();
    if (coupler_starts.size() != coupler_ends.size() ||
        coupler_starts.size() != coupler_weights.size()
    ) {
        throw runtime_error("coupler vectors have mismatched lengths");
    }
    
    // degrees will be a vector of the degrees of each variable
    vector<int> degrees(num_vars, 0);
    // neighbors is a vector of vectors, such that neighbors[i][j] is the jth
    // neighbor of variable i
    vector<vector<int>> neighbors(num_vars);
    // neighbour_couplings is another vector of vectors with the same structure
    // except neighbour_couplings[i][j] is the weight on the coupling between i
    // and its jth neighbor
    vector<vector<double>> neighbour_couplings(num_vars);

    // build the degrees, neighbors, and neighbour_couplings vectors by
    // iterating over the input coupler vectors
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

        // increase the degrees of both variables
        degrees[u]++;
        degrees[v]++;
    }

    // run the steepest descent for `num_samples` times, each time seeded with
    // the initial state from `states`
    for (int sample = 0; sample < num_samples; sample++) {
        // get initial state from states buffer; the solution overwrites the same buffer
        char *state = states + sample * num_vars;

        steepest_gradient_descent_solver(
            state, linear_biases, degrees, neighbors, neighbour_couplings
        );

        // compute the energy of the sample
        energies[sample] = get_state_energy(state, linear_biases, coupler_starts,
                                            coupler_ends, coupler_weights);
    }
}
