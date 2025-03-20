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

#include <cstdint>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <limits>
#include "cpu_rotormc.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// xorshift128+ as defined https://en.wikipedia.org/wiki/Xorshift#xorshift.2B
#define FASTRAND(rand) do {                       \
    uint64_t x = rng_state[0];                    \
    uint64_t const y = rng_state[1];              \
    rng_state[0] = y;                             \
    x ^= x << 23;                                 \
    rng_state[1] = x ^ y ^ (x >> 17) ^ (y >> 26); \
    rand = rng_state[1] + y;                      \
} while (0)

#define RANDMAX ((uint64_t)-1L)

using namespace std;

// this holds the state of our thread-safe/local RNG
thread_local uint64_t rng_state[2];

// Returns the energy delta from flipping variable at index `var`
// @param var the index of the variable to flip
// @param state the current state of all variables
// @param h vector of h or field value on each variable
// @param degrees the degree of each variable
// @param neighbors lists of the neighbors of each variable, such that 
//     neighbors[i][j] is the jth neighbor of variable i.
// @param neighbour_couplings same as neighbors, but instead has the J value.
//     neighbour_couplings[i][j] is the J value or weight on the coupling
//     between variables i and neighbors[i][j]. 
// @return delta energy
double get_marginal_state_fieldC(
    int var,
    std::uint8_t *state,
    const vector<double>& h,
    const vector<int>& degrees,
    const vector<vector<int>>& neighbors,
    const vector<vector<double>>& neighbour_couplings,
    const vector<double>& state_to_costheta
) {
    double energy = h[var];
    // iterate over the neighbors of variable `var`
    for (int n_i = 0; n_i < degrees[var]; n_i++) {
        // increase `energy` by the state of the neighbor variable * the
        // corresponding coupler weight
        energy += state_to_costheta[state[neighbors[var][n_i]]] * neighbour_couplings[var][n_i];
    }
    // the value of the variable `energy` is now equal to the sum of the
    // coefficients of `var`.  we then multiply this by -2 * the state of `var`
    // because the energy delta is given by: (x_i_new - x_i_old) * sum(coefs),
    // and (x_i_new - x_i_old) = -2 * x_i_old
    return energy;
}

inline double get_marginal_state_energyQ(int var,
                                         std::uint8_t *state,
                                         const vector<double>& trans_fields,
                                         const vector<double>& state_to_sintheta) {
  return(trans_fields[var]*state_to_sintheta[state[var]]);
}

// Performs a single run of simulated annealing with the given inputs.
// @param state a int8 array where each int8 holds the state of a
//        variable. Note that this will be used as the initial state of the
//        run.
// @param h vector of h or field value on each variable
// @param degrees the degree of each variable
// @param neighbors lists of the neighbors of each variable, such that 
//        neighbors[i][j] is the jth neighbor of variable i. Note
// @param neighbour_couplings same as neighbors, but instead has the J value.
//        neighbour_couplings[i][j] is the J value or weight on the coupling
//        between variables i and neighbors[i][j]. 
// @param sweeps_per_beta The number of sweeps to perform at each beta value.
//        Total number of sweeps is `sweeps_per_beta` * length of
//        `Hp_field`.
// @param Hp_field A list of the beta values to run `sweeps_per_beta`
//        sweeps at.
// @return Nothing, but `state` now contains the result of the run.
template <bool randomize_order, Proposal proposal>
void simulated_annealing_run(
    std::uint8_t* state,
    const vector<double>& h,
    const vector<int>& degrees,
    const vector<vector<int>>& neighbors,
    const vector<vector<double>>& neighbour_couplings,
    const vector<double>& trans_fields,
    const int sweeps_per_beta,
    const vector<double>& Hp_field,
    const vector<double>& Hd_field,
    const vector<double>& state_to_costheta,
    const vector<double>& state_to_sintheta,
    std::uint8_t* statistics,
    const int schedule_sample_interval
) {
    const int num_vars = h.size();
    // unint8 is interpreted as an angle: theta=uint8*PI/max_abs_size, max_abs_size=256:
    
    // this double array will hold the delta energy for every variable
    // delta_energy[v] is the delta energy for variable `v`
    double *delta_energy = (double*)malloc(num_vars * sizeof(double));
    uint64_t rand; // this will hold the value of the rng

    // build the delta_energy array by getting the delta energy for each
    // variable
    // A flip is a reversal of the angle theta -> pi-theta. This is not ergodic,
    // but obeys detailed balance and can be a valid within a subspace.

    for (int var = 0; var < num_vars; var++) {
      // Local fields
      delta_energy[var] = get_marginal_state_fieldC(var, state, h, degrees,
						    neighbors, neighbour_couplings,
						    state_to_costheta);
    }
    if constexpr (proposal==GibbsNonErgodic || proposal==MetropolisNonErgodic) {
        for (int var = 0; var < num_vars; var++) {
            delta_energy[var] *= -2*state_to_costheta[state[var]];
        }
    }
    for(int i=0; i<32; i++)
      // Burn in .. found to be pathological.
      FASTRAND(rand);
    
    const int discretization = 256;
    const int half_discretization = discretization/2;
    bool flip_spin;
    // perform the sweeps
    for (int beta_idx = 0; beta_idx < (int)Hp_field.size(); beta_idx++) {
        // get the beta value for this sweep
        const double beta = Hp_field[beta_idx];
        const double beta_tf = Hd_field[beta_idx];
	const double rat = beta_tf/beta;
        for (int sweep = 0; sweep < sweeps_per_beta; sweep++) {
            for (int varI = 0; varI < num_vars; varI++) {
                int var;
                if (randomize_order) { // constexpr, get C++17 linking correct! 
                    FASTRAND(rand);
                    var = rand%num_vars;
                } else {
                    var = varI;
                }
                uint8_t proposed_angle;
                double proposed_energy;
                flip_spin = false;

                switch (proposal) {
                case MetropolisNonErgodic:
                    // For uint8_t can use wrap around: anti-polar point
                    // is equivalent to flip in Z-basis in current scheme.
                    proposed_angle = uint8_t(128) + state[var];
                    proposed_energy = - delta_energy[var]; //delta
                    // Metropolis-Hastings acceptance rule
                    if (delta_energy[var] <= 0.0) {
                        // automatically accept any flip that results in a lower
                        // energy
                        flip_spin = true;
                    } else {
                        // get a random number, storing it in rand
                        FASTRAND(rand);
                        // accept the flip if exp(-delta_energy*beta) > random(0, 1)
                        if (exp(-delta_energy[var]*beta) * RANDMAX > rand) {
                            flip_spin = true;
                        }
                    }
                    break;
                case GibbsNonErgodic:
                    // For uint8_t can use wrap around: anti-polar point
                    // is equivalent to flip in Z-basis in current scheme.
                    proposed_angle = uint8_t(128) + state[var];
                    proposed_energy = - delta_energy[var]; //delta
                    // Gibbs update: Sample fairly from the two available
                    // states, independent of the current value
                    FASTRAND(rand);
                    if (RANDMAX > rand * (1+exp(delta_energy[var]*beta))) {
                        flip_spin = true;
                    }
                    break;
                case MetropolisUniform:
                    FASTRAND(rand);
                    proposed_angle = uint8_t(rand%discretization);
                    break;
                case MetropolisTF:
		    // Routine makes no sense, should be tf*rat/classical_gap? Should
		    // depend on {beta_f Hd_i} versus {beta Hp_i}, so makes qualitative
		    // sense only in special O(1) cases of Hd and Hp. 
                    FASTRAND(rand);
                    if( rat >= 1 ) {
                      proposed_angle = uint8_t(rand%discretization);
                    } else {
		      // Match properties of continuous process:
		      // (1) Match expected distance at zero measure (i.e. non-zero move prob.)
		      // (2) uniform p over fully covered intervals.

		      //Number of complete intervals one-side:
		      //r 1: 3/256<, 2: 5/256<, 3: 7/256<, 127: 255/256< 
		      int osi = int(rat*half_discretization-0.5); // NB:Rounds toward 0 small rat.
		      //Mean distance from fair sampling complete intervals:
		      double distanceC = (osi*(osi+1))/double(2*osi+1);
		      // (rat*half_discretization)/2 = p*(osi+1) + (1-p)*distanceC

		      // Probability to sample a boundary state +/- (osi+1)
		      double p = (rat/2*half_discretization - distanceC)/(osi + 1 - distanceC);
		      
		      /*std::assert(p<1);
		      std::assert(p>=0);
		      std::assert(osi<=half_discretization);
		      std::assert(osi>=0);*/
		      auto randL = rand;
		      FASTRAND(rand); // Perhaps safe to reuse high precision bits, but why risk it.. 
		      uint8_t pac = (p*RANDMAX > randL) ?
			(osi + 1)*(2*(rand%2) - 1):
			rand%(2*osi + 1) - osi;
		      proposed_angle = state[var] + pac;
		    }
                    break;
                default:
                    throw std::invalid_argument( "invalid proposal method" );
                    break;
                }
                if constexpr (proposal==MetropolisUniform || proposal==MetropolisTF) {
                    
                    // Rewrite later so this needn't be done every time, i.e. update neighbor local fields
		    proposed_energy = get_marginal_state_fieldC(var, state, h, degrees,
								neighbors, neighbour_couplings,
								state_to_costheta); // Local field
                    double delta_logmeasure = - beta*proposed_energy*(state_to_costheta[proposed_angle] - state_to_costheta[state[var]]) +
                        beta_tf*trans_fields[var]*(state_to_sintheta[proposed_angle] -
						   state_to_sintheta[state[var]]);
                    // Metropolis-Hastings acceptance rule
                    if (delta_logmeasure >= 0.0) {
                        // automatically accept any flip that results in a lower
                        // energy
                        flip_spin = true;
                    } else {
                        // get a random number, storing it in rand
                        FASTRAND(rand);
                        // accept the flip if exp(-delta_energy*beta) > random(0, 1)
                        if (exp(delta_logmeasure) * RANDMAX > rand) {
                            flip_spin = true;
                        }
                    }
                }
                if (flip_spin) {
                    // since we have accepted the spin flip of variable `var`, 
                    // we need to adjust the delta energies of all the 
                    // neighboring variables
		    if constexpr (proposal==MetropolisNonErgodic || proposal==GibbsNonErgodic) {
			// TO DO: use this loop for TF and Uniform routines
			// to reduce complexity for dense graphs. Watch for factor 2!
                        const double multiplier = 2 * (state_to_costheta[state[var]] -
                                                   state_to_costheta[proposed_angle]); 
                        // iterate over the neighbors of `var`
                        for (int n_i = 0; n_i < degrees[var]; n_i++) {
                            int neighbor = neighbors[var][n_i];
			    // Only classical coupling part impacts the neighboring spin:
			    delta_energy[neighbor] += multiplier * 
			        neighbour_couplings[var][n_i] * state_to_costheta[state[neighbor]];
			}
			delta_energy[var] = proposed_energy;
		    }
                    state[var] = proposed_angle;
                }
            }
        }
	if (schedule_sample_interval!=0 && beta_idx%schedule_sample_interval==0) {
	  for(int i=0; i<num_vars; i++) {
	    *statistics++ = state[i];
	  }
	}
    }
    free(delta_energy);
}

// Returns the energy of a given state and problem
// @param state a int8 array containing the spin state to compute the energy of
// @param h vector of h or field value on each variable
// @param coupler_starts an int vector containing the variables of one side of
//        each coupler in the problem
// @param coupler_ends an int vector containing the variables of the other side 
//        of each coupler in the problem
// @param coupler_weights a double vector containing the weights of the 
//        couplers in the same order as coupler_starts and coupler_ends
// @return A double corresponding to the energy for `state` on the problem
//        defined by h and the couplers passed in
double get_state_energyC(
    std::uint8_t* state,
    const vector<double>& h,
    const vector<int>& coupler_starts,
    const vector<int>& coupler_ends,
    const vector<double>& coupler_weights,
    const vector<double>& state_to_costheta
) {
    double energy = 0.0;
    // sum the energy due to local fields on variables
    for (unsigned int var = 0; var < h.size(); var++) {
        energy += state_to_costheta[state[var]] * h[var];
    }
    // sum the energy due to coupling weights
    for (unsigned int c = 0; c < coupler_starts.size(); c++) {
        energy += state_to_costheta[state[coupler_starts[c]]] * coupler_weights[c] *
          state_to_costheta[state[coupler_ends[c]]];
    }
    return energy;
}

double get_state_energyQ(std::uint8_t* state,
                         const vector<double>& trans_fields,
                         const vector<double>& state_to_sintheta) {
  double en = 0;
  for(int var = 0; var < trans_fields.size(); var++) {
    en += trans_fields[var]*state_to_sintheta[state[var]];
  }
  return(en);
}

// Perform simulated annealing on a general problem
// @param states a int8 array of size num_samples * number of variables in the
//        problem. Will be overwritten by this function as samples are filled
//        in. The initial state of the samples are used to seed the simulated
//        annealing runs.
// @param energies a double array of size num_samples. Will be overwritten by
//        this function as energies are filled in.
// @param num_samples the number of samples to get.
// @param h vector of h or field value on each variable
// @param coupler_starts an int vector containing the variables of one side of
//        each coupler in the problem
// @param coupler_ends an int vector containing the variables of the other side 
//        of each coupler in the problem
// @param coupler_weights a double vector containing the weights of the couplers
//        in the same order as coupler_starts and coupler_ends
// @param sweeps_per_beta The number of sweeps to perform at each beta value.
//        Total number of sweeps is `sweeps_per_beta` * length of
//        `Hp_field`.
// @param Hp_field A list of the beta values to run `sweeps_per_beta`
//        sweeps at.
// @param interrupt_callback A function that is invoked between each run of simulated annealing
//        if the function returns True then it will stop running.
// @param interrupt_function A pointer to contents that are passed to interrupt_callback.
// @return the number of samples taken. If no interrupt occured, will equal num_samples.

int general_simulated_annealing(
    std::uint8_t* states,
    double* energies,
    const int num_samples,
    const vector<double> h,
    const vector<int> coupler_starts,
    const vector<int> coupler_ends,
    const vector<double> coupler_weights,
    const vector<double> trans_fields,
    const int sweeps_per_beta,
    const vector<double> Hp_field,
    const vector<double> Hd_field,
    const uint64_t seed,
    const bool randomize_order,
    const Proposal proposal,
    std::uint8_t* statistics,
    const int schedule_sample_interval,
    callback interrupt_callback,
    void * const interrupt_function
) {
    const int discretization = 256;
    std::vector<double> state_to_costheta(discretization),
        state_to_sintheta(discretization);
    // Discretization of uniform measure over [0,pi],
    // NB: measure choice here is canonical, but heuristic for purposes of
    // physical emulation of dynamics/equilibirum
    // (2-angle uniform measure -> 1 angle measure has several options)
    for(auto index = 0; index < discretization; index++){
        state_to_costheta[index] = std::cos(2*index*M_PI/discretization);
        state_to_sintheta[index] = std::abs(std::sin(2*index*M_PI/discretization)); //theta bound on [0,pi], theta in (pi,2pi) interpreted by reflection (this gives uniform density approximation over [0,pi] - inclusive of interval boundaries.
    }
    
    // the number of variables in the problem
    const int num_vars = h.size();
    if (!((coupler_starts.size() == coupler_ends.size()) &&
                (coupler_starts.size() == coupler_weights.size()))) {
        throw runtime_error("coupler vectors have mismatched lengths");
    }
    
    // set the seed of the RNG
    // note that xorshift+ requires a non-zero seed
    rng_state[0] = seed ? seed : RANDMAX;
    rng_state[1] = 0;

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
    // iterating over the inputted coupler vectors
    for (unsigned int cplr = 0; cplr < coupler_starts.size(); cplr++) {
        int u = coupler_starts[cplr];
        int v = coupler_ends[cplr];

        if ((u < 0) || (v < 0) || (u >= num_vars) || (v >= num_vars)) {
            throw runtime_error("coupler indexes contain an invalid variable");
        }

        // add v to u's neighbors list and vice versa
        neighbors[u].push_back(v);
        neighbors[v].push_back(u);
        // add the weights
        neighbour_couplings[u].push_back(coupler_weights[cplr]);
        neighbour_couplings[v].push_back(coupler_weights[cplr]);

        // increase the degrees of both variables
        degrees[u]++;
        degrees[v]++;
    }


    // get the simulated annealing samples
    int sample = 0;
    int stat_block = schedule_sample_interval ? num_vars * (Hp_field.size()/schedule_sample_interval) : 0;
    while (sample < num_samples) {
        // states is a giant spin array that will hold the resulting states for
        // all the samples, so we need to get the location inside that vector
        // where we will store the sample for this sample
        std::uint8_t *state = states + sample*num_vars;
        // then do the actual sample. this function will modify state, storing
        // the sample there
        // Branching here is designed to make expicit compile time optimizations
        if (randomize_order) {
            switch(proposal) {
            case MetropolisNonErgodic:
                simulated_annealing_run<true, MetropolisNonErgodic>(
                    state, h, degrees,
		    neighbors, neighbour_couplings,
		    trans_fields,
		    sweeps_per_beta, Hp_field, Hd_field,
		    state_to_costheta, state_to_sintheta,
		    statistics + stat_block * sample,
		    schedule_sample_interval);
                break;
            case GibbsNonErgodic:
                simulated_annealing_run<true, GibbsNonErgodic>(
	            state, h, degrees,
		    neighbors, neighbour_couplings,
		    trans_fields,
		    sweeps_per_beta, Hp_field, Hd_field,
		    state_to_costheta, state_to_sintheta,
		    statistics + stat_block * sample,
		    schedule_sample_interval);
                break;
            case MetropolisUniform:
                simulated_annealing_run<true, MetropolisUniform>(
		    state, h, degrees,
		    neighbors, neighbour_couplings,
		    trans_fields,
		    sweeps_per_beta, Hp_field, Hd_field,
		    state_to_costheta, state_to_sintheta,
		    statistics + stat_block * sample,
		    schedule_sample_interval);
                break;
            case MetropolisTF:
                simulated_annealing_run<true, MetropolisTF>(
		    state, h, degrees,
		    neighbors, neighbour_couplings,
		    trans_fields,
		    sweeps_per_beta, Hp_field, Hd_field,
		    state_to_costheta, state_to_sintheta,
		    statistics + stat_block * sample,
		    schedule_sample_interval);
                break;
            default:
                throw std::invalid_argument( "invalid proposal method" );
		break;
            }
        }
        else {
            switch(proposal) {
            case MetropolisNonErgodic:
                simulated_annealing_run<false, MetropolisNonErgodic>(
		    state, h, degrees,
		    neighbors, neighbour_couplings,
		    trans_fields,
		    sweeps_per_beta, Hp_field, Hd_field,
		    state_to_costheta, state_to_sintheta,
		    statistics + stat_block * sample,
		    schedule_sample_interval);
                break;
            case GibbsNonErgodic:
                simulated_annealing_run<false, GibbsNonErgodic>(
		    state, h, degrees,
		    neighbors, neighbour_couplings,
		    trans_fields,
		    sweeps_per_beta, Hp_field, Hd_field,
		    state_to_costheta, state_to_sintheta,
		    statistics + stat_block * sample,
		    schedule_sample_interval);
                break;
            case MetropolisUniform:
                simulated_annealing_run<false, MetropolisUniform>(
		    state, h, degrees,
		    neighbors, neighbour_couplings,
		    trans_fields,
		    sweeps_per_beta, Hp_field, Hd_field,
		    state_to_costheta, state_to_sintheta,
		    statistics + stat_block * sample,
		    schedule_sample_interval);
                break;
            case MetropolisTF:
                simulated_annealing_run<false, MetropolisTF>(
		    state, h, degrees,
		    neighbors, neighbour_couplings,
		    trans_fields,
		    sweeps_per_beta, Hp_field, Hd_field,
		    state_to_costheta, state_to_sintheta,
		    statistics + stat_block * sample,
		    schedule_sample_interval);
                break;
            default:
                throw std::invalid_argument( "invalid proposal method" );
                break;
            }
        }
        // compute the energy of the sample and store it in `energies`
        energies[sample] = get_state_energyC(state, h, coupler_starts,
                                             coupler_ends, coupler_weights,
                                             state_to_costheta);

        // Classical only (ambiguous weighting Hd and Hp), interpret as average.
        sample++;

        // if interrupt_function returns true, stop sampling
        if (interrupt_function && interrupt_callback(interrupt_function)) break;
    }

    // return the number of samples we actually took
    return sample;
}
