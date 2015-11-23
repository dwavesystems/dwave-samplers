
#ifndef __ORANG_C_INTERFACE_H__
#define __ORANG_C_INTERFACE_H__

extern "C"{

extern const int MAX_ERROR_LENGTH;

// Heuristics for finding elimination order
enum {
    HEURISTIC_MIN_DEG    = 0,
    HEURISTIC_W_MIN_DEG  = 1,
    HEURISTIC_MIN_FILL   = 2,
    HEURISTIC_W_MIN_FILL = 3
};

struct TableEntry{
    int num_vars;
    // The list of variable indecies. Must be in increasing order.
    int * vars;
    // The domain sizes of varaibles
    int * domSizes;
    // A vector with domSizes[0]*domSizes[1]*...*domSizes[num_vars] elements
    // The first domSizes[0] elements correspond to vars[0] taking all possible
    // domain values and other values staying at the first domain value.
    // Example:
    //   x_5 \in {0, 5, 7}
    //   x_8 \in {1, 2}
    //   f(x_5=0, x_8=1) = 1
    //   f(x_5=5, x_8=1) = 2
    //   f(x_5=7, x_8=1) = 3
    //   f(x_5=0, x_8=2) = 4
    //   f(x_5=5, x_8=2) = 5
    //   f(x_5=7, x_8=2) = 6
    //   Should be passed as :
    //
    //   num_vars = 2
    //   vars     = 5, 8
    //   domSizes = 3, 2
    //   values   = 1, 2, 3, 4, 5, 6
    //
    double * values;
};

// The struct of unary and pairwise marginals.
struct Marginal{
    // The length of vars vector -> either 1 (unary) or 2 (pairwise)
    int vars_len;
    // The list of variables
    int * vars;
    // The marginal values, the order is the same as described in TableEntry struct
    // The length pf values is 2 for unary marginals and 4 for pairwise marginals
    // The order of elements of values is the same as values in TableEntry
    double * values;
};


// greedyVarOrder : Find elimination order
//
//  Inputs:
//
//     tables: problem data as an array of TableEntry.
//
//     tables_len: tables length.
//
//     maxComplexity: maximum tree decomposition complexity of the returned
//       elimination order.  The complexity of a tree decomposition node is
//       log_2 of the product of the node's variables' domain sizes.  The
//       complexity of the whole tree decomposition is the maximum complexity
//       of any node.  Note that if all variables are binary (ie. domain size
//       = 2), then the complexity is the treewidth plus one.
//
//     clampRanks: variable clamping hint.  Clamped variables are those not
//       included in the elimination order.  This parameter must be a vector
//       indexed by variable.  Values indicate relative preferences for
//       clamping variables.  A negative value causes the corresponding
//       variable to be clamped immediately.  If, during the course of the
//       algorithm, a variable needs to be clamped, the choice is confined to
//       those variables with minimum clamp rank.  After clamping a variable,
//       the clamp ranks of all variables with higher clamp rank a decreased
//       by one.
//
//     clampRanks_len: clampRanks length.
//
//     heuristic: variable elimination heuristic to use.
//
//     selectionScale: scaling factor for range size when randomly
//       eliminating or clamping variables.  Normally, the variable to be
//       eliminated is chosen randomly from all variables with lowest cost
//       whose elimination don't violate the maximum complexity bound. The
//       size of that set of variables is multiplied by this value, possibly
//       shrinking it (<1) or enlarging it (>1).  The set is enlarged by
//       adding variables with low elimination cost.  The situation is similar
//       when clamping a variable, though the variable clamp ranks are always
//       respected.
//
//   Returns:
//
//     success(0) or failure(non-zero). Error message is returned in errorMessage.
//
//  Outputs:
//     variableOrder: eliminating order. Missing varaibles are clamped.
//
//     variableOrder_len: variableOrder length.
//
//     errorMessage: error message in case of error. Must have space for MAX_ERROR_LENGTH
//                   characters.

int greedyVarOrder(TableEntry * tables, int tables_len, int maxComplexity,
                   int * clampRanks, int clampRanks_len,
                   int heuristic, float selectionScale,
                   int ** variableOrder, int * variableOrder_len, char * errorMessage);

// optimize: Find lowest energy values/states
//
// Inputs:
//
//     tables: problem data as an array of TableEntry.
//
//     tables_len: tables length.
//
//     varOrder: variable elimination order.  Missing variables are clamped.
//
//     variableOrder_len: variableOrder length.
//
//     maxComplexity: maximum tree decomposition complexity of the returned
//       elimination order.  The complexity of a tree decomposition node is
//       log_2 of the product of the node's variables' domain sizes.  The
//       complexity of the whole tree decomposition is the maximum complexity
//       of any node.  Note that if all variables are binary (ie. domain size
//       = 2), then the complexity is the treewidth plus one.
//
//     maxSolutions: maximum number of solutions and objective values (unless is zero) to
//       return. If zero, the minimum objective value is returned but no solutions.
//
//     initState: initial state.  Values are domain indices (ie. in 0:N-1 for a domain
//       of size N).  Only the values of clamped variables are used.  This
//       parameter is optional.  If NULL, it defaults to all ones.
//
//     initState_len: initState length. Must be zero if initState in NULL.
//
//     minVars: by default, the number of variables is the largest variable
//       index appearing in an input table.  This parameter can be used to
//       increase the number of variables (affecting things such as output
//       sample sizes, required initState size, and variables allowed to appear in
//       varOrder).
//
//   Returns:
//
//     success(0) or failure(non-zero). Error message is returned in errorMessage.
//
//   Outputs:
//
//     energies: array of up to maxSolutions best objective values.  If
//       maxSolutions is zero, the optimal objective value is still returned.
//
//     energies_len: energies length.
//
//     states: best solutions stored as solution 1 followed by solution 2, etc.
//       Length of each solution is returned in state_len.
//       Solutions with equal objective values are sorted lexicographically.
//       Each solution entry is in the range 0:N-1, where N is the size of the corresponding
//       variable domain.
//
//     state_len: length of each state. Length of states is state_len*energies_len.
//
//     errorMessage: error message in case of error. Must have space for MAX_ERROR_LENGTH
//                   characters.

int optimize(TableEntry * tables, int tables_len, int * variableOrder, int variableOrder_len,
             int maxComplexity, int maxSolutions, int * initState, int initState_len, int minVars,
             double ** energies, int * energies_len, int ** states, int * state_len, char * errorMessage);


// sample: Sample from a Boltzmann distribution with p(s)=exp(E(s)), also return the
//         log partition function and marginals
//
// Inputs:
//
//     tables: problem data as an array of TableEntry.
//
//     tables_len: tables length.
//
//     varOrder: variable elimination order.  Missing variables are clamped.
//
//     variableOrder_len: variableOrder length.
//
//     maxComplexity: maximum tree decomposition complexity of the returned
//       elimination order.  The complexity of a tree decomposition node is
//       log_2 of the product of the node's variables' domain sizes.  The
//       complexity of the whole tree decomposition is the maximum complexity
//       of any node.  Note that if all variables are binary (ie. domain size
//       = 2), then the complexity is the treewidth plus one.
//
//     sampleNum: number of samples to be drawn.
//
//     initState: initial state.  Values are domain indices (ie. in 0:N-1 for a domain
//       of size N).  Only the values of clamped variables are used.  This
//       parameter is optional.  If NULL, it defaults to all ones.
//
//     initState_len: initState length. Must be zero if initState in NULL.
//
//     minVars: by default, the number of variables is the largest variable
//       index appearing in an input table.  This parameter can be used to
//       increase the number of variables (affecting things such as output
//       sample sizes, required initState size, and variables allowed to appear in
//       varOrder).
//
//     seed: random number generator seed. If negative, the seed is based on the
//           current date and time.
//
//     returnMarginals: If non-zero, return marginals as an array of Marginal.
//
//   Returns:
//
//     success(0) or failure(non-zero). Error message is returned in errorMessage.
//
//   Outputs:
//
//     logZ: log partition function.
//
//     samples: samples stored as sample 1 followed by sample 2, etc.
//       Each sample entry is in the range 0:N-1, where N is the size of the corresponding
//       variable domain.
//
//     sample_len: length of each sample. Length of samples is sample_len*sampleNum.
//
//     marginals: array of marginals.
//
//     marginals_len: marginals length.
//
//     errorMessage: error message in case of error. Must have space for MAX_ERROR_LENGTH
//                   characters.

int sample(TableEntry * tables, int tables_len, int * variableOrder, int variableOrder_len,
             int maxComplexity, int sampleNum, int * initState, int initState_len, int minVars,
             int seed, int returnMarginals,
             double * logZ, int ** samples, int * sample_len,
             Marginal ** marginals, int * marginals_len,
             char * errorMessage);

// Free array of variable orders returned by greedyVarOrder function.
void free_varOrder(int * varOrder);

// Free array of energies returned by optimize function.
void free_energies(double * energies);

// Free states returned by optimize(states) or sample(samples) functions.
void free_states(int * states);

// Free array of marginals returned by sample. Requires marginals_len returned by sample function.
void free_marginals(Marginal * marginals, int marginals_len);

};

#endif // __ORANG_C_INTERFACE_H__
