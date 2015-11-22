
#ifndef __ORANG_C_INTERFACE_H__
#define __ORANG_C_INTERFACE_H__

extern "C"{

extern const int MAX_ERROR_LENGTH;

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

int greedyVarOrder(TableEntry * tables, int tables_len, int maxComplexity,
                   int * clampRanks, int clampRanks_len,
                   int heuristic, float selectionScale,
                   int ** variableOrder, int * variableOrder_len, char * errorMessage);

int optimize(TableEntry * tables, int tables_len, int * variableOrder, int variableOrder_len,
             int maxComplexity, int maxSolutions, int * initState, int initState_len, int minVars,
             double ** energies, int * energies_len, int ** states, int * varNum, char * errorMessage);

void free_varOrder(int * varOrder);
void free_energies(double * energies);
void free_samples(int * samples);

};

#endif // __ORANG_C_INTERFACE_H__
