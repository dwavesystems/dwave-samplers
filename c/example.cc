
#include <vector>
#include <stdio.h>

#include "interface.h"

using std::vector;

#define CREATE_UNARY_ISING(tb, index, val)          \
    tb.num_vars = 1;                                \
    tb.vars = new int[1]; tb.vars[0] = index;       \
    tb.domSizes = new int[1]; tb.domSizes[0] = 2;   \
    tb.values = new double[2];                      \
    tb.values[0] = -val; tb.values[1] = val;

#define CREATE_PAIRWISE_ISING(tb, index1, index2, val)                  \
    tb.num_vars = 2;                                                    \
    tb.vars = new int[2]; tb.vars[0] = index1; tb.vars[1] = index2;     \
    tb.domSizes = new int[2]; tb.domSizes[0] = tb.domSizes[1] =  2;     \
    tb.values = new double[4];                                          \
    tb.values[0] =  val; tb.values[1] = -val;                           \
    tb.values[2] = -val; tb.values[3] = val;

void cleanup(vector<TableEntry> & tables){
    for (int i = 0; i < tables.size(); i++){
        delete [] tables[i].vars;
        delete [] tables[i].domSizes;
        delete [] tables[i].values;
    }
}

int main(){
    int varNum = 5;
    vector<TableEntry> tables(11);
    vector<int> clampRanks;
    char * errMsg;
    int * varOrder;
    int varOrderLen;

    // Unary terms
    CREATE_UNARY_ISING(tables[0], 0, 0);
    CREATE_UNARY_ISING(tables[1], 1, 0);
    CREATE_UNARY_ISING(tables[2], 2, 0);
    CREATE_UNARY_ISING(tables[3], 3, 0);
    CREATE_UNARY_ISING(tables[4], 4, 0);

    CREATE_PAIRWISE_ISING(tables[5] , 0, 1, 0);
    CREATE_PAIRWISE_ISING(tables[6] , 1, 2, 0);
    CREATE_PAIRWISE_ISING(tables[7] , 1, 3, 0);
    CREATE_PAIRWISE_ISING(tables[8] , 2, 3, 0);
    CREATE_PAIRWISE_ISING(tables[9] , 2, 4, 0);
    CREATE_PAIRWISE_ISING(tables[10], 3, 4, 0);


    int res = greedyVarOrder(tables.data(), tables.size(), 3,
                             clampRanks.data(), clampRanks.size(),
                             HEURISTIC_MIN_FILL, 1, &varOrder, &varOrderLen, &errMsg);
    if (res){
        cleanup(tables);
        printf("%s\n", errMsg);
        free_errorMessage(errMsg);
        return 1;
    }
    printf("order = ");
    for (int i = 0; i < varOrderLen; i++)
        printf("%d ", varOrder[i]);
    printf("\n");
    free_varOrder(varOrder);
    cleanup(tables);
    return 0;
}
