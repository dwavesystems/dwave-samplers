
#ifndef __ORANG_C_HELPER_H__
#define __ORANG_C_HELPER_H__

#include <vector>
#include <orang/orang.h>
#include "interface.h"

using orang::Table;
using orang::Var;
using orang::VarVector;
using orang::DomIndexVector;

using orang::DummyOperations;
using orang::Task;
typedef Task<DummyOperations> task_type;
typedef DummyOperations::CtorArgs ctor_args;

template<typename Tinp, typename Tout>
Tout doNothing(const Tinp & a) {return a;}

template<typename Y>
std::vector<typename Table<Y>::smartptr> createTables(TableEntry * tables, int tableNum){
    std::vector<typename Table<Y>::smartptr> tb;
    tb.reserve(tableNum);
    for (int i = 0; i < tableNum; i++){
        int varNum = tables[i].num_vars;
        VarVector vars;
        DomIndexVector domSizes;
        vars.reserve(varNum);
        domSizes.reserve(varNum);
        for (int j = 0; j < varNum; j++){
            vars.push_back(tables[i].vars[j]);
            domSizes.push_back(tables[i].domSizes[j]);
        }
        typename Table<Y>::smartptr tableptr( new Table<Y>(vars, domSizes) );
        tb.push_back(tableptr);
    }
    return tb;
}

#endif // __ORANG_C_HELPER_H__
