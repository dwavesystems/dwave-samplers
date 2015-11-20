
#include <vector>
#include <boost/random.hpp>
#include <boost/iterator/indirect_iterator.hpp>

#include "orang/orang.h"
#include "orang/varorder.h"
#include "interface.h"
#include "helper.h"
#include "errors.h"


typedef boost::uniform_01<boost::mt19937> Rng;
using orang::Exception;
using orang::Table;
using orang::VarVector;
using orang::Task;
using orang::DummyOperations;
using orang::greedyvarorder::Heuristics;
using orang::greedyVarOrder;


using std::vector;

typedef Task<DummyOperations> task_type;
typedef DummyOperations::CtorArgs ctor_args;

#define MAX_ERROR_SIZE 100

extern "C"{


int greedyVarOrder(TableEntry * tables, int tables_len, int maxComplexity,


                   int * clampRanks, int clampRanks_len,


                   int heuristic, float selectionScale,


                   int ** variableOrder, int * variableOrder_len, char ** errorMessage){
    try{
        if (errorMessage == NULL)
            return 1;
        if (variableOrder == NULL)
            throw Errors("variableOrder is NULL");
        if (variableOrder_len == NULL)
            throw Errors("variableOrder_len is NULL");
        static Rng rng(boost::mt19937(static_cast<unsigned int>(std::time(0))));
        vector<Table<char>::smartptr> tb = createTables<char>(tables, tables_len);
        vector<int> cr;
        cr.reserve(clampRanks_len);
        for (int i = 0; i < clampRanks_len; i++)
            cr.push_back(clampRanks[i]);
        task_type task(make_indirect_iterator(tb.begin()), make_indirect_iterator(tb.end()),
            ctor_args(), static_cast<Var>(cr.size()));
        if (clampRanks_len == 0)
            cr.assign(task.numVars(), 0);
        if (cr.size() != task.numVars())
            throw Errors("'clampRanks' parameter must be empty or have size no less than the largest variable index");
        Heuristics h;
        switch (heuristic){
            case HEURISTIC_MIN_DEG:
                h = orang::greedyvarorder::MIN_DEGREE;
                break;
            case HEURISTIC_W_MIN_DEG:
                h = orang::greedyvarorder::WEIGHTED_MIN_DEGREE;
                break;
            case HEURISTIC_MIN_FILL:
                h = orang::greedyvarorder::MIN_FILL;
                break;
            case HEURISTIC_W_MIN_FILL:
                h = orang::greedyvarorder::WEIGHTED_MIN_FILL;
                break;
            default:
                throw Errors("Invalid heuristic");
        };
        VarVector varOrder = greedyVarOrder(task, maxComplexity, cr, h, rng, selectionScale);
        *variableOrder_len = varOrder.size();
        *variableOrder = new int[varOrder.size()];
        for (int i = 0; i < varOrder.size(); i++)
            (*variableOrder)[i] = varOrder[i];
        return 0;
    }
    catch(Errors & e){
        *errorMessage = new char[MAX_ERROR_SIZE];
        sprintf(*errorMessage, "%s", e.what());
        return 1;
    }
    catch (std::bad_alloc &){
        *errorMessage = new char[MAX_ERROR_SIZE];
        sprintf(*errorMessage, "Out of memory");
        return 1;
    }
    catch (Exception & e){
        *errorMessage = new char[MAX_ERROR_SIZE];
        sprintf(*errorMessage, "%s", e.what().c_str());
        return 1;
    }
}

void free_varOrder(int * varOrder){
    if (varOrder)
        delete [] varOrder;
}

void free_errorMessage(char * errMsg){
    if (errMsg)
        delete [] errMsg;
}


} // extern "C"