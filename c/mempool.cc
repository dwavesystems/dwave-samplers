#include <map>
#include <sstream>
#include <cstdlib>

#include "errors.h"
#include "mempool.h"

void* MemPool::malloc(size_t s)
{
    void * v = ::malloc(s);
    if (!v){
        std::ostringstream st;
        st << "Unable to allocate memory: " << s << " bytes";
        throw Errors(st.str());
    }
    allocated_mem[v] = s;
    return v;
}

void* MemPool::calloc(size_t num, size_t s)
{
    void * v = ::calloc(num, s);
    if (!v){
        std::ostringstream st;
        st << "Unable to allocate memory: " << s * num << " bytes";
        throw Errors(st.str());
    }
    allocated_mem[v] = s;
    return v;
}

void* MemPool::realloc(void* ptr, size_t s)
{
    void * v = ::realloc(ptr, s);
    if (!v && s > 0){
        std::ostringstream st;
        st << "Unable to allocate memory: " << s << " bytes";
        throw Errors(st.str());
    }
    allocated_mem.erase(ptr);
    allocated_mem[v] = s;
    return v;
}

void MemPool::free(void* addr)
{
    allocated_mem.erase(addr);
    ::free(addr);
}

void MemPool::release()
{
    allocated_mem.clear();
}

MemPool::~MemPool()
{
    for (std::map<void *, int>::iterator i = allocated_mem.begin(); i != allocated_mem.end(); i++){
        if (i->second)
            ::free(i->first);
    }
}
