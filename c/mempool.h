#ifndef __MEMPOOL_H__
#define __MEMPOOL_H__

#include <map>

class MemPool{
    public:
        MemPool(){};
        ~MemPool();
        
        void * malloc(size_t s);
        void * calloc(size_t num, size_t s);
        void * realloc(void * ptr, size_t s);
        void free(void * addr);
        void release();
    private:
        std::map<void *, int> allocated_mem;
};

#endif // __MEMPOOL_H__
