#pragma once
#include <cstddef>
#include <cstdint>
#include <exception>
#include <stdexcept>
#include <vector>

namespace initialization {
#ifdef MP_PROC_NUM
constexpr int kDefaultInitializationThreads = MP_PROC_NUM;
#else
constexpr int kDefaultInitializationThreads = 1;
#endif

// Workers write distinct preallocated outputs. Exceptions are rethrown in
// input order on the caller; no exception may escape an OpenMP worker.
template<class Function>
void ForEachInitializationTask(size_t count,int workers,Function work) {
    if(workers<1) throw std::invalid_argument("initialization worker count must be positive");
    if(workers==1) {
        for(size_t i=0;i<count;++i) work(i);
        return;
    }
    std::vector<std::exception_ptr> errors(count);
    #pragma omp parallel for num_threads(workers) schedule(dynamic,1)
    for(int64_t i=0;i<static_cast<int64_t>(count);++i) {
        try {work(static_cast<size_t>(i));}
        catch(...) {errors[i]=std::current_exception();}
    }
    for(const auto& error:errors) if(error) std::rethrow_exception(error);
}
}
