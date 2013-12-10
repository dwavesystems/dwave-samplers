#include <cstdlib>
#include <cerrno>
#include <sys/time.h>
#include <sys/resource.h>
#include <limits>

#include <iostream>

#include "uai.h"

using std::numeric_limits;
using std::getenv;
using std::strtod;

namespace {
const double gb = static_cast<double>(1ul << 30);
}

void limitMemory() {
  char* memStr = getenv("UAI_MEMORY");
  if (memStr) {
    std::cerr << "UAI_MEMORY = '" << memStr << "'\n";
    char* endPtr;
    double memLimitD = strtod(memStr, &endPtr) * gb;
    if (*endPtr == '\0' && memLimitD > 0.0 && memLimitD <= numeric_limits<rlim_t>::max()) {
      rlim_t memLimit = static_cast<rlim_t>(memLimitD);
      rlimit rl;
      getrlimit(RLIMIT_AS, &rl);
      std::cerr << "Previous address space limits: cur=" << rl.rlim_cur << " max=" << rl.rlim_max << "\n";
      if (rl.rlim_max == RLIM_INFINITY || rl.rlim_max >= memLimit) {
        std::cerr << "Setting address space limit to " << memLimit << "\n";
        rl.rlim_cur = memLimit;
        setrlimit(RLIMIT_AS, &rl);
      }
    }
  }
}
