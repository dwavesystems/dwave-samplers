#ifndef INCLUDED_ORANG_BASE_H
#define INCLUDED_ORANG_BASE_H

#include <cstddef>
#include <set>
#include <vector>

#include <boost/cstdint.hpp>

namespace orang {

typedef boost::uint32_t Var;
typedef std::vector<Var> VarVector;
typedef std::set<Var> VarSet;

typedef boost::uint16_t DomIndex;
typedef std::vector<DomIndex> DomIndexVector;

typedef std::vector<std::size_t> SizeVector;

} // namespace orang

#endif
