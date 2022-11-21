#pragma once

#include <boost/serialization/export.hpp>

#include "sgpp/distributedcombigrid/fault_tolerance/FaultCriterion.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/StaticFaults.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/WeibullFaults.hpp"
#include "sgpp/distributedcombigrid/legacy/CombiLinearBasisFunction.hpp"

// this header should be included once for every compilation unit; if there are
// "not registered" or "not exported"-type errors, maybe this header was called before all
// the relevant fixtures were included?
// for serializable classes outside the library libsgppdistributedcombigrid, each compilation unit
// needs to call BOOST_CLASS_EXPORT for those classes.

BOOST_CLASS_EXPORT(combigrid::BasisFunctionBasis)
BOOST_CLASS_EXPORT(combigrid::LinearBasisFunction)
// BOOST_CLASS_EXPORT(combigrid::HierarchicalHatBasisFunction);
// cf. https://www.boost.org/doc/libs/1_71_0/libs/serialization/doc/serialization.html#export
// BOOST_CLASS_EXPORT_GUID(combigrid::LinearBasisFunction, "LinearBasisFunction")
// BOOST_CLASS_EXPORT(combigrid::HierarchicalHatBasisFunction);
BOOST_CLASS_EXPORT_GUID(combigrid::HierarchicalHatBasisFunction, "HierarchicalHatBasisFunction")
BOOST_CLASS_EXPORT_GUID(combigrid::FullWeightingBasisFunction, "FullWeightingBasisFunction")
BOOST_CLASS_EXPORT_GUID(combigrid::FullWeightingPeriodicBasisFunction,
                        "FullWeightingPeriodicBasisFunction")
BOOST_CLASS_EXPORT_GUID(combigrid::BiorthogonalBasisFunction, "BiorthogonalBasisFunction")
BOOST_CLASS_EXPORT_GUID(combigrid::BiorthogonalPeriodicBasisFunction,
                        "BiorthogonalPeriodicBasisFunction")

BOOST_CLASS_EXPORT(combigrid::FaultCriterion)
BOOST_CLASS_EXPORT(combigrid::StaticFaults)
BOOST_CLASS_EXPORT(combigrid::WeibullFaults)

BOOST_CLASS_EXPORT(combigrid::CombiParameters)
