#pragma once

#include <boost/serialization/export.hpp>

#include "discotec/fault_tolerance/FaultCriterion.hpp"
#include "discotec/fault_tolerance/StaticFaults.hpp"
#include "discotec/fault_tolerance/WeibullFaults.hpp"
// this header should be included once for every compilation unit; if there are
// "not registered" or "not exported"-type errors, maybe this header was called before all
// the relevant fixtures were included?
// for serializable classes outside the library libsgppdistributedcombigrid, each compilation unit
// needs to call BOOST_CLASS_EXPORT for those classes.

BOOST_CLASS_EXPORT(combigrid::FaultCriterion)
BOOST_CLASS_EXPORT(combigrid::StaticFaults)
BOOST_CLASS_EXPORT(combigrid::WeibullFaults)

BOOST_CLASS_EXPORT(combigrid::CombiParameters)
