#pragma once

#include <vector>

#include "fullgrid/DistributedFullGrid.hpp"
#include "hierarchization/DistributedHierarchization.hpp"
#include "utils/Stats.hpp"
#include "task/Task.hpp"

namespace combigrid {

template <typename FG_ELEMENT>
static void hierarchizeOverlapping(std::vector<std::unique_ptr<Task>>& tasks,
                                      const std::vector<bool>& hierarchizationDims,
                                      const std::vector<BasisFunctionBasis*>& hierarchicalBases,
                                      const LevelVector& lmin) {
  if (tasks.empty()) return;
  const auto& firstDfg = tasks[0]->getDistributedFullGrid();
  assert(firstDfg.getDimension() > 0);
  assert(firstDfg.getDimension() == hierarchizationDims.size());
  assert(lmin.size() == firstDfg.getDimension());

  std::vector<RemoteDataCollector<FG_ELEMENT>> remoteDataPerTask(tasks.size());

  // get cartesian comm dimensions from mpi cartesian utils
  const auto& cartDims = firstDfg.getCartesianUtils().getCartesianDimensions();
  auto maxCartesianDim = *std::max_element(cartDims.cbegin(), cartDims.cend());

  std::vector<std::pair<std::vector<MPI_Request>, std::vector<MPI_Request>>> requests(
      tasks.size(), std::make_pair(std::vector<MPI_Request>(),
                                   std::vector<MPI_Request>()));

  std::set<DimType> setDimensions;
  for (DimType dim = 0; dim < firstDfg.getDimension(); ++dim) {
    if (hierarchizationDims[dim]) {
      setDimensions.insert(dim);
    }
  }
  // const iterator for setDimensions
  auto setDimIt = setDimensions.cbegin();
  // TODO when doing this all at once, we are duplicating the full grid data of this rank --
  // experimental!!!
  for (size_t taskIndex = 0; taskIndex < tasks.size(); ++taskIndex) {
    auto dim = *(setDimIt);
    if (dynamic_cast<HierarchicalHatBasisFunction*>(hierarchicalBases[dim]) != nullptr ||
        dynamic_cast<HierarchicalHatPeriodicBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
      exchangeData1d(tasks[taskIndex]->getDistributedFullGrid(), dim, remoteDataPerTask[taskIndex], requests[taskIndex].first,
                     requests[taskIndex].second, lmin[dim]);
    } else {
      exchangeAllData1d(tasks[taskIndex]->getDistributedFullGrid(), dim, remoteDataPerTask[taskIndex], requests[taskIndex].first,
                        requests[taskIndex].second);
    }
  }

  // iterate dimensions first, then grids
  for (; setDimIt != setDimensions.cend(); ++setDimIt) {
    auto dim = *(setDimIt);
    for (size_t taskIndex = 0; taskIndex < tasks.size(); ++taskIndex) {
      // wait for requests, then hierarchize
      Stats::startEvent("h wait");
      MPI_Waitall(static_cast<int>(requests[taskIndex].first.size()),
                  requests[taskIndex].first.data(), MPI_STATUSES_IGNORE);
      MPI_Waitall(static_cast<int>(requests[taskIndex].second.size()),
                  requests[taskIndex].second.data(), MPI_STATUSES_IGNORE);
      Stats::stopEvent("h wait");
      DistributedHierarchization::hierarchizeLocalData<FG_ELEMENT>(
          tasks[taskIndex]->getDistributedFullGrid(), remoteDataPerTask[taskIndex],
          hierarchicalBases[dim], dim, lmin[dim]);
      remoteDataPerTask[taskIndex].clear();

      // after hierarchization, start communication in next dimension
      if (setDimIt != std::prev(setDimensions.cend())) {
      Stats::startEvent("h transfer");
        auto nextDim = *(std::next(setDimIt));
        if (dynamic_cast<HierarchicalHatBasisFunction*>(hierarchicalBases[nextDim]) != nullptr ||
            dynamic_cast<HierarchicalHatPeriodicBasisFunction*>(hierarchicalBases[nextDim]) !=
                nullptr) {
          exchangeData1d(tasks[taskIndex]->getDistributedFullGrid(), nextDim, remoteDataPerTask[taskIndex], requests[taskIndex].first,
                         requests[taskIndex].second, lmin[nextDim]);
        } else {
          exchangeAllData1d(tasks[taskIndex]->getDistributedFullGrid(), nextDim, remoteDataPerTask[taskIndex], requests[taskIndex].first,
                            requests[taskIndex].second);
        }
      }
      Stats::stopEvent("h transfer");
    }
  }
}


template <typename FG_ELEMENT>
static void dehierarchizeOverlapping(std::vector<std::unique_ptr<Task>>& tasks,
                                      const std::vector<bool>& hierarchizationDims,
                                      const std::vector<BasisFunctionBasis*>& hierarchicalBases,
                                      const LevelVector& lmin) {
  if (tasks.empty()) return;
  const auto& firstDfg = tasks[0]->getDistributedFullGrid();
  assert(firstDfg.getDimension() > 0);
  assert(firstDfg.getDimension() == hierarchizationDims.size());
  assert(lmin.size() == firstDfg.getDimension());

  std::vector<RemoteDataCollector<FG_ELEMENT>> remoteDataPerTask(tasks.size());

  // get cartesian comm dimensions from mpi cartesian utils
  const auto& cartDims = firstDfg.getCartesianUtils().getCartesianDimensions();
  auto maxCartesianDim = *std::max_element(cartDims.cbegin(), cartDims.cend());

  std::vector<std::pair<std::vector<MPI_Request>, std::vector<MPI_Request>>> requests(
      tasks.size(), std::make_pair(std::vector<MPI_Request>(),
                                   std::vector<MPI_Request>()));

  std::set<DimType> setDimensions;
  for (DimType dim = 0; dim < firstDfg.getDimension(); ++dim) {
    if (hierarchizationDims[dim]) {
      setDimensions.insert(dim);
    }
  }
  // const iterator for setDimensions
  auto setDimIt = setDimensions.cbegin();
  // TODO when doing this all at once, we are duplicating the full grid data of this rank --
  // experimental!!!
  for (size_t taskIndex = 0; taskIndex < tasks.size(); ++taskIndex) {
    auto dim = *(setDimIt);
    if (dynamic_cast<HierarchicalHatBasisFunction*>(hierarchicalBases[dim]) != nullptr ||
        dynamic_cast<HierarchicalHatPeriodicBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
      exchangeData1dDehierarchization(tasks[taskIndex]->getDistributedFullGrid(), dim, remoteDataPerTask[taskIndex], requests[taskIndex].first,
                     requests[taskIndex].second, lmin[dim]);
    } else {
      exchangeAllData1d(tasks[taskIndex]->getDistributedFullGrid(), dim, remoteDataPerTask[taskIndex], requests[taskIndex].first,
                        requests[taskIndex].second);
    }
  }

  // iterate dimensions first, then grids
  for (; setDimIt != setDimensions.cend(); ++setDimIt) {
    auto dim = *(setDimIt);
    for (size_t taskIndex = 0; taskIndex < tasks.size(); ++taskIndex) {
      // wait for requests, then hierarchize
      Stats::startEvent("d wait");
      MPI_Waitall(static_cast<int>(requests[taskIndex].first.size()),
                  requests[taskIndex].first.data(), MPI_STATUSES_IGNORE);
      MPI_Waitall(static_cast<int>(requests[taskIndex].second.size()),
                  requests[taskIndex].second.data(), MPI_STATUSES_IGNORE);
      Stats::stopEvent("d wait");
      DistributedHierarchization::dehierarchizeLocalData<FG_ELEMENT>(
          tasks[taskIndex]->getDistributedFullGrid(), remoteDataPerTask[taskIndex],
          hierarchicalBases[dim], dim, lmin[dim]);
      remoteDataPerTask[taskIndex].clear();

      // after hierarchization, start communication in next dimension
      if (setDimIt != std::prev(setDimensions.cend())) {
      Stats::startEvent("d transfer");
        auto nextDim = *(std::next(setDimIt));
        if (dynamic_cast<HierarchicalHatBasisFunction*>(hierarchicalBases[nextDim]) != nullptr ||
            dynamic_cast<HierarchicalHatPeriodicBasisFunction*>(hierarchicalBases[nextDim]) !=
                nullptr) {
          exchangeData1dDehierarchization(tasks[taskIndex]->getDistributedFullGrid(), nextDim, remoteDataPerTask[taskIndex], requests[taskIndex].first,
                         requests[taskIndex].second, lmin[nextDim]);
        } else {
          exchangeAllData1d(tasks[taskIndex]->getDistributedFullGrid(), nextDim, remoteDataPerTask[taskIndex], requests[taskIndex].first,
                            requests[taskIndex].second);
        }
      }
      Stats::stopEvent("d transfer");
    }
  }
}

}  // namespace combigrid
