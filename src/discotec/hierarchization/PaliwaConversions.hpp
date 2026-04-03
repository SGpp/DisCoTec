#pragma once

#ifdef DISCOTEC_USE_PALIWA

#include <array>
#include <ddc/ddc.hpp>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <utility>

#include "discotec/utils/Types.hpp"
#include "paliwa/paliwa_dimensions.hpp"
#include "paliwa/paliwa_domains.hpp"
#include "paliwa/paliwa_transform.hpp"
#include "paliwa/paliwa_wavelets.hpp"

namespace combigrid {

// --------------------------------------------------------------------------
// Compile-time mapping from DisCoTec DIM to paliwa/DDC dimension types.
//
// DDC requires distinct named types for each dimension slot.  Paliwa
// pre-defines DDimA..DDimU, assigning each dimensionality a disjoint sets
//
// We encode this mapping once in PaliwaDDimList and derive everything else
// from it via index_sequence.
// --------------------------------------------------------------------------

// All paliwa DDim types in a single type list (one per letter, A-U = 21 types).
template <typename... Ts>
struct PaliwaDDimList {};

using AllPaliwaDDims =
    PaliwaDDimList<paliwa::DDimA, paliwa::DDimB, paliwa::DDimC, paliwa::DDimD, paliwa::DDimE,
                   paliwa::DDimF, paliwa::DDimG, paliwa::DDimH, paliwa::DDimI, paliwa::DDimJ,
                   paliwa::DDimK, paliwa::DDimL, paliwa::DDimM, paliwa::DDimN, paliwa::DDimO,
                   paliwa::DDimP, paliwa::DDimQ, paliwa::DDimR, paliwa::DDimS, paliwa::DDimT,
                   paliwa::DDimU>;

// Select the I-th type from a PaliwaDDimList.
template <std::size_t I, typename List>
struct PaliwaDDimAt;

template <std::size_t I, typename Head, typename... Tail>
struct PaliwaDDimAt<I, PaliwaDDimList<Head, Tail...>>
    : PaliwaDDimAt<I - 1, PaliwaDDimList<Tail...>> {};

template <typename Head, typename... Tail>
struct PaliwaDDimAt<0, PaliwaDDimList<Head, Tail...>> {
  using type = Head;
};

template <std::size_t I>
using PaliwaDDim = typename PaliwaDDimAt<I, AllPaliwaDDims>::type;

// Build PaliwaDimTraits<DIM> from an index_sequence selecting the right DDims.
namespace detail {

template <DimType DIM, typename IndexSeq>
struct PaliwaDimTraitsImpl;

template <DimType DIM, std::size_t... Is>
struct PaliwaDimTraitsImpl<DIM, std::index_sequence<Is...>> {
  static constexpr std::size_t offset = static_cast<std::size_t>(DIM) * (DIM - 1) / 2;

  using Domain = ddc::StridedDiscreteDomain<PaliwaDDim<offset + Is>...>;
  using Vector = ddc::DiscreteVector<PaliwaDDim<offset + Is>...>;

  template <typename T>
  static Domain domainFromLevel(const std::array<T, DIM>& level,
                                const std::array<T, DIM>& maxLevel) {
    return paliwa::strided_domain_from_level<PaliwaDDim<offset + Is>...>(level, maxLevel);
  }

  static Vector toVector(const std::array<long int, static_cast<std::size_t>(DIM)>& v) {
    return Vector{v[Is]...};
  }
};

}  // namespace detail

template <DimType DIM>
using PaliwaDimTraits =
    detail::PaliwaDimTraitsImpl<DIM, std::make_index_sequence<static_cast<std::size_t>(DIM)>>;

template <typename FG_ELEMENT, DimType DIM>
auto dfgToPaliwaDomain(const DistributedFullGrid<FG_ELEMENT, DIM>& dfg,
                       const std::array<long int, static_cast<std::size_t>(DIM)>& maxLevel) ->
    typename PaliwaDimTraits<DIM>::Domain {
  constexpr auto dim = static_cast<std::size_t>(DIM);
  const auto& levels = dfg.getLevels();
  std::array<long int, dim> ddcLevel{};
  for (DimType d = 0; d < DIM; ++d) {
    ddcLevel[d] = static_cast<long int>(levels[d]);
  }
  return PaliwaDimTraits<DIM>::domainFromLevel(ddcLevel, maxLevel);
}

template <DimType DIM>
typename PaliwaDimTraits<DIM>::Vector levelToPaliwaVector(const LevelVector& lv) {
  constexpr auto dim = static_cast<std::size_t>(DIM);
  std::array<long int, dim> arr{};
  for (DimType d = 0; d < DIM; ++d) {
    arr[d] = static_cast<long int>(lv[d]);
  }
  return PaliwaDimTraits<DIM>::toVector(arr);
}

template <DimType DIM>
typename PaliwaDimTraits<DIM>::Vector levelToPaliwaVector(const LevelArray<DIM>& la) {
  constexpr auto dim = static_cast<std::size_t>(DIM);
  std::array<long int, dim> arr{};
  for (DimType d = 0; d < DIM; ++d) {
    arr[d] = static_cast<long int>(la[d]);
  }
  return PaliwaDimTraits<DIM>::toVector(arr);
}

inline std::string paliwaWaveletName(BasisFunctionType basis) {
  if (basis == BasisFunctionType::HAT || basis == BasisFunctionType::HAT_PERIODIC) {
    return "hat";
  } else if (basis == BasisFunctionType::BIORTHOGONAL ||
             basis == BasisFunctionType::BIORTHOGONAL_PERIODIC) {
    return "biorthogonal";
  } else if (basis == BasisFunctionType::FULLWEIGHTING ||
             basis == BasisFunctionType::FULLWEIGHTING_PERIODIC) {
    return "fullweighting";
  } else {
    throw std::runtime_error("paliwaWaveletName: unsupported basis function type");
  }
}

}  // namespace combigrid

#endif  // DISCOTEC_USE_PALIWA
