#include "sgpp/distributedcombigrid/combischeme/StaticCombiScheme.hpp"

namespace combigrid {

LevelVector StaticCombiScheme::makeClassicalMin(LevelVector lmin, LevelVector lmax) {
	assert(lmax >= lmin);
	const LevelVector diffs = lmax - lmin;
	LevelType minDiff = 0;

	//search for the smallest max - min difference that is not zero (which would
	//correspond to a dummy dimension)
	//The loop is constructed this way to avoid problems when there are only dimensions d
	//where lmin[d] = lmax[d] (where minDiff will be zero in the end)
	for (size_t i = 0; i < diffs.size(); ++i) {
		const auto currDiff = diffs.at(i);
		if (currDiff != 0) {
			if (minDiff == 0) {
				minDiff = currDiff;
			} else {
				minDiff = currDiff < minDiff ? currDiff : minDiff;
			}
		}
	}

	//set lmin accordingly
	for (size_t i = 0; i < lmin.size(); ++i) {
		//change the lmin only in non dummy dimensions
		if (diffs[i] != 0) {
			lmin[i] = lmax[i] - minDiff;
		}
	}

	return lmin;
}

/*
 * Generates the fault tolerant combination technique with extra
 * grids used in case of faults
 * */
void StaticCombiScheme::makeFaultTolerant() {
	const int extraDiags = 2;
	if (lmin_ == lmax_) return;
	// Add extra combiSpaces to ensure fault tolerance
	for (const LevelVector& level : levels_) {
		for (LevelType p = LevelType(dim()); p < LevelType(dim()) + extraDiags; ++p) {
			if (level >= lmin_ && sum(level) == getMaxLevelSum() - p) {
				combiSpaces_.push_back(level);
				coefficients_.push_back(0.0);
			}
		}
	}
}

LevelType StaticCombiScheme::getMaxLevelSum() {
	const LevelVector diff = lmax_ - lmin_;
	return *std::max_element(diff.begin(), diff.end());
}

void StaticCombiScheme::createLevelsRec(LevelType currMaxSum, size_t currentIndex, LevelVector l) {
		LevelType x = 1;
		if(isDummyDim(currentIndex)){
			x = lmin_.at(currentIndex);
		}
		for (; x <= std::min(lmin_.at(currentIndex) + currMaxSum, lmax_.at(currentIndex)); ++x) {
			l.at(currentIndex) = x;
			assert(l <= lmax_);
			if (currentIndex == 0) {
				levels_.push_back(l);
			} else {
				createLevelsRec(std::min(currMaxSum, currMaxSum - (x - lmin_.at(currentIndex))),	currentIndex - 1, l);
			}
		}
}

void StaticCombiScheme::computeCombiCoeffs() {
	combiSpaces_.clear();
	coefficients_.clear();
	for (const auto& combiSpace : levels_) {
		real coeff = 0;
		LevelVector maxNeighbour {combiSpace + LevelVector(dim(), 1)};

		for (const auto& possibleNeighbour : levels_) {
			if (possibleNeighbour >= combiSpace && possibleNeighbour <= maxNeighbour) {
				LevelVector diff = possibleNeighbour - combiSpace;
				coeff += std::pow(-1.0, std::accumulate(diff.begin(), diff.end(), 0));
			}
		}

		if (coeff != 0) {
			//Grids beyond the border should never have a non-zero index

			combiSpaces_.push_back(combiSpace);
			coefficients_.push_back(coeff);
		}
	}

}

}  // end namespace combigrid
