
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"

namespace combigrid{
    std::string toString(combigrid::LevelVector const& l){
        std::stringstream ss;
        for(size_t i = 0; i < l.size(); ++i)
        {
            if(i != 0)
                ss << ",";
            ss << l[i];
        }
        return ss.str();
    }
}

