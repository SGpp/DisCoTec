//#ifdef USE_VTK
#ifndef VTK_HPP_
#define VTK_HPP_

#include <vtk/vtkSmartPointer.h>
#include <vtk/vtkImageData.h>
#include <vtk/vtkPointData.h>
#include <vtk/vtkDoubleArray.h>
#include <vtk/vtkNew.h>
#include <vtk/vtkXMLPImageDataWriter.h>
#include <vtk/vtkXMLImageDataWriter.h>
#include <string>
#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"

namespace combigrid {

  class DFGPlotFileWriter {
    private:
      static constexpr char outputDir[] = "./vtk/";

      static void writePlotMasterFile(size_t numberOfPieces,
                                      const std::string& filenamePrefix);

      static std::string createFilenamePrefix(DimType dim,
                                              const LevelVector& level,
                                              IndexType speciesId);

    public:
      /** Writes the given distributed full grid into a parallel vtk file */
      template <typename FG_ELEMENT>
      inline static void writePlotFile(
                                     const DistributedFullGrid<FG_ELEMENT>& dfg,
                                     IndexType speciesId) {
        assert(false && "Only implemented for double so far");
      }
  };

}
#endif /* VTK_HPP_ */
//#endif /* USE_VTK */
