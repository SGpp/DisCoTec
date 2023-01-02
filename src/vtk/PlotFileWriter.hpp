#ifdef USE_VTK
#ifndef VTK_HPP_
#define VTK_HPP_

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkNew.h>
#include <vtkXMLPImageDataWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <string>
#include "fullgrid/DistributedFullGrid.hpp"

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
#endif /* USE_VTK */
