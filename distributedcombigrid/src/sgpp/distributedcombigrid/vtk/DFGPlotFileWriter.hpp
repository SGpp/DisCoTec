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
#include <boost/property_tree/xml_parser.hpp>
#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"

namespace combigrid {

  /** Helper class to generate vtk ImageData (.vti/.pvti) plot file(s)
   *  of a given distributed full grid.
   */
  template <typename FG_ELEMENT>
  class DFGPlotFileWriter {
    public:
      DFGPlotFileWriter(DistributedFullGrid<FG_ELEMENT>& dfg,
                        IndexType speciesId,
                        std::string outputDirectory = "./vtk/");

      /** Writes the distributed full grid into vtk file(s) */
      void writePlotFile();

    private:
      /** Output directory of the vtk files*/
      std::string outputDirectory_;

      /** Prefix for all vtk files without the leading directory*/
      std::string filenamePrefixWithoutDirectory_;

      /** The distributed full grid which will be written out */
      DistributedFullGrid<FG_ELEMENT>& dfg_;

      /** Identifies the species of the dfg */
      IndexType speciesId_;

      /** Vtk style min and max indices per dimension of the whole domain */
      std::array<int, 6> globalExtent_; // = {min1, max1, ... }

      /** Vtk style min and max indices per dimension of the local domain part*/
      std::array<int, 6> localExtent_; // = {min1, max1, ... }

      /** Mesh width of the domain */
      std::array<double, 3> spacing_;

      /** Coordinate of the domain's origin */
      std::array<double, 3> origin_; // = {0., 0., 0.}

      /** Sets global/whole extent of the domain */
      void setGlobalExtent();

      /** Sets local extent of domain */
      void setLocalExtent();

      /** Sets mesh width of domain */
      void setSpacing();

      /** Sets filename prefix without leading directory */
      void setFilenamePrefixWithoutDirectory();

      /** Creates piece filename without leading directory for given piece id */
      std::string createPieceFilenameWithoutDirectory(size_t pieceId);

      /** Write the .vti piece file
       *
       * @return vtk image data containing meta information and point data
       */
      vtkSmartPointer<vtkImageData> writePieceFile(
          const std::string& filename);

      /** Write the .pvti master file
       *
       * @param pieceData Used to extract meta information
       */
      void writeMasterFile(vtkSmartPointer<vtkImageData> pieceData);

  };


  template <typename FG_ELEMENT>
  DFGPlotFileWriter<FG_ELEMENT>::DFGPlotFileWriter(
      DistributedFullGrid<FG_ELEMENT>& dfg,
      IndexType speciesId,
      std::string outputDirectory) : dfg_ {dfg}, speciesId_ {speciesId},
                                     outputDirectory_ {outputDirectory},
                                     origin_ {0., 0., 0.}
  {
    setGlobalExtent();
    setLocalExtent();
    setSpacing();
    setFilenamePrefixWithoutDirectory();
  }


  template <typename FG_ELEMENT>
  inline void DFGPlotFileWriter<FG_ELEMENT>::writePlotFile() {
    DimType dim = dfg_.getDimension();
    assert(dim < 4 && "Vtk does not support image domains with more than 3 "
                      "dimensions");

    // create piece file
    RankType myRank = theMPISystem()->getLocalRank();
    std::string myPieceFilename = outputDirectory_ +
                                  createPieceFilenameWithoutDirectory(
                                      static_cast<size_t>(myRank));
    auto myPieceImageData = writePieceFile(myPieceFilename);

    // create master file
    // TODO this requires exchange of boundary points
    /*
    size_t groupSize = theMPISystem()->getNumProcs();

    if (groupSize > 1) {
      writeMasterFile(myPieceImageData);
    }
    */
  }

  template <typename FG_ELEMENT>
  inline vtkSmartPointer<vtkImageData>
    DFGPlotFileWriter<FG_ELEMENT>::writePieceFile(const std::string& filename) {
    assert(false && "Only implemented for double so far");
  }


  /* TODO remove copy of grid values */
  template <>
  inline vtkSmartPointer<vtkImageData> DFGPlotFileWriter<real>::writePieceFile(
                                                  const std::string& filename) {
    // set data (correct order of dimensions in dfg)
    vtkNew<vtkDoubleArray> pointValues;
    pointValues->SetName("GridValues");
    IndexVector localOffsets = dfg_.getLocalOffsets();
    IndexVector nrLocalPoints = dfg_.getLocalSizes();
    DimType dim = dfg_.getDimension();
    switch (dim) {
      case 1:
        for (IndexType i = 0; i < nrLocalPoints[0]; ++i) {
          pointValues->InsertNextTuple1(dfg_.getData()[i]);
        }
        break;
      case 2:
        for (IndexType i = 0; i < nrLocalPoints[1]; ++i) {
          IndexType offset = localOffsets[1] * i;
          for (IndexType j = 0; j < nrLocalPoints[0]; ++j) {
            pointValues->InsertNextTuple1(dfg_.getData()[offset + j]);
          }
        }
        break;
      case 3:
        for (IndexType k = 0; k < nrLocalPoints[2]; ++k) {
          IndexType offsetK = localOffsets[2] * k;
          for (IndexType i = 0; i < nrLocalPoints[1]; ++i) {
            IndexType offset = offsetK + localOffsets[1] * i;
            for (IndexType j = 0; j < nrLocalPoints[0]; ++j) {
              pointValues->InsertNextTuple1(dfg_.getData()[offset + j]);
            }
          }
        }
        break;
      default:
        assert(false && "wrong number of dimensions");
    }

    // set meta information of piece file
    auto imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetExtent(localExtent_.data());
    imageData->SetOrigin(origin_.data());
    imageData->SetSpacing(spacing_.data());
    imageData->GetPointData()->AddArray(pointValues);

    // write piece file
    auto writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(imageData);
    writer->Write();

    return imageData;
  }

  /**
   * Possible manual creation of the master file without using the poorly
   * documented features for creating parallel files directly with vtk.
   *
   * So far the master file will open with an error in paraview because of
   * missing ghost points in the piece files.
   * -> this requires the exchange of boundary points, such that the extents
   *  overlap once: e.g. in 1D with 2 procs and m points the extents/indices
   *  must be [0, n] in piece 1 and [n, m-1] in piece 2.
   */
  template <typename FG_ELEMENT>
  void DFGPlotFileWriter<FG_ELEMENT>::writeMasterFile(
                                      vtkSmartPointer<vtkImageData> pieceData) {
    size_t groupSize = theMPISystem()->getNumProcs();
    MASTER_EXCLUSIVE_SECTION {
      // let vtk generate the master file with the master piece as source
      // this way we can be sure that the meta information matches the pieces.
      auto pwriter = vtkSmartPointer<vtkXMLPImageDataWriter>::New();
      pwriter->SetNumberOfPieces(static_cast<int>(groupSize));
      pwriter->SetInputData(pieceData);
      pwriter->WriteToOutputStringOn();
      pwriter->SetFileName(filenamePrefixWithoutDirectory_.c_str());
      pwriter->Write();

      // store xml output in a property tree for further processing
      std::stringstream masterFileTemplate { pwriter->GetOutputString() };
      boost::property_tree::ptree ptree;
      boost::property_tree::xml_parser::read_xml(masterFileTemplate, ptree);

      // remove already written piece file (wrong file if masterRank != 0)
      ptree.erase("VTKFile.PImageData.Piece");

      // set correct whole extent
      auto& wholeExtendAttr = ptree.get_child("VTKFile.PImageData.<xmlattr>"
                                              ".WholeExtent");
      std::stringstream wholeExtendSS;
      wholeExtendSS << globalExtent_[0] << " " << globalExtent_[1] << " "
                    << globalExtent_[2] << " " << globalExtent_[3] << " "
                    << globalExtent_[4] << " " << globalExtent_[5];
      wholeExtendAttr.put_value(wholeExtendSS.str());

      // gather group extents on master
      std::vector<int> groupExtents(localExtent_.size() * groupSize);
      MPI_Gather(localExtent_.data(), localExtent_.size(), MPI_INT,
                 groupExtents.data(), localExtent_.size(), MPI_INT,
                 theMPISystem()->getMasterRank(),
                 theMPISystem()->getLocalComm());

      // set group pieces
      for (int i = 0; i < static_cast<int>(groupSize); ++i) {
        std::string pieceFilename {
                  createPieceFilenameWithoutDirectory(static_cast<size_t>(i)) };
        std::vector<int> pieceExtent { groupExtents.begin() + i*6,
                                       groupExtents.begin() + i*6+6 };
        std::stringstream pieceExtentSS;
        pieceExtentSS << pieceExtent[0] << " " << pieceExtent[1] << " "
                      << pieceExtent[2] << " " << pieceExtent[3] << " "
                      << pieceExtent[4] << " " << pieceExtent[5];
        boost::property_tree::ptree pieceTree;
        pieceTree.put("<xmlattr>.Extent", pieceExtentSS.str());
        pieceTree.put("<xmlattr>.Source", pieceFilename);
        auto& pImageData = ptree.get_child("VTKFile.PImageData");
        pImageData.add_child("Piece", pieceTree);
      }

      // write file
      std::string masterFilename = outputDirectory_
                                   + filenamePrefixWithoutDirectory_
                                   + ".pvti";
      boost::property_tree::write_xml(masterFilename, ptree,
           std::locale(),
           boost::property_tree::xml_writer_make_settings<std::string>(' ', 2));
    } else {
      // gather group extents on master
      MPI_Gather(localExtent_.data(), localExtent_.size(), MPI_INT,
                 nullptr, 0, MPI_INT,
                 theMPISystem()->getMasterRank(),
                 theMPISystem()->getLocalComm());
    }
  }

  template <typename FG_ELEMENT>
  void DFGPlotFileWriter<FG_ELEMENT>::setFilenamePrefixWithoutDirectory() {
    DimType dim = dfg_.getDimension();
    const LevelVector& level = dfg_.getLevels();

    std::stringstream filenamePrefix;
    filenamePrefix <<  "dfg_l";
    for (DimType d = 0; d < dim; ++d)
      filenamePrefix << level[d] << "_";
    filenamePrefix << "s" << speciesId_;

    filenamePrefixWithoutDirectory_ = filenamePrefix.str();
  }

  template <typename FG_ELEMENT>
  std::string
  DFGPlotFileWriter<FG_ELEMENT>::createPieceFilenameWithoutDirectory(
      size_t pieceId) {
    return filenamePrefixWithoutDirectory_ + "_" + std::to_string(pieceId)
           + ".vti";
  }

  template <typename FG_ELEMENT>
  void DFGPlotFileWriter<FG_ELEMENT>::setGlobalExtent() {
    std::vector<int> dimensionsGlob { dfg_.getGlobalSizes().begin(),
                                      dfg_.getGlobalSizes().end() };
    dimensionsGlob.resize(3, 1);
    globalExtent_ = { 0, dimensionsGlob[0]-1, 0, dimensionsGlob[1]-1, 0,
                      dimensionsGlob[2]-1 };
  }

  template <typename FG_ELEMENT>
  void DFGPlotFileWriter<FG_ELEMENT>::setLocalExtent() {
    DimType dim = dfg_.getDimension();
    IndexVector minExtend(dim);
    IndexVector maxExtend(dim);
    IndexType minIndex = dfg_.getGlobalLinearIndex(0);
    IndexType maxIndex = dfg_.getGlobalLinearIndex(dfg_.getNrLocalElements()-1);
    dfg_.getGlobalVectorIndex(minIndex, minExtend);
    dfg_.getGlobalVectorIndex(maxIndex, maxExtend);
    minExtend.resize(3, 0);
    maxExtend.resize(3, 0);
    localExtent_ = { static_cast<int>(minExtend[0]),
                     static_cast<int>(maxExtend[0]),
                     static_cast<int>(minExtend[1]),
                     static_cast<int>(maxExtend[1]),
                     static_cast<int>(minExtend[2]),
                     static_cast<int>(maxExtend[2]) };
  }

  template <typename FG_ELEMENT>
  void DFGPlotFileWriter<FG_ELEMENT>::setSpacing() {
    std::vector<double> spacing {dfg_.getGlobalSizes().begin(),
                                 dfg_.getGlobalSizes().end() };
    std::for_each(spacing.begin(), spacing.end(),
                  [](double& x){x = 1./(x-1.);});
    spacing.resize(3, 0.);
    spacing_ = {spacing[0], spacing[1], spacing[2]};
  }
}
#endif /* VTK_HPP_ */
#endif /* USE_VTK */
