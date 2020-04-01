#include "sgpp/distributedcombigrid/vtk/PlotFileWriter.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/FTUtils.hpp"
#ifdef USE_VTK

namespace combigrid {


constexpr char DFGPlotFileWriter::outputDir[];

void DFGPlotFileWriter::writePlotMasterFile(size_t numberOfPieces,
                                            const std::string& filenamePrefix) {
  std::string filename = filenamePrefix + ".pvti";
  vtkSmartPointer<vtkXMLPImageDataWriter> writer =
    vtkSmartPointer<vtkXMLPImageDataWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetNumberOfPieces(static_cast<int>(numberOfPieces));
  writer->Update();
}

template <>
void DFGPlotFileWriter::writePlotFile<real>(
                                          const DistributedFullGrid<real>& dfg,
                                          IndexType speciesId) {
  DimType dim = dfg.getDimension();
  assert(dim < 4 && "Vtk does not support domains with more than 3 dimensions");
  // create filename
  const LevelVector& level = dfg.getLevels();
  std::string filenamePrefix {createFilenamePrefix(dim, level, speciesId)};
  RankType myRank = theMPISystem()->getLocalRank();
  std::string filename = filenamePrefix + "_" + std::to_string(myRank) + ".vti";

  // set metadata
  IndexVector dimensions = dfg.getLocalSizes();
  std::vector<double> spacing(dimensions.begin(), dimensions.end());
  std::for_each(spacing.begin(), spacing.end(),
      [](double& x){x = 1./(x-1.);});
  vtkSmartPointer<vtkImageData> imageData =
    vtkSmartPointer<vtkImageData>::New();
  imageData->AllocateScalars(VTK_DOUBLE, 1);
  std::vector<real> origin(dim);
  dfg.getCoordsLocal(0, origin);

  // reverse due to wrong order in dfg
  std::reverse(dimensions.begin(), dimensions.end());
  std::reverse(origin.begin(), origin.end());
  std::reverse(spacing.begin(), spacing.end());

  switch(dim) {
    case 1:
      imageData->SetDimensions(static_cast<int>(dimensions[0]), 1, 1);
      imageData->SetSpacing(spacing[0], 0., 0.);
      imageData->SetOrigin(origin[0], 0., 0.);
      break;
    case 2:
      imageData->SetDimensions(static_cast<int>(dimensions[0]),
                               static_cast<int>(dimensions[1]), 1);
      imageData->SetSpacing(spacing[0], spacing[1], 0.);
      imageData->SetOrigin(origin[0], origin[1], 0.);
      break;
    case 3:
      imageData->SetDimensions(static_cast<int>(dimensions[0]),
                               static_cast<int>(dimensions[1]),
                               static_cast<int>(dimensions[2]));
      imageData->SetSpacing(spacing[0], spacing[1], spacing[2]);
      imageData->SetOrigin(origin[0], origin[1], origin[2]);
      break;
    default:
      assert(false);
      break;
  }

  // loop over rows i -> dim2
  vtkNew<vtkDoubleArray> scalars;
  IndexVector localOffsets = dfg.getLocalOffsets();
  IndexVector nrLocalPoints = dfg.getLocalSizes();
  for (IndexType i = 0; i < nrLocalPoints[1]; ++i) {
    IndexType offset = localOffsets[1] * i;
    for (IndexType j = 0; j < nrLocalPoints[0]; ++j) {
      scalars->InsertNextTuple1(dfg.getData()[offset + j]);
    }
  }
  imageData->GetPointData()->SetScalars(scalars);
  /*
  // also loop over dummy dimensions
  for (int z = 0; z < dimensions[2]; z++) {
    for (int y = 0; y < dimensions[1]; y++) {
      for (int x = 0; x < dimensions[0]; x++) {
        IndexVector locAxisId {x, y, z};
        locAxisId.resize(dim); // remove dummy dims
        IndexType idx = dfg.getLocalLinearIndex(locAxisId);
        real value = dfg.getData()[idx];
        double* pixel = static_cast<double*>(
            imageData->GetScalarPointer(x,y,z));
        pixel[0] = value;
      }
    }
  }
  */

  // write
  vtkSmartPointer<vtkXMLImageDataWriter> writer =
    vtkSmartPointer<vtkXMLImageDataWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData(imageData);
  writer->Write();

  // write master file
  MASTER_EXCLUSIVE_SECTION {
    //size_t numProcs = theMPISystem()->getNumProcs();
    //writePlotMasterFile(numProcs, filenamePrefix);
    //
    size_t numberOfPieces = theMPISystem()->getNumProcs();
    std::string filename = filenamePrefix + ".pvti";
    vtkSmartPointer<vtkXMLPImageDataWriter> writer =
      vtkSmartPointer<vtkXMLPImageDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetNumberOfPieces(static_cast<int>(numberOfPieces));
    writer->SetInputData(imageData);
    writer->Update();
  }
}

std::string DFGPlotFileWriter::createFilenamePrefix(DimType dim,
                                                    const LevelVector& level,
                                                    IndexType speciesId) {
  std::string filenamePrefix = std::string(outputDir) + "dfg_l";
  for (DimType d = 0; d < dim; ++d)
    filenamePrefix += std::to_string(level[d]) + "_";
  filenamePrefix += "s" + std::to_string(speciesId);
  return filenamePrefix;
}

}
#endif /* USE_VTK */
