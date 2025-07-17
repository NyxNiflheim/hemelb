// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

// LocalPropertyHdf5Output.h
#include "extraction/LocalPropertyOutput.h"
#include "hdf5.h"
#include <vector>
#pragma once
#ifdef USE_HDF5
namespace hemelb
{
  namespace extraction
  {
    class LocalPropertyHdf5Output final : public LocalPropertyOutput
    {
    public:
      LocalPropertyHdf5Output(IterableDataSource& dataSource,
                              const PropertyOutputFile& outputSpec,
                              const net::IOCommunicator& ioComms);
      ~LocalPropertyHdf5Output();
      
      void Write(unsigned long timestepNumber, unsigned long totalSteps) override;

    private:
      void WriteXDMFFile(unsigned long timestep);
      void CreateCompoundType();
      
      // HDF5 MPI communicator
      MPI_Comm mpi_comm;
      
      // HDF5 handles
      hid_t file_id = -1;
      hid_t dataset_id = -1;
      hid_t filespace_id = -1;
      hid_t compound_type_id = -1;
      size_t compound_type_size = 0;       
      std::vector<size_t> field_offsets; 
    };
  }
}
#endif // USE_HDF5
