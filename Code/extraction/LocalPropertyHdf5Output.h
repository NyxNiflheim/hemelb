// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

// LocalPropertyHdf5Output.h
// single h5 file + Group + compound + xmf
#pragma once

#ifdef USE_HDF5

#include "extraction/LocalPropertyOutput.h"
#include <hdf5.h>
#include <string>
#include <vector>

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
      void CreateCompoundType();
      void WriteXDMFFile();

      MPI_Comm mpi_comm;
      
      // file handler
      hid_t file_id = -1;
      
      // compound type handler
      hid_t compound_type_id = -1;
      size_t compound_type_size = 0;
      std::vector<size_t> field_offsets;
      
      // vecter for all field names
      std::vector<unsigned long> written_timesteps;
    };
  }
}

#endif // USE_HDF5