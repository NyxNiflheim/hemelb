// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

// LocalPropertyHdf5Output.h
// In extraction/LocalPropertyHdf5Output.h
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
    class VtkPropertyOutput; // Forward declaration for friendship

    class LocalPropertyHdf5Output final : public LocalPropertyOutput
    {
      // Give VtkPropertyOutput friend access to reuse data gathering logic
      friend class VtkPropertyOutput;

    public:
      LocalPropertyHdf5Output(IterableDataSource& dataSource,
                              const PropertyOutputFile& outputSpec,
                              const net::IOCommunicator& ioComms);
      
      ~LocalPropertyHdf5Output() = default;
      
      void Write(unsigned long timestepNumber, unsigned long totalSteps) override;

    private:
      void WriteXDMFFile(const std::string& h5_filename, unsigned long timestep);

      MPI_Comm mpi_comm;
      std::string output_file_pattern;
    };
  }
}

#endif // USE_HDF5