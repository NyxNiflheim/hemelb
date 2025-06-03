// LocalPropertyHdf5Output.h
#pragma once
#ifdef USE_HDF5
#include "LocalPropertyOutputBase.h"
#include <hdf5.h>

namespace hemelb::io {

class LocalPropertyHdf5Output final : public LocalPropertyOutputBase {
public:
  LocalPropertyHdf5Output(IterableDataSource& src,
                          const PropertyOutputFile& spec,
                          const net::IOCommunicator& comms);

  void StartFile(const std::string& fn) override;
  void Write(unsigned long ts, unsigned long total) override;
  void Close() override;

private:
  MPI_Comm comm_;
  hid_t file_      = -1;
  hid_t dset_      = -1;
  hid_t memType_   = -1;
  hid_t fileSpace_ = -1;
};
} // namespace hemelb::io
#endif // USE_HDF5
