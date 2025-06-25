// LocalPropertyHdf5Output.h
#include "extraction/LocalPropertyOutput.h"
#include "hdf5.h"
#pragma once
#ifdef USE_HDF5

namespace hemelb::extraction 
{
class LocalPropertyHdf5Output final : public LocalPropertyHdf5Output
{
  public:
    LocalPropertyHdf5Output(IterableDataSource& dataSource,
                            const PropertyOutputFile& outputSpec,
                            const net::IOCommunicator& ioComms);
    ~LocalPropertyHdf5Output();

    void Write(unsigned long timestepNumber, unsigned long totalStep) override;
    // Close() in ~LocalPropertyHdf5Output

  private:
    void StartFile(const std::string& fn) override;
    void CreateCompoundType();
    MPI_Comm mpi_comm;

    // HDF5 handles
    hid_t file_id      = -1;
    hid_t dataset_id   = -1;
    hid_t memSpace_id  = -1;
    hid_t fileSpace_id = -1;
    hid_t compound_type_id = -1;
};
}
#endif // USE_HDF5
