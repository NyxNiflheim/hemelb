// LocalPropertyHdf5Output.cc
// without hdf5writer
#include "extraction/LocalPropertyHdf5Output.h"
// #include "io/writers/Hdf5Writer.h"// for Hdf5Writer
// #ifdef USE_HDF5
#include <vector>
#include <cassert>
#include <numeric>
#include <stdexcept>

namespace hemelb::extraction 
{
  LocalPropertyHdf5Output::LocalPropertyHdf5Output(IterableDataSource& dataSource,
                                                   const PropertyOutputFile& outputSpec,
                                                   const net::IOCommunicator& comms)
    : LocalPropertyOutput(dataSource, outputSpec, comms), mpi_comm(comms.GetCommunicator())
  {
    // create compound type
    CreateCompoundType();

    // create file and dataset
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl, mpi_comm, MPI_INFO_NULL);
    file_id = H5Fcreate(outputSpec.filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    H5Pclose(fapl);

    // create dataset
    hsize_t global_dims[1] = { global_site_count };
    filespace_id = H5Screate_simple(1, global_dims, NULL);
    dataset_id = H5Dcreate2(file_id, "results", compound_type_id, filespace_id,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }

  LocalPropertyHdf5Output::~LocalPropertyHdf5Output()
  {
    if (dataset_id >= 0) H5Dclose(dataset_id);
    if (filespace_id >= 0) H5Sclose(filespace_id);
    if (compound_type_id >= 0) H5Tclose(compound_type_id);
    if (file_id >= 0) H5Fclose(file_id);
  }

  void LocalPropertyHdf5Output::CreateCompoundType() {
    // traverse the outputSpec to determine the size of the compound type
    compound_type_size = 0;
    compound_type_id = H5Tcreate(H5T_COMPOUND, 1); // init size =1

    // first 3 fields for coordinates (x, y, z)
    H5Tinsert(compound_type_id, "x", compound_type_size, H5T_NATIVE_UINT32);
    field_offsets.push_back(compound_type_size);
    compound_type_size += sizeof(uint32_t);
    H5Tinsert(compound_type_id, "y", compound_type_size, H5T_NATIVE_UINT32);
    field_offsets.push_back(compound_type_size);
    compound_type_size += sizeof(uint32_t);
    H5Tinsert(compound_type_id, "z", compound_type_size, H5T_NATIVE_UINT32);
    field_offsets.push_back(compound_type_size);
    compound_type_size += sizeof(uint32_t);
    
    for (const auto& field_spec : outputSpec.fields) {
      unsigned field_len = GetFieldLength(field_spec.src);
      hid_t field_type = H5T_NATIVE_DOUBLE;// use double for storing all fields

      if (field_len == 1) { // scalar, like pressure, shear stress
        H5Tinsert(compound_type_id, field_spec.name.c_str(), compound_type_size, field_type);
        field_offsets.push_back(compound_type_size);
        compound_type_size += sizeof(double);
      } else { // velocity
        hsize_t dims[1] = {field_len};
        hid_t array_tid = H5Tarray_create(field_type, 1, dims);
        H5Tinsert(compound_type_id, field_spec.name.c_str(), compound_type_size, array_tid);
        field_offsets.push_back(compound_type_size);
        H5Tclose(array_tid);
        compound_type_size += sizeof(double) * field_len;
      }
    }
    // compound type size is set to the total size of all fields
    H5Tset_size(compound_type_id, compound_type_size);
  }

  void LocalPropertyHdf5Output::Write(unsigned long timestepNumber, unsigned long totalSteps)
  {
    if (!ShouldWrite(timestepNumber)) {
        return;
    }

    // create a char buffer for storing the data
    std::vector<char> buffer(local_site_count * compound_type_size);
    char* buffer_ptr = buffer.data();

    // traverse the data source and write to buffer
    dataSource.Reset();
    while(dataSource.ReadNext()) {
      if (outputSpec.geometry->Include(dataSource, dataSource.GetPosition())) {
        // write position
        const auto& pos = dataSource.GetPosition();
        *(reinterpret_cast<uint32_t*>(buffer_ptr + field_offsets[0])) = pos.x();
        *(reinterpret_cast<uint32_t*>(buffer_ptr + field_offsets[1])) = pos.y();
        *(reinterpret_cast<uint32_t*>(buffer_ptr + field_offsets[2])) = pos.z();
        
        // traverse each field and write its data
        for (size_t i = 0; i < outputSpec.fields.size(); ++i) {
          const auto& field_spec = outputSpec.fields[i];
          // +3: former 3 offsets for x,y,z
          char* field_ptr = buffer_ptr + field_offsets[i + 3];

          std::visit([&](auto&& src_type) {
            double* dest = reinterpret_cast<double*>(field_ptr);
            if constexpr (std::is_same_v<std::decay_t<decltype(src_type)>, source::Pressure>) {
              dest[0] = dataSource.GetPressure();
            } else if constexpr (std::is_same_v<std::decay_t<decltype(src_type)>, source::Velocity>) {
              const auto& vel = dataSource.GetVelocity();
              dest[0] = vel.x(); dest[1] = vel.y(); dest[2] = vel.z();
            } else if constexpr (std::is_same_v<std::decay_t<decltype(src_type)>, source::ShearStress>) {
              dest[0] = dataSource.GetShearStress();
            }
          }, field_spec.src);
        }
        
        // next buffer pointer
        buffer_ptr += compound_type_size;
      }
    }

    // memory space for the dataset
    hsize_t count[1] = { local_site_count };
    hid_t memspace_id = H5Screate_simple(1, count, NULL);
    
    hsize_t start[1] = { comms.Scan(local_site_count, MPI_SUM) - local_site_count };
    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, start, NULL, count, NULL);

    // write the data to the dataset, collectively
    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dataset_id, compound_type_id, memspace_id, filespace_id, dxpl, buffer.data());
    
    // clean up
    H5Pclose(dxpl);
    H5Sclose(memspace_id);
  }
}
#endif // USE_HDF5


