// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

// LocalPropertyHdf5Output.cc
#include "extraction/LocalPropertyHdf5Output.h"
#include "net/IOCommunicator.h"
#include <vector>
#include <cassert>
#include <numeric>
#include <stdexcept>
#include <cmath> // check inf
#include <fstream>
#ifdef USE_HDF5

namespace hemelb::extraction 
{
  namespace { 
    struct SiteDataRecord 
    {
      uint32_t x, y, z;
      double pressure;
      double velocity[3];
      double shearstress;
    };
  } // namespace

  LocalPropertyHdf5Output::LocalPropertyHdf5Output(IterableDataSource& dataSource,
                                                   const PropertyOutputFile& outputSpec,
                                                   const net::IOCommunicator& comms)
    : LocalPropertyOutput(dataSource, outputSpec, comms), mpi_comm(comms)
  {
    // create compound type in Write
  }

  LocalPropertyHdf5Output::~LocalPropertyHdf5Output()
  {
    if (dataset_id >= 0) H5Dclose(dataset_id);
    if (filespace_id >= 0) H5Sclose(filespace_id);
    if (compound_type_id >= 0) H5Tclose(compound_type_id);
    if (file_id >= 0) H5Fclose(file_id);
  }

  void LocalPropertyHdf5Output::WriteXDMFFile(unsigned long timestep)
  {
    // rank=0
    if (comms.Rank() != 0) {
      return; // only rank 0 writes the XDMF file
    }
    std::string xmf_filename = outputSpec.filename;
    size_t pos = xmf_filename.rfind(".h5");
    if (pos != std::string::npos){
      xmf_filename.replace(pos,3,".xdmf");
    }else {
      xmf_filename += ".xdmf";
    }
    // Open the XDMF file
    std::ofstream xmf_file(xmf_filename);
    xmf_file << "<?xml version=\"1.0\" ?>\n";
    xmf_file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    xmf_file << "<Xdmf Version=\"3.0\">\n";
    xmf_file << "  <Domain>\n";
    xmf_file << "    <Grid Name=\"HemeLB Data\" GridType=\"Uniform\">\n";

    // Write the topology
    xmf_file << "      <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << global_site_count << "\"/>\n";

    // Write the geometry
    xmf_file << "      <Geometry GeometryType=\"X_Y_Z\">\n";
    xmf_file << "        <DataItem Format=\"HDF\" Dimensions=\"" << global_site_count << " 3\" DataType=\"UInt\">";
    xmf_file << outputSpec.filename << ":/results/i</DataItem>\n";
    xmf_file << "      </Geometry>\n";

    // Write the attributes
    for (const auto& field : outputSpec.fields) {
        unsigned field_len = GetFieldLength(field.src);
        std::string type = (field_len > 1) ? "Vector" : "Scalar";
        xmf_file << "      <Attribute Name=\"" << field.name << "\" AttributeType=\"" << type << "\" Center=\"Node\">\n";
        xmf_file << "        <DataItem Format=\"HDF\" Dimensions=\"" << global_site_count << " " << field_len << "\" DataType=\"Float\" Precision=\"8\">";
        xmf_file << outputSpec.filename << ":/results/" << field.name << "</DataItem>\n";
        xmf_file << "      </Attribute>\n";
    }

    xmf_file << "    </Grid>\n";
    xmf_file << "  </Domain>\n";
    xmf_file << "</Xdmf>\n";
    xmf_file.close();
  }

  void LocalPropertyHdf5Output::CreateCompoundType() {
    // Calculate memory layout 
    size_t current_offset = 0;
    field_offsets.clear();

    // three fields for coordinates (i, j, k)
    field_offsets.push_back(current_offset); // i offset
    current_offset += sizeof(uint32_t);
    field_offsets.push_back(current_offset); // j offset
    current_offset += sizeof(uint32_t);
    field_offsets.push_back(current_offset); // k offset
    current_offset += sizeof(uint32_t);
    
    // Calculate other fields
    for (const auto& field_spec : outputSpec.fields) {
      unsigned field_len = GetFieldLength(field_spec.src);
      field_offsets.push_back(current_offset);
      current_offset += sizeof(double) * field_len;
    }
    // total size
    compound_type_size = current_offset;

    // Create an HDF5 composite type using the computed layout
    compound_type_id = H5Tcreate(H5T_COMPOUND, compound_type_size); 

    // insert i, j, k coordinates
    H5Tinsert(compound_type_id, "i", field_offsets[0], H5T_NATIVE_UINT32);
    H5Tinsert(compound_type_id, "j", field_offsets[1], H5T_NATIVE_UINT32);
    H5Tinsert(compound_type_id, "k", field_offsets[2], H5T_NATIVE_UINT32);
    
    // insert other fields
    for (size_t i = 0; i < outputSpec.fields.size(); ++i) {
      const auto& field_spec = outputSpec.fields[i];
      unsigned field_len = GetFieldLength(field_spec.src);
      hid_t field_type = H5T_NATIVE_DOUBLE;
      
      // +3 to skip the first three offsets for i, j, k
      size_t member_offset = field_offsets[i + 3];

      if (field_len == 1) { // scalar
        H5Tinsert(compound_type_id, field_spec.name.c_str(), member_offset, field_type);
      } else { // array
        hsize_t dims[1] = {field_len};
        hid_t array_tid = H5Tarray_create(field_type, 1, dims);
        H5Tinsert(compound_type_id, field_spec.name.c_str(), member_offset, array_tid);
        H5Tclose(array_tid);
      }
    }
  }

  void LocalPropertyHdf5Output::Write(unsigned long timestepNumber, unsigned long totalSteps)
  {
    if (!ShouldWrite(timestepNumber)) { 
        return;
    }
    // Check if the file is already open
    if (file_id < 0) {
        hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(fapl, mpi_comm, MPI_INFO_NULL);
        file_id = H5Fcreate(outputSpec.filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
        H5Pclose(fapl);
        
        CreateCompoundType();

        hsize_t global_dims[1] = { global_site_count }; // global_site_count from the base class
        filespace_id = H5Screate_simple(1, global_dims, NULL);
        dataset_id = H5Dcreate2(file_id, "results", compound_type_id, filespace_id,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    // create a char buffer for storing the data
    std::vector<char> buffer(local_site_count * compound_type_size); // local_site_count from the base class
    char* buffer_ptr = buffer.data();
    dataSource.Reset(); 
    while(dataSource.ReadNext()) {
      if (outputSpec.geometry->Include(dataSource, dataSource.GetPosition())) {
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
            // check inf and replace it
            if constexpr (std::is_same_v<std::decay_t<decltype(src_type)>, source::Pressure>) {
              double pressure = dataSource.GetPressure();
              dest[0] = std::isinf(pressure) ? -1.0 : pressure; // if inf, set to -1.0
            }
            else if constexpr (std::is_same_v<std::decay_t<decltype(src_type)>, source::Velocity>) {
              const auto& vel = dataSource.GetVelocity();
              dest[0] = std::isinf(vel.x()) ? 0.0 : vel.x();
              dest[1] = std::isinf(vel.y()) ? 0.0 : vel.y();
              dest[2] = std::isinf(vel.z()) ? 0.0 : vel.z();
            }
            else if constexpr (std::is_same_v<std::decay_t<decltype(src_type)>, source::ShearStress>) {
              double stress = dataSource.GetShearStress();
              dest[0] = std::isinf(stress) ? -1.0 : stress;
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
    
    // Write the data to the dataset
    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dataset_id, compound_type_id, memspace_id, filespace_id, dxpl, buffer.data());
    
    //clean up
    H5Pclose(dxpl);
    H5Sclose(memspace_id); 

    // Write XDMF file
    WriteXDMFFile(timestepNumber);   
  }
}
#endif // USE_HDF5


